#ifndef _NS_H_
#define _NS_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "../utility/Clock.h"
#include "../io/mmio.h"

#include <vector>
#include <memory>

//#define __Debug__

#define __x__  -0.4
#define __y__  0.45

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class NS_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef NS_<cvt, vt, Dim> Self;
	typedef NS_<cvt, vt, Dim>& ref_Self;
	typedef NS_<cvt, vt, Dim>* pSelf;
	typedef Domain_<cvt, vt, Dim> Domain;
	typedef Domain_<cvt, vt, Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, Dim>& const_ref_Domain;
	typedef Domain_<cvt, vt, Dim>* pDomain;
	typedef const Domain_<cvt, vt, Dim>* const_pDomain;
	typedef BoundaryCondition_<cvt, vt> BoundaryCondition;
	typedef BoundaryCondition_<cvt, vt>* pBoundaryCondition;
	typedef const BoundaryCondition_<cvt, vt>* const_pBoundaryCondition;

	typedef Grid_<cvt, vt, Dim> Grid;
	typedef Grid_<cvt, vt, Dim> *pGrid;
	typedef const Grid_<cvt, vt, Dim> * const_pGrid;
	typedef Ghost_<cvt, vt, Dim> Ghost;
	typedef const Ghost_<cvt, vt, Dim> const_Ghost;
	typedef Ghost_<cvt, vt, Dim>* pGhost;
	typedef const Ghost_<cvt, vt, Dim>* const_pGhost;
	typedef typename Ghost::GhostNode GhostNode;
	typedef Cell_<cvt, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<vt, Dim> Data;
	typedef Data *pData;
	typedef Node_<cvt, vt, Dim> Node;
	typedef Node_<cvt, vt, Dim> *pNode;
	typedef const Node_<cvt, vt, Dim> *const_pNode;

	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode> *pFace;
	typedef const Face_<Node, pNode> * const_pFace;
	typedef Shape_<cvt, Dim> Shape;
	typedef Shape_<cvt, Dim>* pShape;

	typedef Stencil_<cvt, vt, Dim, Dim> Stencil;

	typedef Interpolate_<cvt, vt, Dim> Interpolate;

	typedef std::function<vt(cvt, cvt, cvt)> Function;
	//typedef std::function<void(pSelf)> Event;

	typedef Expression_<cvt, vt, Dim> Exp;
	typedef Exp* pExp;
	typedef std::shared_ptr<Exp> spExp;
	typedef std::shared_ptr<Face> spFace;

	typedef typename Exp::Term Term;
	typedef typename Exp::iterator ExpIter;
	typedef typename Exp::const_iterator const_ExpIter;

	typedef std::function<vt(vt)> Limiter;
protected:
	typedef MatrixSCR_<vt> Mat;
	typedef ArrayListV<vt> Arr;

	typedef ArrayListT<spExp> Arr_spExp;

	struct comp_spFace {
		bool operator()(const spFace& lhs, const spFace& rhs) {
			// just compare pface
			if (lhs->pori() < rhs->pori()) {
				return true;
			} else if (lhs->pori() == rhs->pori()) {
				if (lhs->pnei() < rhs->pnei()) {
					return true;
				} else if (lhs->pnei() == rhs->pnei()) {
					if (lhs->dir() < rhs->dir()) {
						return true;
					} else {
						return false;
					}
				} else {
					return false;
				}
			} else {
				return false;
			}
		}
	};
	typedef std::map<spFace, Arr_spExp, comp_spFace> Map_FE;
	typedef std::map<spFace, Arr_spExp, comp_spFace>* pMap_FE;
	typedef typename std::map<spFace, Arr_spExp, comp_spFace>::iterator Iter_FE;
	typedef typename std::map<spFace, Arr_spExp, comp_spFace>::const_iterator const_Iter_FE;
	typedef std::pair<spFace, Arr_spExp> Pair_FE;
	typedef std::pair<Iter_FE, bool> Ret_FE;

	void _new_map_fe(pNode pn) {
		utPointer& utp = pn->utp(_utp_fe);
		if (utp == nullptr) {
			utp = new Map_FE();

		} else {
			Map_FE& mfe = CAST_REF(pMap_FE, utp);
			mfe.clear();
		}
	}

	void _delete_map_fe(pNode pn) {
		utPointer& utp = pn->utp(_utp_fe);
		if (utp != nullptr) {
			pMap_FE mfe = CAST(pMap_FE, utp);
			mfe->clear();
			delete mfe;
			utp = nullptr;
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void _face_val_gra(Face& pf, Exp& exp_val, Exp& exp_gra) {
		// face direction
		pNode pori = pf.pori();
		Direction dir = pf.dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);

		// interpolate on face
		exp_val = Interpolate::OnFace(pori, dir, 1);
		exp_gra = Interpolate::OnFace_Gradient(pori, dir, 1);
	}

	void _insert_face_val_gra(pNode pn) {
		utPointer utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		for (st i = 0; i < NumFaces; i++) {
			Direction dir = FaceDirectionInOrder(i);
			pNode pnei = pn->get_neighbor_fast(dir);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face
			FaceType ft = GetFaceType(pn, pnei);
			Arr_spExp arr_exp(2);
			for (st ii = 0; ii < arr_exp.size(); ++ii) {
				spExp sp(new Exp());
				arr_exp[ii] = sp;
			}
			spFace spface(new Face(pn, pnei, dir, ft));
			if ((ft == _Boundary_)				//1
			|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
					|| (ft == _FineCoarse_))				//3
					{
				// work on spface
				//
				_face_val_gra((*spface), (*(arr_exp[0])), (*(arr_exp[1])));
				//
				mfe.insert(Pair_FE(spface, arr_exp));
			}
			// case 2 3
			if ((ft == _Equal_ && (IsFacePDirection(dir)))				//2
			|| (ft == _FineCoarse_))				//3
					{
				// work on pnei
				utPointer utp_nei = pnei->utp(_utp_fe);
				ASSERT(utp_nei != nullptr);
				Map_FE& mfe_nei = CAST_REF(pMap_FE, utp_nei);
				FaceType oft = (ft == _FineCoarse_) ? _CoarseFine_ : _Equal_;
				spFace spface_nei(new Face(pnei, pn, Opposite(dir), oft));
				mfe_nei.insert(Pair_FE(spface_nei, arr_exp));
			}
		}
	}

	void _build_FE_on_leaf() {
		Grid& grid = _pdomain->grid();
		// new map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_new_map_fe(pn);
		}
		// inert face value and gradient
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_insert_face_val_gra(pn);
		}
	}

	void _delete_FE_on_leaf() {
		Grid& grid = _pdomain->grid();
		// delete map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_delete_map_fe(pn);
		}
	}

	/*
	 * Input: 1  face
	 *        2  expression value on face
	 *        3  expression gradient on face
	 *        4  constant index
	 *           if index < 0  the constant is 1
	 *           if index >=0  get the constant from node
	 * Out  : 5  expression of diffusion on face
	 */
	void _face_dif(Face& f, const Exp& exp_val, const Exp& exp_gra,
			const int& idx, Exp& exp_dif) {
		// mu * ▽(phi(x,y))
		vt mu_f = 1.0;
		if (idx >= 0) {
			mu_f = exp_val.substitute(st(idx));
		}
		int sign = IsFacePDirection(f.dir()) ? 1 : -1;
		exp_dif.clear();
		exp_dif = exp_gra;
		exp_dif.times(sign * mu_f);
#ifdef __Debug__
		pNode pn = f.pori();
		if (pn->is_in_on(__x__, __y__)) {
			std::cout << "_face_dif_  ========================\n";
			std::cout << "Direction   " << ToString(f.dir()) << "\n";
			std::cout << " sign        " << sign << "\n";
			std::cout << " mu_f        " << mu_f << "\n";
			std::cout << " coe         " << sign * mu_f << "\n";
			exp_gra.show();
		}
#endif
	}

	void _face_adv(Face& f, const Exp& exp_val, Exp& exp_adv) {
		//const Direction& dir = f.dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(f.dir(), o, a);

		// get velocity on face
		vt veo_f = exp_val.substitute(_v_idx[a]);

		exp_adv.clear();

		// times area
		int sign = IsP(o) ? 1 : -1;
		vt area = f.area();
		// first order up wind -----------------------------
		pNode pC = _find_C(&f, veo_f);
		exp_adv.insert(sign * area * veo_f, pC, 1.0);
		// center difference -------------------------------
		//exp_adv.plus(exp_val);
		//exp_adv.times(sign * area * veo_f);
#ifdef __Debug__
		pNode pn = f.pori();
		if (pn->is_in_on(__x__, __y__)) {
			std::cout << "_face_adv_  ========================\n";
			std::cout << " Direction   " << ToString(dir) << "\n";
			std::cout << " sign        " << sign << "\n";
			std::cout << " area        " << area << "\n";
			std::cout << " coe         " << sign * area * veo_f << "\n";
		}
#endif

	}

	void _face_exp_adv_dif(Iter_FE& iter, Axes ai) {
		const spFace spf = iter->first;
		Arr_spExp& arr_spexp = iter->second;
		spExp spexp_adv(new Exp());
		spExp spexp_dif(new Exp());
		_face_adv(*spf, *(arr_spexp[0]), *spexp_adv);
		_face_dif(*spf, *(arr_spexp[0]), *(arr_spexp[1]), _mu_idx, *spexp_dif);
		//
		arr_spexp.resize(4);
		arr_spexp[2] = spexp_adv;
		arr_spexp[3] = spexp_dif;
	}

	void _node_face_exp_adv_dif(pNode pn, Axes ai) {
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			_face_exp_adv_dif(iter, ai);
		}
	}

	/*
	 * get advection term on node
	 */
	void _node_exp_adv_term(pNode pn, Exp& exp) {
		//assert(pn->d());
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			const spFace spf = iter->first;
			Arr_spExp& arr_spexp = iter->second;
			exp.plus(*(arr_spexp[2]));
		}
	}

	/*
	 * get advection term on node
	 */
	void _node_exp_dif_term(pNode pn, Exp& exp) {
		//assert(pn->d());
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			const spFace spf = iter->first;
			Arr_spExp& arr_spexp = iter->second;
			Exp texp(*(arr_spexp[3]));
			//
			texp.times(spf->area());
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				std::cout << "node_exp_dif_term =======================\n";
				arr_spexp[3]->show();
				texp.show();
			}
#endif
			exp.plus(texp);
		}
	}

	/*
	 * this is the main function of the predict step
	 *
	 */
	void _node_exp_u_star(pNode pn, Exp& exp) {
		// Traversal all the faces
		// cfl number check
		vt cfl = _CFL_number(pn, _dt);
		ASSERT(cfl < 1.0);
		// calculate the term separately
		// 1 ----- advection term
		Exp exp_adv;
		_node_exp_adv_term(pn, exp_adv);
		// divide volume and negative;
		exp_adv.times(-1.0 / pn->volume());

		// 2 ----- diffusion term
		Exp exp_dif;
		_node_exp_dif_term(pn, exp_dif);
#ifdef __Debug__
		if (pn->is_in_on(__x__, __y__)) {
			std::cout << "u star mid =======================\n";
			exp_dif.show();
		}
#endif
		//divide volume
		exp_dif.times(1.0 / pn->volume());

		// 3 ----- source term
		Exp exp_src;
		//_node_exp_src_term(pn, ai, exp_src);
		//   add up diffusion and source
		exp_dif.plus(exp_src);
		exp_dif.times(1.0 / pn->cdva(_rho_idx));
		exp_dif.plus(exp_adv);   // add advection term
		exp_dif.times(_dt);      // time
		// 4 u term
		exp.insert(1.0, pn, 1.0);
		exp.plus(exp_dif);
	}

	void _u_star_on_center(Axes ai, vt dt) {
#ifdef __Debug__
		debug_log << "->_u_star_on_center at line " << __LINE__ << "\n";
#endif
		Clock clock;
		pGrid pgrid = this->_pdomain->p_grid();
		//
		//1 Traverse face
		//  Build up the expression on utp
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_node_face_exp_adv_dif(pn, ai);
		}
		clock.break_point("> face exp adv and dif");
		//2 Traverse face to calculate the results
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			Exp exp;
			pNode pn = it.get_pointer();
			_node_exp_u_star(pn, exp);
			//exp.show();
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				debug_log << "  _u_star_on_center: after _node_exp_u_star\n";
				debug_log << "   pn " << "(" << pn->cp(_X_) << ", "
				<< pn->cp(_Y_) << ")" << " dx = " << pn->d(_X_) << "\n";
				exp.show(debug_log);
				//exp.show_substitute(_v_idx[ai]);
			}
#endif
			_exp_substitute_ghost(pn, exp, _v_idx[ai]);
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				debug_log
				<< "  _u_star_on_center: after _exp_substitute_ghost\n";
				exp.show(debug_log);
			}
#endif
			vt u_star = exp.substitute(_v_idx[ai]);
			pn->cd(_vt_idx[ai]) = u_star;
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				debug_log << "  _u_star_on_center: after substitute\n";
				debug_log << "  Axes = " << ToString(ai) << " ,u_star = "
				<< u_star << "\n";
			}
#endif
		}
		clock.break_point("> 2 loop subsitute");
		//clock.show();
#ifdef __Debug__
		debug_log << "<-_u_star_on_center" << "\n";
#endif

	}

	void _u_star_on_ghost(Axes ai, vt dt) {
#ifdef __Debug__
		debug_log << "->_u_star_on_ghost at line " << __LINE__ << "\n";
#endif
		typedef typename Ghost::GhostNode GhostNode;
		std::function<void(GhostNode&)> _fun = [this, &ai](GhostNode& node) {
			pNode pg = node.second.pghost;
			if (pg != nullptr) {
				// set ustart on boundary
				pg->cd(_vt_idx[ai]) = 0.0;
			}
		};
		this->_pdomain->p_ghost()->for_each_node(_fun);

#ifdef __Debug__
		debug_log << "<-_u_star_on_ghost" << "\n";
#endif

	}

	/*
	 * solve the pressure on center of node
	 */
	void _face_pressure(const Face& f, const Exp& exp_val, const Exp& exp_gra,
			Exp& exp) {
		vt beta = _dt / exp_val.substitute(_rho_idx);
		int sign = IsFacePDirection(f.dir()) ? 1 : -1;
		exp.clear();
		exp = exp_gra;
		vt area = f.area();
		exp.times(sign * beta * area);
	}

	void _face_ustar(Face& f, const Exp& exp_val, Exp& exp) {
		pNode pn = f.pori();
		const Direction& dir = f.dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(f.dir(), o, a);

		// get velocity on face
		//vt veo_f = exp_val.substitute(_vt_idx[a]);
		//exp.insert(veo_f, pn, 0.0);
		exp = exp_val;
		// times area
		int sign = IsFacePDirection(dir) ? 1 : -1;
		vt area = f.area();
		exp.times(sign * area);
	}

	void _face_exp_pressure_ustar(Iter_FE& iter) {
		const spFace spf = iter->first;
		Arr_spExp& arr_spexp = iter->second;
		spExp spexp_ustart(new Exp());
		spExp spexp_pressure(new Exp());
		_face_ustar(*spf, *(arr_spexp[0]), *spexp_ustart);
		_face_pressure(*spf, *(arr_spexp[0]), *(arr_spexp[1]), *spexp_pressure);
		//
		//arr_spexp.resize(4);
		arr_spexp[2] = spexp_ustart;
		arr_spexp[3] = spexp_pressure;
	}

	void _node_face_exp_pressure_ustar(pNode pn) {
		//assert(pn->d());
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			_face_exp_pressure_ustar(iter);
		}
	}

	void _node_exp_pressure_equation(pNode pn, Exp& exp) {
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		vt sum_ustar = 0;
		vt sumd = 0;
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			const spFace& spf = iter->first;
			Arr_spExp& arr_spexp = iter->second;
			Exp texp(*(arr_spexp[3]));
			//texp.times(spf->area());
			exp.plus(texp);
			//
			Orientation o;
			Axes a;
			FaceDirectionToOrientationAndAxes(spf->dir(), o, a);
			sum_ustar += arr_spexp[2]->substitute(_vt_idx[a]);

#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				std::cout << ToString(spf->dir()) << "======\n";
				std::cout << "star f  " << arr_spexp[2]->substitute(_vt_idx[a])
				<< "\n";
				//std::cout<< "p   f  \n";
				//texp.show();
				//arr_spexp[2]->show();
				//arr_spexp[2]->show_substitute(_vt_idx[a]);
			}
#endif
		}
		exp.plus(-sum_ustar, pn, 0);
		//exp.times(1.0 / pn->volume());

#ifdef __Debug__
		if (pn->is_in_on(__x__, __y__)) {
			std::cout << "sum    " << sum_ustar << "\n";
		}
#endif
	}

	void _build_pressure_equation(Mat& mat, Arr& b) {
#ifdef __Debug__
		debug_log << "->_build_pressure_equation at line " << __LINE__ << "\n";
#endif
		//
		//1 Traverse face
		//  Build up the expression on utp
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_node_face_exp_pressure_ustar(pn);
		}
		//2 Treverse face
		//  to calculate the expression for solver
		typedef std::list<st> ListST;
		typedef std::list<vt> ListVT;
		ListST l_rowptr;
		l_rowptr.push_back(0);
		ListST l_colid;
		ListVT l_val;
		ListVT l_b;
		int countnz = 0;

		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			Exp exp;
			pNode pn = it.get_pointer();
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {

			}
#endif
			_node_exp_pressure_equation(pn, exp);
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				debug_log << "  _build_pressure_equation: after _node_exp_pressure_equation\n";
				exp.show(debug_log);
				//_show_exp(exp, *_pdomain, pn);
			}
#endif
			_exp_substitute_ghost_pressure(pn, exp);
			exp.concise();
#ifdef __Debug__
			if (pn->is_in_on(__x__, __y__)) {
				std::cout << "pressure sub =======================\n";
				exp.show_substitute(_p_idx);
				//_show_exp(exp, *_pdomain, pn);
			}
#endif

			int fconst = 0;
			for (typename Exp::iterator ite = exp.begin(); ite != exp.end();
					++ite) {
				if (!Exp::IsConstant(ite)) {
					const_pNode pn = Exp::GetpNode(ite);
					vt val = Exp::GetCoe(ite);
					l_colid.push_back(pn->d_idx());
					l_val.push_back(val);
					countnz++;
				} else {
					fconst = 1;
					vt val = Exp::GetCoe(ite);
					l_b.push_back(-val);    //!!!!! negative added here
				}
			}
			if (fconst == 0) {
				l_b.push_back(0);
			}
			l_rowptr.push_back(countnz);
			//exp.show();
		}
		//copy list to array  ======================
		st nr = l_rowptr.size() - 1;
		st nz = l_val.size();
		ASSERT(nz == l_colid.size());
		ASSERT(nr <= nz);
		mat.newsize(nr, nr, nz);
		b.reconstruct(l_b.size());
		int i = 0;
		for (typename ListST::iterator it = l_colid.begin();
				it != l_colid.end(); ++it) {
			mat.col_ind(i) = (*it);
			i++;
		}
		i = 0;
		for (typename ListST::iterator it = l_rowptr.begin();
				it != l_rowptr.end(); ++it) {
			mat.row_ptr(i) = (*it);
			i++;
		}
		i = 0;
		for (typename ListVT::iterator it = l_val.begin(); it != l_val.end();
				++it) {
			mat.val(i) = (*it);
			i++;
		}
		i = 0;
		for (typename ListVT::iterator it = l_b.begin(); it != l_b.end();
				++it) {
			b[i] = (*it);
			i++;
		}
#ifdef __Debug__
		debug_log << "<-_build_pressure_equation " << "\n";
#endif
	}

	void _node_correct_veo_on_center(pNode pn) {
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		Arr sum_p(Dim);
		Arr sum_a(Dim);
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			const spFace& spf = iter->first;
			Arr_spExp& arr_spexp = iter->second;
			//
			Orientation o;
			Axes a;
			FaceDirectionToOrientationAndAxes(spf->dir(), o, a);
			// pressure gradient on face
			vt pressure_gra_f = arr_spexp[1]->substitute(_p_idx);
			vt rho_f = arr_spexp[0]->substitute(_rho_idx);
			sum_p[a] += rho_f * pressure_gra_f * _dt * spf->area();
			sum_a[a] += spf->area();
		}
		for (st di = 0; di < Dim; ++di) {
			vt av_p = sum_p[di] / sum_a[di];
			pn->cd(_v_idx[di]) = pn->cd(_vt_idx[di]) - av_p;
		}
	}

	void _node_correct_veo_on_center2(pNode pn) {
		utPointer& utp = pn->utp(_utp_fe);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);

		for (st di = 0; di < Dim; ++di) {
			vt sum_pp = 0, sum_pm = 0;
			vt sum_ap = 0, sum_am = 0;
			for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
				const spFace& spf = iter->first;
				Arr_spExp& arr_spexp = iter->second;
				//
				Orientation o;
				Axes a;
				FaceDirectionToOrientationAndAxes(spf->dir(), o, a);
				if (a == ToAxes(di)) {
					if (o == _M_) {
						vt pressure_f = arr_spexp[0]->substitute(_p_idx);
						sum_pm += pressure_f * spf->area();
						sum_am += spf->area();
					} else {
						vt pressure_f = arr_spexp[0]->substitute(_p_idx);
						sum_pp += pressure_f * spf->area();
						sum_ap += spf->area();
					}
				}
			}
			vt av_dp = (sum_pp / sum_ap - sum_pm / sum_am) / pn->d(ToAxes(di));
			pn->cd(_v_idx[di]) = pn->cd(_vt_idx[di])
					- av_dp * _dt / pn->cd(_rho_idx);
		}

	}

	void _correct_veo_on_center() {
		//
		//1 Traverse face
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_node_correct_veo_on_center2(pn);
		}
	}
protected:
	//matrx
	typedef ArrayListV<st> Arr_st;

	pDomain _pdomain;

	vt _dt;
	vt _max_step;

	st _scheme;
	std::vector<Limiter> _v_limiter;

	//function list
	//std::list<Event> _levent_start;
	//std::list<Event> _levent_step;
	//std::list<Event> _levent_end;


	st _vt_idx[Dim]; //tmp velo term
	st _v_idx[Dim];  //velo     term
	st _s_idx[Dim];  //source   term
	st _p_idx;
	st _rho_idx;
	st _mu_idx;

	//
	st _utp_fe;

#ifdef __Debug__
	std::fstream debug_log;
#endif
public:
	/*
	 * constructor
	 *   2D
	 *
	 * pressure       0
	 * velocity  u    1
	 * velocity  v    2
	 * velocity  u*   3
	 * velocity  v*   4
	 * source    sx   5
	 * source    sy   6
	 * Density   rho  7
	 * viscocity mu   8
	 *
	 * utp
	 * utp map          0
	 */
	NS_(pDomain pf, Float dt, st scheme = 0) :
			_pdomain(pf), _scheme(scheme) {
		_push_limiter();
		_p_idx = 0;
		for (st i = 0; i < Dim; ++i) {
			_v_idx[i] = i + 1;
			_vt_idx[i] = i + 1 + Dim;
			_s_idx[i] = i + 1 + 2 * Dim;
		}
		_rho_idx = 3 * Dim + 1;
		_mu_idx = 3 * Dim + 2;
		//
		_utp_fe = 0;
		_pdomain->new_data(this->max_idx() + 1, 0, 0, 1);
		// delta t
		_dt = dt;
#ifdef __Debug__
		std::cout << "  Construct NS Class === \n";
		debug_log.open("debug.log", std::fstream::out);
		debug_log << "  Construct NS Class === \n";
		debug_log << "  Debug point " << __x__ << " , " << __y__ << "\n";
#endif
	}

	~NS_() {
#ifdef __Debug__
		debug_log << "  Deconstruct NS Class ~~~ \n";
		debug_log.close();
#endif
	}

	st idx_u() const {
		return _v_idx[0];
	}
	st idx_v() const {
		ASSERT(Dim >= 2);
		return _v_idx[1];
	}
	st idx_w() const {
		ASSERT(Dim >= 3);
		return _v_idx[2];
	}
	st idx_ut() const {
		return _vt_idx[0];
	}
	st idx_vt() const {
		ASSERT(Dim >= 2);
		return _vt_idx[1];
	}
	st idx_wt() const {
		ASSERT(Dim >= 3);
		return _vt_idx[2];
	}
	st idx_p() const {
		return _p_idx;
	}
	st idx_mu() const {
		return _mu_idx;
	}
	st idx_rho() const {
		return _rho_idx;
	}

public:

	/*
	 * build face exp to utp
	 */
	int solve(vt tol = 1e-6, int max_iter = 1000, int info = 0) {
		std::list<vt> lr;
		int rcode = this->solve(tol, max_iter, lr, info);
		//std::cout<< "max iter "<<max_iter<< " tol "<< tol <<"\n";
		return rcode;
	}

	int solve(vt& tol, int& max_iter, std::list<vt>& lr, int info = 0) {
		Clock clock;
#ifdef __Debug__
		debug_log << "->solve at line " << __LINE__ << "\n";
#endif

#ifdef __Debug__
		debug_log << "  solve: 1 build u star on center\n";
#endif
		for (st i = 0; i < Dim; ++i) {
#ifdef __Debug__
			debug_log << "  solve: on axe" << ToString(ToAxes(i)) << "\n";
#endif
			_u_star_on_center(ToAxes(i), _dt);
		}
		clock.break_point("u* on center");
#ifdef __Debug__
		debug_log << "  solve: 2 build u star on boundary\n";
#endif
		// set u star bc
		for (st i = 0; i < Dim; ++i) {
#ifdef __Debug__
			debug_log << "  solve: on axe "<< ToString(ToAxes(i)) <<"\n";
#endif
			_u_star_on_ghost(ToAxes(i), _dt);
		}
		clock.break_point("u* on ghost");

		Mat mat;
		Arr b;
		_build_pressure_equation(mat, b);
		clock.break_point("build pressure equation");
		Arr x(b.size());
		// initial x
		_grid_to_arr(x, _p_idx);
		clock.break_point("copy array to x");
#ifdef __Debug__

//		b.show();
//		std::list<Gnuplot_actor> lga;
//		Gnuplot_actor ga;
//		GnuplotActor_MatrixSCR(ga, mat);
//		lga.push_back(ga);
//		Gnuplot gp;
//		gp.set_equal_ratio();
//		gp.plot(lga);

#endif
		//std::cout<<" diagnally dominant " << mat.is_diagonally_dominant() <<"\n";
		//mm_write_mtx_sparse("mat", mat);
		int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);

		clock.break_point("solve");
		//int sf = Jacobi(mat, x, b, max_iter, tol, lr);
		if (info == 1) {
			std::cout << "iter n " << lr.size() << " residual " << tol << "\n";
		}
		if (sf != 0) {
			std::cerr << " >! Poisson solve failed \n";
		}
		//x.show();
		_arr_to_grid(x, _p_idx);
		clock.break_point("copy x to grid");
		// set Pressure star bc
		for (st i = 0; i < Dim; ++i) {
			this->_pdomain->set_val_ghost_by_bc(_p_idx);
		}
		clock.break_point("set pressure on ghost");
		_correct_veo_on_center();
		clock.break_point("correct veo");
		if (info == 1) {
			clock.show();
		}
#ifdef __Debug__
		debug_log << "<-solve " << "\n";
#endif
		return 1;
	}

	int advance(st max_step, vt tol = 1e-6, int max_iter = 1000, int info = 0) {
#ifdef __Debug__
		debug_log << "->advance at line " << __LINE__ << "\n";
		debug_log << "  tol " << tol << ", max_iter " << max_iter << ", info "
		<< info << "\n";
#endif
		_build_FE_on_leaf();
		tick_t tick_cpu = Clock::SystemTime();
		tick_t tick_wall = Clock::Tick();
		for (st i = 0; i < max_step; ++i) {
			if (info == 1) {
				fmt::print("step: {:>8} dt: {:>8.5f} cpu: {:>8.5f} wall: {:>8.5f}\n",
						i, i*_dt,
						Clock::TimespanToSecondsD(tick_cpu, Clock::SystemTime()),
						Clock::TimespanToSecondsD(tick_wall, Clock::Tick()));
				tick_cpu = Clock::SystemTime();
				//tick_wall = Clock::Tick();
			}
#ifdef __Debug__
			debug_log << "  advance step = " << i << "\n";
#endif
			solve(tol, max_iter, info);
		}
		_delete_FE_on_leaf();
		return 1;
	}

	/*
	 * Boundary condition
	 */
	void set_boundary_condition(st si, st segi, st vi, pBoundaryCondition pbc) {
		_pdomain->_pbindex->insert(si, segi, vi, pbc);
	}
	void set_BC_velocity(st si, st segi, Axes a, pBoundaryCondition pbc) {
		this->set_boundary_condition(si, segi, _v_idx[a], pbc);
		// the boundary condition of velo star is the same as the velo
		//this->set_boundary_condition(si, segi, _vt_idx[a], pbc);
	}
	void set_BC_velocity_t(st si, st segi, Axes a, pBoundaryCondition pbc) {
		this->set_boundary_condition(si, segi, _vt_idx[a], pbc);
		// the boundary condition of velo star is the same as the velo
		//this->set_boundary_condition(si, segi, _vt_idx[a], pbc);
	}

	/*
	 * set
	 */
	void set_v_grid(Function pfun, Axes a) {
		ASSERT(a < Dim);
		_pdomain->set_val_grid(_v_idx[a], pfun);
	}

	void set_v_ghost(Axes a) {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_v_idx[a]);
	}
	void set_velo(Function pfun, Axes a) {
		this->set_v_grid(pfun, a);
		this->set_v_ghost(a);
	}
	void set_source_grid(Function pfun, Axes a) {
		ASSERT(a < Dim);
		_pdomain->set_val_grid(_v_idx[a], pfun);
	}

	void set_source_ghost(Function pfun, Axes a) {
		_pdomain->set_val_ghost(_v_idx[a], pfun);
	}
	void set_source(Function pfun, Axes a) {
		this->set_source_grid(pfun, a);
		this->set_source_ghost(pfun, a);
	}

	void set_rho_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_rho_idx);
	}
	void set_rho_grid(Function pfun) {
		_pdomain->set_val(_rho_idx, pfun);
	}
	void set_rho(Function pfun) {
		this->set_rho_grid(pfun);
		this->set_rho_ghost();
	}
	void set_mu_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_mu_idx);
	}
	// can be used to set initial boundary condition
	void set_mu_grid(Function pfun) {
		_pdomain->set_val(_mu_idx, pfun);
	}
	void set_mu(Function pfun) {
		this->set_mu_grid(pfun);
		this->set_mu_ghost();
	}

	/*
	 * this function returns the max index of the valuables in
	 * ns, the max index makes sure the array size is ok.
	 */
	st max_idx() {
		st max = 0;
		st m1 = Max(_v_idx, Dim);
		st m2 = Max(_vt_idx, Dim);
		st m3 = Max(_s_idx, Dim);
		max = Max(m1, m2, m3);
		max = Max(max, _rho_idx);
		max = Max(max, _mu_idx);
		max = Max(max, _p_idx);
		return max;
	}

protected:

	void _push_limiter() {
		Limiter fou = [](vt r) {return 0;};

		Limiter minmod = [](vt r) {
			return Max(0.0, Min(1.0, r));
		};

		Limiter superbee = [](vt r) {
			return Max(0.0, Min(2.0 * r, 1.0), Min(r, 2.0));
		};

		Limiter vanLeer = [](vt r) {
			return (r + Abs(r)) / (1 + r);
		};

		_v_limiter.push_back(fou);
		_v_limiter.push_back(minmod);
		_v_limiter.push_back(superbee);
		_v_limiter.push_back(vanLeer);
	}

	/*
	 *
	 */
	const_pBoundaryCondition _find_bc(const_pNode pg, st vali) {
		//input a ghost node
		BoundaryCondition res;
		typename Ghost::GhostID gid = Ghost::ToGhostID(pg);
		// find in map ghost
		typename Ghost::iterator iter = this->_pdomain->p_ghost()->find(gid);
		ASSERT(iter != this->_pdomain->p_ghost()->end());
		// find in BoundaryIndex
		st si = iter->second.shape_idx;
		st segi = iter->second.seg_idx;
		return this->_pdomain->find_bc(si, segi, vali);
	}

	int _exp_substitute_ghost(pNode pn, Exp& exp, st vali) {
		for (ExpIter iter = exp.begin(); iter != exp.end();) {
			if (Exp::IsGhostNode(iter)) {
				ExpIter itere = iter;
				++iter;
				// find boundary condition
				const_pNode pg = Exp::GetpNode(itere);
				if (!Exp::IsZeroCoe(itere)) { // if coe = 0 the term will be deleted
					const_pBoundaryCondition pbc = _find_bc(pg, vali);
					if (pbc->get_type() == BoundaryCondition::_BC1_) {
						// Boundary condition 1
						vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
								pg->cp(_Z_));
						//vt vol = pg->volume();
						vt coe = Exp::GetCoe(itere);
						exp.insert(coe * val, pn, 0);
					} else { // bc2
						//
						//SHOULD_NOT_REACH;
						typename Ghost::GhostID gid = Ghost::ToGhostID(pg);
						ASSERT(gid.step == 0);
						const_pNode po = pg->father;
						Axes a = FaceDirectionToAxes(gid.direction);
						vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
								pg->cp(_Z_), a);
						cvt dl = (pg->cp(a) - po->cp(a));
						vt coe = Exp::GetCoe(itere);
						// pg = po - val*dl
						exp.insert(coe, po, 1);
						exp.insert(-val * dl * coe, po, 0);			//constant
						//SHOULD_NOT_REACH;
						//vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
						//		pg->cp(_Z_));
						//exp.insert(Term(val, pn, 0));
					}
				}
				exp.erase(itere);
			} else {
				++iter;
			}
		}
		exp.concise();
		return 1;
	}
	int _exp_substitute_ghost_pressure(pNode pn, Exp& exp) {
		for (ExpIter iter = exp.begin(); iter != exp.end();) {
			if (Exp::IsGhostNode(iter)) {
				ExpIter itere = iter;
				++iter;
				// find boundary condition
				const_pNode pg = Exp::GetpNode(itere);
				if (!Exp::IsZeroCoe(itere)) { // if coe = 0 the term will be deleted
					typename Ghost::GhostID gid = Ghost::ToGhostID(pg);
					ASSERT(gid.step == 0);
					const_pNode po = pg->father;
					Axes a;
					Orientation o;
					FaceDirectionToOrientationAndAxes(gid.direction, o, a);
					// dp = u* - u;
					// vt sign = IsP(o)? 1:0;
					//vt val = (pg->cd(_vt_idx[a]) - pg->cd(_v_idx[a]));
					vt val = (pg->cd(_vt_idx[a]));
					cvt dl = (po->cp(a) - pg->cp(a));
					vt coe = Exp::GetCoe(itere);
					// pg = po - val*dl
					exp.insert(coe, po, 1);
					exp.insert(val * dl * coe, po, 0);			//constant
					//SHOULD_NOT_REACH;
					//vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
					//		pg->cp(_Z_));
					//exp.insert(Term(val, pn, 0));
				}
				exp.erase(itere);
			} else {
				++iter;
			}
		}
		exp.concise();
		return 1;
	}

	vt _CFL_number(pNode pn, vt dt) {
		vt veo[Dim];
		vt cfl[Dim];
		for (st i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(_v_idx[i]);
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
	}

	pNode _find_C(pFace pface, Float veo_f) {
		pNode pori = pface->pori();
		pNode pnei = pface->pnei();
		Direction dir = pface->dir();
		if (IsFacePDirection(dir)) {
			if (veo_f > 0) {
				return pori;  //o
			} else {
				return pnei;
			}
		} else {
			if (veo_f > 0) {
				return pnei;
			} else {
				return pori;
			}
		}
		SHOULD_NOT_REACH;
		return nullptr;
	}
	void _arr_to_grid(const Arr& x, st idx) {
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			it->cd(idx) = x[it->d_idx()];
		}
	}
	void _grid_to_arr(Arr& x, st idx) {
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			x[it->d_idx()] = it->cd(idx);
		}
	}

#ifdef __Debug__
	void _show_exp(Exp& e, Domain& domain, pNode pn) {
		Exp exp(e);
		exp.concise();
		exp.show();
		std::list<Gnuplot_actor> lga;
		Gnuplot_actor ga;
		GnuplotActor_LeafNodes(ga, domain.grid());
		lga.push_back(ga);
		GnuplotActor_Expression(ga, exp);
		lga.push_back(ga);
		GnuplotActor_Node(ga, *pn);
		lga.push_back(ga);
		Gnuplot gp;
		gp.set_equal_ratio();
		//gp.set_xrange(2.0,3.0);
		//gp.set_yrange(1.5,2.5);
		//gp.set_cbrange(-2.0, 2.0);
		gp.plot(lga);
		//delete shape
	}
#endif

}
;

}

#endif
