#ifndef _NS_H_
#define _NS_H_

#include "calculation_define.hpp"
#include "expression.hpp"

#include <vector>

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

	typedef NS_<COO_VALUE, VALUE, DIM> Self;
	typedef NS_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef Domain_<COO_VALUE, VALUE, DIM> Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>& ref_Domain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>& const_ref_Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>* pDomain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>* const_pDomain;
	typedef typename Domain::BoundaryIndex::BoundaryCondition BoundaryCondition;
	typedef typename Domain::BoundaryIndex::pBoundaryCondition pBoundaryCondition;
	typedef typename Domain::BoundaryIndex::const_pBoundaryCondition const_pBoundaryCondition;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Ghost_<COO_VALUE, VALUE, DIM> Ghost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM> const_Ghost;
	typedef Ghost_<COO_VALUE, VALUE, DIM>* pGhost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>* const_pGhost;
	typedef typename Ghost::GhostNode GhostNode;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> *const_pNode;

	typedef typename Domain::Face Face;
	typedef typename Domain::pFace pFace;
	typedef Shape_<COO_VALUE, DIM> Shape;
	typedef Shape_<COO_VALUE, DIM>* pShape;

	typedef Stencil_<COO_VALUE, VALUE, DIM, DIM> Stencil;

	typedef Interpolate_<COO_VALUE, VALUE, DIM> Interpolate;

	typedef std::function<vt(cvt, cvt, cvt)> Function;

	typedef Expression_<COO_VALUE, VALUE, DIM> Exp;
	typedef Exp* pExp;
	typedef typename Exp::Term Term;
	typedef typename Exp::iterator ExpIter;
	typedef typename Exp::const_iterator const_ExpIter;

	typedef std::function<vt(vt)> Limiter;

protected:
	//matrix
	typedef MatrixSCR_<vt> Mat;
	typedef ArrayListV<vt> Arr;
	typedef ArrayListV<st> Arr_st;

	pDomain _pdomain;

	vt _dt;
	vt _max_step;

	st _scheme;
	std::vector<Limiter> _v_limiter;

	st _vt_idx[Dim]; //tmp velo term
	st _v_idx[Dim];  //velo     term
	st _s_idx[Dim];  //source   term
	st _p_idx;
	st _rho_idx;
	st _mu_idx;

	st _adv_utp_idx[Dim];
	st _dif_utp_idx[Dim];

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
	 * adv utp   x u   0
	 * adv utp   y v   1
	 * dif utp   x u   2
	 * dif utp   y u   3
	 */
	NS_(pDomain pf, st scheme = 0) :
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
		for (st i = 0; i < Dim; ++i) {
			_adv_utp_idx[i] = i;
			_dif_utp_idx[i] = i + Dim;
		}
		_pdomain->new_data(this->max_idx() + 1, 0, 0, 2 * Dim);
	}
	/*
	 * build face exp to utp
	 */
protected:
	typedef std::list<std::pair<pFace, pExp> > ListFE;
	typedef std::pair<pFace, pExp> PairFE;

public:
	/*
	 * predict step
	 */
	int _face_scheme_adv_equal(pFace pface, Exp& exp) {
		// face direction
		pNode pori = pface->pori();
		Direction dir = pface->dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);
		// interpolate on face
		Exp expv = Interpolate::OnFace(pori, dir, 1);
		Float veo_f = expv.subsitute(_v_idx[a]);

		//get U C D ------------------------------
		pNode pC = _find_C(pface, veo_f);
		//
		// exp should be empty before insert
		exp.insert(1.0, pC, 1.0);
		return 1;
	}

	int _face_scheme_dif_equal_gradient(pFace pface, Exp& exp) {
		// face direction
		pNode pori = pface->pori();
		pNode pnei = pface->pnei();
		Direction dir = pface->dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);
		int sign = IsFacePDirection(dir) ? 1 : -1;

		// interpolate beta f
		Exp expb = Interpolate::OnFace(pori, dir, 1);
		Float beta_f = expb.subsitute(_mu_idx);

		Exp expg = Interpolate::OnFace_Gradient(pori, dir, 2);
		expg.times(sign * beta_f);
		//if (pori->is_in_on(0.4875, 0.328)
		//		&& IsFacePDirection(pface->direction)) {
		//	_show_exp(expg, *_pdomain);
		//	std::cout << "stop\n";
		//}
		exp.plus(expg);
		return 1;
	}
	/*
	 * face scheme advection and diffusion
	 * -A  and  D
	 */
	void _face_scheme_adv_dif(pFace pf, Exp& exp_adv, Exp& exp_dif) {
		// face direction
		pNode pori = pf->pori();
		Direction dir = pf->dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);
		int sign = IsFacePDirection(dir) ? 1 : -1;

		// interpolate on face
		Exp expf = Interpolate::OnFace(pori, dir, 1);
		vt veo_f = expf.subsitute(_v_idx[a]);
		vt mu_f = expf.subsitute(_mu_idx);

		// advection =============================
		pNode pC = _find_C(pf, veo_f);
		exp_adv.insert(1.0, pC, 1.0);

		// diffusion =============================
		Exp expg = Interpolate::OnFace_Gradient(pori, dir, 2);
		expg.times(sign * mu_f);
		//if (pori->is_in_on(0.4875, 0.328)
		//		&& IsFacePDirection(pface->direction)) {
		//	_show_exp(expg, *_pdomain);
		//	std::cout << "stop\n";
		//}
		exp_dif.plus(expg);
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

	void __push_listFE(pFace pf, pExp pe, utPointer utp) {
		// utp ---> ListFE
		ListFE& lpexp = CAST_REF(ListFE*, utp);
		PairFE pfe(pf, pe);
		lpexp.push_back(pfe);
	}
	/*
	 * face exp adv dif
	 */
	void _face_exp_adv_dif(pNode pn, Axes ai) {
		st aui = _adv_utp_idx[ai];
		st dui = _dif_utp_idx[ai];
		_new_list_face_exp(pn, aui);
		_new_list_face_exp(pn, dui);
		utPointer& utp_adv = pn->utp(aui);
		utPointer& utp_dif = pn->utp(dui);
		for (st i = 0; i < NumFaces; i++) {
			Direction dir = FaceDirectionInOrder(i);
			pNode pnei = pn->get_neighbor_fast(dir);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face
			FaceType ft = GetFaceType(pn, pnei);
			pExp pexp_adv = new Exp();
			pExp pexp_dif = new Exp();
			bool flag = true;
			if ((ft == _Boundary_)				//1
			|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
					|| (ft == _FineCoarse_))				//3
					{
				// work on pn
				pFace pf = new Face(pn, pnei, dir, ft);
				//
				_face_scheme_adv_dif(pf, *pexp_adv, *pexp_dif);
				//
				__push_listFE(pf, pexp_adv, utp_adv);
				__push_listFE(pf, pexp_dif, utp_dif);
				flag = false;
			}
			// case 2 3
			if ((ft == _Equal_ && (IsFacePDirection(dir)))				//2
			|| (ft == _FineCoarse_))				//3
					{
				// work on pnei
				_new_list_face_exp(pnei, aui);
				_new_list_face_exp(pnei, dui);
				utPointer& utpn_adv = pnei->utp(aui);
				utPointer& utpn_dif = pnei->utp(dui);
				FaceType oft = (ft == _FineCoarse_) ? _CoarseFine_ : _Equal_;
				pFace pfn = new Face(pnei, pn, Opposite(dir), oft);
				pExp pexpn_adv = new Exp(*pexp_adv);
				pExp pexpn_dif = new Exp(*pexp_dif);
				pexpn_dif->times(-1.0);
				__push_listFE(pfn, pexpn_adv, utpn_adv);
				__push_listFE(pfn, pexpn_dif, utpn_dif);
			}
			if (flag) {
				delete pexp_adv;
				delete pexp_dif;
			}
		}
	}
	/*
	 * get advection term on node
	 */
	void _node_exp_adv_term(pNode pn, Axes ai, Exp& exp) {
		//assert(pn->d());
		st aui = _adv_utp_idx[ai];
		ListFE& lpexp = CAST_REF(ListFE*, pn->utp(aui));
		Exp FF[NumFaces];
		Arr_st countCF(NumFaces);
		countCF.assign(0);
		for (typename ListFE::iterator iter = lpexp.begin();
				iter != lpexp.end(); ++iter) {
			Face* pface = iter->first;
			Exp* pexp = iter->second;
			ASSERT(pface->pori() == pn); //
			Direction dir = pface->dir();
			st face_idx = FaceDirectionInOrder(dir);
			// times area
			int sign = IsFacePDirection(dir) ? 1 : -1;
			vt area = pface->area();
			if (pface->ft() == _CoarseFine_) {
				area = pface->pnei()->face_area(Opposite(dir));
			}
			pexp->times(sign * area);
			FF[face_idx].plus(*pexp);
			countCF[face_idx]++;
		}
		for (st i = 0; i < NumFaces; i++) {
			exp.plus(FF[i]);
		}

		_delete_list_face_exp(pn, aui);
	}

	void _node_exp_dif_term(pNode pn, Axes ai, Exp& exp) {
		st utp_idx = _dif_utp_idx[ai];
		ListFE& lpexp = CAST_REF(ListFE*, pn->utp(utp_idx));
		Exp sumCF[NumFaces];
		ArrayListV<st> countCF(NumFaces);
		countCF.assign(0);
		for (typename ListFE::iterator iter = lpexp.begin();
				iter != lpexp.end(); ++iter) {
			Face* pface = iter->first;
			Exp* pexp = iter->second;
			ASSERT(pface->pori() == pn); //
			if (pface->ft() == _Equal_ || pface->ft() == _Boundary_
					|| pface->ft() == _FineCoarse_) {
				cvt a_f = pface->pori()->face_area(pface->dir());
				pexp->times(a_f);
				//if (pn->is_in_on(0.4875, 0.328) && IsFacePDirection(pface->direction)) {
				//	_show_exp(*pexp, *_pdomain);
				//}
				exp.plus((*pexp));
			} else if (pface->ft() == _CoarseFine_) {
				sumCF[FaceDirectionInOrder(pface->dir())].plus(*pexp);
				countCF[FaceDirectionInOrder(pface->dir())]++;
			}
		}
		for (st i = 0; i < NumFaces; i++) {
			if (countCF[i] > 0) {
				// get average face gradient;
				sumCF[i].times(1.0 / Float(countCF[i]));
				// get face area;
				Direction dir = FaceDirectionInOrder(i);
				cvt a_f = pn->face_area(dir);
				sumCF[i].times(a_f);
				exp.plus(sumCF[i]);
			}
		}
		//
		_delete_list_face_exp(pn, utp_idx);

	}

	/*
	 *  cal source term
	 */
	void _node_exp_src_term(pNode pn, Axes ai, Exp& exp) {

	}

	/*
	 * this is the main function of the predict step
	 *
	 */
	void _node_exp_u_star(pNode pn, Axes ai, Exp& exp, vt dt) {
		// Traversal all the faces
		// cfl number check
		vt cfl = _CFL_number(pn, dt);
		ASSERT(cfl < 1.0);
		// calculate the term separately
		// 1 advection term
		Exp exp_adv;
		_node_exp_adv_term(pn, ai, exp_adv);
		exp_adv.times(-1.0 / pn->volume()); //negative;
				// 2 diffusion term
		Exp exp_dif;
		_node_exp_dif_term(pn, ai, exp_dif);
		exp_dif.times(1.0 / pn->volume()); //divide volume
				// 3 source term
		Exp exp_src;
		_node_exp_src_term(pn, ai, exp_src);
		//   add up diffusion and source
		exp_dif.plus(exp_src);
		exp_dif.times(1.0 / pn->cdva(_rho_idx));
		exp_dif.plus(exp_adv);  // add advection term
		exp_dif.times(dt);      // time
		// 4 u term
		exp.plus(exp_dif);
		exp.insert(1.0, pn, 1.0);
	}

	void u_star_on_center(Axes ai, vt dt) {
		//
		//1 Traverse face
		//  Build up the expression on utp
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_face_exp_adv_dif(pn, ai);
		}
		//2 Treverse face to calculate the results
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			Exp exp;
			pNode pn = it.get_pointer();
			_node_exp_u_star(pn, ai, exp, dt);
			vt u_star = exp.subsitute(_v_idx[ai]);
			pn->cd(_vt_idx[ai]) = u_star;
		}
	}

	/*
	 * Boundary condition
	 */
	void set_boundary_condition(st si, st segi, st vi, pBoundaryCondition pbc) {
		_pdomain->_pbindex->insert(si, segi, vi, pbc);
	}

	/*
	 * set
	 */
	void set_v_grid(Function pfun, Axes a) {
		ASSERT(a < Dim);
		_pdomain->set_val_grid(_v_idx[a], pfun);
	}

	void set_v_ghost(Function pfun, Axes a) {
		_pdomain->set_val_ghost(_v_idx[a], pfun);
	}
	void set_v(Function pfun, Axes a) {
		this->set_v_grid(pfun, a);
		this->set_v_ghost(pfun, a);
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
	void set_rho(Function pfun) {
		_pdomain->set_val(_rho_idx, pfun);
	}
	void set_mu_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_mu_idx);
	}
	// can be used to set initial boundary condition
	void set_mu(Function pfun) {
		_pdomain->set_val(_mu_idx, pfun);
	}

	/*
	 * this function return the max index of valuable used in
	 * poisson, the max index make sure the array size is ok.
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

	void run() {
		// 1 predict step
		for (st ia = 0; ia < Dim; ia++) {
			u_star_on_center(ToAxes(ia), _dt);
		}
		// 2 slove pressure equation

		// 3 correct veolity
	}

protected:
	int _new_list_face_exp(pNode& pn, st utp_idx) {
		utPointer& utp = pn->utp(utp_idx);
		if (utp == nullptr) {
			utp = new ListFE();   //new !!!!!
			return 1;
		}
		return 0;
	}
	int _delete_list_face_exp(pNode pn, st utp_idx) {
		utPointer& utp = pn->utp(utp_idx);
		if (utp != nullptr) {
			ListFE& lpexp = CAST_REF(ListFE*, utp);
			for (typename ListFE::iterator iter = lpexp.begin();
					iter != lpexp.end(); ++iter) {
				delete iter->first;
				delete iter->second;
			}
			lpexp.clear();
			delete &lpexp;
			utp = nullptr;
			return 1;
		}
		return 0;
	}

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

	vt _CFL_number(pNode pn, vt dt) {
		vt veo[Dim];
		vt cfl[Dim];
		for (st i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(_v_idx[i]);
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
	}

}
;

}

#endif
