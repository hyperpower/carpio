#ifndef _POISSON_H_
#define _POISSON_H_

#include "calculation_define.hpp"
#include "expression.hpp"

namespace carpio {
//This file use to solve poisson equation
//
//   ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//

/*
 * the Poisson class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Poisson_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Poisson_<COO_VALUE, VALUE, DIM> Self;
	typedef Poisson_<COO_VALUE, VALUE, DIM>& ref_Self;
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
	/*
	 struct pNode_less: public std::binary_function<pNode, pNode, bool> {
	 bool operator()(const pNode& __x, const pNode& __y) const {
	 return __x->d_idx() < __y->d_idx();
	 }
	 };
	 typedef Polynomial_<vt, pNode, int, IsZero_<vt>, IsZero_<int>, pNode_less,
	 std::less<int> > Exp;
	 typedef Exp* pExp;
	 typedef typename Exp::Term Term;
	 typedef typename Exp::iterator ExpIter;
	 typedef typename Exp::const_iterator const_ExpIter;

	 bool _is_ghost_node(const ExpIter& iter) const {
	 pNode pn = iter->first.first;
	 ASSERT(pn != nullptr);
	 return (pn->get_type() == _Ghost_);
	 }

	 const pNode& _get_pnode(ExpIter& iter) const {
	 // return the reference
	 return iter->first.first;
	 }

	 const int& _get_exp(ExpIter& iter) const {
	 return iter->first.second;
	 }

	 vt& _get_coe(ExpIter& iter) const {
	 return iter->second;
	 }
	 */
protected:
	pDomain _pdomain;

	st _beta_idx;
	st _phi_idx;
	st _f_idx;

	st _utp_idx;

	//matrix
	typedef MatrixSCR_<vt> Mat;
	typedef ArrayListV<vt> Arr;
public:
	/*
	 * constructor
	 */
	Poisson_(pDomain pf, st bi, st pi, st fi) :
			_pdomain(pf), _beta_idx(bi), _phi_idx(pi), _f_idx(fi) {
		_utp_idx = 0;
		_pdomain->new_data(this->max_idx() + 1, 0, 0, 1);
	}
	/*
	 * set
	 */
	void set_f(Function pfun) {
		_pdomain->set_val_grid(_f_idx, pfun);
	}

	void set_beta(Function pfun) {
		_pdomain->set_val_grid(_beta_idx, pfun);
	}
	void set_f_ghost(Function pfun) {
		_pdomain->set_val_ghost(_f_idx, pfun);
	}
	void set_beta_ghost(Function pfun) {
		_pdomain->set_val_ghost(_beta_idx, pfun);
	}
	void set_f_all(Function pfun) {
		this->set_f(pfun);
		this->set_f_ghost(pfun);
	}
	void set_beta_all(Function pfun) {
		this->set_beta(pfun);
		this->set_beta_ghost(pfun);
	}
	void set_phi_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_phi_idx);
	}
	// can be used to set initial boundary condition
	void set_phi(Function pfun) {
		_pdomain->set_val(_phi_idx, pfun);
	}

	/*
	 * this function return the max index of valuable used in
	 * poisson, the max index make sure the array size is ok.
	 */
	st max_idx() {
		st max = 0;
		max = Max(max, _beta_idx);
		max = Max(max, _phi_idx);
		max = Max(max, _f_idx);
		return max;
	}
	/*
	 * Boundary condition
	 */
	void set_boundary_condition(st si, st segi, st vi, pBoundaryCondition pbc) {
		_pdomain->_pbindex->insert(si, segi, vi, pbc);
	}
	/*
	 * this function can be deleted
	 */
	int _face_scheme_gradient(pFace pface, Exp& exp) {
		//face type
		switch (pface->ft()) {
		//_Null_ = -1, _Boundary_ = 0, _Equal_ = 1, _FineCoarse_ = 2, _CoarseFine_ = 3,
		case _Null_:
			SHOULD_NOT_REACH;
			break;
		case _Boundary_: {
			this->_face_scheme_equal_gradient(pface, exp);
			break;
		}
		case _Equal_: {
			this->_face_scheme_equal_gradient(pface, exp);
			break;
		}
		case _FineCoarse_: {
			this->_face_scheme_equal_gradient(pface, exp);
			break;
		}
		case _CoarseFine_: {
			SHOULD_NOT_REACH;
			std::cout << "Corase Fine function unfinish\n";
			break;
		}
		default:
			return -1;
		}
		return pface->ft();
	}

	int _face_scheme_equal_gradient(pFace pface, Exp& exp) {
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
		Float beta_f = expb.subsitute(_beta_idx);

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
	int _face_scheme_boundary_gradient(pFace pface, Exp& exp) {
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
		Float beta_f = expb.subsitute(_beta_idx);
		Exp expg = Interpolate::OnFace_Gradient(pori, dir, 2);
		//_exp_substitute_ghost(pori, expg);
		expg.times(sign * beta_f);
		exp.plus(expg);
		return 1;
	}
	int _face_scheme_fine_coarse_gradient(pFace pface, Exp& exp) {
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
		Float beta_f = expb.subsitute(_beta_idx);
		// gradient on face
		Exp expg = Interpolate::OnFace_Gradient(pori, dir, 2);
		expg.times(sign * beta_f);
		exp.plus(expg);
		return 1;
	}
	/*
	 * build face exp to utp
	 */
protected:
	typedef std::list<std::pair<pFace, pExp> > ListFE;
	typedef std::pair<pFace, pExp> PairFE;
public:
	void _node_face_exp(pNode& pn) {
		_new_list_face_exp(pn, _utp_idx);
		utPointer& utp = pn->utp(_utp_idx);
		for (st i = 0; i < NumFaces; i++) {
			Direction dir = FaceDirectionInOrder(i);
			pNode pnei = pn->get_neighbor_fast(dir);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face
			FaceType ft = GetFaceType(pn, pnei);
			pExp pexp = new Exp();
			bool flag = true;
			if ((ft == _Boundary_)				//1
			|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
					|| (ft == _FineCoarse_))				//3
					{
				// work on pn
				pFace pf = new Face(pn, pnei, dir, ft);
				_face_scheme_gradient(pf, *pexp);
				ListFE& lpexp = CAST_REF(ListFE*, utp);
				PairFE pfe(pf, pexp);
				lpexp.push_back(pfe);
				flag = false;
			}
			// case 2 3
			if ((ft == _Equal_ && (IsFacePDirection(dir)))				//2
			|| (ft == _FineCoarse_))				//3
					{
				// work on pnei
				_new_list_face_exp(pnei, _utp_idx);
				utPointer& utpn = pnei->utp(_utp_idx);
				FaceType oft = (ft == _FineCoarse_) ? _CoarseFine_ : _Equal_;
				pFace fn = new Face(pnei, pn, Opposite(dir), oft);
				pExp pexpn = new Exp(*pexp);
				pexpn->times(-1.0);
				ListFE& lpexp = (*CAST(ListFE*, utpn));
				PairFE pairn(fn, pexpn);
				lpexp.push_back(pairn);
			}
			if(flag){
				delete pexp;
			}
		}
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
	int _exp_substitute_ghost(pNode pn, Exp& exp) {
		for (ExpIter iter = exp.begin(); iter != exp.end();) {
			if (Exp::IsGhostNode(iter)) {
				ExpIter itere = iter;
				++iter;
				// find boundary condition
				const_pNode pg = Exp::GetpNode(itere);
				if (!Exp::IsZeroCoe(itere)) { // if coe = 0 the term will be deleted
					const_pBoundaryCondition pbc = _find_bc(pg, _phi_idx);
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
	void _node_exp(pNode pn, Exp& exp) {
		ListFE& lpexp = CAST_REF(ListFE*, pn->utp(_utp_idx));
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
		_f_term(pn, exp);
		_delete_list_face_exp(pn, _utp_idx);
	}

	void _f_term(pNode& pn, Exp& exp) {
		// f term is const;
		// debug
		//if (pn->is_in_on(0.41, 0.3)) {
		//	std::cout << " f " << -(pn->cda(_f_idx) * pn->volume()) << '\n';
		//}
		Term f_(-(pn->cda(_f_idx) * pn->volume()), pn, 0);
		exp.insert(f_);
	}

	void _show_exp(Exp& exp, Domain& domain) {
		exp.show();
		std::list<Gnuplot_actor> lga;
		Gnuplot_actor ga;
		GnuplotActor_LeafNodes(ga, domain.grid());
		lga.push_back(ga);
		GnuplotActor_Expression(ga, exp);
		lga.push_back(ga);
		Gnuplot gp;
		gp.set_equal_ratio();
		//gp.set_xrange(2.0,3.0);
		//gp.set_yrange(1.5,2.5);
		//gp.set_cbrange(-2.0, 2.0);
		gp.plot(lga);
		//delete shape
	}

	int _bulid_matrix(Mat& mat, Arr& b) { //
        //1 Traverse face
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_node_face_exp(pn);
		}

		typedef std::list<st> ListST;
		typedef std::list<vt> ListVT;
		ListST l_rowptr;
		l_rowptr.push_back(0);
		ListST l_colid;
		ListVT l_val;
		ListVT l_b;
		int countnz = 0;
		//2 build matrix
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			Exp exp;
			pNode pn = it.get_pointer();
			_node_exp(pn, exp);
			_exp_substitute_ghost(pn, exp);
			// debug
			//if (pn->is_in_on(-0.4, -0.4)) {
			//	_show_exp(exp, *_pdomain);
			//}
			exp.concise();
			// debug

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
		return 1;
	}
	int solve() {
		vt tol = 1e-6;
		int max_iter = 1000;
		std::list<vt> lr;
		int rcode = this->solve(tol, max_iter, lr);
		//std::cout<< "max iter "<<max_iter<< " tol "<< tol <<"\n";
		return rcode;
	}

	int solve(vt& tol, int& max_iter, std::list<vt>& lr) {
		Mat mat;
		Arr b;
		this->_bulid_matrix(mat, b);
		Arr x(b.size());
		// initial x
		_grid_to_arr(x);
		//x.show();
		//defalt set up ========
		//std::list<Float> lr;	//list residual
		//solver =======================
#ifdef VIENNACL_WITH_OPENCL
		Solver_<vt>::Solve(mat, x, b, tol, max_iter);
#else
		int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);
		if (sf != 0) {
			std::cerr << " >! Poisson solve failed \n";
			return -1;
		}
#endif
		//gnuplot_show_ylog(lr);
		//put the value back
		x.show();
		_arr_to_grid(x);
		return 1;
	}

protected:
	void _arr_to_grid(const Arr& x) {
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			it->cd(_phi_idx) = x[it->d_idx()];
		}
	}
	void _grid_to_arr(Arr& x) {
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			x[it->d_idx()] = it->cd(_phi_idx);
		}
	}

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

}
;
}

#endif
