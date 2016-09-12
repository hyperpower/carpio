#ifndef _EQUATION_H_
#define _EQUATION_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "event.hpp"
#include "../utility/Clock.h"
#include "../io/mmio.h"

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

//#define __Debug__

#define __x__  -0.4
#define __y__  0.45

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class Equation_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Equation_<cvt, vt, Dim> Self;
	typedef Equation_<cvt, vt, Dim>& ref_Self;
	typedef Equation_<cvt, vt, Dim>* pSelf;
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

	// function define ------------------------------
	typedef std::function<vt(cvt, cvt, cvt)> Function;

	// Events ----------------------------------------
	typedef Event_<cvt, vt, Dim> Event;
	typedef std::shared_ptr<Event> spEvent;
	typedef std::unordered_map<std::string, spEvent> Events;

	// Variables ------------------------------------
	typedef std::map<std::string, st> Variables;

	typedef Expression_<cvt, vt, Dim> Exp;
	typedef Exp* pExp;
	typedef std::shared_ptr<Exp> spExp;
	typedef std::shared_ptr<Face> spFace;

	typedef typename Exp::Term Term;
	typedef typename Exp::iterator ExpIter;
	typedef typename Exp::const_iterator const_ExpIter;

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

	void _new_map_fe(pNode pn, st idx) {
		utPointer& utp = pn->utp(idx);
		if (utp == nullptr) {
			utp = new Map_FE();
		} else {
			Map_FE& mfe = CAST_REF(pMap_FE, utp);
			mfe.clear();
		}
	}

	void _delete_map_fe(pNode pn, st idx) {
		utPointer& utp = pn->utp(idx);
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
		exp_gra = Interpolate::OnFace_Gradient(pori, dir, 2);
	}

	void _insert_face_val_gra(pNode pn, st ut_idx) {
		utPointer utp = pn->utp(ut_idx);
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
				utPointer utp_nei = pnei->utp(ut_idx);
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
			_new_map_fe(pn, this->_vars_ut["FE"]);
		}
		// inert face value and gradient
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_insert_face_val_gra(pn, this->_vars_ut["FE"]);
		}
	}

	void _delete_FE_on_leaf() {
		Grid& grid = _pdomain->grid();
		// delete map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_delete_map_fe(pn, this->_vars_ut["FE"]);
		}
	}

protected:
	//matrx
	typedef ArrayListV<st> Arr_st;

	pDomain _pdomain;

	// time relates variables
	vt _dt;
	vt _max_step;

	// Flag
	Events _events;

	// Variables
	Variables _vars_c;  //!< variables on the center of node
	Variables _vars_ut; //!< untype variables on the node

#ifdef __Debug__
	std::fstream debug_log;
#endif
public:
	/*
	 * default constructor
	 *   2D
	 *
	 * utp
	 * utp map          0
	 */
	Equation_(pDomain pf, vt dt = -1, st maxstep = 1) :
			_pdomain(pf) {
		_dt = dt;
		_max_step = maxstep;
	}

	virtual int run_one_step(st step){
		std::cout<<step<<"  Equation: run one step \n";
		return -1;
	}

	virtual int solve(){
		std::cout<<"  Equation: solve \n";
		return -1;
	}

	void run(){
		// the equation don't have time
		if(!this->has_time()){
			_dt = 0;
			_max_step = 1;
			// start events
			run_events(0, 0.0, -1);

			solve();

			run_events(_max_step, -1, 1);
			// start events
		}
		vt t = 0.0;
		st step = 0;
		// events before calculation
		run_events(0 ,0.0 , -1);
		// loop
		for(;step < _max_step; ++step){
			//
			// events before each step
			run_events(step ,t , -1);

			// run one step =================
			run_one_step(step);
			// ==============================

			// events after each step
			run_events(step ,t , 1);
			//
			t += this->_dt;
		}
		// events after calculation
		run_events( _max_step, t , 1);
	}

	virtual ~Equation_() {

	}

	void run_events(st step, vt t, int fob){
		for(auto& event : this->_events){
			if(event.second->_do_execute(step, t,fob)){
				event.second->execute(step, t, fob);
			}
		}
	}

public:

	/**
	 * @brief this function returns the max index of the valuables in
	 *        ns, the max index makes sure the array size is ok.
	 */
	st max_idx(const Variables& vars) {
		st max = 0;
		for (auto& item : vars) {
			if (item.second > max) {
				max = item.second;
			}
		}
		return max;
	}

	/**
	 * @brief this function check the events
	 *        if flags contain key and value == val return true
	 *        the default val == 1
	 */
	bool has(const std::string& key, const int& val = 1) const {
		auto it = this->_events.find(key);
		if (it != this->_events.end()) {
			if ((*it).second->flag() == val) {
				return true;
			}
		}
		return false;
	}

	void set_output_time(int is = -1, int ie = -1, int istep = -1, std::FILE * f = nullptr) {
		// output to screen
spEvent pse(new EventOutputTime_<cvt,vt, Dim>(is, ie, istep, f));
		this->_events["OutputTime"] = pse ;
	}

	/*
	 * Does the equation has term of time
	 */
	bool has_time() const {
		return this->has("time");
	}

protected:
	void _output_time(st step, vt dt, vt t, vt cpu, vt wall) {

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

	vt _CFL_number(pNode pn, vt dt, st vidx) {
		vt veo[Dim];
		vt cfl[Dim];
		for (st i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(vidx);
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
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

}
;

}

#endif
