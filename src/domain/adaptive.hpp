#ifndef ADAPTIVE_HPP_
#define ADAPTIVE_HPP_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"
#include "shape.hpp"

namespace carpio {
template<typename COO_VALUE, typename VALUE, st DIM>
void _visit_current_info(Node_<COO_VALUE, VALUE, DIM> *pn, utPointer utp) {
	ArrayListV<st> &arr = CAST_REF(ArrayListV<st>*, utp);
	st &min_l = arr[0];
	st &max_l = arr[1];
	st &num_n = arr[2];
	st &num_leaf = arr[3];
	num_n++;

	st cl = pn->get_level();
	if (cl > max_l) {
		max_l = cl;
	}
	if (pn->is_leaf()) {
		min_l = (cl < min_l) ? pn->get_level() : min_l;
		num_leaf++;
	}
}

template<typename COO_VALUE, typename VALUE, st DIM>
class Adaptive_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef typename SpaceT<pNode, Dim>::reference reference;
	typedef typename SpaceT<pNode, Dim>::const_reference const_reference;
	typedef typename SpaceT<pNode, Dim>::size_type size_type;

	typedef void (*pfunction)(pNode, utPointer);

	typedef void (*pfunction_conditional)(arrayList &, pNode, utPointer);

protected:

	st _c_min_l;    //current min level
	st _c_max_l;    //current max level
	st _c_num_n;    //current num of node
	st _c_num_leaf; //current num of leaf

	st _min_l;
	st _max_l;
	//st _min_num_n;
	//st _max_num_n;

	pGrid _grid;

	/**
	 *  protected function
	 */
	void _update_current_info() {
		ArrayListV<st> arr_info(4);
		for (st i = 0; i < _grid->size(); i++) {
			pNode pn = _grid->nodes.at_1d(i);
			if (pn != nullptr) {
				//
				pn->traversal(_visit_current_info, &arr_info);
			}
		}
		_c_min_l = arr_info[0];
		_c_max_l = arr_info[1];
		_c_num_n = arr_info[2];
		_c_num_leaf = arr_info[3];
	}

public:
	/*
	 * constructor
	 */
	Adaptive_(Grid_<COO_VALUE, VALUE, DIM> *pg,  //the pointer of grid
			st minl = 1,  //min level
			st maxl = 3   //max level
			) {
		_min_l = minl;
		_max_l = maxl;
		_grid = pg;
		ASSERT(minl <= maxl);
		_update_current_info();
	}

	/*
	 * adapt
	 */
	void adapt_full() {
		std::function<void(pNode&, st)> fun = [](pNode pn, st min_l) {
			if(pn->is_leaf()&& pn->get_level()<min_l) {
				pn->new_full_child();
			}
		};
		for (st i = 0; i < _grid->size(); i++) {
			pNode pn = _grid->nodes.at_1d(i);
			if (pn != nullptr) {
				pn->traversal(fun, _min_l);
			}
		}
		_update_current_info();
	}
	void adapt_inner_solid(const Shape_<VALUE, DIM>& shape, vt vr = 1.0) {
		// 1 adapt to min level
		// adapt_full();
		// 2 all out   -> adapt to min level
		// 3 all in    -> do not initial
		// 4 intersect -> adapt to max level
		std::function<void(pNode&, st)> fun =
				[this, &shape, &vr](pNode pn, st max_l) {
					//if(pn->is_leaf()) {
						Shape2D sn,res;
						CreatCube(sn,
								pn->p(_M_,_X_), pn->p(_M_,_Y_),
								pn->p(_P_,_X_), pn->p(_P_,_Y_));
						Intersect(sn,shape,res);
						vt vsn = sn.volume();
						vt vres = res.volume();
						if(res.empty() || Abs((vsn-vres)/vsn)>=vr) { //node is all out
							if(pn->get_level() < this->_min_l) {
								pn->new_full_child();
							}
						} else if(Abs((vres-vsn)/vsn)< 1e-5) { // all in
							if(pn->is_root()) {
								delete pn;
								pn=nullptr;
							} else {
								pNode f = pn->father;
								f->child[pn->get_idx()] = nullptr;
								delete pn;
								pn=nullptr;
							}
						} else {  //intersect
							if(pn->get_level() < max_l) {
								pn->new_full_child();
							}
						}
					//}
				};

		for (st i = 0; i < _grid->size(); i++) {
			pNode pn = _grid->nodes.at_1d(i);
			if (pn != nullptr) {
				pn->traversal(fun, _max_l);
			}
		}
		_update_current_info();
	}
	void adapt_bound_solid(const Shape_<VALUE, DIM>& shape, vt vr = 1.0) {
		// 1 adapt to min level
		// adapt_full();
		// 2 all out   -> adapt to min level
		// 3 all in    -> do not initial
		// 4 intersect -> adapt to max level
		std::function<void(pNode&, st)> fun =
				[this, &shape, &vr](pNode& pn, st max_l) {
					ASSERT(pn!=nullptr);

					//if(pn->is_leaf()) {
					Shape2D sn,res;
					CreatCube(sn,
							pn->p(_M_,_X_), pn->p(_M_,_Y_),
							pn->p(_P_,_X_), pn->p(_P_,_Y_));
					Intersect(sn,shape,res);
					vt vsn = sn.volume();
					vt vres = res.volume();
					if(res.empty() ||//
							Abs((vsn-vres)/vsn)>=vr) { //node is all out
						if(pn->is_root() ) {
							delete pn;
							pn=nullptr;
						} else {
							pNode f = pn->father;
							f->child[pn->get_idx()] = nullptr;
							delete pn;
							pn=nullptr;
						}
					} else if(Abs((vsn-vres)/vsn) <= 1e-5) { // all in
						if(pn->get_level() < this->_min_l) {
							pn->new_full_child();
						}
					} else {  //intersect
						if(pn->get_level() < max_l) {
							pn->new_full_child();
						}
					}
					//}
				};

		for (st i = 0; i < _grid->size(); i++) {
			pNode& pn = _grid->nodes.at_1d(i);
			if (pn != nullptr) {
				Traversal(pn, fun, _max_l);
			}
		}

		_update_current_info();
	}
	void adapt_shape_boundary(const Shape_<VALUE, DIM>& shape) {
		// 1 adapt to min level
		// adapt_full();
		// 2 all out   -> adapt to min level
		// 3 all in    -> adapt to min level
		// 4 intersect -> adapt to max level
		std::function<void(pNode&, st)> fun =
				[this, &shape](pNode& pn, st max_l) {
					if(pn->is_leaf()) {
						Shape2D sn,res;
						CreatCube(sn,
								pn->p(_M_,_X_), pn->p(_M_,_Y_),
								pn->p(_P_,_X_), pn->p(_P_,_Y_));
						Intersect(sn,shape,res);
						if(res.empty()) { //node is all out
							if(pn->get_level() < this->_min_l) {
								pn->new_full_child();
							}
						} else if(Abs(res.volume() - sn.volume())<1e-8) { // all in
							if(pn->get_level() < this->_min_l) {
								pn->new_full_child();
							}
						} else {  //intersect
							if(pn->get_level() < max_l) {
								pn->new_full_child();
							}
						}
					}
				};

		for (typename Grid::iterator_leaf iter = _grid->begin_leaf();
				iter != _grid->end_leaf(); ++iter) {
			pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				pn->traversal(fun, _max_l);
			}
		}
		_update_current_info();
	}
	void adapt_shape_inner(const Shape_<VALUE, DIM>& shape) {
			// 1 adapt to min level
			// adapt_full();
			// 2 all out   -> adapt to min level
			// 3 all in    -> adapt to min level
			// 4 intersect -> adapt to max level
			std::function<void(pNode&, st)> fun =
					[this, &shape](pNode& pn, st max_l) {
						if(pn->is_leaf()) {
							Shape2D sn,res;
							CreatCube(sn,
									pn->p(_M_,_X_), pn->p(_M_,_Y_),
									pn->p(_P_,_X_), pn->p(_P_,_Y_));
							Intersect(sn,shape,res);
							if(res.empty()) { //node is all out
								if(pn->get_level() < this->_min_l) {
									pn->new_full_child();
								}
							} else if(Abs(res.volume() - sn.volume())<1e-8) { // all in
								if(pn->get_level() < max_l) {
									pn->new_full_child();
								}
							} else {  //intersect
								if(pn->get_level() < max_l) {
									pn->new_full_child();
								}
							}
						}
					};

			for (typename Grid::iterator_leaf iter = _grid->begin_leaf();
					iter != _grid->end_leaf(); ++iter) {
				pNode pn = iter.get_pointer();
				if (pn != nullptr) {
					pn->traversal(fun, _max_l);
				}
			}
			_update_current_info();
		}

	/*
	 * show information
	 */
	void show_i() const {

	}
};

}

#endif
