#ifndef STENCIL_HPP_
#define STENCIL_HPP_

#include "../carpio_define.hpp"
#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "../algebra/space.hpp"
#include "../algebra/array_list.hpp"

#include <functional>
#include <math.h>

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM, int DIMST>
class Stencil_ {
public:
	static const st Dim = DIMST;  //Dimension of stencil

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> const_Node;
	typedef const Node_<COO_VALUE, VALUE, DIM>* const_pNode;
	typedef Cell_<COO_VALUE, DIM> Cell;
	typedef Cell_<COO_VALUE, DIM> *pCell;
	typedef Data_<VALUE, DIM> Data;
	typedef Data *pData;
protected:
	/*
	 *  data
	 */
	template<typename CV, typename V, int D>
	struct sNode_ { //Stencil node
		Node_<CV, V, D>* pnode;
		int type;  //the type indicates that the node is new or already exist.
	};
	typedef sNode_<COO_VALUE, VALUE, DIM> sNode;
	typedef sNode_<COO_VALUE, VALUE, DIM>* psNode;

	SpaceT<sNode, Dim> _pnodes;
	ArrayListT<Axes> _axes;
	ArrayListT<st> _steps_f;
	ArrayListT<st> _steps_b;
protected:
	/*
	 * function
	 */
	pNode _neighbor_forward(pNode pn, Axes a) {
		ASSERT(pn != nullptr);
		pNode res = nullptr;
		if (a == _X_) {
			res = pn->get_neighbor(_XP_);
		}
		if (a == _Y_) {
			res = pn->get_neighbor(_YP_);
		}
		if (a == _Z_) {
			res = pn->get_neighbor(_ZP_);
		}
		return res;
	}
	pNode _neighbor_backward(pNode pn, Axes a) {
		ASSERT(pn != nullptr);
		pNode res = nullptr;
		if (a == _X_) {
			res = pn->get_neighbor(_XM_);
		}
		if (a == _Y_) {
			res = pn->get_neighbor(_YM_);
		}
		if (a == _Z_) {
			res = pn->get_neighbor(_ZM_);
		}
		return res;
	}
	st _get_idx_c(st sf, st sb) const {
		return sb;
	}
	void _construct_1d(pNode pnc, Axes a, st sf, st sb) {
		ASSERT(pnc != nullptr);
		ASSERT(Dim == 1);
		_axes[0] = a;
		_steps_f[0] = sf;
		_steps_b[0] = sb;
		_pnodes.reconstruct(sf + sb + 1);
		// set null
		for (st i = 0; i < _pnodes.size(); ++i) {
			_pnodes.at_1d(i).pnode = nullptr;
			_pnodes.at_1d(i).type = 0;   //whish is not created by this class
		}
		// set center node
		_pnodes.at_1d(sb).pnode = pnc;
		//find neighbor forward
		pNode pc = pnc;
		for (st i = 0; i < sf; ++i) {
			pNode pnt = nullptr;
			pnt = _neighbor_forward(pc, a);
			if (pnt == nullptr) {
				break;
			} else {
				_pnodes.at_1d(sb + i + 1).pnode = pnt;
				pc = pnt;
			}
		}
		//find neighbor backward
		pc = pnc;
		for (st i = 0; i < sb; ++i) {
			pNode pnt = nullptr;
			pnt = _neighbor_backward(pc, a);
			if (pnt == nullptr) {
				break;
			} else {
				_pnodes.at_1d(sb - 1 - i).pnode = pnt;
				pc = pnt;
			}
		}
	}
public:
	/*
	 *  1d constructor
	 */
	Stencil_(pNode pnc, Axes a1, st sf1, st sb1) :
			_pnodes(), _axes(Dim), _steps_f(Dim), _steps_b(Dim) {
		_construct_1d(pnc, a1, sf1, sb1);
	}

protected:
	/*
	 * iterator
	 */

public:
	pNode operator()(st i, st j = 0, st k = 0) {
		return _pnodes(i, j, k).pnode;
	}
	const_pNode operator()(st i, st j = 0, st k = 0) const {
		return _pnodes(i, j, k).pnode;
	}
	pNode at_1d(st i) {
		return _pnodes.at_1d(i).pnode;
	}
	const_pNode at_1d(st i) const {
		return _pnodes.at_1d(i).pnode;
	}
	st size() const {
		return _pnodes.size();
	}
	st null_nodes() const {
		st res = 0;
		for (st i = 0; i < _pnodes.size(); ++i) {
			if (at_1d(i) == nullptr) {
				res++;
			}
		}
		return res;
	}
	st non_null_nodes() const {
		st res = 0;
		for (st i = 0; i < _pnodes.size(); ++i) {
			if (at_1d(i) != nullptr) {
				res++;
			}
		}
		return res;
	}
	st dim() const {
		return Dim;
	}
	/*
	 *  get
	 */
	pNode center_pnode() {
		if (Dim == 1) {
			return _pnodes(_steps_b[0], 0, 0).pnode;
		} else if (Dim == 2) {
			return _pnodes(_steps_b[0], _steps_b[1], 0).pnode;
		} else {
			return _pnodes(_steps_b[0], _steps_b[1], _steps_b[2]).pnode;
		}
	}
	const_pNode center_pnode() const {
		if (Dim == 1) {
			return _pnodes(_steps_b[0], 0, 0).pnode;
		} else if (Dim == 2) {
			return _pnodes(_steps_b[0], _steps_b[1], 0).pnode;
		} else {
			return _pnodes(_steps_b[0], _steps_b[1], _steps_b[2]).pnode;
		}
	}
	Axes axes(st i = 0) const {
		return _axes[i];
	}
protected:
	bool _is_valid_axes(Axes a) {
		for (st i = 0; i < _axes.size(); ++i) {
			if (_axes[i] == a) {
				return true;
			}
		}
		return false;
	}
	st _to_arraylist_idx(Axes a) const {
		st i = 0;
		for (; i < _axes.size(); ++i) {
			if (_axes[i] == a) {
				return i;
			}
		}
		ASSERT_MSG(false, "Not valid axes");
		return i;
	}
public:
	pNode forward_pnode(st step, Axes a = _X_) {
		st ai = _to_arraylist_idx(a);
		if (!(step > 0 && step < _steps_f[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] + step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] + step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] + step).pnode;
		}
	}
	const_pNode forward_pnode(st step, Axes a = _X_) const {
		st ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_f[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] + step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] + step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] + step).pnode;
		}
	}
	pNode backward_pnode(st step, Axes a = _X_) {
		st ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_b[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] - step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] - step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] - step).pnode;
		}
	}
	const_pNode backward_pnode(st step, Axes a = _X_) const {
		st ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_b[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] - step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] - step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] - step).pnode;
		}
	}
	/*
	 *  show
	 */
	void show() const {
		std::cout << " stencil dim = " << Dim << "\n";
		std::cout << " size        = " << size() << "\n";
		std::cout << " null        = " << null_nodes() << "\n";
		std::cout << " no null     = " << size() - null_nodes() << "\n";
		if (Dim == 1) {
			std::cout << " ---";
			for (st i = 0; i < _pnodes.size(); ++i) {
				if (i == _steps_b(0) && at_1d(i) != nullptr) {
					std::cout << "C";
				} else if (at_1d(i) != nullptr) {
					std::cout << "O";
				} else {
					std::cout << "X";
				}
				std::cout << "---";
			}
			std::cout << "--> " << ToString(_axes(0)) << "\n";
		}
	}

};

}

#endif
