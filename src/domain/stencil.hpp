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
	typedef const Node const_Node;
	typedef const Node* const_pNode;
	typedef Cell_<COO_VALUE, DIM> Cell;
	typedef Cell_<COO_VALUE, DIM> *pCell;
	typedef Data_<VALUE, DIM> Data;
	typedef Data *pData;
protected:
	/*
	 *  data
	 */
	SpaceT<pNode, Dim> _pnodes;
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
	st _get_idx_c(st sf, st sb) const{
		return sb;
	}
	void _construct(pNode pnc, Axes a, st sf, st sb){
		ASSERT(pnc != nullptr);

	}
public:
	Stencil_() :
			_pnodes(), _axes(Dim), _steps_f(Dim), _steps_b(Dim) {

	}

};

}
