#include "interpolate.h"

namespace carpio {
/*
 *  There is only one node in the stencil
 *  The point (x,y) is in that node, and
 *  the value is equal to the center point.
 *  Stencil:  ---X---C---X---->
 */
int _1Node(Float& res, const st& idx, const Stencil_2D1& stc, Float x,
		Float y) {
	ASSERT(stc.non_null_nodes() == 1);
	// the non null nodes must be center node
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	ASSERT(pc->is_in_on(x, y));
	// ---
	res = pc->cd(idx);
	return _SUCCESS;
}
int _1Node(PData_2D& res, const Stencil_2D1& stc) {
	//assert
	ASSERT(stc.non_null_nodes() == 1);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	ASSERT(pc->is_in_on(res.x(), res.y()));
	// ---
	for (st i = 0; i < res.size(); ++i) {
		res.set_val(i, pc->cd(res.idx(i)), res.idx(i), 1);
	}
	return _SUCCESS;
}

/*
 *  There are two nodes in the stencil, x is
 *  the location, somewhere on the stencil axes.
 *
 *  Stencil: ---X---C---O---->  or
 *           ---O---C---X---->
 *                ^
 *                x
 */
int _2NodeOnAxes(Float& res, const st& idx, const Stencil_2D1& stc, Float x) {
	// assert
	ASSERT(stc.non_null_nodes() == 2);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	// another pnode is close to center
	int flag = 0;
	typename Stencil_2D1::const_pNode po = stc.forward_pnode(1, stc.axes());
	if (po == nullptr) {
		flag = 1;
		po = stc.backward_pnode(1, stc.axes());
	}
	// ---
	typename Stencil_2D1::const_pNode ps = nullptr, pb = nullptr;
	//                                 s small       b big
	if (flag == 0) {
		ps = pc;
		pb = po;
	} else {
		ps = po;
		pb = pc;
	}

}

/*
 *  the idx in res is defined
 */
int _GetPDataFromPNodeCenter(PData_2D& res, const_pNode_2D pc) {
	ASSERT(pc != nullptr);
	res.set_point(pc->cp(_X_), pc->cp(_Y_));  //set point as the center of pNode
	for (st i = 0; i < res.size(); ++i) {
		res.flag(i) = PData_2D::Flag_Center;  //set flag as Center
		res.val(i) = pc->cd(res.idx(i));      //get val
	}
	return _SUCCESS;
}

void _AverangeValueFromCenterLeaf(        //
		PData_2D& res,        //
		const_pNode_2D pn) {
//improve for special case
	if (pn->is_leaf()) {
		_GetPDataFromPNodeCenter(res, pn);
		return;
	}
//========================
	res.set_all_center();
	st num = 0;
	std::function<void(const_pNode_2D, int)> fun =
			[&res, &num](const_pNode_2D pnode, int dummy) {
				if (pnode->is_leaf()) {
					for(st i = 0; i<res.size();++i) {
						res.val(i) = res.val(i) + pnode->cd(res.idx(i));
					}
					num++;
				}
			};
	int dummy = 1;
	pn->traversal(fun, dummy);
	for (st i = 0; i < res.size(); ++i) {
		res.val(i) /= num;
	}
}
template<class Node>
inline int __2NodeRelation(const Node* pc, const Node* p) {
	st pcl = pc->get_level();
	st pl = p->get_level();
	if (pcl == pl && p->has_child() == false) {
		//case 1  on the same level
		return _E_;
	}
	if (pcl > pl) {
		//case 1  neighbor is coarse
		return _F_C_;
	}
	if (pcl == pl) {
		//case 3  neighbor is fine
		return _C_F_;
	}
	SHOULD_NOT_REACH;
	return -1;
}

/*
 * Linear interpolation
 */
template<typename TYPE>
inline Float __LinearInterpolation( //
		const TYPE& x, //
		const TYPE& x0, const TYPE& y0, //
		const TYPE& x1, const TYPE& y1  //
		) {
	return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}

int _2NodeOnAxes(PData_2D& res, const Stencil_2D1& stc) {
	// assert
	ASSERT(stc.non_null_nodes() == 2);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	PData_2D pdc(res);
	_AverangeValueFromCenterLeaf(pdc, pc);
	// another pnode is close to center
	typename Stencil_2D1::const_pNode pn = stc.forward_pnode(1, stc.axes());
	if (pn == nullptr) {
		pn = stc.backward_pnode(1, stc.axes());
	}
	PData_2D pdn(res);
	switch (__2NodeRelation(pc, pn)) {
	case _E_: {
		_AverangeValueFromCenterLeaf(pdn, pn);
		break;
	}
	case _F_C_: {
		//unfinish ==================
		SHOULD_NOT_REACH;
		break;
	}
	case _C_F_: {
		_AverangeValueFromCenterLeaf(pdn, pn);
		break;
	}
	}
	// ---
	for (st i = 0; i < res.size(); ++i) {
		Float x = res.p(stc.axes(0));
		Float x1 = pdc.p(stc.axes(0));
		Float y1 = pdc.val(i);
		Float x2 = pdn.p(stc.axes(0));
		Float y2 = pdn.val(i);
		res.val(i) = __LinearInterpolation(x, x1, y1, x2, y2);
	}
	return _SUCCESS;
}
/*
 * This function construct a stencil
 * the point must in the Node
 */
int _Stencil1f(Stencil_2D1& stc, pNode_2D pn, Cvt x, Cvt y){
	ASSERT(pn->is_in_on(x, y));

}

}
