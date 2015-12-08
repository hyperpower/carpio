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

}
