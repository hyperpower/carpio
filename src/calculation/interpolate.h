#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include "../carpio_define.hpp"
#include "../domain/domain.hpp"

#include <functional>
#include <math.h>

namespace carpio {
/*
 *  There is only one node in the stencil
 *  The point (x,y) is in that node, and
 *  the value is equal to the center point.
 *  Stencil:  ---X---C---X---->
 */
int _1Node(Float& res,const st& idx, const Stencil_2D1& stc, Float x, Float y);
int _1Node(PData_2D& res, const Stencil_2D1& stc);

/*
 *  There are two nodes in the stencil, x is
 *  the location, somewhere on the stencil axes.
 *
 *  Stencil: ---X---C---O---->  or
 *           ---O---C---X---->
 *                ^
 *                x
 */
int _2NodeOnAxes(PData_2D& res, const Stencil_2D1& stc);

/*
 *  There are two nodes in the stencil, x is
 *  the location, somewhere on the stencil axes.
 *
 *  Stencil: ---C---O---O---->  or
 *           ---O---C---O---->  or
 *           ---O---O---C---->
 *                ^
 *                x
 */
int _3NodeOnAxes(Float& res, Stencil_2D1& stc, Float x);

void InterpolateOnFace_1Order(        // 2D QuadTree Node
		pNode_2D pn,                  //node
		const Direction& dir,                //face
		const ArrayListV<st>& arridx,            //data index
		ArrayListV<Vt>& arrres                //data res
		);
void InterpolateOnFace_1Order( // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,                //face
		const st& idx,            //data index
		Vt& res                //data res
		);
void InterpolateOnFace_2Order(        // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,                //face
		const ArrayListV<st>& arridx,            //data index
		ArrayListV<Vt>& arrres                //data res
		);
void InterpolateOnFace_2Order( // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,                //face
		const st& idx,            //data index
		Vt& res                //data res
		);
void Interpolate_1Order(        // 2D QuadTree Node
		pNode_2D pn,                  //node
		const Cvt& x, const Cvt& y,              //point
		const ArrayListV<st>& arridx,            //data index
		ArrayListV<Vt>& arrres                //data res
		);
}

#endif
