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

/*
 *  There are two nodes in the stencil, x is
 *  the location, somewhere on the stencil axes.
 *
 *  Stencil: ---X---C---O---->  or
 *           ---O---C---X---->
 *                ^
 *                x
 */
int _2NodeOnAxes(Float& res, Stencil_2D1& stc, Float x);

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

}

#endif
