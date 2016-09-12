#ifndef __TEST_EQUATION_H_
#define __TEST_EQUATION_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/equation.hpp"
#include "../algebra/arithmetic.hpp"
#include "gtest/gtest.h"
#include "test_define.hpp"
#include <math.h>
#include "../io/matplot.h"
#include "../io/matplot_actor.h"
#include "../calculation/event.hpp"

namespace carpio {

/*
 * this set get from
 * Shashank et al. Journal of Computational Physics 229 (2010) 4425-4430
 */
inline Float set_u(Float x, Float y, Float z) {
	return -cos(PI * x) * sin(PI * y);
}

inline Float set_v(Float x, Float y, Float z) {
	return sin(PI * x) * cos(PI * y);
}

TEST(equation, unigrid) {
	std::cout << "test equatioon   \n";
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = -0.5, y0 = -0.5, x1 = 0.5, y1 = 0.5;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 4, 5);
	domain.build();
	//
	Equation_<Float, Float, dim> ns(&domain, 3);
	domain.grid().show_info();


}







}

#endif
