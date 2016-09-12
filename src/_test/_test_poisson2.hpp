#ifndef __TEST_POISSON_H_
#define __TEST_POISSON_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/poisson.hpp"
#include "gtest/gtest.h"
#include "../io/io_gnuplot_domain.h"
#include <math.h>

namespace carpio {
Float coe_set_b(Float x, Float y, Float z) {
	return 1;
}
Float coe_set_f(Float x, Float y, Float z) {
	if (IsInRange(-0.25, x, 0.0, _oo_)) {
		return 1;
	}
	return 0;
}

TEST(Poisson, unigrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	Float x1 = -0.5, y1 = -0.5, x2 = 0.5, y2 = 0.5;
	CreatCube(shape, x1, y1, x2, y2);
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 2, 10);
	domain.build();
	//Equation_<Float, Float, dim> (&domain, 0.1, 10);
	//equation.set_output_time(0,10,1);
	//equation.run();
	Poisson_<Float, Float, dim> poisson(&domain, 0.2, 10);
	poisson.set_output_time(0,10, 2);
	poisson.run();

}

}
#endif
