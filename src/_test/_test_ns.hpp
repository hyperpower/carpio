#ifndef __TEST_NS_H_
#define __TEST_NS_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/ns.hpp"
#include "gtest/gtest.h"
#include <math.h>

namespace carpio {

TEST(ns, unigrid) {
	std::cout << "test ns   \n";
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = 0, y0 = 0, x1 = 1, y1 = 1;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 4, 5);
	domain.build();
	NS_<Float, Float, dim> ns(&domain, 3);
	domain.grid().show_info();
	ns.run();

}

}

#endif
