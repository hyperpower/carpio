#ifndef __TEST_STENCIL_H_
#define __TEST_STENCIL_H_

#include "../io/gnuplot.h"

#include "gtest/gtest.h"
#include <math.h>

#include "test_define.hpp"

namespace carpio {

TEST(stencil, index) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape, cir;
	Float x1 = 1.5, y1 = 1.5, r = 0.3;
	Float x0 = 1.5, y0 = 1.5, x3 = 3.5, y3 = 3.5;
	CreatCircle(cir, x1, y1, r, 359);
	CreatCube(shape, x0, y0, x3, y3);
	// define unit length
	Float UL = 0.5;
	// build grid ------------------
	typedef Domain_<Float, Float, dim> Domain;
	Domain_<Float, Float, dim> domain(&shape, UL, 3, 5);
	domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	domain.new_data(1, 0, 0, 0);
	domain.set_val(0, set_sin);
	//
	Domain::pNode pn = domain.grid().get_pnode(1.729, 1.661);

	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	//GnuplotActor_LeafNodesContours(ga, domain.grid(),0);
	//lga.push_back(ga);
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_GhostNodes(ga, domain.ghost());
	lga.push_back(ga);

	GnuplotActor_Node(ga, *pn);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);

	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	gp.set_palette_blue_red();
	gp.set_xrange(1.4, 2.0);
	gp.set_yrange(1.4, 2.0);
	gp.plot(lga);

}

}

#endif
