#ifndef _TEST_DOMAIN_HPP_
#define _TEST_DOMAIN_HPP_

#include "../carpio_define.hpp"
//#include "../geometry/geometry.hpp"
//#include "../utility/clipper.hpp"
#include "../domain/domain.hpp"
#include "../algebra/matrix_SparCompRow.hpp"
#include "../algebra/matrix_SparCompCol.hpp"
#include <iostream>
#include <cmath>

namespace carpio {
void test_domain_1() {
	//
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	//CreatCube(shape, 1.5, 1.5, 2.5, 2.5);
	// shape is out bound
	//
	// define unit length
	Float UL = 0.499;
	// build grid ------------------
	Float max_x = shape.max_x();
	Float max_y = shape.max_y();
	Float min_x = shape.min_x();
	Float min_y = shape.min_y();
	st n_x = std::ceil((max_x - min_x) / UL);
	st n_y = std::ceil((max_y - min_y) / UL);
	std::cout<<n_x << " "<<n_y<<"\n";
	Grid_<Float, Float, dim> g(n_x, min_x, UL, //
			n_y, min_y, UL);
	// g.show_info();
	// build adaptive
	Adaptive_<Float, Float, dim> adp(&g, 2, 2);
	std::cout<<" here0 ----\n";
	adp.adapt_bound_solid(shape);
	std::cout<<" here3 ----\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_RootNodes(ga, g);
	lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, g);
	//lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
	//delete shape
}
void test_domain_2() {
	//
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
    //CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	CreatCube(shape, 1.5, 1.5, 3.5, 3.5);
	// shape is out bound
	//
	// define unit length
	Float UL = 0.25;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 2, 3);
	domain.build();
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	//GnuplotActor_LeafNodes(ga, domain.grid());
	//lga.push_back(ga);
	GnuplotActor_GhostNodes(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_GhostNodesContour_BoundaryIndex(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);

	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
	//delete shape
}
}

#endif
