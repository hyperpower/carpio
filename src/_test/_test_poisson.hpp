#ifndef __TEST_POISSON_H_
#define __TEST_POISSON_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/poisson.hpp"
#include "gtest/gtest.h"
#include <math.h>

namespace carpio {
Float coe_set_b(Float x, Float y, Float z) {
	return 1;
}
Float coe_set_f(Float x, Float y, Float z) {
	return 0;
}

TEST(Poisson, unigrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	CreatCube(shape, -0.5, -0.5, 0.5, 0.5);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 4, 4);
	domain.build();
	Poisson_<Float, Float, dim> poisson(&domain, 0, 1, 2);
	poisson.set_beta_all(coe_set_b);
	poisson.set_f_all(coe_set_f);
	// boundary condition
	Poisson_<Float, Float, dim>::BoundaryCondition bc;
	bc.set_default_1_bc(10);
	Poisson_<Float, Float, dim>::BoundaryCondition bc2;
	bc2.set_default_1_bc(0);
	poisson.set_boundary_condition(0, 0, 1, &bc);
	poisson.set_boundary_condition(0, 1, 1, &bc2);
	poisson.set_boundary_condition(0, 2, 1, &bc);
	poisson.set_boundary_condition(0, 3, 1, &bc2);
	poisson.set_phi_ghost();
	cout << "solve -----------\n";
	poisson.solve();
	cout << "end solve -------\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 1);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(), 1);
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
	//delete shape
}

TEST(DISABLED_Poisson, adpgrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	Shape2D cir;
	CreatCircle(cir, 0.00, 0.0, 0.3, 359);
	CreatCube(shape, -0.5, -0.5, 0.5, 0.5);
	// define unit length
	Float UL = 0.5;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 3, 5);
	domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	Poisson_<Float, Float, dim> poisson(&domain, 0, 1, 2);
	poisson.set_beta_all(coe_set_b);
	poisson.set_f_all(coe_set_f);
	// boundary condition
	Poisson_<Float, Float, dim>::BoundaryCondition bc;
	bc.set_default_1_bc(10);
	Poisson_<Float, Float, dim>::BoundaryCondition bc2;
	bc2.set_default_1_bc(0);
	poisson.set_boundary_condition(0, 0, 1, &bc);
	poisson.set_boundary_condition(0, 1, 1, &bc2);
	poisson.set_boundary_condition(0, 2, 1, &bc);
	poisson.set_boundary_condition(0, 3, 1, &bc2);
	poisson.set_phi_ghost();
	std::cout << "solve -----------\n";
	poisson.solve();
	cout << "end solve -------\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 1);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(),1);
	lga.push_back(ga);
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	gp.plot(lga);
	//delete shape
}

Float f_fun(Float x, Float y, Float z) {
	Float k = 3;
	Float l = 3;
	Float pi = _PI_;
	//return 0;
	return -pi * pi * (k * k + l * l) * sin(pi * k * x) * sin(pi * l * y);
}

TEST(DISABLED_Poisson, test2_unigrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	Shape2D cir;
	//CreatCircle(cir, 2.1, 2.1, 0.8, 359);
	CreatCube(shape, -0.5, -0.5, 0.5, 0.5);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 4, 5);
	//domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	Poisson_<Float, Float, dim> poisson(&domain, 0, 1, 2);
	poisson.set_beta_all(coe_set_b);
	poisson.set_f_all(f_fun);
	// boundary condition
	Poisson_<Float, Float, dim>::BoundaryCondition bc;
	bc.set_default_1_bc(0);
	poisson.set_boundary_condition(0, 0, 1, &bc);
	poisson.set_boundary_condition(0, 1, 1, &bc);
	poisson.set_boundary_condition(0, 2, 1, &bc);
	poisson.set_boundary_condition(0, 3, 1, &bc);
	poisson.set_phi_ghost();
	std::cout << "solve -----------\n";
	poisson.solve();
	cout << "end solve -------\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 1);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(),1);
	lga.push_back(ga);
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	gp.plot(lga);
	//delete shape
}

}
#endif
