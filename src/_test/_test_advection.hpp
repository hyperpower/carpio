#ifndef __TEST_ADVECTION_H_
#define __TEST_ADVECTION_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/advection.hpp"
#include "gtest/gtest.h"
#include <math.h>

namespace carpio {
Float v_half1half0(Float x, Float y, Float z) {
	if (x < y) {
		return 1;
	} else {
		return 0;
	}
}
Float set_v(Float x, Float y, Float z) {
	return 1.0;
}

void show_x05(Domain_<Float, Float, 2>& domain) {
	typedef Domain_<Float, Float, 2> Domain;
	typedef Domain_<Float, Float, 2>::pNode pNode;
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 4);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, domain.grid(), 4);
	//lga.push_back(ga);
	std::list<pNode> lpn = domain.grid().get_leaf(_Y_, 0.5);
	GnuplotActor_NodesContours(ga, lpn, 4);
	lga.push_back(ga);
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	//GnuplotActor_Shape2D(ga, shape, 0);
	//lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
}

void show_x05_line(Domain_<Float, Float, 2>& domain) {
	typedef Domain_<Float, Float, 2> Domain;
	typedef Domain_<Float, Float, 2>::pNode pNode;
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;

	//GnuplotActor_GhostNodesContours(ga, domain.ghost(), 4);
	//lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, domain.grid(), 4);
	//lga.push_back(ga);
	std::list<pNode> lpn = domain.grid().get_leaf(_Y_, 0.5);
	GnuplotActor_NodesValues(ga,lpn, 4 , _X_);
	//GnuplotActor_NodesContours(ga, lpn, 4);
	lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, domain.grid());
	//lga.push_back(ga);
	//GnuplotActor_Shape2D(ga, shape, 0);
	//lga.push_back(ga);
	Gnuplot gp;
	//gp.set_equal_ratio();
	gp.plot(lga);
}

// The test case is based on
// M.S. Darwish, F. Moukalled , TVD schemes for unstructured grids
TEST(DISABLED_Advection, unigrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = 0, y0 = 0, x1 = 1, y1 = 1;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 5, 8);
	domain.build();
	Advection_<Float, Float, dim> advection(&domain, 3);
	advection.set_v(set_v, _X_);
	advection.set_v(set_v, _Y_);
	domain.grid().show_info();
	// boundary condition
	Advection_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(0);
	Advection_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(1);
	advection.set_boundary_condition(0, 0, 4, &bc0);
	advection.set_boundary_condition(0, 1, 4, &bc0);
	advection.set_boundary_condition(0, 2, 4, &bc1);
	advection.set_boundary_condition(0, 3, 4, &bc1);
	advection.set_phi_ghost();
	cout << "solve -----------\n";
	advection.solve_tvd(1e-6, 1000, 1);
	cout << "end solve -------\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 4);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(), 4);
	lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, domain.grid());
	//lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.plot(lga);
	//delete shape
	// splot show ================================
	std::list<Gnuplot_actor> slga;
	Gnuplot_actor sga;
	GnuplotActor_LeafNodesSurface(sga, domain.grid(), 4);
	//GnuplotActor_NodesSurface(sga, pn, 0);
	slga.push_back(sga);
	GnuplotActor_GhostNodesSurface(sga, domain.ghost(), 4);
	slga.push_back(sga);
	//sga.show_data();
	Gnuplot sgp;
	//sgp.set_equal_ratio();
	sgp.set_view(45, 10, 1, 1);
	sgp.set_palette_blue_red();
	sgp.set("ticslevel 0");
	//sgp.set_xrange(1.4, 2.0);
	//sgp.set_yrange(1.4, 2.0);
	sgp.splot(slga);
	//
	show_x05_line(domain);
}

TEST(Advection, adpgrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = 0, y0 = 0, x1 = 1, y1 = 1;
	CreatCube(shape, x0, y0, x1, y1);
	Float xr = 0.3, yr =0.3, r=0.3;
	Shape2D cir;
	CreatCircle(cir, xr, yr, r, 100);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 5, 8);
	domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	Advection_<Float, Float, dim> advection(&domain, 3);
	advection.set_v(set_v, _X_);
	advection.set_v(set_v, _Y_);
	domain.grid().show_info();
	// boundary condition
	Advection_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(0);
	Advection_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(1);
	advection.set_boundary_condition(0, 0, 4, &bc0);
	advection.set_boundary_condition(0, 1, 4, &bc0);
	advection.set_boundary_condition(0, 2, 4, &bc1);
	advection.set_boundary_condition(0, 3, 4, &bc1);
	advection.set_phi_ghost();
	cout << "solve -----------\n";
	advection.solve_tvd(1e-6, 1000, 1);
	cout << "end solve -------\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 4);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(), 4);
	lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, domain.grid());
	//lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.plot(lga);
	//delete shape
	// splot show ================================
	std::list<Gnuplot_actor> slga;
	Gnuplot_actor sga;
	GnuplotActor_LeafNodesSurface(sga, domain.grid(), 4);
	//GnuplotActor_NodesSurface(sga, pn, 0);
	slga.push_back(sga);
	GnuplotActor_GhostNodesSurface(sga, domain.ghost(), 4);
	slga.push_back(sga);
	//sga.show_data();
	Gnuplot sgp;
	//sgp.set_equal_ratio();
	sgp.set_view(45, 10, 1, 1);
	sgp.set_palette_blue_red();
	sgp.set("ticslevel 0");
	//sgp.set_xrange(1.4, 2.0);
	//sgp.set_yrange(1.4, 2.0);
	sgp.splot(slga);
	//
	show_x05_line(domain);
}

}
#endif

