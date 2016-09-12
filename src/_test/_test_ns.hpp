#ifndef __TEST_NS_H_
#define __TEST_NS_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/ns.hpp"
#include "../algebra/arithmetic.hpp"
#include "gtest/gtest.h"
#include "test_define.hpp"
#include <math.h>
#include "../io/matplot.h"
#include "../io/matplot_actor.h"

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

TEST(ns, DISABLED_unigrid) {
	std::cout << "test ns   \n";
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
	NS_<Float, Float, dim> ns(&domain, 3);
	domain.grid().show_info();
	// boundary condition
	NS_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(0);
	NS_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(1);
	NS_<Float, Float, dim>::BoundaryCondition bc2_0;
	// u
	ns.set_BC_velocity(0, 0, _X_, &bc0);
	ns.set_BC_velocity(0, 1, _X_, &bc2_0);
	ns.set_BC_velocity(0, 2, _X_, &bc0);
	ns.set_BC_velocity(0, 3, _X_, &bc1);
	// v
	ns.set_BC_velocity(0, 0, _Y_, &bc0);
	ns.set_BC_velocity(0, 1, _Y_, &bc0);
	ns.set_BC_velocity(0, 2, _Y_, &bc0);
	ns.set_BC_velocity(0, 3, _Y_, &bc0);
	// u t
	ns.set_BC_velocity_t(0, 0, _X_, &bc0);
	ns.set_BC_velocity_t(0, 1, _X_, &bc2_0);
	ns.set_BC_velocity_t(0, 2, _X_, &bc0);
	ns.set_BC_velocity_t(0, 3, _X_, &bc1);
	// v t
	ns.set_BC_velocity_t(0, 0, _Y_, &bc0);
	ns.set_BC_velocity_t(0, 1, _Y_, &bc0);
	ns.set_BC_velocity_t(0, 2, _Y_, &bc0);
	ns.set_BC_velocity_t(0, 3, _Y_, &bc0);
	// pressure
	ns.set_boundary_condition(0, 3, 0, &bc0);

	// initial condition
	ns.set_velo(set_v_0, _X_);
	ns.set_velo(set_v_0, _Y_);
	ns.set_mu(set_v_1);
	ns.set_rho(set_v_1);

	//show_plot_contour(domain, ns.idx_u(), "u field");
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	ns.advance(4, 1e-6, 1000, 1);
	show_plot_contour(domain, ns.idx_ut(), "u star field");
	// show
	// show_plot_contour(domain, ns.idx_ut());

	//show_plot_contour(domain, ns.idx_mu(), "mu field");
	//show_plot_contour(domain, ns.idx_rho(), "rho field");
	//show_plot_contour(domain, ns.idx_u(), "u field");
	//show_plot_contour(domain, ns.idx_v(), "v field");
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	//show_plot_contour(domain, ns.idx_vt(), "v star field");
	show_plot_contour(domain, ns.idx_p(), "p field");

	//show_value_on_line(domain, _X_, 0, ns.idx_u());
	show_veo_field(domain, ns.idx_u(), ns.idx_v());
}

TEST(ns2, unigrid) {
	std::cout << "test ns 2  cavity lid driven flow\n";
	int Re = 100;
	int mesh = 6;
	int maxstep = 500;
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = -0.5, y0 = -0.5, x1 = 0.5, y1 = 0.5;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, mesh, 10);
	domain.build();
	NS_<Float, Float, dim> ns(&domain, 1e-4, 3);
	// domain.grid().show_info();
	// boundary condition
	NS_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(0);
	NS_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(1);
	NS_<Float, Float, dim>::BoundaryCondition bc2_0;
	// u
	ns.set_BC_velocity(0, 0, _X_, &bc0);
	ns.set_BC_velocity(0, 1, _X_, &bc0);
	ns.set_BC_velocity(0, 2, _X_, &bc1);
	ns.set_BC_velocity(0, 3, _X_, &bc0);
	// v
	ns.set_BC_velocity(0, 0, _Y_, &bc0);
	ns.set_BC_velocity(0, 1, _Y_, &bc0);
	ns.set_BC_velocity(0, 2, _Y_, &bc0);
	ns.set_BC_velocity(0, 3, _Y_, &bc0);
	// u t
	ns.set_BC_velocity_t(0, 0, _X_, &bc0);
	ns.set_BC_velocity_t(0, 1, _X_, &bc0);
	ns.set_BC_velocity_t(0, 2, _X_, &bc0);
	ns.set_BC_velocity_t(0, 3, _X_, &bc0);
	// v t
	ns.set_BC_velocity_t(0, 0, _Y_, &bc0);
	ns.set_BC_velocity_t(0, 1, _Y_, &bc0);
	ns.set_BC_velocity_t(0, 2, _Y_, &bc0);
	ns.set_BC_velocity_t(0, 3, _Y_, &bc0);
	// initial condition
	ns.set_velo(set_v_0, _X_);
	ns.set_velo(set_v_0, _Y_);
	ns.set_rho(set_v_1);
	// =======================
	NS_<Float, Float, dim>::Function set_val =
			[Re](Float x, Float y, Float z) {return 1.0/float(Re);};
	ns.set_mu(set_val);

	//show_plot_contour(domain, ns.idx_u(), "u field");
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	ns.advance(maxstep, 1e-2, 30, 1);
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	show_value_on_line(domain, _Y_, 0.0, ns.idx_v());
	show_value_on_line(domain, _X_, 0.0, ns.idx_u());
	string dir = "/home/zhou/Dropbox/Paper/NS_test/1_Lid_Driven_Cavity/carpio/";
	string fxu = dir + "xu_re" + ToString(Re, mesh, "_");
	string fyv = dir + "yv_re" + ToString(Re, mesh, "_");
	//string fxu = "/home/zhou/Dropbox/Paper/NS_test/1_Lid_Driven_Cavity/carpio/xu_re100_4";
	//string fyv = "/home/zhou/Dropbox/Paper/NS_test/1_Lid_Driven_Cavity/carpio/yv_re100_4";
	//string fxu = "xu";
	//string fyv = "yv";
	output_value_on_line(fxu, domain, _X_, 0.0, ns.idx_u());
	output_value_on_line(fyv, domain, _Y_, 0.0, ns.idx_v());
	// show
	// show_plot_contour(domain, ns.idx_ut());

	//show_plot_contour(domain, ns.idx_mu(), "mu field");
	//show_plot_contour(domain, ns.idx_rho(), "rho field");
	//show_plot_contour(domain, ns.idx_u(), "u field");
	//show_plot_contour(domain, ns.idx_v(), "v field");
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	//show_plot_contour(domain, ns.idx_vt(), "v star field");
	//show_plot_contour(domain, ns.idx_p(), "p field");

	//show_value_on_line(domain, _X_, 0, ns.idx_u());
	show_veo_field(domain, ns.idx_u(), ns.idx_v());
	//Matplot p;
	//p.import_matplot();
	//Matplot_actors la;
	//MatplotActors_LeafNodes(la, domain.grid(), "color = \"blue\"");
	//MatplotActors_GhostNodes(la, domain.ghost(), "color = \"green\"");
	//Matplot_actor ma;
	//MatplotActor_GridCenterValue(ma, domain.grid(), ns.idx_ut(), "");
	//la.push_back(ma);
	//MatplotActor_GhostCenterValue(ma, domain.ghost(), ns.idx_vt(), "");
	//la.push_back(ma);
	//p.plot(la);
	//p.set_equal_ratio();
	//p.show();
}

Float U_field(Float x, Float y, Float z, Float t) {
	return -cos(x) * sin(y) * pow(EXP, -2.0 * t);
}

Float V_field(Float x, Float y, Float z, Float t) {
	return sin(x) * cos(y) * pow(EXP, -2.0 * t);
}

Float P_field(Float x, Float y, Float z, Float t) {
	return -0.25 * (cos(2.0 * x) + cos(2.0 * y)) * pow(EXP, -4.0 * t);
}

Float U_field0(Float x, Float y, Float z) {
	return U_field(x, y, z, 0.0);
}

Float V_field0(Float x, Float y, Float z) {
	return V_field(x, y, z, 0.0);
}

Float P_field0(Float x, Float y, Float z) {
	return P_field(x, y, z, 0.0);
}

TEST(ns3, DISABLED_unigrid) {
	std::cout << "test ns 3  two-dimensional unsteady flow\n";
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = -0.5, y0 = -0.5, x1 = 0.5, y1 = 0.5;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 4, 4);
	domain.build();
	NS_<Float, Float, dim> ns(&domain, 1e-4, 3);
	domain.grid().show_info();
	// boundary condition
	NS_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(U_field0);
	NS_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(V_field0);
	NS_<Float, Float, dim>::BoundaryCondition bc2_0;
	// u
	ns.set_BC_velocity(0, 0, _X_, &bc0);
	ns.set_BC_velocity(0, 1, _X_, &bc0);
	ns.set_BC_velocity(0, 2, _X_, &bc0);
	ns.set_BC_velocity(0, 3, _X_, &bc0);
	// v
	ns.set_BC_velocity(0, 0, _Y_, &bc1);
	ns.set_BC_velocity(0, 1, _Y_, &bc1);
	ns.set_BC_velocity(0, 2, _Y_, &bc1);
	ns.set_BC_velocity(0, 3, _Y_, &bc1);

	// initial condition
	ns.set_velo(U_field0, _X_);
	ns.set_velo(V_field0, _Y_);
	ns.set_rho(set_v_1);
	ns.set_mu(set_v_1);

	//show_plot_contour(domain, ns.idx_u(), "u field");
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	ns.advance(3500, 1e-2, 1000, 1);
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	show_value_on_line(domain, _Y_, 0.0, ns.idx_v());
	show_value_on_line(domain, _X_, 0.0, ns.idx_u());
	//output_value_on_line("xu", domain, _X_, 0.0, ns.idx_u());
	//output_value_on_line("yv", domain, _Y_, 0.0, ns.idx_v());
	// show
	// show_plot_contour(domain, ns.idx_ut());

	//show_plot_contour(domain, ns.idx_mu(), "mu field");
	//show_plot_contour(domain, ns.idx_rho(), "rho field");
	//show_plot_contour(domain, ns.idx_u(), "u field");
	//show_plot_contour(domain, ns.idx_v(), "v field");
	//show_plot_contour(domain, ns.idx_ut(), "u star field");
	//show_plot_contour(domain, ns.idx_vt(), "v star field");
	//show_plot_contour(domain, ns.idx_p(), "p field");

	//show_value_on_line(domain, _X_, 0, ns.idx_u());
	show_veo_field(domain, ns.idx_u(), ns.idx_v());
	//Matplot p;
	//p.import_matplot();
	//Matplot_actors la;
	//MatplotActors_LeafNodes(la, domain.grid(), "color = \"blue\"");
	//MatplotActors_GhostNodes(la, domain.ghost(), "color = \"green\"");
	//Matplot_actor ma;
	//MatplotActor_GridCenterValue(ma, domain.grid(), ns.idx_ut(), "");
	//la.push_back(ma);
	//MatplotActor_GhostCenterValue(ma, domain.ghost(), ns.idx_vt(), "");
	//la.push_back(ma);
	//p.plot(la);
	//p.set_equal_ratio();
	//p.show();
}

TEST(ns, DISABLED_unigrid2) {
	Matplot p;
	p.import_matplot();
	Matplot_actors la;
	Matplot_actor ma;
	MatplotActor_Segment2D(ma, 0.0, 0.0, 1.0, 1.0, "-r");
	la.push_back(ma);
	p.plot(la);
	p.show();
}

}

#endif
