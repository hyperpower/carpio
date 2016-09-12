#ifndef _TEST_POISSON_H_
#define _TEST_POISSON_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "test_define.hpp"
#include "gtest/gtest.h"
#include <math.h>
using namespace std;
namespace carpio {


Float coe_set_b(Float x, Float y, Float z) {
	return 1;
}

/*
 * this is for test case 2
 */
Float f_fun_0(Float x, Float y, Float z) {
	return 0;
}
Float f_fun_2(Float x, Float y, Float z) {
	Float k = 3;
	Float l = 3;
	Float pi = PI;
	//return 0;
	return -pi * pi * (k * k + l * l) * sin(pi * k * x) * sin(pi * l * y);
}

Float exact_fun_2(Float x, Float y, Float z) {
	Float k = 3;
	Float l = 3;
	Float pi = PI;
	//return 0;
	return sin(pi * k * x) * sin(pi * l * y);
}





void test_2_run(int level, Float& e1, Float& e2, Float& e3) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	Shape2D cir;
	//CreatCircle(cir, 2.1, 2.1, 0.8, 359);
	Float x1=-0.5, y1=-0.5, x2=0.5, y2=0.5;
	CreatCube(shape, x1, y1, x2, y2);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, level, level + 1);
	//domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	//domain.new_data(4,0,0,0); // idx = 3 used for exact
	//domain.set_val(3,exact_fun_2);
	Poisson_<Float, Float, dim> poisson(&domain, 1, 2, 3);
	poisson.set_beta_all(coe_set_b);
	poisson.set_f_all(f_fun_2);
	// set exact
	domain.set_val(0, exact_fun_2);
	// boundary condition
	Poisson_<Float, Float, dim>::BoundaryCondition bc;
	//bc.set_default_1_bc(exact_fun_2);
	poisson.set_boundary_condition(0, 0, 1, &bc);
	poisson.set_boundary_condition(0, 1, 1, &bc);
	poisson.set_boundary_condition(0, 2, 1, &bc);
	poisson.set_boundary_condition(0, 3, 1, &bc);
	poisson.set_phi_ghost();
	std::cout << "solve -----------\n";
	poisson.solve();
	cout << "end solve -------\n";
	e1 = error_1(domain, 2, 0);
	e2 = error_2(domain, 2, 0);
	e3 = error_i(domain, 2, 0);
	//cout << "error 1  " << e1 << "\n";
	//cout << "error 2  " << e2 << "\n";
	//cout << "error 3  " << e3 << "\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 0);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(), 2);
	lga.push_back(ga);
	//GnuplotActor_LeafNodes(ga, domain.grid());
	//lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	//gp.set_cbrange(-2.0, 2.0);
	//gp.plot(lga);
	//delete shape
}

void test_2_run2(int level, Float& e1, Float& e2, Float& e3) {
	const st dim = 2;
	// new shape--------------------
	Shape_<Float, 2> shape;
	Shape_<Float, 2> cir;
	Float xo=0.0, yo=0.0, r =0.3;
	CreatCircle(cir, xo, yo, r, 5);
	Float x1=-0.5, y1=-0.5, x2=0.5, y2=0.5;
	CreatCube(shape, x1, y1, x2, y2);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, level, level + 2);
	domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	//domain.new_data(4,0,0,0); // idx = 3 used for exact
	//domain.set_val(3,exact_fun_2);
	Poisson_<Float, Float, dim> poisson(&domain, 1, 2, 3);
	poisson.set_beta_all(coe_set_b);
	poisson.set_f_all(f_fun_2);
	// set exact
	domain.set_val(0, exact_fun_2);
	// boundary condition
	Poisson_<Float, Float, dim>::BoundaryCondition bc,bc2;
	bc2.set_default_2_bc(1);
	poisson.set_boundary_condition(0, 0, 1, &bc);
	poisson.set_boundary_condition(0, 1, 1, &bc);
	poisson.set_boundary_condition(0, 2, 1, &bc);
	poisson.set_boundary_condition(0, 3, 1, &bc);
	poisson.set_phi_ghost();
	std::cout << "solve -----------\n";
	poisson.solve();
	cout << "end solve -------\n";
	e1 = error_1(domain, 2, 0);
	e2 = error_2(domain, 2, 0);
	e3 = error_i(domain, 2, 0);
	//cout << "error 1  " << e1 << "\n";
	//cout << "error 2  " << e2 << "\n";
	//cout << "error 3  " << e3 << "\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_GhostNodesContours(ga, domain.ghost(), 0);
	lga.push_back(ga);
	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
	//lga.push_back(ga);
	GnuplotActor_LeafNodesContours(ga, domain.grid(), 2);
	lga.push_back(ga);
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	//gp.set_cbrange(-2.0, 2.0);
	gp.plot(lga);
	//delete shape
}

TEST(Poisson, DISABLED_unigrid) {
	int n = 5;
	int bl = 3;
	arrayList_st arr_level(n);
	arrayList arr_norm1(n);
	arrayList arr_norm2(n);
	arrayList arr_normi(n);
	for (int i = 0; i < n; i++) {
		arr_level[i] = bl + i;
	}
	// put the function here
	for (int i = 0; i < n; i++) {
		test_2_run(arr_level[i], arr_norm1[i], arr_norm2[i],
				arr_normi[i]);
	}
	//output all result
	cout << "error output =======================" << endl;
	cout << "  i         n1         n2         ni" << endl;
	for (int ii = 0; ii < n; ii++) {
		cout << std::scientific;
		cout << bl + ii << "  ";
		cout << arr_norm1[ii] << "  ";
		cout << arr_norm2[ii] << "  ";
		cout << arr_normi[ii] << endl;
	}
	cout << "order output =======================" << endl;
	cout << "  i         n1o         n2o         nio" << endl;
	Float sumn1 = 0.0;
	Float sumn2 = 0.0;
	Float sumni = 0.0;
	for (int ii = 1; ii < n; ii++) {
		cout << " " << bl + ii - 1 << "-" << bl + ii << "  ";
		cout << cal_order(arr_norm1[ii - 1], arr_norm1[ii]) << "  ";
		cout << cal_order(arr_norm2[ii - 1], arr_norm2[ii]) << "  ";
		cout << cal_order(arr_normi[ii - 1], arr_normi[ii]) << "  " << endl;
		sumn1 += cal_order(arr_norm1[ii - 1], arr_norm1[ii]);
		sumn2 += cal_order(arr_norm2[ii - 1], arr_norm2[ii]);
		sumni += cal_order(arr_normi[ii - 1], arr_normi[ii]);
	}
	cout << " ================================ " << endl;
	cout << "aver  ";
	cout << sumn1 / (n - 1) << "  ";
	cout << sumn2 / (n - 1) << "  ";
	cout << sumni / (n - 1) << "  " << endl;
}

TEST(Poisson, adpgrid) {
	int n = 5;
	int bl = 3;
	arrayList_st arr_level(n);
	arrayList arr_norm1(n);
	arrayList arr_norm2(n);
	arrayList arr_normi(n);
	for (int i = 0; i < n; i++) {
		arr_level[i] = bl + i;
	}
	// put the function here
	for (int i = 0; i < n; i++) {
		test_2_run2(arr_level[i], arr_norm1[i], arr_norm2[i],
				arr_normi[i]);
	}
	//output all result
	cout << "error output =======================" << endl;
	cout << "  i         n1         n2         ni" << endl;
	for (int ii = 0; ii < n; ii++) {
		cout << std::scientific;
		cout << bl + ii << "  ";
		cout << arr_norm1[ii] << "  ";
		cout << arr_norm2[ii] << "  ";
		cout << arr_normi[ii] << endl;
	}
	cout << "order output =======================" << endl;
	cout << "  i         n1o         n2o         nio" << endl;
	Float sumn1 = 0.0;
	Float sumn2 = 0.0;
	Float sumni = 0.0;
	for (int ii = 1; ii < n; ii++) {
		cout << " " << bl + ii - 1 << "-" << bl + ii << "  ";
		cout << cal_order(arr_norm1[ii - 1], arr_norm1[ii]) << "  ";
		cout << cal_order(arr_norm2[ii - 1], arr_norm2[ii]) << "  ";
		cout << cal_order(arr_normi[ii - 1], arr_normi[ii]) << "  " << endl;
		sumn1 += cal_order(arr_norm1[ii - 1], arr_norm1[ii]);
		sumn2 += cal_order(arr_norm2[ii - 1], arr_norm2[ii]);
		sumni += cal_order(arr_normi[ii - 1], arr_normi[ii]);
	}
	cout << " ================================ " << endl;
	cout << "aver  ";
	cout << sumn1 / (n - 1) << "  ";
	cout << sumn2 / (n - 1) << "  ";
	cout << sumni / (n - 1) << "  " << endl;
}

}
#endif

