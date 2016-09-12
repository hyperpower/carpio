#ifndef _TEST_ADVECTION_H_
#define _TEST_ADVECTION_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/advection.hpp"
#include "test_define.hpp"
#include "gtest/gtest.h"
#include <math.h>
using namespace std;
namespace carpio {



Float exact_case1(Float x, Float y, Float z) {
	if (y >= x) {
		return 1;
	} else {
		return 0;
	}
}

Float exact_case2(Float x, Float y, Float z) {
	if (y > x && y < (x + 0.3)) {
		return 1;
	} else {
		return 0;
	}
}

void test_1_run(int level, Float& e1, Float& e2, Float& e3) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = 0, y0 = 0, x1 = 1, y1 = 1;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, level, level + 1);
	domain.build();
	Advection_<Float, Float, dim> advection(&domain, 3);
	advection.set_v(set_v_1, _X_);
	advection.set_v(set_v_1, _Y_);
	// domain.grid().show_info();
	// boundary condition
	Advection_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(0);
	Advection_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(1);
	advection.set_boundary_condition(0, 0, 2, &bc0);
	advection.set_boundary_condition(0, 1, 2, &bc0);
	advection.set_boundary_condition(0, 2, 2, &bc1);
	advection.set_boundary_condition(0, 3, 2, &bc1);
	advection.set_phi_ghost();
	cout << "solve -----------\n";
	advection.solve_tvd(1e-6, 1000, 1);
	cout << "end solve -------\n";
	// exact ===============================
	domain.resize_data(4, 0, 0, 1);
	domain.set_val(3, exact_case1);
	// error ===============================
	e1 = error_1(domain, 2, 3);
	e2 = error_2(domain, 2, 3);
	e3 = error_i(domain, 2, 3);
	// splot show ================================
	//show_splot_surface(domain, 2);
}

void test_2_run(int level, Float& e1, Float& e2, Float& e3) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x0 = 0, y0 = 0, x1 = 1, y1 = 1;
	CreatCube(shape, x0, y0, x1, y1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, level, level);
	domain.build();
	Advection_<Float, Float, dim> advection(&domain, 2);
	advection.set_v(set_v_1, _X_);
	advection.set_v(set_v_1, _Y_);
	// domain.grid().show_info();
	// boundary condition
	Advection_<Float, Float, dim>::BoundaryCondition bc_out;
	Advection_<Float, Float, dim>::BoundaryCondition bc0;
	bc0.set_default_1_bc(0);
	Advection_<Float, Float, dim>::BoundaryCondition bc1;
	std::function<Float(Float, Float, Float)> fun =
			[](Float x, Float y,Float z) {
				if(y<0.3) {
					return 1;
				} else {
					return 0;
				}
			};
	bc1.set_default_1_bc(fun);
	advection.set_boundary_condition(0, 0, 2, &bc0);
	advection.set_boundary_condition(0, 1, 2, &bc_out);
	advection.set_boundary_condition(0, 2, 2, &bc_out);
	advection.set_boundary_condition(0, 3, 2, &bc1);
	advection.set_phi_ghost();
	cout << "solve -----------\n";
	advection.solve_tvd(1e-6, 1000, 1);
	cout << "end solve -------\n";
	// exact ===============================
	domain.resize_data(4, 0, 0, 1);
	domain.set_val(3, exact_case2);
	// show ================================
	// error ===============================
	e1 = error_1(domain, 2, 3);
	e2 = error_2(domain, 2, 3);
	e3 = error_i(domain, 2, 3);
	// splot show ================================
	//show_splot_surface(domain, 2);
}

void output_result( //
		Domain& d, // Domain
		string& outputf, //the output folder
		int c,   //test case number
		int s1,  //scheme major version
		int s2   //scheme second version
		){
	//

}

TEST(Advection, unigrid) {
	int n = 5;
	int bl = 3;
	arrayList_st arr_level(n);
	arrayList arr_norm1(n);
	arrayList arr_norm2(n);
	arrayList arr_normi(n);
	for (int i = 0; i < n; i++) {
		arr_level[i] = bl + i;
	}
	// put the function here ==============================
	for (int i = 0; i < n; i++) {
		std::cout << "Mesh level" << arr_level[i] << "\n";
		test_2_run(arr_level[i], arr_norm1[i], arr_norm2[i], arr_normi[i]);
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
