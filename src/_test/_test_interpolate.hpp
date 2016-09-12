#ifndef _TEST_INTERPOLATE_H_
#define _TEST_INTERPOLATE_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../io/io_gnuplot_domain.h"
#include "../calculation/interpolate.h"
#include "../utility/format.h"
#include "test_define.hpp"
#include "gtest/gtest.h"

namespace carpio {

Float set_sin(Float x, Float y, Float z) {
	//center
	Float ox = 0.5;
	Float oy = 0.5;
	Float dis = Distance(ox, oy, x, y);
	//return 1;
	return sin(3 * dis);
}

TEST(DISABLED_Interpolate, direction) {
	EXPECT_GE(true, is_on_direction(5, _ZP_));
	EXPECT_GE(true, is_on_direction(1, _XP_));
	EXPECT_GE(true, is_on_direction(3, _XP_));
	Direction d = 56;
	std::cout << ToString(d) << "\n";
	EXPECT_GE(true, is_on_direction(0, d));
	EXPECT_GE(false, is_on_direction(1, d));
	EXPECT_GE(false, is_on_direction(5, d));
	EXPECT_GE(false, is_on_direction(7, d));

	Direction df = ToFaceDirection(_M_, _X_);
	std::cout << ToString(df) << "\n";
	df = ToFaceDirection(_P_, _X_);
	std::cout << ToString(df) << "\n";
	df = ToFaceDirection(_M_, _Y_);
	std::cout << ToString(df) << "\n";
	df = ToFaceDirection(_M_, _Z_);
	std::cout << ToString(df) << "\n";
	df = ToFaceDirection(_P_, _Z_);
	std::cout << ToString(df) << "\n";

}

TEST(Interpolate, DISABLED_unigrid) {
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
	Domain::pNode pn = domain.grid().get_pnode(1.79, 1.5);
	typedef Interpolate_<Float, Float, dim> Interpolate;
	typedef Expression_<Float, Float, dim> Exp;
	Axes a = _X_;
	Orientation o = _M_;
	Float dis = 0.001;
	Direction d = ToFaceDirection(o, a);
	Exp exp = Interpolate::OnAxes(pn, d, dis, 2);
	exp.show();
	Float x = pn->cp(_X_) + ((a == _X_) ? dis : 0.0);
	Float y = pn->cp(_Y_) + ((a == _Y_) ? dis : 0.0);
	std::cout << " x " << x << " y " << y << "\n";
	std::cout << "res   = " << exp.substitute(0);
	std::cout << "   acc   = " << set_sin(x, y, 0.0) << "\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	//GnuplotActor_LeafNodesContours(ga, domain.grid(),0);
	//lga.push_back(ga);
	GnuplotActor_Expression(ga, exp);
	lga.push_back(ga);
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
	gp.set_xrange(1.4, 2.0);
	gp.set_yrange(1.4, 2.0);
	gp.plot(lga);
}

TEST(DISABLED_Interpolate, corner) {
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
	typedef Interpolate_<Float, Float, dim> Interpolate;
	typedef Expression_<Float, Float, dim> Exp;
	Axes a1 = _X_;
	Orientation o1 = _M_;
	Axes a2 = _Y_;
	Orientation o2 = _M_;
	Direction d = ToCornerDirection(o1, a1, o2, a2);
	Exp exp = Interpolate::OnCorner_1Order(pn, d);
	exp.show();
	Float x = pn->p(o1, _X_);
	Float y = pn->p(o2, _Y_);
	std::cout << "   x " << x << " y " << y << "\n";
	std::cout << "   res   = " << exp.substitute(0);
	std::cout << "   acc   = " << set_sin(x, y, 0.0) << "\n";
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	//GnuplotActor_LeafNodesContours(ga, domain.grid(),0);
	//lga.push_back(ga);
	GnuplotActor_Expression(ga, exp);
	lga.push_back(ga);
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

	// show ================================
	std::list<Gnuplot_actor> slga;
	Gnuplot_actor sga;
	GnuplotActor_LeafNodesSurface(sga, domain.grid(), 0);
	//GnuplotActor_NodesSurface(sga, pn, 0);
	slga.push_back(sga);
	GnuplotActor_GhostNodesSurface(sga, domain.ghost(), 0);
	slga.push_back(sga);
	//sga.show_data();
	Gnuplot sgp;
	//sgp.set_equal_ratio();
	sgp.set_view(45, 0, 1, 1);
	sgp.set_palette_blue_red();
	//sgp.set_xrange(1.4, 2.0);
	//sgp.set_yrange(1.4, 2.0);
	sgp.splot(slga);
}

TEST(Interpolate, plane) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape, cir;
	Float x1 = 1.6, y1 = 1.6, r = 0.3;
	Float x0 = 1.5, y0 = 1.5, x3 = 3.5, y3 = 3.5;
	CreatCircle(cir, x1, y1, r, 359);
	CreatCube(shape, x0, y0, x3, y3);
	// define unit length
	Float UL = 0.5;
	// build grid ------------------
	typedef Domain_<Float, Float, dim> Domain;
	Domain_<Float, Float, dim> domain(&shape, UL, 4, 6);
	domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	domain.new_data(1, 0, 0, 0);
	domain.set_val(0, set_sin);
	// generate points to get error
	int nmax = 1000000;
	ArrayListV<Float> arrx(nmax);
	ArrayListV<Float> arry(nmax);
	ArrayListV<Float> err(nmax);
	srand((unsigned) time(0));
	for (int i = 0; i < nmax; i++) {
		Float x, y;
		x = arrx[i] = generate_random_number(x0, x3);
		y = arry[i] = generate_random_number(y0, y3);
		//fmt::print("i   = {:>9d} x   = {:>9.5f} y   = {:>9.5f}\n", i, arrx[i], arry[i]);
		Domain::pNode pn = domain.grid().get_pnode(x, y);
		typedef Interpolate_<Float, Float, dim> Interpolate;
		typedef Expression_<Float, Float, dim> Exp;
		Exp&& exp = Interpolate::OnPlane(pn, _X_, _Y_, arrx[i], arry[i], 1);
		Float res = exp.substitute(0);
		Float acc = set_sin(x, y, 0.0);
		//fmt::print("res = {:>9.5f} acc = {:>9.5f} err = {:>9.5f}\n", res, acc, acc-res);
		err[i] = acc - res;
	}

	Float e1 = nrm1(err);
	Float e2 = nrm2(err);
	Float e3 = nrminf(err);

	fmt::print("e1   = {:>9.5e} \ne2   = {:>9.5e} \neinf = {:>9.5e}\n", e1, e2,
			e3);
	ASSERT_LE(e2, 1e-1);

}

void test_order(int minlevel, int dl, Float& e1, Float& e2, Float &einf) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape, cir;
	Float x1 = 1.6, y1 = 1.6, r = 0.3;
	Float x0 = 1.5, y0 = 1.5, x3 = 3.5, y3 = 3.5;
	CreatCircle(cir, x1, y1, r, 359);
	CreatCube(shape, x0, y0, x3, y3);
	// define unit length
	Float UL = 0.5;
	// build grid ------------------
	typedef Domain_<Float, Float, dim> Domain;
	Domain_<Float, Float, dim> domain(&shape, UL, minlevel, minlevel + dl);
	domain.adaptive().adapt_shape_boundary(cir);
	domain.build();
	domain.new_data(1, 0, 0, 0);
	domain.set_val(0, set_sin);
	// generate points to get error
	int nmax = 1000;
	ArrayListV<Float> arrx(nmax);
	ArrayListV<Float> arry(nmax);
	ArrayListV<Float> err(nmax);
	srand((unsigned) time(0));
	for (int i = 0; i < nmax; i++) {
		Float x, y;
		x = arrx[i] = generate_random_number(x0, x3);
		y = arry[i] = generate_random_number(y0, y3);
		//fmt::print("i   = {:>9d} x   = {:>9.5f} y   = {:>9.5f}\n", i, arrx[i], arry[i]);
		Domain::pNode pn = domain.grid().get_pnode(x, y);
		typedef Interpolate_<Float, Float, dim> Interpolate;
		typedef Expression_<Float, Float, dim> Exp;
		// get two face direction
		Orientation o1 = ((x - pn->cp(_X_)) >= 0) ? _P_ : _M_;
		Direction d1 = ToFaceDirection(o1, _X_);
		// the funciton is here ----------------------------------------
        Exp exp = Interpolate::OnPlane(pn, _X_, _Y_, arrx[i], arry[i], 1);
		//Exp exp = Interpolate::OnCorner_1Order();
		//--------------------------------------------------------------
		Float res = exp.substitute(0);
		Float acc = set_sin(x, y, 0.0);
		//fmt::print("res = {:>9.5f} acc = {:>9.5f} err = {:>9.5f}\n", res, acc, acc-res);
		err[i] = acc - res;
	}
	e1 = nrm1(err);
	e2 = nrm2(err);
	einf = nrminf(err);
}

TEST(DISABLED_Interpolate, plane_order) {
	int n = 5;
	int bl = 3;
	arrayList_st arr_level(n);
	arrayList arr_norm1(n);
	arrayList arr_norm2(n);
	arrayList arr_normi(n);
	for (int i = 0; i < n; i++) {
		arr_level[i] = bl + i;
	}
	int dl = 5;
	// put the function here
	for (int i = 0; i < n; i++) {
		fmt::print("Working on level: {:>5d} --> {:>5d}\n", arr_level[i],
				arr_level[i] + dl);
		test_order(arr_level[i], dl, arr_norm1[i], arr_norm2[i], arr_normi[i]);
		fmt::print("e1   = {:>9.5e} e2   = {:>9.5e} einf = {:>9.5e}\n",
				arr_norm1[i], arr_norm2[i], arr_normi[i]);
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
