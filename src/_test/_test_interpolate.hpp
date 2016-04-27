#ifndef _TEST_INTERPOLATE_H_
#define _TEST_INTERPOLATE_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/interpolate.h"
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

TEST(Interpolate, direction) {
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

TEST(Interpolate, unigrid) {
	const st dim = 2;
	// new shape--------------------
	Shape2D shape, cir;
	CreatCircle(cir, 1.5, 1.5, 0.3, 359);
	CreatCube(shape, 1.5, 1.5, 3.5, 3.5);
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
	std::cout << "res   = " << exp.subsitute(0);
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

}

#endif
