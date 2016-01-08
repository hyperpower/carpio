#ifndef _TEST_GEOMETRY_HPP_
#define _TEST_GEOMETRY_HPP_

#include "../carpio_define.hpp"
#include "../geometry/geometry.hpp"
#include "../utility/clipper.hpp"
#include "../domain/shape.hpp"

#include <iostream>

namespace carpio {
void test_geo() {
	Point_<Float, 2> p1(1, 3);
	Point_<Float, 2> p2(0, 2);
	Point_<Float, 2> p3(5, 5);

	p1.show();
	p2.show();
	Segment_<Float, 2> seg1(p1, p2);
	Segment_<Float, 2> seg2(p2, p3);
	std::cout << "Is intersect : " << IsIntersect(seg1, seg2) << "\n";
	std::cout << " =============End test==================\n";
}

void test_clipper() {
	Shape2D cir, cube;
	CreatCircle(cir, 0.0, 0.0, 1.5, 359);
	CreatCube(cube, 0.0,0.0,1.0, 1.0);
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Shape2D(ga, cir);
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, cube);
	lga.push_back(ga);
	//clip
	Shape2D res;
	Intersect(cir,cube, res);
	GnuplotActor_Shape2D(ga, res);
	lga.push_back(ga);
	GnuplotShow(lga);
	std::cout<<std::scientific<<res.volume()<<"\n";
}


}

#endif
