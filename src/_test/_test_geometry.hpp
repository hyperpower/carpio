#ifndef _TEST_GEOMETRY_HPP_
#define _TEST_GEOMETRY_HPP_

#include "../carpio_define.hpp"
#include "../geometry/geometry.hpp"

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
}

#endif
