
#ifndef _TEST_GEOMETRY_HPP_
#define _TEST_GEOMETRY_HPP_


#include "../carpio_define.hpp"
#include "../geometry/_point.hpp"
#include "../geometry/_line.hpp"
#include "../geometry/_plane.hpp"
#include "../geometry/_segment.hpp"


#include <iostream>

namespace carpio {
void test_geo(){
	Point_<int, 3> p(1,2,3);
	p.show();
	Segment_<Float, 3> seg;
	std::cout<<" =============End test==================\n";
}
}

#endif
