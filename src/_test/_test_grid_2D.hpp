#ifndef _TEST_GRID_2D_HPP_
#define _TEST_GRID_2D_HPP_

#include "../domain/domain.hpp"
#include "../carpio_define.hpp"

namespace carpio {

inline void test_grid_2d() {
	Grid_<Float, Float, 2> g(
			1, 0, 1, //
			1, 0, 1 );
	g.connect_root();

	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt();


}

}

#endif
