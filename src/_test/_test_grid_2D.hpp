#ifndef _TEST_GRID_2D_HPP_
#define _TEST_GRID_2D_HPP_

#include "../domain/domain.hpp"
#include "../carpio_define.hpp"
#include "../io/io_gnuplot_domain.h"

namespace carpio {

inline void test_grid_2d() {
	Grid_<Float, Float, 2> g(
			2, 0, 1, //
			2, 0, 1 );
	g.connect_root();
	std::cout<< "num root : "<<g.get_num_root()<<std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt();

	GnuplotShow_RootNodes(g);


}

}

#endif
