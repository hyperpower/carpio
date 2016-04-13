#ifndef _TEST_PATH_HPP_
#define _TEST_PATH_HPP_

#include "../domain/domain.hpp"
#include "../io/io_gnuplot_domain.h"
#include "../carpio_define.hpp"

namespace carpio {
void test_path_1() {
	Path_<3> path(3);
	std::cout << "size  :  " << path.size() << std::endl;
	path.show();

	Path_<3> path1(7);
	path1.set();
	path1.show();
	std::cout << "is valid " << path1.is_valid() << std::endl;
	path1.copy(path, 0, 3);
	path1.show();

	Path_<3> pathn = path1;
	path1.show();
	path1.append(path);
	path1.show();

	std::cout << " --------------------- \n";
	path.set();
	path.show();
	path1.show();
	std::cout<< (path < path) <<"\n";

}

inline void test_path_2() {
	Grid_<Float, Float, 2> g(2, 0, 1, //
			2, 0, 1);
	g.connect_root();
	std::cout << "num root : " << g.get_num_root() << std::endl;
	Adaptive_<Float, Float, 2> adp(&g, 3, 9);
	adp.adapt_full();

	Grid_<Float, Float, 2>::pNode p = g.get_pnode(0.8, 0.35);
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_LeafNodes(ga, g);
	lga.push_back(ga);
	GnuplotActor_Node(ga, (*p));
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	gp.set_xrange(0, 2);
	gp.set_yrange(0, 2);
	GnuplotShow(gp, lga);
	p->show("p");

}

}

#endif
