#ifndef _TEST_GRID_2D_HPP_
#define _TEST_GRID_2D_HPP_

#include "../domain/domain.hpp"
#include "../carpio_define.hpp"
#include "../io/io_gnuplot_domain.h"
#include "../domain/stencil.hpp"
#include "../calculation/scalar.h"
#include "../calculation/interpolate.h"

namespace carpio {

inline void test_grid_2d() {
	Grid_<Float, Float, 2> g(2, 0, 1, //
			2, 0, 1);
	g.connect_root();
	std::cout << "num root : " << g.get_num_root() << std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt_full();

	//GnuplotShow_RootNodes(g);
	// test iterator
	int i = 0;
	Grid_2D::iterator_leaf iter = g.begin_leaf();
	for (; iter != g.end_leaf(); ++iter) {
		i++;
	}
	std::cout << "all     : " << g(0, 0)->count_all() << std::endl;
	std::cout << "leaf    : " << g(0, 0)->count_leaf() << std::endl;
	std::cout << "level 1 : " << g(0, 0)->count_level(1) << std::endl;
	std::cout << "level 2 : " << g(0, 0)->count_level(2) << std::endl;
}

inline void test_stencil_1d() {
	Grid_<Float, Float, 2> g(2, 0, 1, //
			2, 0, 1);
	g.connect_root();
	std::cout << "num root : " << g.get_num_root() << std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt_full();
	//stencil
	Grid_<Float, Float, 2>::iterator_leaf il = g.begin_leaf();
	++il;
	++il;
	++il;
	Stencil_2D1 st(il.get_pointer(), _Y_, 2, 1);
	st.show();
	std::cout << "Stencil \n";

	//GnuplotShow_LeafNodes(g);
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga_leaf, ga_stencil;
	GnuplotActor_LeafNodes(ga_leaf, g);
	GnuplotActor_Stencil(ga_stencil, st);
	lga.push_back(ga_leaf);
	lga.push_back(ga_stencil);
	GnuplotShow(lga);
	// test iterator
}
Vt _set(Cvt x, Cvt y, Cvt z) {
	return (sin(x) - 0.5) * (x - 0.5) + y * y;
}
inline void test_interpolate_1d() {
	Grid_<Float, Float, 2> g(2, 0, 1, //
			2, 0, 1);
	g.connect_root();
	std::cout << "num root : " << g.get_num_root() << std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt_full();
	// new value on leaf
	NewCenterDataOnLeaf(g, 2);
	// set value on leaf
	SetScalarOnCenterLeaf(g, 0, _set);
	//
	//stencil
	Grid_<Float, Float, 2>::iterator_leaf il = g.begin_leaf();
	++il;
	++il;
	++il;
	st idx = 0;
	Float res;
	il->show();
	InterpolateOnFace_1Order(il.get_pointer(), _XM_, idx, res);
	std::cout << "pnode c  : " << il->cd(0) << std::endl;
	std::cout << "res      : " << res << std::endl;
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga_leaf, ga_stencil;
	GnuplotActor_LeafNodesContours(ga_leaf, g, 0);
	lga.push_back(ga_leaf);
	GnuplotShow(lga);
}

void test_adapt_solid() {
	Grid_<Float, Float, 2> g(2, 0.0, 1, //
			2, 0.0, 1);
	g.connect_root();
	std::cout << "num root : " << g.get_num_root() << std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 7);
	//adp.adapt_full();
	// shape
	Shape2D cir;
	CreatCircle(cir, 0.0, 0.0, 1.0, 359);
	adp.adapt_solid(cir);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_LeafNodes(ga, g);
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, cir);
	lga.push_back(ga);
	GnuplotShow(lga);
}

void test_stencil_2D() { //jan 13, 2016
	std::cout << "test stencil 2D   ================\n";
	Grid_<Float, Float, 2> g(2, 0.0, 1, //
			2, 0.0, 1);
	g.connect_root();
	std::cout << "num root : " << g.get_num_root() << std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	//adp.adapt_full();
	// shape
	Shape2D cir;
	CreatCircle(cir, 0.0, 0.0, 1.0, 359);
	adp.adapt_solid(cir);
	//stencil ==============================
	Grid_<Float, Float, 2>::iterator_leaf il = g.begin_leaf();
	++il;
	++il;
	++il;
	++il;
	Stencil_2D1 st(il.get_pointer(), _Y_, 2, 1);
	st.show();
	std::cout << "Stencil \n";

	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_LeafNodes(ga, g);
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, cir);
	lga.push_back(ga);
	GnuplotActor_Stencil(ga, st);
	lga.push_back(ga);
	GnuplotShow(lga);
}

}

#endif
