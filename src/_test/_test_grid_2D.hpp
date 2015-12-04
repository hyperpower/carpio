#ifndef _TEST_GRID_2D_HPP_
#define _TEST_GRID_2D_HPP_

#include "../domain/domain.hpp"
#include "../carpio_define.hpp"
#include "../io/io_gnuplot_domain.h"
#include "../domain/stencil.hpp"

namespace carpio {

inline void test_grid_2d() {
	Grid_<Float, Float, 2> g(
			2, 0, 1, //
			2, 0, 1 );
	g.connect_root();
	std::cout<< "num root : "<<g.get_num_root()<<std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt();

	//GnuplotShow_RootNodes(g);
	// test iterator
	int i=0;
	Grid_2D::iterator_leaf iter=g.begin_leaf();
	for(; iter!=g.end_leaf(); ++iter){
		i++;
	}
	std::cout<< "all     : "<<g(0,0)->count_all()<< std::endl;
	std::cout<< "leaf    : "<<g(0,0)->count_leaf()<< std::endl;
	std::cout<< "level 1 : "<<g(0,0)->count_level(1)<< std::endl;
	std::cout<< "level 2 : "<<g(0,0)->count_level(2)<< std::endl;
}

inline void test_stencil_1d() {
	Grid_<Float, Float, 2> g(
			2, 0, 1, //
			2, 0, 1 );
	g.connect_root();
	std::cout<< "num root : "<<g.get_num_root()<<std::endl;
	Adaptive<Float, Float, 2> adp(&g, 2, 5);
	adp.adapt();
	//stencil
	Grid_<Float, Float, 2>::iterator_leaf il = g.begin_leaf();
	Stencil_2D1 st(il.get_pointer(),_X_, 1, 1);
	st.show();
	std::cout<<"Stencil \n";

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

}

#endif
