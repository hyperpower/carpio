#ifndef _TEST_DEFINE_H_
#define _TEST_DEFINE_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/poisson.hpp"
#include "gtest/gtest.h"
#include <math.h>
using namespace std;
namespace carpio {

typedef Domain_<Float, Float, 2> Domain;

inline Float error_1(Domain& d, st ires, st iexact) {
	Float norm1 = 0;
	Float svol = 0;
	for (typename Domain::Grid::iterator_leaf iterf = d.grid().begin_leaf();
			iterf != d.grid().end_leaf(); ++iterf) {
		Float res = iterf->cdva(ires);
		Float exa = iterf->cdva(iexact);
		Float vol = iterf->volume();
		Float err = res - exa;
		norm1 += (Abs(err) * vol);
		svol += vol;
	}
	return norm1 / svol;
}

inline Float error_2(Domain& d, st ires, st iexact) {
	Float norm2 = 0;
	Float svol = 0;
	for (typename Domain::Grid::iterator_leaf iterf = d.grid().begin_leaf();
			iterf != d.grid().end_leaf(); ++iterf) {
		Float res = iterf->cdva(ires);
		Float exa = iterf->cdva(iexact);
		Float vol = iterf->volume();
		Float err = res - exa;
		norm2 += (err * err * vol);
		svol += vol;
	}
	return sqrt(norm2) / svol;
}

inline Float error_i(Domain& d, st ires, st iexact) {
	Float normi = 0;
	for (typename Domain::Grid::iterator_leaf iterf = d.grid().begin_leaf();
			iterf != d.grid().end_leaf(); ++iterf) {
		Float res = iterf->cdva(ires);
		Float exa = iterf->cdva(iexact);
		Float err = res - exa;
		if (iterf == d.grid().begin_leaf()) {
			normi = Abs(err);
		} else {
			if (normi < Abs(err)) {
				normi = Abs(err);
			}
		}
	}
	return normi;
}

inline Float cal_order(Float ec, Float ef) {
	return log(Abs(ec) / Abs(ef)) / log(2);
}

inline void show_splot_surface(Domain& d, st idx) {
	std::list<Gnuplot_actor> slga;
	Gnuplot_actor sga;
	GnuplotActor_LeafNodesSurface(sga, d.grid(), 2);
	//GnuplotActor_NodesSurface(sga, pn, 0);
	slga.push_back(sga);
	GnuplotActor_GhostNodesSurface(sga, d.ghost(), 2);
	slga.push_back(sga);
	//sga.show_data();
	Gnuplot sgp;
	//sgp.set_equal_ratio();
	sgp.set_view(45, 10, 1, 1);
	sgp.set_palette_blue_red();
	sgp.set("ticslevel 0");
	//sgp.set_xrange(1.4, 2.0);
	//sgp.set_yrange(1.4, 2.0);
	sgp.splot(slga);
}

}

#endif
