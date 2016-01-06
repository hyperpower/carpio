#ifndef IO_GNUPLOT_DOMAIN_H_
#define IO_GNUPLOT_DOMAIN_H_

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "../carpio_define.hpp"
#include "../domain/domain.hpp"
#include "gnuplot.h"
#include <string>
#include <fstream>

namespace carpio {

int GnuplotActor_Cell(Gnuplot_actor& actor, const Cell_2D& c);
int GnuplotActor_Node(Gnuplot_actor& actor, const Node_2D& node);
int GnuplotActor_RootNodes(Gnuplot_actor& actor, const Grid_2D& g);
int GnuplotActor_LeafNodes(Gnuplot_actor& actor, const Grid_2D& g);
int GnuplotActor_Nodes(Gnuplot_actor& actor, const std::list<pNode_2D>& lpn);
int GnuplotActor_Stencil(Gnuplot_actor& actor, const Stencil_2D1& s);
int GnuplotActor_LeafNodesContours(Gnuplot_actor& actor, const Grid_2D& g,
		st idx);

int GnuplotShow_RootNodes(const Grid_2D& grid);
int GnuplotShow_LeafNodes(const Grid_2D& grid);
int GnuplotShow(const std::list<Gnuplot_actor> lga);
}

#endif /* IO_H_ */