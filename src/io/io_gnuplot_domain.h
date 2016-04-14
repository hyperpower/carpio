#ifndef IO_GNUPLOT_DOMAIN_H_
#define IO_GNUPLOT_DOMAIN_H_

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "../carpio_define.hpp"
#include "../domain/domain.hpp"
#include "../calculation/vof.h"
#include "../geometry/geometry.hpp"
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
int GnuplotActor_Stencil(Gnuplot_actor& actor, const Stencil_2D2& s);
int GnuplotActor_StencilContour(Gnuplot_actor& actor, const Stencil_2D2& s, st idx);
int GnuplotActor_LeafNodesContours(Gnuplot_actor& actor, const Grid_2D& g,
		st idx);
int GnuplotActor_GhostNodes(Gnuplot_actor& actor, const Ghost_2D& g);
int GnuplotActor_GhostNodesContour(Gnuplot_actor& actor, const Ghost_2D& g);
int GnuplotActor_Shape2D(Gnuplot_actor& actor, const Shape2D& g);
int GnuplotActor_Shape2D(Gnuplot_actor& actor, const Shape2D& g, st base_idx);
int GnuplotActor_Vof2D(Gnuplot_actor& actor, const Vof_<Float, Float, 2>& vof);

// this part is creat gnuplot actor for geometry
int GnuplotActor_Segment2D(Gnuplot_actor& actor, const Segment_2D& g);
int GnuplotActor_Polygon(Gnuplot_actor& actor, const Polygon& g);
int GnuplotActor_Polygon_vector(Gnuplot_actor& actor, const Polygon& g);


int GnuplotShow_RootNodes(const Grid_2D& grid);
int GnuplotShow_LeafNodes(const Grid_2D& grid);

}

#endif /* IO_H_ */
