#ifndef IO_GNUPLOT_DOMAIN_H_
#define IO_GNUPLOT_DOMAIN_H_

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "../carpio_define.hpp"
#include "../domain/domain.hpp"
#include "../calculation/vof.h"
#include "../geometry/geometry.hpp"
#include "../calculation/poisson.hpp"
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
int GnuplotActor_StencilContour(Gnuplot_actor& actor, const Stencil_2D2& s,
		st idx);
int GnuplotActor_LeafNodesContours(Gnuplot_actor& actor, const Grid_2D& g,
		st idx);
int GnuplotActor_NodesContours(Gnuplot_actor& actor,
		const std::list<pNode_2D>& pn, st idx);
int GnuplotActor_NodesValues(Gnuplot_actor&, std::list<pNode_2D>&, st, Axes);
int GnuplotActor_LeafNodesSurface(Gnuplot_actor& actor, Grid_2D& g, st idx);
int GnuplotActor_GhostNodesSurface(Gnuplot_actor& actor, Ghost_2D& g, st idx);
int GnuplotActor_NodesSurface(Gnuplot_actor& actor, pNode_2D pn, st idx);
int GnuplotActor_LeafNodesDataIndex(Gnuplot_actor& actor, const Grid_2D& g);
int GnuplotActor_GhostNodes(Gnuplot_actor& actor, const Ghost_2D& g);
int GnuplotActor_GhostNodesContours(Gnuplot_actor& actor, const Ghost_2D& g,
		st vi);
int GnuplotActor_GhostNodesContour_BoundaryIndex(Gnuplot_actor& actor,
		const Ghost_2D& g);
int GnuplotActor_GhostNodesDataIndex(Gnuplot_actor& actor, const Ghost_2D& g);
int GnuplotActor_Shape2D(Gnuplot_actor& actor, const Shape2D& g);
int GnuplotActor_Shape2D(Gnuplot_actor& actor, const Shape2D& g, st base_idx);
int GnuplotActor_Vof2D(Gnuplot_actor& actor, const Vof_<Float, Float, 2>& vof);
// this part is create gnuplot actor for calculation
//int GnuplotActor_Expression(Gnuplot_actor& actor, const Poisson_<Float, Float, 2>::Exp& exp);
int GnuplotActor_Expression(Gnuplot_actor& actor,
		const Expression_<Float, Float, 2>& exp);

// this part is creat gnuplot actor for geometry
int GnuplotActor_Segment2D(Gnuplot_actor& actor, const Segment_2D& g);
int GnuplotActor_Polygon(Gnuplot_actor& actor, const Polygon& g);
int GnuplotActor_Polygon_vector(Gnuplot_actor& actor, const Polygon& g);

int GnuplotActor_MatrixSCR(Gnuplot_actor&, const MatrixSCR_<Float>&);
int GnuplotActor_ArrayList(Gnuplot_actor&, const ArrayListV<Float>&);

int GnuplotShow_RootNodes(const Grid_2D& grid);
int GnuplotShow_LeafNodes(const Grid_2D& grid);

template<class VALUE>
void gnuplot_show_ylog(const std::list<VALUE>& list) {
	if (list.size() == 0) {
		return;
	}
	Gnuplot gp;
	gp.set_ylogscale(10);
	ArrayListV<VALUE> arr(list.size());
	int i = 0;
	for (typename std::list<VALUE>::const_iterator it = list.begin();
			it != list.end(); ++it) {
		arr[i] = (*it);
		i++;
	}
	gp.plot_1(arr, "");
}


}

#endif /* IO_H_ */
