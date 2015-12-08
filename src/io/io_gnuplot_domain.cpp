#include "io_gnuplot_domain.h"
#include "io_define.hpp"
#include "../carpio_define.hpp"
#include "../domain/domain.hpp"
#include "gnuplot.h"

namespace carpio {

int GnuplotActorDataPushBack(std::list<std::string>& ldata, const Cell_2D& c) {
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back("");
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back("");
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back("");
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back("");
	return _SUCCESS;
}

int GnuplotActor_Cell(Gnuplot_actor& actor, const Cell_2D& c) {
	actor.command() = "using 1:2 ";
	GnuplotActorDataPushBack(actor.data(), c);
	return _SUCCESS;
}

int GnuplotActor_Node(Gnuplot_actor& actor, const Node_2D& node) {
	actor.command() = "using 1:2 ";
	GnuplotActorDataPushBack(actor.data(), *(node.cell));
	return _SUCCESS;
}

int GnuplotActor_Nodes(Gnuplot_actor& actor, const std::list<pNode_2D>& lpn) {
	actor.command() = "using 1:2 ";
	for (std::list<pNode_2D>::const_iterator iter = lpn.begin();
			iter != lpn.end(); ++iter) {
		pNode_2D pn = (*iter);
		GnuplotActorDataPushBack(actor.data(), *(pn->cell));
	}
	return _SUCCESS;
}

int GnuplotActor_RootNodes(Gnuplot_actor& actor, const Grid_2D& g) {
	actor.command() = "using 1:2 title \"\" ";
	for (Grid_2D::const_iterator iter = g.begin(); iter != g.end(); ++iter) {
		Grid_2D::pNode proot = (*iter);
		if (proot != nullptr) {
			GnuplotActorDataPushBack(actor.data(), *(proot->cell));
		}
	}
	return _SUCCESS;
}

int GnuplotActor_LeafNodes(Gnuplot_actor& actor, const Grid_2D& g) {
	actor.command() = "using 1:2 title \"\" ";
	for (Grid_2D::const_iterator_leaf iter = g.begin_leaf();
			iter != g.end_leaf(); ++iter) {
		GnuplotActorDataPushBack(actor.data(), *(iter->cell));
	}
	return _SUCCESS;
}

int GnuplotActor_Stencil(Gnuplot_actor& actor, const Stencil_2D1& s) {
	actor.command() = "using 1:2 title \"\" ";
	for (st i = 0; i < s.size(); ++i) {
		Stencil_2D1::const_pNode pn = s.at_1d(i);
		if (pn != nullptr) {
			GnuplotActorDataPushBack(actor.data(), *(pn->cell));
		}
	}
	return _SUCCESS;
}

int GnuplotShow_RootNodes(const Grid_2D& grid) {
	// creat an actor
	Gnuplot_actor actor;
	GnuplotActor_RootNodes(actor, grid);
	//
	Gnuplot gp;
	gp.set_equal_ratio();
	std::ostringstream ss;
	ss << "plot \"-\" " << actor.command() << "with lines lw 1";
	gp.cmd(ss.str() + "\n");
	ss.str("");
	gp.output_inline_data(actor);

	return _SUCCESS;
}

int GnuplotShow_LeafNodes(const Grid_2D& grid) {
	// creat an actor
	Gnuplot_actor actor;
	GnuplotActor_LeafNodes(actor, grid);
	//
	Gnuplot gp;
	gp.set_equal_ratio();
	std::ostringstream ss;
	ss << "plot \"-\" " << actor.command() << "with lines lw 1";
	gp.cmd(ss.str() + "\n");
	ss.str("");
	gp.output_inline_data(actor);

	return _SUCCESS;
}

int GnuplotShow(const std::list<Gnuplot_actor> lga) {
	//
	Gnuplot gp;
	gp.set_equal_ratio();
	std::ostringstream ss;
	ss << "plot ";
	for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
			iter != lga.end(); ++iter) {
		ss << "\"-\" " << iter->command() << "with lines lw 1";
		if (lga.size() >= 2 && (iter != (--lga.end()))) {
			ss << ",\\\n";
		}
	}
	gp.cmd(ss.str() + "\n");
	ss.str("");
	for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
			iter != lga.end(); ++iter) {
		gp.output_inline_data((*iter));
	}
	return _SUCCESS;
}

}
