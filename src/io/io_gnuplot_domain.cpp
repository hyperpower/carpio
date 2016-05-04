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

int GnuplotActorDataPushBack(std::list<std::string>& ldata,
		const Segment_2D& c) {
	ldata.push_back(ToString(c.psx(), c.psy(), " "));
	ldata.push_back(ToString(c.pex(), c.pey(), " "));
	ldata.push_back("");
	return _SUCCESS;
}

int GnuplotActorDataPushBack_Contour(std::list<std::string>& ldata,
		const_pNode_2D& pn, st idx) {
	ldata.push_back(
			ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
					pn->p(_M_, _Y_), pn->p(_P_, _Y_), pn->cda(idx), " ")); //point" "));
	return _SUCCESS;
}

int GnuplotActorDataPushBack_Contour(std::list<std::string>& ldata,
		const_pNode_2D& pn, Float vt) {
	ldata.push_back(
			ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
					pn->p(_M_, _Y_), pn->p(_P_, _Y_), vt, " ")); //point" "));
	return _SUCCESS;
}

int GnuplotActorDataIndexPushBack_Contour(std::list<std::string>& ldata,
		const_pNode_2D& pn) {
	ldata.push_back(
			ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
					pn->p(_M_, _Y_), pn->p(_P_, _Y_), pn->d_idx(), " ")); //point" "));
	return _SUCCESS;
}

int GnuplotActor_Cell(Gnuplot_actor& actor, const Cell_2D& c) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	GnuplotActorDataPushBack(actor.data(), c);
	return _SUCCESS;
}

int GnuplotActor_Node(Gnuplot_actor& actor, const Node_2D& node) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	GnuplotActorDataPushBack(actor.data(), *(node.cell));
	return _SUCCESS;
}

int GnuplotActor_Nodes(Gnuplot_actor& actor, const std::list<pNode_2D>& lpn) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	for (std::list<pNode_2D>::const_iterator iter = lpn.begin();
			iter != lpn.end(); ++iter) {
		pNode_2D pn = (*iter);
		GnuplotActorDataPushBack(actor.data(), *(pn->cell));
	}
	return _SUCCESS;
}

int GnuplotActor_RootNodes(Gnuplot_actor& actor, const Grid_2D& g) {
	actor.clear();
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
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	for (Grid_2D::const_iterator_leaf iter = g.begin_leaf();
			iter != g.end_leaf(); ++iter) {
		GnuplotActorDataPushBack(actor.data(), *(iter->cell));
	}
	return _SUCCESS;
}
int GnuplotActor_GhostNodes(Gnuplot_actor& actor, const Ghost_2D& g) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	for (typename Ghost_2D::const_iterator iter = g.begin(); iter != g.end();
			++iter) {
		GnuplotActorDataPushBack(actor.data(), *(iter->second.pghost->cell));
	}
	return _SUCCESS;
}

int GnuplotActor_GhostNodesContour_BoundaryIndex(Gnuplot_actor& actor,
		const Ghost_2D& g) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	typedef typename Ghost_2D::GhostNode Node;
	std::function<void(const Node&)> fun = [&](const Node& n) {
		//
			typename Ghost_2D::GhostVal gval=n.second;
			typename Ghost_2D::pNode pn = gval.pghost;
			typename Ghost_2D::GhostID gid = n.first;
			//assume segments for each shape less than 10000
			Float v = gval.shape_idx * 10 + gval.seg_idx;
			//std::cout<< v << " " <<gval.seg_idx << " " << gval.shape_idx<<"\n";
			if (pn != nullptr) {
				actor.data().push_back(
						ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
								pn->p(_M_, _Y_), pn->p(_P_, _Y_), v, " "));
			}
		};
	g.for_each_ghost_node(fun);

	return _SUCCESS;
}

int GnuplotActor_GhostNodesContours(Gnuplot_actor& actor, const Ghost_2D& g,
		st vi) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	typedef typename Ghost_2D::GhostNode Node;
	std::function<void(const Node&)> fun = [&](const Node& n) {
		//
			typename Ghost_2D::GhostVal gval=n.second;
			typename Ghost_2D::pNode pn = gval.pghost;
			typename Ghost_2D::GhostID gid = n.first;
			//assume segments for each shape less than 10000
			Float v = pn->cda(vi);
			//std::cout<< v << " " <<gval.seg_idx << " " << gval.shape_idx<<"\n";
			if (pn != nullptr) {
				actor.data().push_back(
						ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
								pn->p(_M_, _Y_), pn->p(_P_, _Y_), v, " "));
			}
		};
	g.for_each_ghost_node(fun);

	return _SUCCESS;
}

int GnuplotActor_GhostNodesDataIndex(Gnuplot_actor& actor, const Ghost_2D& g) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	typedef typename Ghost_2D::GhostNode Node;
	std::function<void(const Node&)> fun = [&](const Node& n) {
		//
			typename Ghost_2D::GhostVal gval=n.second;
			typename Ghost_2D::pNode pn = gval.pghost;
			typename Ghost_2D::GhostID gid = n.first;
			//assume segments for each shape less than 10000
			Float v = pn->d_idx();
			//std::cout<< v << " " <<gval.seg_idx << " " << gval.shape_idx<<"\n";
			if (pn != nullptr) {
				actor.data().push_back(
						ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
								pn->p(_M_, _Y_), pn->p(_P_, _Y_), v, " "));
			}
		};
	g.for_each_ghost_node(fun);

	return _SUCCESS;
}
int GnuplotActor_StencilContour(Gnuplot_actor& actor, const Stencil_2D2& s,
		st idx) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	for (st i = 0; i < s.size(); ++i) {
		Stencil_2D1::const_pNode pn = s.at_1d(i);
		if (pn != nullptr) {
			GnuplotActorDataPushBack_Contour(actor.data(), pn, idx);
		}
	}
	return _SUCCESS;
}

int GnuplotActor_LeafNodesContours(Gnuplot_actor& actor, const Grid_2D& g,
		st idx) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	for (Grid_2D::const_iterator_leaf iter = g.begin_leaf();
			iter != g.end_leaf(); ++iter) {
		const_pNode_2D pn = iter.get_pointer();
		GnuplotActorDataPushBack_Contour(actor.data(), pn, idx);
	}
	return _SUCCESS;
}

int GnuplotActor_NodesContours(Gnuplot_actor& actor,
		const std::list<pNode_2D>& lp, st idx) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	for (std::list<pNode_2D>::const_iterator iter = lp.begin();
			iter != lp.end(); ++iter) {
		const_pNode_2D pn = *iter;
		GnuplotActorDataPushBack_Contour(actor.data(), pn, idx);
	}
	return _SUCCESS;
}
int GnuplotActor_NodesValues(Gnuplot_actor& actor, std::list<pNode_2D>& lp,
		st idx, Axes aix) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	actor.style() = "with linespoints";
	// sort
	// define compare function
	struct compare_pNode_x {
		bool operator()(pNode_2D first, pNode_2D second) {
			return first->cp(_X_) < second->cp(_X_);
		}
	};
	struct compare_pNode_y {
		bool operator()(pNode_2D first, pNode_2D second) {
			return first->cp(_Y_) < second->cp(_Y_);
		}
	};

	if (aix == _X_) {
		lp.sort(compare_pNode_x());
	} else {
		lp.sort(compare_pNode_y());
	}
	//
	for (std::list<pNode_2D>::const_iterator iter = lp.begin();
			iter != lp.end(); ++iter) {
		const_pNode_2D pn = *iter;
		actor.data().push_back(ToString(pn->cp(aix), pn->cda(idx), " ")); //point" "));
	}
	return _SUCCESS;
}

// ! this actor need splot
int GnuplotActor_NodesSurface(Gnuplot_actor& actor, pNode_2D pn, st idx) {
	actor.clear();
	actor.command() = "using 1:2:3 title \"\" ";
	actor.style() = "with lines lw 1 lc palette";
	std::function<void(Orientation, Orientation)> _fun =
			[&actor, pn, idx](Orientation o1, Orientation o2) {
				typedef Interpolate_<Grid_2D::cvt, Grid_2D::vt, Grid_2D::Dim> Inter;
				Direction dir = ToCornerDirection(o1,_X_,o2,_Y_);
				Inter::Exp exp = Inter::OnCorner_1Order(pn,dir);
				Inter::vt x = pn->p(dir,_X_);
				Inter::vt y = pn->p(dir,_Y_);
				Inter::vt v = exp.subsitute(idx);
				actor.data().push_back(ToString(x, y, v, " "));
				std::cout<< ToString(x, y, v, " ")<<"\n";
			};
	_fun(_M_, _M_);
	_fun(_P_, _M_);
	actor.data().push_back("");
	_fun(_P_, _M_);
	_fun(_P_, _P_);
	actor.data().push_back("");
	_fun(_P_, _P_);
	_fun(_M_, _P_);
	actor.data().push_back("");
	_fun(_M_, _P_);
	_fun(_M_, _M_);
	actor.data().push_back("");
	return _SUCCESS;
}
// ! this actor need splot
int GnuplotActor_LeafNodesSurface(Gnuplot_actor& actor, Grid_2D& g, st idx) {
	actor.clear();
	actor.command() = "using 1:2:3 title \"\" ";
	actor.style() = "with lines lw 1 lc palette";
	for (Grid_2D::iterator_leaf iter = g.begin_leaf(); iter != g.end_leaf();
			++iter) {
		pNode_2D pn = iter.get_pointer();
		std::function<void(Orientation, Orientation)> _fun =
				[&actor, pn, idx](Orientation o1, Orientation o2) {
					typedef Interpolate_<Grid_2D::cvt, Grid_2D::vt, Grid_2D::Dim> Inter;
					Direction dir = ToCornerDirection(o1,_X_,o2,_Y_);
					Inter::Exp exp = Inter::OnCorner_1Order(pn,dir);
					Inter::vt x = pn->p(dir,_X_);
					Inter::vt y = pn->p(dir,_Y_);
					Inter::vt v = exp.subsitute(idx);
					actor.data().push_back(ToString(x, y, v, " "));
				};
		_fun(_M_, _M_);
		_fun(_P_, _M_);
		_fun(_P_, _P_);
		_fun(_M_, _P_);
		_fun(_M_, _M_);
		actor.data().push_back("");
		actor.data().push_back(""); //may have problem
	}
	return _SUCCESS;
}

int GnuplotActor_GhostNodesSurface(Gnuplot_actor& actor, Ghost_2D& g, st idx) {
	actor.clear();
	actor.command() = "using 1:2:3 title \"\" ";
	actor.style() = "with lines lw 1 lc palette";
	for (typename Ghost_2D::iterator iter = g.begin(); iter != g.end();
			++iter) {
		pNode_2D pn = (*iter).second.pghost;
		std::function<void(Orientation, Orientation)> _fun =
				[&actor, pn, idx](Orientation o1, Orientation o2) {
					typedef Interpolate_<Grid_2D::cvt, Grid_2D::vt, Grid_2D::Dim> Inter;
					Direction dir = ToCornerDirection(o1,_X_,o2,_Y_);
					Inter::Exp exp = Inter::OnCorner_1Order(pn,dir);
					Inter::vt x = pn->p(dir,_X_);
					Inter::vt y = pn->p(dir,_Y_);
					Inter::vt v = exp.subsitute(idx);
					actor.data().push_back(ToString(x, y, v, " "));
				};
		_fun(_M_, _M_);
		_fun(_P_, _M_);
		_fun(_P_, _P_);
		_fun(_M_, _P_);
		_fun(_M_, _M_);
		actor.data().push_back("");
		actor.data().push_back(""); //may have problem
	}
	return _SUCCESS;
}
int GnuplotActor_LeafNodesDataIndex(Gnuplot_actor& actor, const Grid_2D& g) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	for (Grid_2D::const_iterator_leaf iter = g.begin_leaf();
			iter != g.end_leaf(); ++iter) {
		const_pNode_2D pn = iter.get_pointer();
		GnuplotActorDataIndexPushBack_Contour(actor.data(), pn);
	}
	return _SUCCESS;
}

/*
 * shape
 */
int GnuplotActor_Shape2D(Gnuplot_actor& actor, const Shape2D& g) {
	actor.clear();
	actor.command() = "using 1:2 title \"\"";
	if (g.empty()) {
		actor.data().push_back("");
		return _ERROR;
	}
	typedef typename Shape2D::S2D::Point Poi;
	for (st i = 0; i < g.size_vertexs(); ++i) {
		const Poi& p = g.v(i);
		actor.data().push_back(ToString(p.x(), p.y(), " "));
	}
	const Poi& pstart = g.v(0);
	actor.data().push_back(ToString(pstart.x(), pstart.y(), " "));
	actor.data().push_back("");
	return _SUCCESS;
}
int GnuplotActor_Shape2D(Gnuplot_actor& actor, const Shape2D& g, st base_idx) {
	actor.clear();
	actor.command() = "using 1:2:3 title \"\" ";
	actor.style() = "with lines lc variable";
	if (g.empty()) {
		actor.data().push_back("");
		return _ERROR;
	}
	typedef typename Shape2D::S2D::Point Poi;
	for (st i = 0; i < g.size_vertexs(); ++i) {
		const Poi& p = g.v(i);
		actor.data().push_back(ToString(p.x(), p.y(), i + base_idx, " "));
	}
	const Poi& pstart = g.v(0);
	actor.data().push_back(ToString(pstart.x(), pstart.y(), base_idx, " "));
	actor.data().push_back("");
	return _SUCCESS;
}

int GnuplotActor_Stencil(Gnuplot_actor& actor, const Stencil_2D1& s) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	for (st i = 0; i < s.size(); ++i) {
		Stencil_2D1::const_pNode pn = s.at_1d(i);
		if (pn != nullptr) {
			GnuplotActorDataPushBack(actor.data(), *(pn->cell));
		}
	}
	return _SUCCESS;
}
int GnuplotActor_Stencil(Gnuplot_actor& actor, const Stencil_2D2& s) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	for (st i = 0; i < s.size(); ++i) {
		Stencil_2D2::const_pNode pn = s.at_1d(i);
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
	// Create an actor
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

int GnuplotActor_Vof2D(Gnuplot_actor& actor, const Vof_<Float, Float, 2>& vof) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	typename Vof_<Float, Float, 2>::const_pGrid pg = vof.get_pgrid();
	for (Grid_2D::const_iterator_leaf iter = pg->begin_leaf();
			iter != pg->end_leaf(); ++iter) {
		const_pNode_2D pn = iter.get_pointer();
		Vof_<Float, Float, 2>::const_pVoff pvf = vof.get_pface(pn);
		if (pvf->has_seg()) {
			GnuplotActorDataPushBack(actor.data(), (*pvf->get_pSeg()));
		}
	}
	return _SUCCESS;
}

int GnuplotActor_Expression(Gnuplot_actor& actor,
		const Expression_<Float, Float, 2>& exp) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	for (auto iter = exp.begin(); iter != exp.end(); ++iter) {
		const_pNode_2D pn = iter->first.first;
		Float val = iter->second;
		GnuplotActorDataPushBack_Contour(actor.data(), pn, val);
	}
	return _SUCCESS;
}
int GnuplotActor_Segment2D(Gnuplot_actor& actor, const Segment_2D& g) {
	actor.clear();
	actor.command() = "using 1:2 title \"\"";
	if (g.empty()) {
		actor.data().push_back("");
		return _ERROR;
	}
	typedef typename Segment_2D::Point Poi;
	const Poi& ps = g.ps();
	actor.data().push_back(ToString(ps.x(), ps.y(), " "));
	const Poi& pe = g.pe();
	actor.data().push_back(ToString(pe.x(), pe.y(), " "));
	actor.data().push_back("");
	return _SUCCESS;
}

int GnuplotActor_MatrixSCR(Gnuplot_actor& actor, const MatrixSCR_<Float>& m) {
	actor.clear();
	actor.command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor.style() = "with boxxy fs solid palette";
	st k = 0;
	for (st i = 1; i <= m.iLen(); i++) {
		for (st j = k; j < m.row_ptr(i); j++) {
			Float xc = i - 1 + 0.5;
			Float yc = m.col_ind(k) + 0.5;
			Float xm = i - 1;
			Float xp = i;
			Float ym = m.col_ind(k);
			Float yp = m.col_ind(k) + 1;
			Float val = m.val(k);
			actor.data().push_back(ToString(xc, yc, xm, xp, ym, yp, val, " "));
			k++;
		}
	}
	return _SUCCESS;
}

int GnuplotActor_ArrayList(Gnuplot_actor& actor, const ArrayListV<Float>& arr) {
	actor.clear();
	actor.command() = "using 1:2 title \"\" ";
	for (st i = 0; i < arr.size(); i++) {
		actor.data().push_back(ToString(i, arr[i], " "));
	}
	return _SUCCESS;
}

int GnuplotActor_Polygon(Gnuplot_actor& actor, const Polygon& g) {
	actor.clear();
	actor.command() = "using 1:2 title \"\"";
	if (g.empty()) {
		actor.data().push_back("");
		return _ERROR;
	}
	typedef typename Polygon::Point Poi;
	for (st i = 0; i < g.size_vertexs(); ++i) {
		const Poi& p = g.v(i);
		actor.data().push_back(ToString(p.x(), p.y(), " "));
	}
	const Poi& pstart = g.v(0);
	actor.data().push_back(ToString(pstart.x(), pstart.y(), " "));
	actor.data().push_back("");
	return _SUCCESS;
}
int GnuplotActor_Polygon_vector(Gnuplot_actor& actor, const Polygon& g) {
	actor.clear();
	actor.command() = "using 1:2:3:4 title \"\"";
	actor.style() = "with vectors";
	if (g.empty()) {
		actor.data().push_back("");
		return _ERROR;
	}
	typedef typename Polygon::Segment Segment;
	for (st i = 0; i < g.size_segments(); ++i) {
		Segment p = g.get_segment(i);
		std::stringstream sstr;
		sstr << p.psx() << " " << p.psy() << " " << p.dx() << " " << p.dy();
		actor.data().push_back(sstr.str());
		actor.data().push_back("");
	}
	return _SUCCESS;
}

}
