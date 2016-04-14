#ifndef SHAPE_H_
#define SHAPE_H_

#include "../carpio_define.hpp"
#include "domain_define.hpp"

#include "../geometry/geometry.hpp"

#include "../utility/clipper.hpp"
#include <functional>
#include <math.h>

namespace carpio {
/*
 *  if the shape is 2D, the shape is polygon
 *  if the shape is 3D, the shape is surface
 */

template<typename VALUE, st DIM>
class Shape_ {
public:
	static const st Dim = DIM;

	typedef VALUE vt;

	typedef Shape_<VALUE, DIM> Self;

	typedef Polygon_<VALUE> S2D;
	typedef Polygon_<VALUE>* pS2D;
	typedef Polygon_<VALUE>& ref_S2D;
	typedef const Polygon_<VALUE>& const_ref_S2D;
	typedef typename S2D::Segment Seg2D;
	typedef typename S2D::Point Poi2D;
protected:
	pS2D _ps2d;
	int _type;
public:
	/*
	 * constructor
	 */
	Shape_() {
		if (Dim == 2) {
			_ps2d = nullptr;
		}
		_type = 0;
	}

	/*
	 * destructor
	 */
	~Shape_() {
		if (_ps2d != nullptr) {
			delete _ps2d;
		}
		_ps2d = nullptr;
	}
	/*
	 * new
	 */
	void _new() {
		if (Dim == 2) {
			if (_ps2d == nullptr) {
				return;
			} else {
				_ps2d = new S2D();
			}
		}
		_type = 0;
	}
	/*
	 * edge
	 */
	inline st size_segments() const{
		ASSERT(Dim == 2);
		return _ps2d->size_segments();
	}
	typename S2D::Segment seg(st i) const{
		ASSERT(Dim==2);
		return _ps2d->get_segment(i);
	}
	/*
	 * vertex
	 */
	inline st size_vertexs() const {
		ASSERT(Dim == 2);
		return _ps2d->size_vertexs();
	}
	typename Polygon_<VALUE>::const_ref_Point v(st i) const {
		ASSERT(Dim == 2);
		return _ps2d->v(i);
	}
	typename Polygon_<VALUE>::ref_Point v(st i) {
		ASSERT(Dim == 2);
		return _ps2d->v(i);
	}
	/*
	 * set
	 */
	void set(typename Polygon_<VALUE>::ArrP& arrp) {
		ASSERT(Dim == 2);
		if (_ps2d != nullptr) {
			delete _ps2d;
			_ps2d = nullptr;
		}
		_ps2d = new S2D(arrp);
	}
	void set(const Polygon_<VALUE>& poly) {
		ASSERT(Dim == 2);
		if (_ps2d != nullptr) {
			delete _ps2d;
			_ps2d = nullptr;
		}
		_ps2d = new S2D(poly);
	}
	/*
	 * clear
	 */
	void clear() {
		if (_ps2d == nullptr) {
			return;
		}
		_ps2d->clear();
	}
	bool empty() const {
		if (Dim == 2) {
			if (_ps2d == nullptr) {
				return true;
			}
			return _ps2d->empty();
		}
		return true;
	}
	/*
	 *  max and min
	 */
	vt max_x() const {
		if (Dim == 2) {
			return _ps2d->max_x();
		}
		return 0;
	}
	vt max_y() const {
		if (Dim == 2) {
			return _ps2d->max_y();
		}
		return 0;
	}
	vt min_x() const {
		if (Dim == 2) {
			return _ps2d->min_x();
		}
		return 0;
	}
	vt min_y() const {
		if (Dim == 2) {
			return _ps2d->min_y();
		}
		return 0;
	}
	/*
	 * volume
	 */
	vt volume() const {
		if (Dim == 2) {
			if (!this->empty()) {
				return _ps2d->area();
			}
		}
		return 0.0;
	}
	/*
	 * special function
	 */
	void find_all_seg_across(std::list<st>& l_seg_idx, Axes aix, vt v) const{
		l_seg_idx.clear();
		st i=0;
		int flag = GEL(v, _ps2d->v(i).val(aix));
		int nf;
		for(++i;i<_ps2d->size_vertexs(); i++){
			vt pv = _ps2d->v(i).val(aix);
			nf = GEL(v, pv);
			if( flag != nf){
				l_seg_idx.push_back(i-1);
				flag = nf;
			}
		}
		nf = GEL(v, _ps2d->v(0).val(aix));
		if(nf!=flag){
			l_seg_idx.push_back(_ps2d->size_vertexs()-1);
		}
	}

}
;

template<typename VALUE>
void CreatCircle(Shape_<VALUE, 2>& s, VALUE x0, VALUE y0, VALUE r, int n) {
	Polygon_<VALUE> ply;
	CreatCircle(ply, x0, y0, r, n);
	s.set(ply);
}

template<typename VALUE>
void CreatCube(Shape_<VALUE, 2>& s, VALUE x0, VALUE y0, VALUE x1, VALUE y1) {
	Polygon_<VALUE> ply;
	CreatCube(ply, x0, y0, x1, y1);
	s.set(ply);
}

template<typename TYPE, st DIM>
bool IsBoxCross(const Shape_<TYPE, DIM> &s1, const Shape_<TYPE, DIM> &s2) {
	ASSERT(DIM == 2); //temp
	typedef TYPE tF;
	tF s1maxx = s1.max_x();
	tF s1maxy = s1.max_y();
	tF s2maxx = s2.max_x();
	tF s2maxy = s2.max_y();

	tF s1minx = s1.min_x();
	tF s1miny = s1.min_y();
	tF s2minx = s2.min_x();
	tF s2miny = s2.min_y();

	tF s1dx = s1maxx - s1minx;
	tF s2dx = s2maxx - s2minx;
	tF s1dy = s1maxy - s1miny;
	tF s2dy = s2maxy - s2miny;
	if (true
			== (IsInBox(s1minx, s1maxx, s1miny, s1maxy, s2minx, s2miny)
					|| IsInBox(s1minx, s1maxx, s1miny, s1maxy, (s2minx + s2dx),
							s2miny)
					|| IsInBox(s1minx, s1maxx, s1miny, s1maxy, s2maxx, s2maxy)
					|| IsInBox(s1minx, s1maxx, s1miny, s1maxy, s2minx,
							s2miny + s2dy))) {
		return true;
	} else {
		return IsInBox(s2minx, s2maxx, s2miny, s2maxy, s1minx, s1miny)
				|| IsInBox(s2minx, s2maxx, s2miny, s2maxy, (s1minx + s1dx),
						s1miny)
				|| IsInBox(s2minx, s2maxx, s2miny, s2maxy, s1maxx, s1maxy)
				|| IsInBox(s2minx, s2maxx, s2miny, s2maxy, s1minx,
						s1miny + s1dy);
	}

}
template<typename VALUE>
void ToPath(const Shape_<VALUE, 2>& shape, ClipperLib::Path& path,
		ClipperLib::cInt coe) {
	path.resize(shape.size_vertexs() + 1);
	st i = 0;
	for (st i = 0; i < shape.size_vertexs(); i++) {
		path[i].X = ClipperLib::cInt(shape.v(i).x() * coe);
		path[i].Y = ClipperLib::cInt(shape.v(i).y() * coe);
	}
	path[shape.size_vertexs()].X = ClipperLib::cInt(shape.v(i).x() * coe);
	path[shape.size_vertexs()].Y = ClipperLib::cInt(shape.v(i).y() * coe);
}
template<typename VALUE>
void ToShape(const ClipperLib::Path& path, Shape_<VALUE, 2>& shape,
		ClipperLib::cInt coe) {
	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
	ArrPoint arrp(path.size());
	for (st i = 0; i < arrp.size(); ++i) {
		arrp[i].x() = VALUE(path[i].X / VALUE(coe));
		arrp[i].y() = VALUE(path[i].Y / VALUE(coe));
	}
	shape.set(arrp);
}

template<typename VALUE, st DIM>
void Intersect(const Shape_<VALUE, DIM>& sub, const Shape_<VALUE, DIM>& clip,
		Shape_<VALUE, DIM>& res) {
	ASSERT(DIM == 2); //temp
	/*
	 * 2D clip
	 */
	res.clear();
	// 1 box cross check
	if (!IsBoxCross(sub, clip)) {
		return;
	}
	typedef VALUE tF;
	// 2 change the coordinate type
	tF maxx = Max(sub.max_x(), clip.max_x());
	tF minx = Min(sub.min_x(), clip.min_x());
	tF maxy = Max(sub.max_y(), clip.max_y());
	tF miny = Min(sub.min_y(), clip.min_y());
	tF max = Max(maxx, maxy);
	tF min = Min(minx, miny);
	tF coe = ClipperLib::GetCoe(max, min);
	ClipperLib::Path psub;
	ToPath(sub, psub, coe);
	ClipperLib::Paths pssub(1);
	pssub[0] = psub;
	//
	ClipperLib::Path pclip;
	ToPath(clip, pclip, coe);
	ClipperLib::Paths psclip(1);
	psclip[0] = pclip;

	// 3 using clipper lib
	ClipperLib::Clipper clpr;
	clpr.AddPaths(pssub, ClipperLib::ptSubject, true);
	clpr.AddPaths(psclip, ClipperLib::ptClip, true);
	ClipperLib::Paths solution;
	clpr.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftEvenOdd,
			ClipperLib::pftEvenOdd);
	// 4 to shape
	if (!solution.empty()) {
		if (solution[0].size() > 3) {
			ToShape(solution[0], res, coe);
		}
	}
}

template<typename VALUE, st DIM>
bool IsIntersect(const Shape_<VALUE, DIM>& sub,
		const Shape_<VALUE, DIM>& clip) {
	Shape_<VALUE, DIM> res;
	Intersect(sub, clip, res);
	if (res.empty()) {
		return false;
	} else {
		return true;
	}
}

template<typename VALUE, st DIM>
void Difference(const Shape_<VALUE, DIM>& sub, const Shape_<VALUE, DIM>& clip,
		Shape_<VALUE, DIM>& res) {
	ASSERT(DIM == 2); //temp
	/*
	 * 2D clip
	 */
	res.clear();
	typedef VALUE tF;
	// 2 change the coordinate type
	tF maxx = Max(sub.max_x(), clip.max_x());
	tF minx = Min(sub.min_x(), clip.min_x());
	tF maxy = Max(sub.max_y(), clip.max_y());
	tF miny = Min(sub.min_y(), clip.min_y());
	tF max = Max(maxx, maxy);
	tF min = Min(minx, miny);
	tF coe = ClipperLib::GetCoe(max, min);
	ClipperLib::Path psub;
	ToPath(sub, psub, coe);
	ClipperLib::Paths pssub(1);
	pssub[0] = psub;
	//
	ClipperLib::Path pclip;
	ToPath(clip, pclip, coe);
	ClipperLib::Paths psclip(1);
	psclip[0] = pclip;

	// 3 using clipper lib
	ClipperLib::Clipper clpr;
	clpr.AddPaths(pssub, ClipperLib::ptSubject, true);
	clpr.AddPaths(psclip, ClipperLib::ptClip, true);
	ClipperLib::Paths solution;
	clpr.Execute(ClipperLib::ctDifference, solution, ClipperLib::pftEvenOdd,
			ClipperLib::pftEvenOdd);
	// 4 to shape
	if (!solution.empty()) {
		if (solution[0].size() > 3) {
			ToShape(solution[0], res, coe);
		}
	}
}

typedef Shape_<Float, 2> Shape2D;
typedef Shape_<Float, 3> Shape3D;

}
#endif
