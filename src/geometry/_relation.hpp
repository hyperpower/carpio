#ifndef _RELATION_HPP_
#define _RELATION_HPP_

#include "../carpio_define.hpp"
#include <array>
#include "_point.hpp"

namespace carpio {

/*
 *  Point -------- Point
 */
/*===============================================
 * Calculate the distance between 2 Points
 // for Point2D
 // @param   p1 Point1
 // @param   p2 Point1
 // @return     Distance
 */
template<class T>
double Distance(const Point_<T, 2> &p1, const Point_<T, 2> &p2) {
	return sqrt(
			double(
					(p1.x() - p2.x()) * (p1.x() - p2.x())
							+ (p1.y() - p2.y()) * (p1.y() - p2.y())));
}
template<class T>
double Distance(const Point_<T, 3> &p1, const Point_<T, 3> &p2) {
	return sqrt(
			double(
					(p1.x() - p2.x()) * (p1.x() - p2.x())
							+ (p1.y() - p2.y()) * (p1.y() - p2.y())
							+ (p1.z() - p2.z()) * (p1.z() - p2.z())));
}
//===============================================
// Dot multiply (sp-op).(ep-op)
// for Point2D
// @param    p1 Point1
// @param    p2 Point1
// @return      the resualt of dot multiply
//-----------------------------------------------
template<class T>
double Dot(const Point_<T, 2> &sp, const Point_<T, 2> &ep,
		const Point_<T, 2> &op) {
	return ((sp.x() - op.x()) * (ep.x() - op.x())
			+ (sp.y() - op.y()) * (ep.y() - op.y()));
}
template<class T>
Float Cro(const Point_<T, 2> &sp, const Point_<T, 2> &ep,
		const Point_<T, 2> &op) {
	return ((sp.x() - op.x()) * (ep.y() - op.y())
			- (ep.x() - op.x()) * (sp.y() - op.y()));
}
template<class T>
Float Cro(const Point_<T, 3> &v1, const Point_<T, 3> &v2,
		const Point_<T, 3> &v3, const Point_<T, 3> &v4) {
	Float a[3][3];
	for (short i = 0; i != 3; ++i) {
		a[0][i] = v1[i] - v4[i];
		a[1][i] = v2[i] - v4[i];
		a[2][i] = v3[i] - v4[i];
	}

	return a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0]
			+ a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0]
			- a[0][1] * a[1][0] * a[2][2] - a[0][0] * a[1][2] * a[2][1];
}

/*
 * Point ----------------  Segment
 */
/*
 *  function out of class
 */
enum IntersectType {
	NO_INTERSECT = 0x10000,
	INTERSECT = 0x1,
	START_1 = 0x10,
	END_1 = 0x20,
	START_2 = 0x100,
	END_2 = 0x200,
};

template<typename TYPE>
bool IsInBox(const TYPE& xmin, const TYPE& xmax, const TYPE& ymin,
		const TYPE& ymax, const TYPE& x, const TYPE& y) {
	ASSERT(ymin <= ymax);
	ASSERT(xmin <= xmax);
	if (ymin == ymax) {
		return (xmin <= x) && (x <= xmax);
	}
	if (xmin == xmax) {
		return (ymin <= y) && (y <= ymax);
	}
	return ((xmin <= x) && (x <= xmax)) && ((ymin <= y) && (y <= ymax));
}

template<typename TYPE, st DIM>
bool IsInBox(const Segment_<TYPE, DIM> &s, const Point_<TYPE, DIM> &pt) {
	ASSERT(!s.empty());
	if (s.is_horizontal()) {
		return (((s.psx() <= pt.x()) && (pt.x() <= s.pex()))
				|| ((s.pex() <= pt.x()) && (pt.x() <= s.psx())));
	}
	if (s.is_vertical()) {
		return (((s.psy() <= pt.y()) && (pt.y() <= s.pey()))
				|| ((s.pey() <= pt.y()) && (pt.y() <= s.psy())));
	}
	return (((s.psx() <= pt.x()) && (pt.x() <= s.pex()))
			|| ((s.pex() <= pt.x()) && (pt.x() <= s.psx())))
			&& (((s.psy() <= pt.y()) && (pt.y() <= s.pey()))
					|| ((s.pey() <= pt.y()) && (pt.y() <= s.psy())));
}

template<typename TYPE, st DIM>
bool IsBoxCross(const Segment_<TYPE, DIM> &s1, const Segment_<TYPE, DIM> &s2) {
	return IsInBox(s1, s2.ps()) || IsInBox(s1, s2.pe()) || IsInBox(s2, s1.ps())
			|| IsInBox(s2, s1.pe());
}
template<typename TYPE>
int OnWhichSide3(const Segment_<TYPE, 2> &s, const Point_<TYPE, 2> &pt) {
	Float rcro = Cro(s.pe(), pt, s.ps());
	if (rcro == 0.0) {
		return 0;
	} else if (rcro < 0) {
		return -1;
	} else {
		return 1;
	}
}

template<typename TYPE>
int IntersectType(const Segment_<TYPE, 2> &s1, const Segment_<TYPE, 2> &s2) {
	if (!IsBoxCross(s1, s2)) {
		return NO_INTERSECT;
	} else {
		//step 1
		int s12s = OnWhichSide3(s1, s2.ps());
		int s12e = OnWhichSide3(s1, s2.pe());
		if (s12s == s12e) { //ignore the both equal to 0, overlap is not intersect
			return NO_INTERSECT;
		}
		int s21s = OnWhichSide3(s2, s1.ps());
		int s21e = OnWhichSide3(s2, s1.pe());
		if (s21s == s21e) { //ignore the both equal to 0, overlap is not intersect
			return NO_INTERSECT;
		}
		if ((s12s + s12e) == 0 && (s21s + s21e) == 0) {
			return INTERSECT;
		}
		int res = INTERSECT;
		if (s12s == 0)
			res = res | START_2;
		if (s12e == 0)
			res = res | END_2;
		if (s21s == 0)
			res = res | START_1;
		if (s21e == 0)
			res = res | END_1;
		return res;
	}
}

template<typename TYPE, st DIM>
bool IsIntersect(const Segment_<TYPE, DIM> &s1, const Segment_<TYPE, DIM> &s2) {
	int type = IntersectType(s1, s2);
	return (type | INTERSECT) == type ? true : false;
}
template<typename TYPE, st DIM>
bool IsIntersect(const Point_<TYPE, 2>& s1s, const Point_<TYPE, 2>& s1e,
		const Point_<TYPE, 2>& s2s, const Point_<TYPE, 2>& s2e) {
	Segment_<TYPE, DIM> s1(s1s, s1e);
	Segment_<TYPE, DIM> s2(s2s, s2e);
	return IsIntersect(s1, s2);
}
/*
 * Point   -------- Polygon
 */
template<typename TYPE>
double _WindingNumber(const Point_<TYPE, 2>& ref, const Point_<TYPE, 2>& vi,
		const Point_<TYPE, 2>& vip) {
	Point_<TYPE, 2> refh(ref.x() + 1.0, ref.y());
	Float wn = 0;
	Float a = Cro(refh, vi, ref);
	Float b = Cro(refh, vip, ref);
	if (a == 0 && b == 0) {
		return wn;
	}
	if (a * b < 0) {
		//vi vi+1 crosses the x
		Float c = Cro(vip, ref, vi);
		if ((c > 0 && a < 0) || (c < 0 && a > 0)) {
			//vi vi+1 crosses the positive x
			if (a < 0) {
				wn++;
			} else {
				wn--;
			}
		}
	} else if (a == 0 && (vi.x() > ref.x())) {
		if (b > 0) {
			wn = wn + 0.5;
		} else {
			wn = wn - 0.5;
		}
	} else if (b == 0 && (vip.x() > ref.x())) {
		if (a < 0) {
			wn = wn + 0.5;
		} else {
			wn = wn - 0.5;
		}
	}
	return wn;
}
template<typename TYPE>
Float WindingNumber(const Polygon_<TYPE>& poly, const Point_<TYPE, 2>& ref) {
	Float wn = 0;    // the  winding number counter
	// loop through all edges of the polygon
	for (int i = 0; i < poly.size_vertexs() - 1; i++) { // edge from V[i] to  V[i+1]
		wn += WindingNum(ref, poly.v(i), poly.v(i + 1));
	}
	return wn += WindingNum(ref, poly.v(poly.size_vertexs() - 1), poly(0));
}

template<typename TYPE>
bool IsOut(const Polygon_<TYPE>& poly, const Point_<TYPE, 2>& ref) {
	return (0 == WindingNumber(poly, ref)) ? true : false;
}

/*
 * Segment -------- Polygon
 */
template<typename TYPE>
ArrayListT<Segment_<TYPE, 2> > ToArraySegment(const Polygon_<TYPE>& p) {
	st n = p.size_vertexs();
	ArrayListT<Segment_<TYPE, 2> > as(n);
	for (st i = 0; i < n - 1; i++) {
		as[i].reconstruct(p.v(i), p.v(i + 1));
	}
	as[n - 1].reconstruct(p.v(n - 1), p.v(0));
	return as;
}
template<typename TYPE>
bool IsSimple(const Polygon_<TYPE>& ap) {
	int nap = ap.size_vertexs();
	//assert(nap >= 3);
	if (nap == 3) {
		return true;
	}
	int i = 0;
	for (int j = i + 2; j < nap - 1; j++) {
		if (IsIntersect(ap.v(i), ap.v(i + 1), ap.v(j), ap.v(j + 1))) {
			return false;
		}
	}
	for (i = 1; i < nap - 2; i++) {
		for (int j = i + 2; j < nap; j++) {
			if (IsIntersect(ap.v(i), ap.v(i + 1), ap.v(j),
					ap.v((j + 1) % nap))) {
				return false;
			}
		}
	}
	return true;
}

template<typename TYPE>
Float Perimeter(const Polygon_<TYPE>& ap) {
	ASSERT(!ap.empty());
	Float res = 0;
	for (int i = 1; i < ap.size_() - 1; i++) {
		res += Distance(ap.v(i), ap.v(i + 1));
	}
	res += Distance(ap.v(ap.size() - 1), ap.v(0));
	return res;
}
}

#endif
