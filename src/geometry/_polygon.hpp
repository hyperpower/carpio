#ifndef _POLYGON_HPP_
#define _POLYGON_HPP_

#include "../carpio_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include "../algebra/array_list.hpp"
#include <array>

namespace carpio {

template<typename TYPE>
class Polygon_ {
public:
	typedef Point_<TYPE, 2> Point;
	typedef Point_<TYPE, 2>& ref_Point;
	typedef const Point_<TYPE, 2>& const_ref_Point;
	typedef Segment_<TYPE, 2> Segment;
	typedef TYPE vt;
	typedef ArrayListT<Point> ArrP;

public:
	/*
	 *  Contructor
	 */
	Polygon_() :
			_arrp() {
	}
	Polygon_(const Polygon_ &a) {
		this->_arrp = a._arrp;
	}
	Polygon_(const ArrP &a) {
		assert(a.size() >= 3);
		this->_arrp = a;
		_trim_same_points();
		assert(_is_simple(_arrp));
	}

	Polygon_& operator=(const Polygon_ &a) {
		if (this == &a) {
			return *this;
		} else {
			this->_arrp = a._arrp;
		}
		return *this;
	}
	Float area();
	bool empty() const {
		if (_arrp.size() == 0) {
			return true;
		} else {
			return false;
		}
	}
	void show(const std::string& name = "") {
		std::cout << "Polygon  ---- " << name << "\n";
		for (st i = 0; i < _arrp.size(); i++) {
			std::cout << "> AP[ " << i << " ]=( " << _arrp[i].x() << " , "
					<< _arrp[i].y() << " )" << "\n";
		}
	}
	void reverse() {
		if (empty()) {
			return;
		} else {
			_arrp.reverse();
		}
	}

	inline st size_vertexs() const {
		return _arrp.size();
	}
	inline int size_segments() const {  //
		return _arrp.size();
	}
	const_ref_Point v(st i) const {
		assert(i < _arrp.size());
		return _arrp[i];
	}
	ref_Point v(st i) {
		assert(i < _arrp.size());
		return _arrp[i];
	}
	Segment get_segment(st i) const {
		return Segment(_arrp[i], _arrp[(i + 1) % _arrp.size()]);
	}

	vt max_x() const { //
		assert(_arrp.size() > 0);
		Float max = _arrp[0].x();
		for (st i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].x() > max) {
				max = _arrp[i].x();
			}
		}
		return max;
	}
	vt min_x() const { //
		assert(_arrp.size() > 0);
		Float min = _arrp[0].x();
		for (int i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].x() < min) {
				min = _arrp[i].x();
			}
		}
		return min;
	}
	vt max_y() const {
		assert(_arrp.size() > 0);
		Float max = _arrp[0].y();
		for (st i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].y() > max) {
				max = _arrp[i].y();
			}
		}
		return max;
	}
	vt min_y() const { //
		assert(_arrp.size() > 0);
		Float min = _arrp[0].y();
		for (int i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].y() < min) {
				min = _arrp[i].y();
			}
		}
		return min;
	}

	vt perimeter() const {
		ASSERT(empty());
		vt res = 0;
		for (int i = 1; i < _arrp.size() - 1; i++) {
			res += Distance(_arrp[i], _arrp[i + 1]);
		}
		res += Distance(_arrp[_arrp.size() - 1], _arrp[0]);
		return res;
	}

	Float winding_number(const Point_2D& ref) const {
		Float wn = 0;    // the  winding number counter
		// loop through all edges of the polygon
		for (int i = 0; i < _arrp.size() - 1; i++) {   // edge from V[i] to  V[i+1]
			wn += _winding_number(ref, _arrp[i], _arrp[i + 1]);
		}
		return wn += _winding_number(ref, _arrp[_arrp.size() - 1], _arrp[0]);
	}
	bool is_out(const Point_2D& ref) const{
		return (0 == winding_number(ref)) ? true : false;
	}
protected:
	bool _is_simple(const ArrP &ap) {
		int nap = ap.size();
		//assert(nap >= 3);
		if (nap == 3) {
			return true;
		}
		int i = 0;
		for (int j = i + 2; j < nap - 1; j++) {
			if (isIntersect(ap[i], ap[i + 1], ap[j], ap[j + 1])) {
				return false;
			}
		}
		for (i = 1; i < nap - 2; i++) {
			for (int j = i + 2; j < nap; j++) {
				if (isIntersect(ap[i], ap[i + 1], ap[j], ap[(j + 1) % nap])) {
					return false;
				}
			}
		}
		return true;
	}
	void _trim_same_points() {
		for (int i = 0; i < _arrp.size() - 1; i++) {
			if (_arrp[i] == _arrp[i + 1]) {
				_arrp.erase(i);
				i--;
			}
		}
		if (_arrp[0] == _arrp[_arrp.size() - 1]) {
			_arrp.pop_back();
		}
	}
	vt _winding_number(const Point& ref, const Point& vi,
			const Point& vip) const {
		Point_2D refh(ref.x() + 1.0, ref.y());
		Float wn = 0;
		Float a = CROSS(refh, vi, ref);
		Float b = CROSS(refh, vip, ref);
		if (a == 0 && b == 0) {
			return wn;
		}
		if (a * b < 0) {
			//vi vi+1 crosses the x
			Float c = CROSS(vip, ref, vi);
			if ((c > 0 && a < 0) || (c < 0 && a > 0)) {
				//vi vi+1 crosses the positive x
				if (a < 0) {
					wn++;
				} else {
					wn--;
				}
			}
		} else if (a == 0 && (vi.x > ref.x)) {
			if (b > 0) {
				wn = wn + 0.5;
			} else {
				wn = wn - 0.5;
			}
		} else if (b == 0 && (vip.x > ref.x)) {
			if (a < 0) {
				wn = wn + 0.5;
			} else {
				wn = wn - 0.5;
			}
		}
		return wn;
	}
	/*
	 * data
	 */
	ArrP _arrp;
};
/*
 * Function out of class
 */
template<typename TYPE>
ArrayListT<Segment_<TYPE, 2> > ToArraySegment(const Polygon_<TYPE>& p) const {
	st n = p.size_vertexs();
	ArrayListT<Segment_<TYPE, 2> > as(n);
	for (st i = 0; i < n - 1; i++) {
		as[i].reconstruct(p.v(i), p.v(i + 1));
	}
	as[n - 1].reconstruct(p.v(n - 1), p.v(0));
	return as;
}

}

#endif
