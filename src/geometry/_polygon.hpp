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
	typedef Segment_<TYPE, 2>& ref_Segment;
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
		ASSERT(a.size() >= 3);
		this->_arrp = a;
		_trim_same_points();
	}
	void reconstruct(const ArrP & a) {
		ASSERT(a.size() >= 3);
		this->_arrp = a;
		_trim_same_points();
	}

	Polygon_& operator=(const Polygon_ &a) {
		if (this == &a) {
			return *this;
		} else {
			this->_arrp = a._arrp;
		}
		return *this;
	}
	vt area() const {
		if (empty()) {
			return 0.0;
		}
		vt s = 0.0;
		for (st i = 1; i < _arrp.size() - 1; i++) {
			s = s + Cro(_arrp[i + 1], _arrp[i], _arrp[0]); // det to cro
		}
		return Abs(s) / 2.0;
	}
	void clear() {
		_arrp.resize(0);
	}
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
	inline st size_segments() const {  //
		return _arrp.size();
	}
	const_ref_Point v(st i) const {
		ASSERT(i < _arrp.size());
		return _arrp[i];
	}
	ref_Point v(st i) {
		ASSERT(i < _arrp.size());
		return _arrp[i];
	}
	Segment get_segment(st i) const {
		return Segment(_arrp[i], _arrp[(i + 1) % _arrp.size()]);
	}

	vt max_x() const { //
		ASSERT(_arrp.size() > 0);
		Float max = _arrp[0].x();
		for (st i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].x() > max) {
				max = _arrp[i].x();
			}
		}
		return max;
	}
	vt min_x() const { //
		ASSERT(_arrp.size() > 0);
		Float min = _arrp[0].x();
		for (st i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].x() < min) {
				min = _arrp[i].x();
			}
		}
		return min;
	}
	vt max_y() const {
		ASSERT(_arrp.size() > 0);
		Float max = _arrp[0].y();
		for (st i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].y() > max) {
				max = _arrp[i].y();
			}
		}
		return max;
	}
	vt min_y() const { //
		ASSERT(_arrp.size() > 0);
		Float min = _arrp[0].y();
		for (st i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].y() < min) {
				min = _arrp[i].y();
			}
		}
		return min;
	}

	st find_closest_vertex(const Point& p) const {
		ASSERT(_arrp.size() > 0);
		st idx = 0;
		vt mindis = Distance(_arrp[0], p);
		for (st i = 1; i < _arrp.size(); i++) {
			vt dis = Distance(_arrp[i], p);
			if (dis < mindis) {
				idx = i;
				mindis = dis;
			}
		}
		return idx;
	}

	st find_closest_vertex(const vt& x, const vt& y) const {
		Point p(x, y);
		return this->find_closest_vertex(p);
	}
	/*
	 * special function
	 * aix = x, y
	 * v   = (a value on coordinate)
	 * l_seg_idx (return)
	 *
	 * Find all the segments across x=v, y=v
	 */
	void find_seg_across(std::list<st>& l_seg_idx, Axes aix, vt val) {
		l_seg_idx.clear();
		st i = 0;
		int flag = GEL(val, v(i).val(aix));
		int nf;
		for (++i; i < size_vertexs(); i++) {
			vt pv = v(i).val(aix);
			nf = GEL(val, pv);
			if (flag != nf) {
				l_seg_idx.push_back(i - 1);
				flag = nf;
			}
		}
		nf = GEL(val, v(0).val(aix));
		if (nf != flag) {
			l_seg_idx.push_back(size_vertexs() - 1);
		}
	}

	void find_seg_connect_to_vertex(std::list<st>& l_seg_idx,
			st ver_idx) const {
		st n = this->size_vertexs();
		ASSERT(ver_idx < n);
		l_seg_idx.clear();
		st prev = ver_idx - 1;
		if (ver_idx == 0) {
			prev = n - 1;
		}
		l_seg_idx.push_back(prev);
		l_seg_idx.push_back(ver_idx);
	}
protected:
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

	/*
	 * data
	 */
	ArrP _arrp;
};
/*
 * Function out of class
 */
template<typename VALUE>
void CreatCircle(Polygon_<VALUE>& s, VALUE x0, VALUE y0, VALUE r, int n) {
	ASSERT(n >= 3);
	Float pi = 3.141592653589793238;
	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
	typedef typename Polygon_<VALUE>::Point Poi;
	ArrPoint arrp;
	for (int i = 0; i < n; i++) {
		Float x = x0 + r * cos(2. * pi / float(n) * i);
		Float y = y0 + r * sin(2. * pi / float(n) * i);
		arrp.push_back(Poi(x, y));
	}
	s.reconstruct(arrp);
}

template<typename VALUE>
void CreatCube(Polygon_<VALUE>& s, VALUE x0, VALUE y0, VALUE x1, VALUE y1) {
	ASSERT(x1 > x0);
	ASSERT(y1 > y0);
	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
	typedef typename Polygon_<VALUE>::Point Poi;
	ArrPoint arrp;
	VALUE dx = x1 - x0;
	VALUE dy = y1 - y0;
	arrp.push_back(Poi(x0, y0));
	arrp.push_back(Poi(x0 + dx, y0));
	arrp.push_back(Poi(x1, y1));
	arrp.push_back(Poi(x0, y0 + dy));
	s.reconstruct(arrp);
}

}

#endif
