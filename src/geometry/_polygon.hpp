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


}

#endif
