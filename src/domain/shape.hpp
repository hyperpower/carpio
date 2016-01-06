#ifndef SHAPE_H_
#define SHAPE_H_

#include "../carpio_define.hpp"
#include "domain_define.hpp"

#include "../geometry/geometry.hpp"

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
	Shape_()

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
	 * clear
	 */
	void clear() {
		if (_ps2d == nullptr) {
			return;
		}
		_ps2d->clear();
	}
	bool empty() {
		if (Dim == 2) {
			return _ps2d->empty();
		}
		return true;
	}
}
;

template<typename VALUE>
void CreatCircle(Shape_<VALUE, 2>& s, VALUE x, VALUE y, VALUE r, int n) {
	Float pi = 3.141592653589793238;
	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
	typedef typename Polygon_<VALUE>::Point Poi;
	ArrPoint arrp;
	for (int i = 0; i < 360; i++) {
		Float x = x + r * cos(pi / 180 * i);
		Float y = y + r * sin(pi / 180 * i);
		arrp.push_back(Poi(x , y));
	}
	return vcir;
}

}
#endif
