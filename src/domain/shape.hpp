#ifndef SHAPE_H_
#define SHAPE_H_

#include "../carpio_define.hpp"
#include "domain_define.hpp"

#include "../geometry/geometry.hpp"

#include <functional>
#include <math.h>

namespace carpio {

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
	/*
	 * destructor
	 */
	~Shape_() {
		if (_ps2d != nullptr) {
			delete _ps2d;
		}
		_ps2d = nullptr;
	}

};



}
#endif
