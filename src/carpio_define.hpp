#ifndef _DEFINE_HPP_
#define _DEFINE_HPP_

#include <array>
#include <deque>
#include <forward_list>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <assert.h>

//#define NDEBUG  //disable assert()

namespace carpio {

#define ASSERT(expr) assert(expr)
#define ASSERT_MSG(expr, msg) assert((expr)&&(msg))

#define SHOULD_NOT_REACH assert((false)&&(" >! Should not reach"))
#define CAST(type, p)           ((type)p)
#define CAST_REF(type, p)       (*((type)p))
#define _IF_TRUE_RETRUN(expr)   if(expr){return;};
#define _IF_FALSE_RETRUN(expr)  if(false==(expr)){return;};
#define SMALL                   1.0e-11
// value type
typedef std::size_t st; //size type
typedef double Float;
typedef void* utPointer;
typedef const void* const_utPointer;
//

// std  container
template<typename T>
using Vector = std::vector<T>;

//return code
#define _SUCCESS   0
#define _ERROR     1
#define _WARNING   2

template<class TYPE>
inline TYPE Max(const TYPE& a, const TYPE& b) {
	return a >= b ? a : b;
}
template<class TYPE>
inline TYPE Max(const TYPE& a, const TYPE& b, const TYPE& c) {
	return Max(Max(a, b), c);
}
template<class TYPE>
inline TYPE Max(const TYPE a[], st len) {
	TYPE max = a[0];
	for (st i = 1; i < len; ++i) {
		if (max < a[i]) {
			max = a[i];
		}
	}
	return max;
}
template<class TYPE>
inline TYPE Min(const TYPE& a, const TYPE& b) {
	return a <= b ? a : b;
}
template<class TYPE>
inline TYPE Abs(const TYPE& s) {
	return s < 0 ? -s : s;
}
template<class TYPE>
inline bool IsEqual(const TYPE& a, const TYPE& b) {
	return Abs(a - b) < SMALL;
}
template<class TYPE>
inline bool IsZero(const TYPE& a) {
	return Abs(a - 0.0) < SMALL;
}
// Greater, equal or less
template<class TYPE>
inline int GEL(const TYPE& a, const TYPE& v) {
	//greater equal or less
	if (v < a) {
		return -1;
	} else if (v == a) {
		return 0;
	} else {
		return 1;
	}
}
/*
 * geometry
 */
enum Orientation {
	_M_ = 0, //
	_P_ = 1, //
	_C_ = 2, //
};

enum Axes {
	_X_ = 0, //
	_Y_ = 1, //
	_Z_ = 2, //
};
inline Axes ToAxes(const st& i){
	ASSERT(i>=0 && i<3);
	switch(i){
	case 0:
		return _X_;
	case 1:
		return _Y_;
	case 2:
		return _Z_;
	default:
		ASSERT_MSG(false, "Error input i");
	}
	SHOULD_NOT_REACH;
	return _X_;
}
inline Axes VerticalAxes2D(const Axes& a) {
	ASSERT(a != _Z_);
	return a == _X_ ? _Y_ : _X_;
}
inline Axes VerticalAxes1(const Axes& a) {
	switch (a) {
	case _X_: {
		return _Y_;
	}
	case _Y_: {
		return _Z_;
	}
	case _Z_: {
		return _X_;
	}
	default: {
		SHOULD_NOT_REACH;
		return _X_;
	}
	}
}
inline Axes VerticalAxes2(const Axes& a) {
	switch (a) {
	case _X_: {
		return _Z_;
	}
	case _Y_: {
		return _X_;
	}
	case _Z_: {
		return _Y_;
	}
	default: {
		SHOULD_NOT_REACH;
		return _X_;
	}
	}
}

enum Plane {
	_XY_ = 24, _YZ_ = 48, _ZX_ = 40,
};

}

#endif
