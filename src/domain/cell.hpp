#ifndef CELL_H_
#define CELL_H_

#include "domain_define.hpp"

namespace carpio {

template<typename VALUE, size_t DIM>
class Cell {
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3)?8:(DIM + DIM);

	typedef VALUE value_t;

	typedef Cell<VALUE, DIM> self;

	typedef void (*pfunction)(self*, utPointer);

protected:
	value_t _center[Dim];
	value_t _hd[Dim];
public:
	/*
	 *  constructor
	 */
	Cell() {
		for (size_t i = 0; i < Dim; ++i) {
			_center[i] = 0.0;
			_hd[i] = 0.0;
		}
	}
	Cell(const value_t& x, const value_t& dhx, //
			const value_t& y = 0.0, const value_t& dhy = 0.0, //
			const value_t& z = 0.0, const value_t& dhz = 0.0) {
		for (size_t i = 0; i < Dim; ++i) {
			if (i == 0) {
				_center[i] = x;
				ASSERT(dhx > 0.0);
				_hd[i] = dhx;
			} else if (i == 1) {
				_center[i] = y;
				ASSERT(dhy > 0.0);
				_hd[i] = dhy;
			} else if (i == 2) {
				_center[i] = z;
				ASSERT(dhz > 0.0);
				_hd[i] = dhz;
			}
		}
	}
	/*
	 *  get
	 */
	inline value_t Get(const Orientation& ori, const Axes& axes) const {
		value_t res = 0.0;
		if(axes>=Dim){
			return 0.0;
		}
		switch (ori) {
		case _M_: {
			res = _center[axes] - _hd[axes];
			break;
		}
		case _C_: {
			res = _center[axes];
			break;
		}
		case _P_: {
			res = _center[axes] + _hd[axes];
			break;
		}
		default: {
			break;
		}
		}
		return res;
	}

	inline value_t GetD(const Axes& axes) const {
		return 2.0 * _hd[axes];
	}
	inline value_t GetHD(const Axes& axes) const {
		return _hd[axes];
	}

	inline value_t Volume() const {
		value_t res = 1.0;
		for (size_t i = 0; i < Dim; ++i) {
			res *= 2.0 * _hd[i];
		}
		return res;
	}

	void Show(pfunction fun = nullptr, utPointer utp = nullptr) const {
		fun(this, utp);
	}
};


}

#endif /* CELL_H_ */
