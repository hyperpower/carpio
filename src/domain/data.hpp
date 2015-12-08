#ifndef DATA_H_
#define DATA_H_

#include "domain_define.hpp"
#include "../algebra/array_list.hpp"

#include <array>

namespace carpio {

template<typename VALUE, st DIM>
class Data_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM); //

	typedef VALUE vt;
	typedef Data_<VALUE, DIM> Self;

	typedef void (*pfunction)(Self*, utPointer);

protected:
	int _idx;
	ArrayListV<vt> _center;
	ArrayListV<vt> _face[NumFaces];
	ArrayListV<vt> _vertex[NumVertexes];
	utPointer untype;
public:
	Data_() {
		_idx = 0;
		untype = nullptr;
	}
	Data_(const st& nc, const st& nf, const st& nv, utPointer utp) :
			_center(nc) {
		_idx = 0;
		for (int i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			_face[i].reconstruct(nv);
		}
		untype = utp;
	}

	inline vt& center(st i) {
		ASSERT(i < _center.size());
		return _center[i];
	}

	inline const vt& center(st i) const {
		ASSERT(i < _center.size());
		return _center[i];
	}

	inline vt& face(Direction d, st i) {
		ASSERT(i < this->NumFaces);
		return _face[i];
	}

	inline const vt& face(Direction d, st i) const {
		ASSERT(i < this->NumFaces);
		return _face[i];
	}

	inline vt& vertex(Direction d, st i) {
		ASSERT(i < this->NumVertexes);
		return _vertex[i];
	}

	inline const vt& vertex(Direction d, st i) const {
		ASSERT(i < this->NumVertexes);
		return _vertex[i];
	}

	bool empty() const {
		bool res = true;
		res = res && (_center.size() == 0);
		for (int i = 0; i < NumFaces; ++i) {
			res = res && (_face[i].size() == 0);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			res = res && (_vertex[i].size() == 0);
		}
		return res;
	}

	void show_info() const {
		std::cout << "center data:" << this->_center.size() << "\n";
		std::cout << "face data   :" << "\n";
		for (int i = 0; i < NumFaces; ++i) {
			std::cout << "    " << i << "       :" << _face[i].size() << "\n";
		}
		for (int i = 0; i < NumVertexes; ++i) {
			std::cout << "    " << i << "       :" << _vertex[i].size() << "\n";
		}
	}

};

template<typename COO_VALUE, typename VALUE, int DIM>
class PData_ {  //Point Data
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM); //

	static const int Flag_Invalid = 0;
	static const int Flag_Center = 1;
	static const int Flag_Face = 2;
	static const int Flag_Vertex = 3;

	typedef VALUE vt;
	typedef COO_VALUE cvt;
	typedef PData_<COO_VALUE, VALUE, DIM> Self;
	typedef PData_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef const PData_<COO_VALUE, VALUE, DIM>& const_ref_Self;
protected:
	cvt _p[Dim];
	ArrayListV<int> _flag;
	ArrayListV<vt> _val;
	ArrayListV<st> _idx;
	/*
	 * assign p
	 */
	void _assign_p(cvt x, cvt y = 0, cvt z = 0) {
		_p[0] = x;
		if (Dim >= 2) {
			_p[1] = y;
		}
		if (Dim >= 3) {
			_p[2] = z;
		}
	}
	void _assign_p_origin() {
		for (st i = 0; i < Dim; ++i) {
			_p[i] = 0.0;
		}
	}
	void _copy_p(const cvt _lp[]) {
		for (st i = 0; i < Dim; ++i) {
			_p[i] = _lp[i];
		}
	}
	bool is_valid_val(st i) const {
		if (_flag[i] != 0) {
			return true;
		} else {
			return false;
		}
	}
	bool _is_same_size() const {
		return (_val.size() == _flag.size() && _idx.size() == _flag.size());
	}
public:
	/*
	 *  constructor
	 */
	PData_() :
			_flag(), _val(), _idx() {
		_assign_p_origin();
	}
	PData_(st n) :
			_flag(n), _val(n), _idx(n) {
		_flag.assign(0);
		_assign_p_origin();
	}
	PData_(st n, cvt x, cvt y = 0, cvt z = 0) :
			_flag(n), _val(n), _idx(n) {
		_flag.assign(0);
		_assign_p(x, y, z);
	}
	PData_(const Self& _pd) :
			_flag(_pd._flag.size()), _val(_pd._val.size()), _idx(
					_pd._idx.size()) {
		_flag = _pd._flag;
		_val = _pd._val;
		_idx = _pd._idx;
		_copy_p(_pd._p);
	}
	ref_Self& operator=(const const_ref_Self& _pd) {
		if (this == &_pd) {
			return *this;
		}
		_flag = _pd._flag;
		_val = _pd._val;
		_idx = _pd._idx;
		_copy_p(_pd._p);
		return *this;
	}
	/*
	 *  I/O
	 */
	st size() const {
		return _flag.size();
	}

	inline cvt& x() {
		return _p[0];
	}

	inline const cvt& x() const {
		return _p[0];
	}

	inline cvt& y() {
		ASSERT(Dim >= 2);
		return _p[1];
	}

	inline const cvt& y() const {
		ASSERT(Dim >= 2);
		return _p[1];
	}

	inline cvt& z() {
		ASSERT(Dim >= 3);
		return _p[2];
	}

	inline const cvt& z() const {
		ASSERT(Dim >= 3);
		return _p[2];
	}

	inline cvt& p(Axes axes) {  //point
		switch (axes) {
		case _X_: {
			return x();
			break;
		}
		case _Y_: {
			return y();
			break;
		}
		case _Z_: {
			return z();
			break;
		}
		default:
			SHOULD_NOT_REACH;
			return x();
		}
	}

	inline void set_point(cvt x, cvt y = 0, cvt z = 0) {
		_p[0] = x;
		if (Dim >= 2) {
			_p[1] = y;
		}
		if (Dim >= 3) {
			_p[2] = z;
		}
	}

	inline vt& val(st n) {
		ASSERT(is_valid_val(n));
		return _val[n];
	}

	inline const vt& val(st n) const {
		ASSERT(is_valid_val(n));
		return _val[n];
	}

	inline st& idx(st n) {
		return _idx[n];
	}

	inline const st& idx(st n) const {
		return _idx[n];
	}

	inline int& flag(st n) {
		return _flag[n];
	}

	inline const int& flag(st n) const {
		return _flag[n];
	}

	inline void set_val(st i, const vt& val, st idx, int flag) {
		_val[i] = val;
		_idx[i] = idx;
		_flag[i] = flag;
	}

	inline bool has_valid_val() const {
		for (st i = 0; i < _flag.size(); ++i) {
			if (_flag[i] != Flag_Invalid) {
				return true;
			}
		}
		return false;
	}

	inline st count_valid_val() const {
		st res = 0;
		for (st i = 0; i < _flag.size(); ++i) {
			if (_flag[i] != Flag_Invalid) {
				res++;
			}
		}
		return res;
	}

	inline st count_invalid_val() const {
		return size() - count_valid_val();
	}

	void set_all_center(){
		for (st i = 0; i < _flag.size(); ++i) {
			_flag[i] =  Flag_Center;
		}
	}

};

}

#endif /* DATA_H_ */
