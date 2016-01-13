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
	ArrayListT<utPointer> _untype;
public:
	Data_() {
		_idx = 0;
	}
	Data_(const st& nc, const st& nf, const st& nv, utPointer utp) :
			_center(nc), _untype(1) {
		_idx = 0;
		for (int i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			_face[i].reconstruct(nv);
		}
		_untype[0] = utp;
	}

	inline vt& center(const st& i) {
		ASSERT(i < _center.size());
		return _center[i];
	}

	inline const vt& center(const st& i) const {
		ASSERT(i < _center.size());
		return _center[i];
	}

	inline vt& face(const Direction& d, const st& i) {
		ASSERT(i < this->NumFaces);
		return _face[i];
	}

	inline const vt& face(const Direction& d, const st& i) const {
		ASSERT(i < this->NumFaces);
		return _face[i];
	}

	inline vt& vertex(const Direction d, const st& i) {
		ASSERT(i < this->NumVertexes);
		return _vertex[i];
	}

	inline const vt& vertex(const Direction d, const st& i) const {
		ASSERT(i < this->NumVertexes);
		return _vertex[i];
	}

	inline utPointer& utp(const st& i) {
		ASSERT(i < _untype.size());
		return _untype[i];
	}

	inline const utPointer& utp(const st& i) const{
		ASSERT(i < _untype.size());
		return _untype[i];
	}

	/*
	 *  resize
	 */
	void resize_center(st len) {
		ASSERT(len >= 0);
		this->_center.resize(len);
	}
	void resize_face(st idx, st len) {
		ASSERT(idx >= 0 && idx < NumFaces);
		this->_face[idx].resize(len);
	}
	void resize_vertex(st idx, st len) {
		ASSERT(idx >= 0 && idx < NumVertexes);
		this->_vertex[idx].resize(len);
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
	PData_(ArrayListV<st> arridx, cvt x, cvt y = 0, cvt z = 0) :
			_flag(arridx.size()), _val(arridx.size()), _idx(arridx.size()) {
		_flag.assign(0);
		_idx = arridx;
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
		ASSERT(is_valid(n));
		return _val[n];
	}

	inline const vt& val(st n) const {
		ASSERT(is_valid(n));
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

	inline ArrayListV<vt>& arr_val() {
		return _val;
	}

	inline const ArrayListV<vt>& arr_val() const {
		return _val;
	}

	inline ArrayListV<st>& arr_idx() {
		return _idx;
	}

	inline const ArrayListV<st>& arr_idx() const {
		return _idx;
	}

	inline void set_val(st i, const vt& val, st idx, int flag) {
		_val[i] = val;
		_idx[i] = idx;
		_flag[i] = flag;
	}

	inline bool is_valid(st i) const {
		if (_flag[i] == Flag_Invalid) {
			return false;
		} else {
			return true;
		}
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

	void set_all_center() {
		for (st i = 0; i < _flag.size(); ++i) {
			_flag[i] = Flag_Center;
		}
	}

	/*
	 *  show
	 */

	void show_point(std::string name = "") const {
		if (name != "") {
			std::cout << name << "\n";
		}
		for (st i = 0; i < Dim; ++i) {
			std::cout << "p" << i << " = " << _p[i] << "\n";
		}
	}
};

}

#endif /* DATA_H_ */
