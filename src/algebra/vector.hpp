#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "../carpio_define.hpp"
#include <array>
#include <iostream>
#include <math.h>

namespace carpio {
//Vector T ====================================
template<typename TYPE, st DIM>
class Vector_: public std::array<TYPE, DIM> {
public:
	static const st Dim = DIM;
	//typedef Vector__tag self_tag;
	typedef st size_type;
	typedef TYPE vt;
	typedef TYPE& reference;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef const TYPE& const_reference;

	//constructor
	Vector_() :
			std::array<TYPE, Dim>() {
	}

	Vector_(const vt& a, const vt& b, const vt& c = 0, const vt& d = 0) :
			std::array<vt, Dim>() {
		this->at(0) = a;
		if (Dim >= 2) {
			this->at(1) = b;
		}
		if (Dim >= 3) {
			this->at(2) = c;
		}
		if (Dim >= 4) {
			this->at(3) = d;
		}
	}

	/*
	 * the order is x -> y -> z -> w
	 */
	const_reference x() const {
		return this->at(0);
	}

	reference x() {
		return this->at(0);
	}

	const_reference y() const {
		return this->at(1);
	}

	reference y() {
		return this->at(1);
	}

	const_reference z() const {
		ASSERT(Dim == 3);
		return this->at(2);
	}

	reference z() {
		ASSERT(Dim == 3);
		return this->at(2);
	}

	const_reference w() const {
		ASSERT(Dim == 4);
		return this->at(3);
	}

	reference w() {
		ASSERT(Dim == 4);
		return this->at(3);
	}

	void reconstruct(const vt& a, const vt& b, const vt& c = 0,
			const vt& d = 0) {
		this->at(0) = a;
		if (Dim >= 2) {
			this->at(1) = b;
		}
		if (Dim >= 3) {
			this->at(2) = c;
		}
		if (Dim >= 4) {
			this->at(3) = d;
		}

	}

	bool operator==(const Vector_<vt, Dim> &a) const {
		bool res = true;
		for (st i = 0; i < Dim; ++i) {
			res = res && (this->at(i) == a.at(i));
		}
		return res;
	}
	bool operator!=(const Vector_<vt, Dim> &a) const {
		return !(this->operator==(a));
	}

	Vector_& operator+=(const Vector_& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) += r.at(i);
		}
		return *this;
	}

	Vector_& operator-=(const Vector_& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) -= r.at(i);
		}
		return *this;
	}

	Vector_& operator*=(const Vector_& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) *= r.at(i);
		}
		return *this;
	}

	Vector_& operator/=(const Vector_& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) /= r.at(i);
		}
		return *this;
	}

	Vector_& operator+=(const vt& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) += r;
		}
		return *this;
	}
	Vector_& operator-=(const vt& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) -= r;
		}
		return *this;
	}

	Vector_& operator*=(const vt& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) *= r;
		}
		return *this;
	}
	Vector_& operator/=(const vt& r) {
		for (st i = 0; i < Dim; ++i) {
			this->at(i) /= r;
		}
		return *this;
	}

	void show() const {
		std::cout << std::scientific << "( ";
		for (st i = 0; i < Dim; ++i) {
			std::cout << this->at(i);
			if (i == Dim - 1) {
				std::cout << " )\n";
			} else {
				std::cout << " , ";
			}
		}
	}

	inline size_type size() const {
		return size_type(Dim);
	}
	inline pointer ptr() {
		return this->data();
	}
	inline const_pointer ptr() const {
		return this->data();
	}

	Vector_ cross(const Vector_& v) const {
		ASSERT(Dim == 3);
		vt _x = this->y() * v.z() - this->z() * v.y();
		vt _y = this->z() * v.x() - this->x() * v.z();
		vt _z = this->x() * v.y() - this->y() * v.x();

		return Vector_(_x, _y, _z);
	}

	Vector_ normalize() {
		vt sum = 0;
		for (st i = 0; i < Dim; ++i) {
			sum += (this->at(i) * this->at(i));
		}
		vt l = sqrt(sum);
		this->operator /=(l);

		return *this;
	}
};

template<typename TYPE, st DIM>
inline Vector_<TYPE, DIM> operator+(const Vector_<TYPE, DIM>& l,
		const Vector_<TYPE, DIM>& r) {
	Vector_<TYPE, DIM> Ret(l);
	Ret += r;

	return Ret;
}

template<typename TYPE, st DIM>
inline Vector_<TYPE, DIM> operator-(const Vector_<TYPE, DIM>& l,
		const Vector_<TYPE, DIM>& r) {
	Vector_<TYPE, DIM> Ret(l);
	Ret -= r;

	return Ret;
}

template<typename TYPE, st DIM>
inline Vector_<TYPE, DIM> operator*(const Vector_<TYPE, DIM>& l,
		const Vector_<TYPE, DIM>& r) {
	Vector_<TYPE, DIM> Ret(l);
	Ret *= r;

	return Ret;
}

template<typename TYPE, st DIM>
inline Vector_<TYPE, DIM> operator/(const Vector_<TYPE, DIM>& l,
		const Vector_<TYPE, DIM>& r) {
	Vector_<TYPE, DIM> Ret(l);
	Ret /= r;

	return Ret;
}

} //end namespace

#endif /* Vector_H_ */
