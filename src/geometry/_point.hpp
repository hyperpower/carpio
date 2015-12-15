#ifndef _POINT_HPP_
#define _POINT_HPP_

#include "../carpio_define.hpp"
#include <array>

namespace carpio {
//Point T ====================================
template<typename TYPE, st DIM>
class Point_: public std::array<TYPE, DIM> {
public:
	static const st Dim = DIM;
	//typedef point__tag self_tag;
	typedef st size_type;
	typedef TYPE vt;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;

	//constructor
	Point_() :
			std::array<TYPE, Dim>() {
	}

	Point_(const vt& a, const vt& b, const vt& c = 0) :
			std::array<vt, Dim>() {
		this->at(0) = a;
		this->at(1) = b;
		if (Dim == 3) {
			this->at(2) = c;
		}
	}

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

	void reconstruct(const vt& a, const vt& b, const vt& c = 0) {
		this->at(0) = a;
		this->at(1) = b;
		if (Dim == 3) {
			this->at(2) = c;
		}
	}

	bool operator==(const Point_<vt, Dim> &a) const {
		if (Dim == 2) {
			return (this->at(0) == a[0] && this->at(1) == a[1]) ? true : false;
		} else {
			return (this->at(0) == a[0] && this->at(1) == a[1]
					&& this->at(2) == a[2]) ? true : false;
		}
	}
	bool operator!=(const Point_<vt, Dim> &a) const {
		if (Dim == 2) {
			return !((this->at(0) == a[0] && this->at(1) == a[1]) ? true : false);
		} else {
			return !(
					(this->at(0) == a[0] && this->at(1) == a[1]
							&& this->at(2) == a[2]) ? true : false);
		}
	}
	void show() const {
		std::cout << std::scientific << "( " << this->at(0) << " , "
				<< this->at(1);
		if (Dim == 3) {
			std::cout << " , " << this->at(2) << " )\n";
		} else {
			std::cout << " )\n";
		}
	}
	template<typename T>
	void transfer(const T&dx, const T&dy, const T&dz) {
		this->at(0) = this->at(0) + vt(dx);
		this->at(1) = this->at(1) + vt(dy);
		if (Dim == 3) {
			this->at(2) = this->at(2) + vt(dz);
		}
	}
	template<typename T>
	void scale(const T&dx, const T&dy, const T&dz) {
		this->at(0) = this->at(0) * vt(dx);
		this->at(1) = this->at(1) * vt(dy);
		if (Dim == 3) {
			this->at(2) = this->at(2) * vt(dz);
		}
	}
	inline size_type size() const {
		return size_type(Dim);
	}
};


} //end namespace

#endif /* POINT_H_ */
