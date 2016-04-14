#ifndef _POINT_HPP_
#define _POINT_HPP_

#include "../carpio_define.hpp"
#include <array>
#include <sstream>

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
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
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

	const_reference val(Axes axi) const {
		switch (axi) {
		case _X_: {
			return this->at(0);
		}
		case _Y_: {
			ASSERT(Dim>=2);
			return this->at(1);
		}
		case _Z_: {
			ASSERT(Dim>=3);
			return this->at(2);
		}
		default:{
			SHOULD_NOT_REACH;
		}
		}
		return this->at(0); //make compile happy;
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
	std::string to_string() const {
		std::stringstream sstr;
		sstr.precision(4);
		sstr << "( " << this->at(0) << " , " << this->at(1);
		if (Dim == 3) {
			sstr << " , " << this->at(2) << " )";
		} else {
			sstr << " )";
		}
		return sstr.str();
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
	inline pointer ptr() {
		return this->data();
	}
	inline const_pointer ptr() const {
		return this->data();
	}
};

} //end namespace

#endif /* POINT_H_ */
