#ifndef _LINE_HPP_
#define _LINE_HPP_

#include "../carpio_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include <array>

namespace carpio {
template<typename TYPE>
class Line_: public std::array<TYPE, 3> {
	//The Line function defined as ax+by=alpha
public:
	typedef st size_type;
	typedef TYPE vt;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
public:
	Line_() :
			std::array<vt, 3>() {
	}
	Line_(const vt& a, const vt& b, const vt& c) :
			std::array<vt, 3>(a, b, c) {
	}
	Line_(vt ax, vt ay, vt bx, vt by) {
		//assert(!isEqual(ax, bx) || !isEqual(ay,by));
		Point_<vt, 2> p1(ax, ay), p2(bx, by);
		if (p1.x() == p2.x()) {
			this->a() = 1;
			this->b() = 0;
			this->alpha() = p1.x();
		} else if (p1.y() == p2.y()) {
			this->a() = 0;
			this->b() = 1;
			this->alpha() = p1.y();
		} else {
			this->a() = 1.0 / (p1.x() - p2.x());
			this->b() = -1.0 / (p1.y() - p2.y());
			this->alpha() = p2.x() / (p1.x() - p2.x())
					- p2.y() / (p1.y() - p2.y());
		}
	}
	Line_(const Point_<vt, 2> &p1, const Point_<vt, 2> &p2) {
		if (p1.x() == p2.x()) {
			this->a() = 1;
			this->b() = 0;
			this->alpha() = p1.x();
		} else if (p1.y() == p2.y()) {
			this->a() = 0;
			this->b() = 1;
			this->alpha() = p1.y();
		} else {
			this->a() = 1.0 / (p1.x() - p2.x());
			this->b() = -1.0 / (p1.y() - p2.y());
			this->alpha() = p2.x() / (p1.x() - p2.x())
					- p2.y() / (p1.y() - p2.y());
		}
	}
	void reconstruct(vt a, vt b, vt c) {
		if (a == 0.0 && b == 0.0) {
			a = SMALL;
			b = SMALL;
		}
		this->a() = a;
		this->b() = b;
		this->alpha() = c;
	}
	inline reference a() {
		return this->at(0);
	}
	inline reference b() {
		return this->at(1);
	}
	inline reference alpha() {
		return this->at(2);
	}
	inline const_reference a() const {
		return this->at(0);
	}
	inline const_reference b() const {
		return this->at(1);
	}
	inline const_reference alpha() const {
		return this->at(2);
	}
	vt cal_x(vt y) const {
		return (this->alpha() - this->b() * y)
				/ ((this->a() == 0.0) ? SMALL : this->a());
	}
	vt cal_y(vt x) const {
		return (this->alpha() - this->a() * x)
				/ ((this->b() == 0.0) ? SMALL : this->b());
	}
	vt slope() const {
		return -this->a() / (this->b() + SMALL);
	}
	vt intersept_x() const {
		return this->alpha() / (this->a() + SMALL);
	}
	vt intersept_y() const {
		return this->alpha() / (this->b() + SMALL);
	}
	vt norm_x() const {
		return this->a();
	}
	vt norm_y() const {
		return this->b();
	}
	bool empty() const {
		if (this->a() != 0.0 || this->b() != 0.0) {
			return true;
		} else {
			return false;
		}
	}
	void show() const {
		std::cout << this->a() << " X + " << this->b() << " Y= "
				<< this->alpha() << "\n";
	}
};

}

#endif
