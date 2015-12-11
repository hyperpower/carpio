#ifndef _SEGMENT_HPP_
#define _SEGMENT_HPP_

#include "../carpio_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include <array>

namespace carpio {

template<typename TYPE, st DIM>
class Segment_: public std::array<Point_<TYPE, DIM>, 2> {
public:
	static const st Dim = DIM;
	typedef TYPE vt;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef Segment_<TYPE, DIM> Self;
	typedef Segment_<TYPE, DIM>& ref_Self;
	typedef const Segment_<TYPE, DIM>& const_ref_Self;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>* pPoint;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
public:
	Segment_() :
			std::array<Point_<TYPE, DIM>, 2>() {
		_set_empty();
	}
	Segment_(const Point& s, const Point& e) {
		ASSERT(s != e);
		this->ps() = s;
		this->pe() = e;
	}
	Segment_(const vt& ax, const vt& bx, //
			const vt& ay, const vt& by,  //
			const vt& az = 0, const vt& bz = 0) {
		Point s(ax, ay, az);
		Point e(ax, ay, az);
		reconstruct(s, e);
	}
	Segment_(const_ref_Self rhs) {
		this->ps() = rhs.ps();
		this->pe() = rhs.pe();
	}
	void reconstruct(const Point& s, const Point& e) {
		ASSERT(s != e);
		this->ps() = s;
		this->pe() = e;
	}
	void reconstruct(const vt& ax, const vt& bx, //
			const vt& ay, const vt& by,  //
			const vt& az = 0, const vt& bz = 0) {
		Point s(ax, ay, az);
		Point e(ax, ay, az);
		reconstruct(s, e);
	}

	bool operator==(const_ref_Self rhs) const {
		return (this->ps() == rhs.ps() && this->pe() == rhs.pe()) ? true : false;
	}

	ref_Point ps() {
		return this->at(0);
	}
	const_ref_Point ps() const {
		return this->at(0);
	}
	ref_Point pe() {
		return this->at(1);
	}
	const_ref_Point pe() const {
		return this->at(1);
	}
	Point pc() const {
		return Point((pex() + psx()) * 0.5, (pey() + psy()) * 0.5,
				(Dim == 3) ? ((pez() + psz()) * 0.5) : 0);
	}

	vt psx() const {
		return this->ps().x();
	}
	vt pex() const {
		return this->pe().x();
	}
	vt psy() const {
		ASSERT(Dim >= 2);
		return this->ps().y();
	}
	vt pey() const {
		ASSERT(Dim >= 2);
		return this->pe().y();
	}
	vt psz() const {
		ASSERT(Dim >= 3);
		return this->ps().z();
	}
	vt pez() const {
		ASSERT(Dim >= 3);
		return this->pe().z();
	}

	vt length() const {
		vt len = 0.0;
		len = sqrt(
				double(
						(psx() - pex()) * (psx() - pex())
								+ (psy() - pey()) * (psy() - pey())));
		return len;
	}
	vt slope() const {
		ASSERT(Dim == 2);
		return (pey() - psy()) / (pex() - psx() + SMALL);
	}

	void scale(vt xfactor, vt yfactor, vt zfactor = 1) {
		psx() = psx() * xfactor;
		psy() = psy() * yfactor;
		if (Dim == 3) {
			psz() = psz() * zfactor;
		}
		pex() = pex() * xfactor;
		pey() = pey() * yfactor;
		if (Dim == 3) {
			pez() = pez() * zfactor;
		}
		if (pe() == ps()) {
			_set_empty();
		}
	}
	void transfer(vt dx, vt dy, vt dz) {
		if (!empty()) {
			psx() = psx() + dx;
			psy() = psy() + dy;
			pex() = pex() + dx;
			pey() = pey() + dy;
		}
	}

	bool empty() const {
		if (psx() == 0.0 && psy() == 0.0 && pex() == 0.0 && pey() == 0.0
				&& ((Dim == 3) ? (psz() == 0.0 && pez() == 0.0) : true)) {
			return true;
		} else {
			return false;
		}
	}
	void show() const {
		std::cout.precision(4);
		std::cout << "( " << psx() << ", " << psy() << ", "
				<< ((Dim == 3) ? psz() : "") << " )--->(" << pex() << ", "
				<< pey() << ", " << ((Dim == 3) ? psz() : "") << " ) \n";
	}

	/*
	 *  compare
	 */

	bool is_gt_x(const vt& v) const {    //>
		return (this->pex() > v && this->psx() > v);
	}
	bool is_gt_y(const vt& v) const {    //>
		return (this->pey() > v && this->psy() > v);
	}
	bool is_ge_x(const vt& v) const {    //>=
		return (this->pex() >= v && this->psx() >= v);
	}
	bool is_ge_y(const vt& v) const {    //>=
		return (this->pey() >= v && this->psy() >= v);
	}

	bool is_lt_x(const vt& v) const {    //<
		return (this->pex() < v && this->psx() < v);
	}
	bool is_lt_y(const vt& v) const {    //<
		return (this->pey() < v && this->psy() < v);
	}
	bool is_lt_z(const vt& v) const {    //<
		ASSERT(Dim == 3);
		return (this->pez() < v && this->psz() < v);
	}
	bool is_le_x(const vt& v) const {    //<=
		return (this->pex() <= v && this->psx() <= v);
	}
	bool is_le_y(const vt& v) const {    //<=
		return (this->pey() <= v && this->psy() <= v);
	}

protected:
	void _set_empty() {
		ps().x() = 0.0;
		pe().x() = 0.0;
		ps().y() = 0.0;
		pe().y() = 0.0;
		if (Dim == 3) {
			ps().z() = 0.0;
			pe().z() = 0.0;
		}
	}
};
/*
 *  function out of class
 */

}
#endif
