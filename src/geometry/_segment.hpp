#ifndef _SEGMENT_HPP_
#define _SEGMENT_HPP_

#include "../carpio_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include <array>
#include "math.h"
namespace carpio {

template<typename TYPE, st DIM>
class Segment_: public std::array<Point_<TYPE, DIM>, 2> {
public:
	static const st Dim = DIM;
	typedef TYPE vt;
	typedef TYPE& ref_vt;
	typedef const TYPE& const_ref_vt;
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
		Point e(bx, by, bz);
		reconstruct(s, e);
	}
	Segment_(const_ref_Self rhs) {
		this->ps() = rhs.ps();
		this->pe() = rhs.pe();
	}
	void reconstruct(const_ref_Self rhs) {
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
		Point e(bx, by, bz);
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

	const_ref_vt psx() const {
		return this->ps().x();
	}
	const_ref_vt pex() const {
		return this->pe().x();
	}
	const_ref_vt psy() const {
		ASSERT(Dim >= 2);
		return this->ps().y();
	}
	const_ref_vt pey() const {
		ASSERT(Dim >= 2);
		return this->pe().y();
	}
	const_ref_vt psz() const {
		ASSERT(Dim >= 3);
		return this->ps().z();
	}
	const_ref_vt pez() const {
		ASSERT(Dim >= 3);
		return this->pe().z();
	}
	ref_vt psx() {
		return this->ps().x();
	}
	ref_vt pex() {
		return this->pe().x();
	}
	ref_vt psy() {
		ASSERT(Dim >= 2);
		return this->ps().y();
	}
	ref_vt pey() {
		ASSERT(Dim >= 2);
		return this->pe().y();
	}
	ref_vt psz() {
		ASSERT(Dim >= 3);
		return this->ps().z();
	}
	ref_vt pez() {
		ASSERT(Dim >= 3);
		return this->pe().z();
	}
	vt dx() const {
		return pex() - psx();
	}
	vt dy() const {
		return pey() - psy();
	}
	vt dz() const {
		return pez() - psz();
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
		this->ps().x() = psx() * xfactor;
		this->ps().y() = psy() * yfactor;
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
	void transfer(vt dx, vt dy, vt dz = 0.0) {
		if (!empty()) {
			psx() = psx() + dx;
			psy() = psy() + dy;
			pex() = pex() + dx;
			pey() = pey() + dy;
		}
		if (Dim == 3) {
			psz() = psz() + dz;
			pez() = pez() + dz;
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
		std::cout << "( " << psx() << ", " << psy();
		if (Dim == 3) {
			std::cout << ", " << psz();
		} else {
			std::cout << "";
		}
		std::cout << " )--->(" << this->pex() << ", " << pey();
		if (Dim == 3) {
			std::cout << ", " << pez();
		} else {
			std::cout << "";
		}
		std::cout << " )\n";
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

	bool is_vertical() const {
		ASSERT(!empty());
		return psx() == pex();
	}

	bool is_horizontal() const {
		ASSERT(!empty());
		return psy() == pey();
	}
	bool is_in_box(const Point_<TYPE, DIM> &pt) const {
		ASSERT(!empty());
		if (is_horizontal()) {
			return (((psx() <= pt.x) && (pt.x <= pex()))
					|| ((pex() <= pt.x) && (pt.x <= psx())));
		}
		if (is_vertical()) {
			return (((psy() <= pt.y) && (pt.y <= pey()))
					|| ((pey() <= pt.y) && (pt.y <= psy())));
		}
		return (((psx() <= pt.x) && (pt.x <= pex()))
				|| ((pex() <= pt.x) && (pt.x <= psx())))
				&& (((psy() <= pt.y) && (pt.y <= pey()))
						|| ((pey() <= pt.y) && (pt.y <= psy())));
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

}
#endif
