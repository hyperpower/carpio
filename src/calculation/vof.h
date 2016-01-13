#ifndef SCALAR_H_
#define SCALAR_H_

#include "../carpio_define.hpp"
#include "../domain/domain.hpp"
#include "../geometry/geometry.hpp"

#include <functional>
#include <math.h>

namespace carpio {
/*
 * The vof in the node
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Vof_face_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Vof_face_<COO_VALUE, VALUE, DIM> Self;
	typedef Vof_face_<COO_VALUE, VALUE, DIM>& ref_Self;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;

	typedef Line_<vt> Line;
	typedef Line_<vt>* pLine;
	typedef Plane_<vt> Plane;
	typedef Plane_<vt>* pPlane;

	typedef Segment_<cvt, 2> Segment;
	typedef Segment_<cvt, 3>* pSegment;

protected:
	/*
	 * Data
	 */
	pLine _pe2; //equation  2D
	pPlane _pe3; //          3D
	//
	pSegment _ps2; //the segment 2D
	//       _ps3;               3D
public:
	Vof_face_() {
		_set_all_null();
	}
	Vof_face_(vt alpha, vt A, vt B, vt C = 0) {
		if (Dim == 2) {
			_pe2 = new Line(A, B, alpha);
			_pe3 = nullptr;
			_ps2 = nullptr;
		} else {
			_pe2 = nullptr;
			_pe3 = new Plane(A, B, C, alpha);
			_ps2 = nullptr;
		}
	}
	~Vof_face_() {
		if (_pe2 != nullptr) {
			delete _pe2;
			_pe2 = nullptr;
		}
		if (_pe3 != nullptr) {
			delete _pe3;
			_pe3 = nullptr;
		}
		if (_ps2 != nullptr) {
			delete _ps2;
			_ps2 = nullptr;
		}
	}
	bool has_segment() const {
		if (Dim == 2) {
			if (_ps2 == nullptr) {
				return false;
			} else {
				if (_ps2->empty()) {
					return false;
				} else {
					return true;
				}
			}
		} else {
			return false; //unfinish
		}

	}

	bool empty() {
		if (_pe2 == nullptr && _pe3 == nullptr) {
			return true;
		} else {
			if (Dim == 2) {
				return (_pe2 == nullptr) ? true : false;
			}
		}
	}
	/*
	 * calculate face
	 */
	void cal_face() {

	}
protected:
	void _set_all_null() {
		_pe2 = nullptr;
		_pe3 = nullptr;
		_ps2 = nullptr;
	}

};
/*
 * the vof class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Vof_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Vof_<COO_VALUE, VALUE, DIM> Self;
	typedef Vof_<COO_VALUE, VALUE, DIM>& ref_Self;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef Vof_face_<COO_VALUE, VALUE, DIM> Voff;
	typedef Vof_face_<COO_VALUE, VALUE, DIM>* pVoff;
	typedef Vof_face_<COO_VALUE, VALUE, DIM>& ref_Voff;

protected:
	/*
	 * Data
	 */
	st _c_idx;  //idx c on data
	st _f_idx;  //idx vof_face on utPoint

	pGrid _pg;

public:
	/*
	 * constructor
	 */
	Vof_(pGrid pg, st ci, st fi) :
			_pg(pg), _c_idx(ci), _f_idx(fi) {
	}
	/*
	 *  calculate a and b
	 */

	/*
	 *  construct face
	 *  C get line equation
	 */
	pVoff _get_pface(pNode pn) {
		if (_has_face(pn)) {
			return CAST(pVoff, pn->data->utp(_f_idx));
		} else {
			return nullptr;
		}
	}
	bool _has_face(pNode pn) {
		if (pn->data->utp(_f_idx) == nullptr) {
			return false;
		} else {
			return true;
		}
	}
	void _new_face(pNode pn) {
		ASSERT(pn != nullptr);
		if (!_has_face(pn)) {
			pn->data->utp(_f_idx) = new Voff();
		}
	}
	void _delete_face(pNode pn) {
		ASSERT(pn != nullptr);
		if (_has_face(pn)) {
			pVoff pvf = CAST(pVoff, pn->data->utp(_f_idx));
			delete pvf;
			pn->data->utp(_f_idx) = nullptr;
		}
	}
	void _new_faces() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			_new_face(iter->get_pointer());
		}
	}
	void _delete_faces() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			_delete_face(iter->get_pointer());
		}
	}
	void _construct_face(pNode pn) {
		// 1 find stencil

		// 2 calculate normal vector (a, b)

		// 3 get eqution, save to vof_face

	}

protected:
	/**
	 * \brief   known a,b in ax+by=alpha and C, calculate alpha \n
	 *          no matter what a and b are, they will be change to abs(a) and abs(b)
	 *          return alpha, ax+by=alpha, a>b>0;
	 * \param   Float a a in ax+by=alpha
	 * \param   Float b b in ax+by=alpha
	 * \param   Float C the color function
	 * \return  alpha
	 */
	vt calAlpha(vt a, vt b, vt C) {
		vt c1, c2, alpha;
		vt absa = (a < 0) ? (-a) : a;
		vt absb = (b < 0) ? (-b) : b;
		vt m, n;
		n = (absa >= absb) ? (absa) : absb;
		m = (absa <= absb) ? (absa) : absb;
		c1 = m / 2 / n;
		c2 = 1 - c1;
		if (C >= 0 && C <= c1) {
			alpha = sqrt(2 * C * m * n);
		} else if (C > c1 && C < c2) {
			alpha = (2 * C * n + m) / 2;
		} else { //(C>=c2 && C<=1)
			alpha = m + n - sqrt(2 * (1 - C) * m * n);
		}
		return alpha;
	}

};

}
