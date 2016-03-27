/*
 * ts_triangle.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_TRIANGLE_H_
#define TS_TRIANGLE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_edge.h"

#include "ts_tri_moller.h"

namespace LarusTS {

template<class TYPE, st DIM> class Surface;
template<class TYPE, st DIM> class Face;

template<class TYPE, st DIM>
class Triangle {
public:
	typedef Triangle<TYPE, DIM> self_class;
	typedef st size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Seg* pSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef Edg* pEdg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Ver* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Tri* pTri;
	typedef Face<TYPE, DIM> Fac;
	typedef Fac* pFac;
	typedef List<pSeg> list_pSeg;
	typedef List<pVer> list_pVer;
	typedef List<pTri> list_pTri;
public:
	pEdg e1;
	pEdg e2;
	pEdg e3;
protected:
	Int edges_check(pEdg a, pEdg b, pEdg c);
public:

	Triangle(pEdg a, pEdg b, pEdg c, pFac f) {
		e1 = a;
		e2 = b;
		e3 = c;
		assert(edges_check(a, b, c) == NO_ERROR);
		e1->faces.push_back(f);
		e2->faces.push_back(f);
		e3->faces.push_back(f);
	}

	pVer get_vertex1() const {
		return e1->v1;
	}
	pVer get_vertex2() const {
		return e1->v2;
	}
	pVer get_vertex3() const {
		return (e1->v1 == e2->v1) || (e1->v2 == e2->v1) ? e2->v2 : e2->v1;
	}
	pVer get_vertex(st i) const {
		assert(i < 3);
		if (i == 0) {
			return get_vertex1();
		}
		if (i == 1) {
			return get_vertex2();
		}
		if (i == 2) {
			return get_vertex3();
		}
		return nullptr;
	}
	/**
	 * triangle_normal:
	 * @t: a #GtsTriangle.
	 * @x: the x coordinate of the normal.
	 * @y: the y coordinate of the normal.
	 * @z: the z coordinate of the normal.
	 *
	 * Computes the coordinates of the oriented normal of @t as the
	 * cross-product of two edges, using the left-hand rule. The normal is
	 * not normalized.  If this triangle is part of a closed and oriented
	 * surface, the normal points to the outside of the surface.
	 */
	Poi normal() const {
		pVer v1 = get_vertex1();
		pVer v2 = get_vertex2();
		pVer v3 = get_vertex3();

		TYPE x1 = v2->x() - v1->x();
		TYPE y1 = v2->y() - v1->y();
		TYPE z1 = v2->z() - v1->z();

		TYPE x2 = v3->x() - v1->x();
		TYPE y2 = v3->y() - v1->y();
		TYPE z2 = v3->z() - v1->z();

		TYPE x = y1 * z2 - z1 * y2;
		TYPE y = z1 * x2 - x1 * z2;
		TYPE z = x1 * y2 - y1 * x2;

		return Poi(x, y, z);
	}
	/**
	 * revert:
	 * @t: a #GtsTriangle.
	 *
	 * Changes the orientation of triangle @t, turning it inside out.
	 */
	void revert() {
		pEdg e;

		e = this->e1;
		this->e1 = this->e2;
		this->e2 = e;
	}

	/**
	 * gts_triangle_area:
	 * @t: a #GtsTriangle.
	 *
	 * Returns: the area of the triangle @t.
	 */
	TYPE area() const {
		Poi n = this->normal();
		return sqrt(
				n.x() * n.x() + n.y() * n.y()
						+ ((DIM == 3) ? n.z() * n.z() : 0.0)) / 2.;
	}

	Poi centroid() const {
		pVer v1 = get_vertex1();
		pVer v2 = get_vertex2();
		pVer v3 = get_vertex3();
		TYPE x = (v1->x() + v2->x() + v3->x()) / 3.0;
		TYPE y = (v1->y() + v2->y() + v3->y()) / 3.0;
		TYPE z = (DIM == 3) ? (v1->z() + v2->z() + v3->z()) / 3.0 : 0.0;
		return Poi(x, y, z);
	}

	// output ==================
	void output_vtk(const String& fn) const;

};
template<class TYPE, st DIM>
Int Triangle<TYPE, DIM>::edges_check(  //
		pEdg e1,   //
		pEdg e2,   //
		pEdg e3) { //
	_return_val_if_fail(e1 != nullptr, ERR_NULL_POINTER);
	_return_val_if_fail(e2 != nullptr, ERR_NULL_POINTER);
	_return_val_if_fail(e3 != nullptr, ERR_NULL_POINTER);
	_return_val_if_fail(e1 != e2 && e1 != e3 && e2 != e3, ERR_DEGERATE);

	if (e1->v1 == e2->v1) {
		_return_val_if_fail(segment_connect((e3), (e1->v2), (e2->v2)),
				ERR_OTHER);
	} else if (e1->v2 == e2->v1) {
		_return_val_if_fail(segment_connect((e3), (e1->v1), (e2->v2)),
				ERR_OTHER);
	} else if (e1->v2 == e2->v2) {
		_return_val_if_fail(segment_connect((e3), (e1->v1), (e2->v1)),
				ERR_OTHER);
	} else if (e1->v1 == e2->v2) {
		_return_val_if_fail(segment_connect((e3), (e1->v2), (e2->v1)),
				ERR_OTHER);
	}
	return NO_ERROR;
}
template<class TYPE, st DIM>
void Triangle<TYPE, DIM>::output_vtk(const String& fn) const {
	FILE* fptr = fopen(fn.c_str(), "w"); //write
	if (fptr == NULL) {
		std::cerr << "!> Open file error! " << fn << " \n";
		exit(-1);
	}
	fprintf(fptr, "# vtk DataFile Version 2.0\n"
			"Generated by LarusTS\n"
			"ASCII\n"
			"DATASET POLYDATA\n"
			"POINTS %lu float\n", 3);
	Map<pVer, uInt> m_veridx;
	pVer ver1 = this->get_vertex1();
	pVer ver2 = this->get_vertex2();
	pVer ver3 = this->get_vertex3();

	fprintf(fptr, "%f %f %f \n", ver1->x(), ver1->y(), ver1->z());
	fprintf(fptr, "%f %f %f \n", ver2->x(), ver2->y(), ver2->z());
	fprintf(fptr, "%f %f %f \n", ver3->x(), ver3->y(), ver3->z());

	fprintf(fptr, "POLYGONS %lu %lu\n", 1, 4);
	fprintf(fptr, "3 %u %u %u\n", 0, 1, 2);
	fclose(fptr);
}

/*
 *  function out of class
 */
/**
 * gts_triangles_are_compatible:
 * @t1: a #GtsTriangle.
 * @t2: a #GtsTriangle.
 * @e: a #GtsEdge used by both @t1 and @t2.
 *
 * Checks if @t1 and @t2 have compatible orientations i.e. if @t1 and
 * @t2 can be part of the same surface without conflict in the surface
 * normal orientation.
 *
 * Returns: %TRUE if @t1 and @t2 are compatible, %FALSE otherwise.
 */
template<class TYPE, st DIM>
bool AreCompatible( //
		const Triangle<TYPE, DIM> * t1,  //
		const Triangle<TYPE, DIM> * t2,  //
		const Edge<TYPE, DIM> * e) {
	typedef Edge<TYPE, DIM> Edg;
	typedef Edg* pEdg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Ver* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Tri* pTri;
	pEdg e1 = nullptr, e2 = nullptr;

	_return_val_if_fail(t1 != nullptr, false);
	_return_val_if_fail(t2 != nullptr, false);
	_return_val_if_fail(e != nullptr, false);

	if (t1->e1 == e)
		e1 = t1->e2;
	else if (t1->e2 == e)
		e1 = t1->e3;
	else if (t1->e3 == e)
		e1 = t1->e1;
	else
		assert(false);
	if (t2->e1 == e)
		e2 = t2->e2;
	else if (t2->e2 == e)
		e2 = t2->e3;
	else if (t2->e3 == e)
		e2 = t2->e1;
	else
		assert(false);
	if (e1->v1 == e2->v1
			|| e1->v1 == e2->v2
			|| e1->v2 == e2->v1
			|| e1->v2 == e2->v2)
		return false;
	return true;
}

/**
 * triangles_common_edge:
 * @t1: a #GtsTriangle.
 * @t2: a #GtsTriangle.
 *
 * Returns: a #GtsEdge common to both @t1 and @t2 or %NULL if @t1 and @t2
 * do not share any edge.
 */
template<class TYPE, st DIM>
Edge<TYPE,DIM>* GetCommonpEdge ( //
		const Triangle<TYPE, DIM> * t1,  //
		const Triangle<TYPE, DIM> * t2)
{
  _return_val_if_fail (t1 != nullptr, nullptr);
  _return_val_if_fail (t2 != nullptr, nullptr);

  if (t1->e1 == t2->e1 || t1->e1 == t2->e2 || t1->e1 == t2->e3)
    return t1->e1;
  if (t1->e2 == t2->e1 || t1->e2 == t2->e2 || t1->e2 == t2->e3)
    return t1->e2;
  if (t1->e3 == t2->e1 || t1->e3 == t2->e2 || t1->e3 == t2->e3)
    return t1->e3;
  return nullptr;
}

}
#endif /* _TS_TRIANGLE_H_ */
