/************************
 //  \file   ts_surface_constructor.h
 //  \brief
 // 
 //  \author czhou
 //  \date   26 juin 2015 
 ***********************/
#ifndef TS_SURFACE_CONSTRUCTOR_H_
#define TS_SURFACE_CONSTRUCTOR_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_edge.h"
#include "ts_face.h"
#include <fstream>
#include <sstream>
#include <math.h>

namespace LarusTS {
template<class TYPE, st DIM>
void ConstructCircle( //
		Surface<TYPE, DIM>& sur,          // the surface
		uInt n,                           //the number of triangle
		const TYPE& r) {
	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef Edge<TYPE, DIM>* pEdg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Triangle<TYPE, DIM>* pTri;
	typedef Face<TYPE, DIM> Fac;
	typedef Face<TYPE, DIM>* pFac;

	Vector<pVer> v_vertex;
	Vector<pEdg> v_edge;
	Vector<pFac> v_face;
	Vertex<TYPE, DIM>* pverc = new Vertex<TYPE, DIM>(0, 0, 0);
	double da = 2 * PI / n;
	v_vertex.push_back(pverc);
	for (uInt i = 0; i < n; i++) {
		TYPE x = r * cos(i * da);
		TYPE y = r * sin(i * da);
		TYPE z = 0;
		Vertex<TYPE, DIM>* pv = new Vertex<TYPE, DIM>(x, y, z);
		v_vertex.push_back(pv);
	}
	//edge
	for (uInt i = 0; i < n; i++) {
		pEdg pe = new Edg(v_vertex[0], v_vertex[i + 1]);
		v_edge.push_back(pe);
	}
	for (uInt i = 1; i < n; i++) {
		pEdg pe = new Edg(v_vertex[i], v_vertex[i + 1]);
		v_edge.push_back(pe);
	}
	pEdg pe = new Edg(v_vertex[n], v_vertex[1]);
	v_edge.push_back(pe);
	//surface
	sur.clear();
	for (uInt i = 0; i < n; i++) {
		pFac pfac = new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i], &sur);
		sur.faces.insert(pfac);
	}
	for (pVer pv : v_vertex) {
		sur.c_vertex.insert(pv);
	}
	for (pEdg pe : v_edge) {
		sur.c_edge.insert(pe);
	}
}
template<class TYPE, st DIM>
void ConstructTriangle(
		Surface<TYPE, DIM>& sur,          // the surface
		const Point<TYPE, DIM> p1, const Point<TYPE, DIM> p2,
		const Point<TYPE, DIM> p3) {
	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef Edge<TYPE, DIM>* pEdg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Triangle<TYPE, DIM>* pTri;
	typedef Face<TYPE, DIM> Fac;
	typedef Face<TYPE, DIM>* pFac;

	Vertex<TYPE, DIM>* pv1 = new Vertex<TYPE, DIM>(p1.x(), p1.y(), p1.z());
	Vertex<TYPE, DIM>* pv2 = new Vertex<TYPE, DIM>(p2.x(), p2.y(), p2.z());
	Vertex<TYPE, DIM>* pv3 = new Vertex<TYPE, DIM>(p3.x(), p3.y(), p3.z());

	pEdg pe1 = new Edg(pv1, pv2);
	pEdg pe2 = new Edg(pv2, pv3);
	pEdg pe3 = new Edg(pv3, pv1);

	sur.clear();
	pFac pfac = new Fac(pe1, pe2, pe3, &sur);
	sur.faces.insert(pfac);
	sur.c_vertex.insert(pv1);
	sur.c_vertex.insert(pv2);
	sur.c_vertex.insert(pv3);
	sur.c_edge.insert(pe1);
	sur.c_edge.insert(pe2);
	sur.c_edge.insert(pe3);
}


}

#endif /* TS_SURFACE_CONSTRUCTOR_H_ */
