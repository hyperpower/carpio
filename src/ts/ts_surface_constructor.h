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
void construct_circle( //
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

}

#endif /* TS_SURFACE_CONSTRUCTOR_H_ */
