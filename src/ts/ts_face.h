/************************
 //  \file   ts_face.h
 //  \brief
 // 
 //  \author czhou
 //  \date   15 juin 2015 
 ***********************/
#ifndef TS_FACE_H_
#define TS_FACE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_triangle.h"
#include "ts_edge.h"

namespace LarusTS {

template<class TYPE, st DIM> class Surface;

template<class TYPE, st DIM>
class Face: public Triangle<TYPE, DIM> {
public:
	typedef Triangle<TYPE, DIM> base_class;
	typedef Face<TYPE, DIM> self_class;
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
	typedef Surface<TYPE, DIM> Sur;
	typedef Sur* pSur;
	typedef List<pSeg> list_pSeg;
	typedef List<pVer> list_pVer;
	typedef List<pTri> list_pTri;
	typedef List<pSur> list_pSur;
public:
	list_pSur surfaces;
public:
	Face(pEdg a, pEdg b, pEdg c, pSur sur) :
			base_class(a, b, c) {
		surfaces.push_back(sur);
	}

	void show() const {
		std::ios::fmtflags f(std::cout.flags());
		std::cout.setf(std::ios::right);
		this->e1->show();
		this->e2->show();
		this->e3->show();
		std::cout << "    -> sur = " << surfaces.size() << "\n";
		std::cout.setf(f);
	}

};

}

#endif /* TS_FACE_H_ */
