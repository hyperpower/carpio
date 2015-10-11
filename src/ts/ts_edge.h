/*
 * ts_edge.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_EDGE_H_
#define TS_EDGE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"

namespace LarusTS {

template<class TYPE, st DIM> class Triangle;
template<class TYPE, st DIM> class Surface;

template<class TYPE, st DIM>
class Edge: public Segment<TYPE, DIM> {
public:
	typedef Segment<TYPE, DIM> base_class;
	typedef Edge<TYPE, DIM> self_class;
	typedef st size_type;

	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Triangle<TYPE, DIM>* pTri;
	typedef List<pSeg> list_pSeg;
	typedef List<pVer> list_pVer;
	typedef List<pTri> list_pTri;
public:
	list_pTri triangles;

	Edge(pVer a, pVer b) :
		base_class(a, b) {
	}

	// show =====================================
	void show() const;

};

template<class TYPE, st DIM>
void Edge<TYPE, DIM>::show() const{
	std::cout<<"edg ("<< this->v1->at(0)<<" "<<  this->v1->at(1) << " " << this->v1->at(2)<<")";
	std::cout<<"->("<< this->v2->at(0)<<" "<<  this->v2->at(1) << " " << this->v2->at(2)<<") ";
	std::cout<< " -> tri "<<triangles.size()<<"\n";
}

}

#endif /* TS_TS_EDGE_H_ */
