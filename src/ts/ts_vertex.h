/*
 * ts_vertex.h
 *
 *  Created on: May 30, 2015
 *      Author: zhou
 */

#ifndef TS_VERTEX_H_
#define TS_VERTEX_H_

#include "ts_define.h"
#include "ts_point.h"

namespace LarusTS
{

template<class TYPE, st DIM> class Segment;
template<class TYPE, st DIM> class Surface;
template<class TYPE, st DIM> class Edge;

template<class TYPE, st DIM>
class Vertex: public Point<TYPE, DIM>
{
public:
	typedef Point<TYPE, DIM> base_class;
	typedef Vertex<TYPE, DIM> self_class;
	typedef st size_type;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef List<pSeg> list_pSeg;
	typedef List<pVer> list_pVer;
public:
	list_pSeg segments;
public:
	Vertex(const base_class& poi) :
			base_class(poi), segments()
	{
	}
	Vertex(const TYPE&x, const TYPE& y, const TYPE& z) :
			base_class(x, y, z), segments()
	{
	}

	inline bool vertex_is_unattached()
	{
		if (segments.empty())
			return true;
		return false;
	}

	/**
	 * vertex_replace:
	 * @with: another #GtsVertex.
	 *
	 * Replaces vertex this with vertex @with. this and @with must be
	 * different.  All the #GtsSegment which have @v has one of their
	 * vertices are updated.  The segments list of vertex @v is freed and
	 * @v->segments is set to %NULL.
	 */
	void vertex_replace(const Vertex<TYPE, DIM>& with)
	{
		return_if_fail(this != &with);

		for (auto iter = segments.begin(); iter != segments.end(); ++iter) {
			pSeg& s = (*iter);
			if (s->v1 != with && s->v2 != with)
				with->segments.push_front(s);
			if (s->v1 == this)
				s->v1 = with;
			if (s->v2 == this)
				s->v2 = with;
		}
		this->segments.clear();
	}
	/**
	 * vertices_are_connected:
	 * @v2: another #GtsVertex.
	 *
	 * Returns: if @v1 and @v2 are the vertices of the same #GtsSegment
	 * this segment else %NULL.
	 */
	pSeg vertices_are_connected(pVer v2)
	{
		for (auto iter = this->segments->begin(); iter != this->segments->end();
				++iter) {
			pSeg& s = (*iter);
			if (s->v1 == v2 || s->v2 == v2)
				return s;
		}
		return nullptr;
	}
	/**
	 * vertices_from_segments:
	 * @segments: a list of Segment.
	 *
	 * Returns: a list of Vertex, vertices of a Segment in list segments.
	 * Each element in the list is unique (no duplicates).
	 */
	void vertices_from_segments(list_pVer& lpver) const
	{
		lpver.clear();
		std::set<pVer> spver;
		for (auto iter = segments.begin(); iter != segments.end(); ++iter) {
			pSeg& s = (*iter);
			std::pair<typename std::set<pVer>::iterator, bool> ret;
			ret = spver.insert((*iter)->v1);
			if (ret.second == true) {
				lpver.push_back((*iter)->v1);
			}
			ret = spver.insert((*iter)->v2);
			if (ret.second == true) {
				lpver.push_back((*iter)->v2);
			}
		}
	}

	void show() const
	{
		std::ios::fmtflags f(std::cout.flags());
		std::cout.setf(std::ios::right);
		std::cout << "ver ( " << this->at(0) << " , " << this->at(1);
		if (DIM == 3) {
			std::cout << " , " << this->at(2) << " )";
		} else {
			std::cout << " )";
		}
		std::cout << " -> seg = " << segments.size() << "\n";
		std::cout.setf(f);
	}
}
;

}

#endif /* TS_TS_VERTEX_H_ */
