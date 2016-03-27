#ifndef _POLYGON_BOOLEAN_HPP_
#define _POLYGON_BOOLEAN_HPP_

#include "../carpio_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include "_polygon.hpp"
#include "../algebra/array_list.hpp"
#include <list>

namespace carpio {

template<typename T>
class _WANode_ {
public:
	typedef Point_<T, 2> Point;
	typedef Point_<T, 2>& ref_Point;
	typedef const Point_<T, 2>& const_ref_Point;

	typedef std::list<_WANode_<T> > WAList;
	typedef std::list<_WANode_<T> >* pWAList;

	typedef _WANode_<T> Self;
	typedef _WANode_<T>& ref_Self;
	typedef const _WANode_<T>& const_ref_Self;

protected:
	Point m_point;
	int m_flag;
	pWAList m_plist;

public:
	_WANode_() :
			m_point(), m_flag(-1), m_plist(nullptr) {
	}
	_WANode_(const Point& p, int f, pWAList pl) :
			m_point(p), m_flag(f), m_plist(pl) {
	}
	_WANode_(const_ref_Self self) :
			m_point(self.m_point), m_flag(self.m_flag), m_plist(self.m_plist) {
	}
	ref_Self operator=(const_ref_Self self){
		if(&self == this){
			return *this;
		}else{
			this->m_point = self.m_point;
			this->m_flag = self.m_flag;
			this->m_plist = self.m_plist;
			return *this;
		}
	}

	Point get_point() const {
		return m_point;
	}
	ref_Point point() {
		return m_point;
	}
	const_ref_Point point() const{
		return m_point;
	}
	void set_flag(int f){
		m_flag = f;
	}
	int get_flag() const{
		return m_flag;
	}

};

}

#endif
