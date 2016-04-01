#ifndef _POLYGON_BOOLEAN_HPP_
#define _POLYGON_BOOLEAN_HPP_

#include "../carpio_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include "_polygon.hpp"
#include "_relation.hpp"
#include "../algebra/array_list.hpp"
#include <list>
#include <iterator>
#include <sstream>

#include "../io/gnuplot.h"

namespace carpio {

enum WAType {
	WA_ERR = 0x100000, WA_ORIGNAL = 0x10000,

	WA_IN = 0x1, WA_OUT = 0x2, WA_USED_IN = 0x4, WA_USED_OUT = 0x8,

	WA_IN_R = 0x10, WA_OUT_R = 0x20, WA_USED_IN_R = 0x40, WA_USED_OUT_R = 0x80,

	WA_HEAD_C = 0x100, WA_TAIL_C = 0x200, WA_CLIP = 0x400, WA_SUBJ = 0x800,

	WA_HEAD_S = 0x1000, WA_TAIL_S = 0x2000,
};

inline bool HasType(int type, WAType t) {
	return (type | t) == type ? true : false;
}

inline bool HasType(int type, int typecombine) {
	return (type | typecombine) == type ? true : false;
}
inline std::string ParseWAType(int type) {
	std::ostringstream osstr;
	int flag = 0;
	if (HasType(type, WA_ERR)) {
		osstr << " WA_ERR ";
		flag = 1;
	}
	if (HasType(type, WA_ORIGNAL)) {
		osstr << " WA_ORIGNAL ";
		flag = 1;
	}
	if (HasType(type, WA_IN)) {
		osstr << " WA_IN ";
		flag = 1;
	}
	if (HasType(type, WA_OUT)) {
		osstr << " WA_OUT ";
		flag = 1;
	}
	if (HasType(type, WA_USED_IN)) {
		osstr << " WA_USED_IN ";
		flag = 1;
	}
	if (HasType(type, WA_USED_OUT)) {
		osstr << " WA_USED_OUT ";
		flag = 1;
	}
	if (HasType(type, WA_IN_R)) {
		osstr << " WA_IN_R ";
		flag = 1;
	}
	if (HasType(type, WA_OUT_R)) {
		osstr << " WA_OUT_R ";
		flag = 1;
	}
	if (HasType(type, WA_USED_IN_R)) {
		osstr << " WA_USED_IN_R ";
		flag = 1;
	}
	if (HasType(type, WA_USED_OUT_R)) {
		osstr << " WA_USED_OUT_R ";
		flag = 1;
	}
	if (HasType(type, WA_HEAD_C)) {
		osstr << " WA_HEAD_C ";
		flag = 1;
	}
	if (HasType(type, WA_TAIL_C)) {
		osstr << " WA_TAIL_C ";
		flag = 1;
	}
	if (HasType(type, WA_HEAD_S)) {
		osstr << " WA_HEAD_S ";
		flag = 1;
	}
	if (HasType(type, WA_TAIL_S)) {
		osstr << " WA_TAIL_S ";
		flag = 1;
	}
	if (HasType(type, WA_CLIP)) {
		osstr << " WA_CLIP ";
		flag = 1;
	}
	if (HasType(type, WA_SUBJ)) {
		osstr << " WA_SUBJ ";
		flag = 1;
	}
	if (flag == 0) {
		osstr << " Error input ";
	}
//	osstr << "\n";
	return osstr.str();
}
inline void ShowWAType(int type) {
	std::cout << ParseWAType(type) << "\n";
}

template<typename T>
class _WANode_ {
public:
	typedef Point_<T, 2> Point;
	typedef Point_<T, 2>& ref_Point;
	typedef const Point_<T, 2>& const_ref_Point;

	typedef std::list<_WANode_<T> > List;
	typedef std::list<_WANode_<T> >* pList;

	typedef typename std::list<_WANode_<T> >::iterator iterator;
	typedef typename std::list<_WANode_<T> >::const_iterator cosnt_iterator;

	typedef _WANode_<T> Self;
	typedef _WANode_<T>& ref_Self;
	typedef const _WANode_<T>& const_ref_Self;

protected:
	Point m_point;
	int m_flag;
	iterator m_iter;

public:
	_WANode_() :
			m_point(), m_flag(-1), m_iter() {
	}
	_WANode_(const Point& p, int f, iterator pl) :
			m_point(p), m_flag(f), m_iter(pl) {
	}
	_WANode_(const_ref_Self self) :
			m_point(self.m_point), m_flag(self.m_flag), m_iter(self.m_iter) {
	}
	ref_Self operator=(const_ref_Self self) {
		if (&self == this) {
			return *this;
		} else {
			this->m_point = self.m_point;
			this->m_flag = self.m_flag;
			this->m_iter = self.m_iter;
			return *this;
		}
	}
	bool operator==(const_ref_Self self) const {
		if (this->m_flag == self.m_flag && this->m_point == self.m_point
				&& this->m_iter == self.m_iter) {
			return true;
		} else {
			return false;
		}
	}

	Point get_point() const {
		return m_point;
	}
	ref_Point point() {
		return m_point;
	}
	const_ref_Point point() const {
		return m_point;
	}
	void set_flag(int f) {
		m_flag = f;
	}
	int get_flag() const {
		return m_flag;
	}
	bool has_flag(int att) const {
		return HasType(m_flag, att);
	}
	/*
	 * substract some attribute from flag
	 */
	void substract_flag(int att) {
		m_flag = m_flag & ~att;
	}
	void add_flag(int att) {
		m_flag |= att;
	}
	void set_iter(iterator iter) {
		m_iter = iter;
	}
	iterator get_iter() {
		return m_iter;
	}
};

template<typename T>
class WAClipping_ {
	typedef Point_<T, 2> Point;
	typedef Point_<T, 2>& ref_Point;
	typedef const Point_<T, 2>& const_ref_Point;
	typedef Polygon_<T> Polygon;
	typedef Polygon_<T>& ref_Polygon;
	typedef const Polygon_<T>&const_ref_Polygon;
	typedef typename Polygon_<T>::Segment Segment;
	typedef _WANode_<T> Node;
	typedef std::list<_WANode_<T> > List;
	typedef typename List::iterator iterator;
	typedef typename List::const_iterator const_iterator;

	typedef std::list<Polygon> List_Polygon;

protected:
	List l_cli;
	List l_sub;
public:
	WAClipping_(const Polygon& clip, const Polygon& sub) :
			l_cli(), l_sub() {
		_copy_polygon_to_list(clip, l_cli);
		_copy_polygon_to_list(sub, l_sub);
	}

protected:
	int intersect_type(const Segment &c_seg, const Segment &s_seg) {
		if (!IsBoxCross(c_seg, s_seg)) {
			return WA_ERR; //no intersect (overlap is not intersect)
		} else {
			//clip segment
			int s12s = OnWhichSide3(c_seg, s_seg.ps());  //tail
			int s12e = OnWhichSide3(c_seg, s_seg.pe());  //head
			//subj segment
			int s21s = OnWhichSide3(s_seg, c_seg.ps());  //tail
			int s21e = OnWhichSide3(s_seg, c_seg.pe());  //head

			if (s12s != s12e && (s12s + s12e) != -2 && s21s != s21e
					&& (s21s + s21e) != -2) {
				//intersect
				int res = 0;
				if (s12s == 0)
					res = res | WA_TAIL_S | WA_SUBJ;
				if (s12e == 0)
					res = res | WA_HEAD_S | WA_SUBJ;
				if (s21s == 0)
					res = res | WA_TAIL_C | WA_CLIP;
				if (s21e == 0)
					res = res | WA_HEAD_C | WA_CLIP;
				if (s12e == 1 || s12s == -1)
					res = res | WA_IN;
				else
					res = res | WA_OUT;
				if (s21e == 1 || s21s == -1)
					res = res | WA_IN_R;
				else
					res = res | WA_OUT_R;

				return res;
			} else {
				return WA_ERR; //no intersect
			}
		}
	}
	int _copy_polygon_to_list(const Polygon& p, List& list) const {
		ASSERT(!p.empty());
		list.clear();
		for (st i = 0; i < p.size_vertexs(); i++) {
			Node node(p.v(i), WA_ORIGNAL, list.end());
			list.push_back(node);
		}
		return _SUCCESS;
	}
	int copy_list_to_polygon(const List& list, Polygon& p) const {
		ASSERT_MSG(!list.empty(), " >! Empty List");
		typename Polygon_<T>::ArrP arrp(list.size());
		st i = 0;
		const_iterator iter = list.begin();
		for (; iter != list.end();) {
			arrp[i] = iter->point();
			++iter;
			++i;
		}
		p.reconstruct(arrp);
		return _SUCCESS;
	}
	iterator _next(const List& list, const iterator& iter) const {
		// no stop, loop as a circle
		iterator res = std::next(iter);
		if (res != list.end()) {
			return res;
		} else {
			return ++res;
		}
	}
	iterator _next_ori(const List& list, const iterator& iter) const {
		iterator res = _next(list, iter);
		while (!res->has_flag(WA_ORIGNAL)) {
			res = _next(list, res);
		}
		return res;
	}
	iterator _next_ori_with_end(const List& list, const iterator& iter) const {
		iterator res = std::next(iter);
		while (!(res == list.end() || res->has_flag(WA_ORIGNAL))) {
			res = std::next(res);
		}
		return res;
	}
	Segment _to_segment(const iterator& iter_s, const iterator& iter_e) {
		return Segment(iter_s->point(), iter_e->point());
	}
public:
	int clipping(List_Polygon& resp) {
		// outter loop --------------------------
		iterator iter_sub = l_sub.begin();
		int io = 0;
		for (; iter_sub != l_sub.end(); io++) {
			iterator iter_sub_nxt = _next_ori(l_sub, iter_sub);
			Segment seg_sub = _to_segment(iter_sub, iter_sub_nxt);
			// inner loop -----------------------
			// clip
			int ino = 0;
			iterator iter_cli = l_cli.begin();
			for (; iter_cli != l_cli.end(); ino++) {
				std::cout << "out " << io << "in " << ino << "\n";
				iterator iter_cli_nxt = _next_ori(l_cli, iter_cli);
				Segment seg_cli = _to_segment(iter_cli, iter_cli_nxt);

				Point interp;

				int type = intersect_type(seg_cli, seg_sub);
				//seg_cli.show();
				//seg_sub.show();
				if (type != WA_ERR) { //intersect
					//ShowWAType(type);

					if (!HasType(type, WA_SUBJ) && HasType(type, WA_CLIP)) {
						// One point on clip
						if (HasType(type, WA_HEAD_C)) {
							iter_sub_nxt = l_sub.insert(iter_sub_nxt,
									Node(iter_cli_nxt->point(), type,
											iter_cli_nxt));
							iter_cli_nxt->add_flag(type);
							iter_cli_nxt->set_iter(iter_sub_nxt);
						} else {
							//iter_sub_nxt->add_flag(type);
							iter_cli->add_flag(type);
							//iter_cli->set_iter(iter_sub_nxt);
						}
					} else if (HasType(type, WA_SUBJ)
							&& !HasType(type, WA_CLIP)) {
						// One point on subject
						if (HasType(type, WA_HEAD_S)) {
							iter_cli_nxt = l_cli.insert(iter_cli_nxt,
									Node(iter_sub_nxt->point(), type,
											iter_sub_nxt));
							iter_sub_nxt->add_flag(type);
							iter_sub_nxt->set_iter(iter_cli_nxt);
						} else {
							//	iter_cli_nxt = l_sub.insert(iter_cli_nxt,
							//			Node(iter_sub->point(), type, iter_sub));
							iter_sub->add_flag(type);
							//	iter_sub->set_iter(iter_cli_nxt);
						}
					} else if (HasType(type, WA_SUBJ)
							&& HasType(type, WA_CLIP)) {
						// two point
						if (HasType(type, WA_HEAD_S)
								&& HasType(type, WA_HEAD_C)) {
							iter_sub_nxt->add_flag(type);
							iter_cli_nxt->add_flag(type);
							iter_sub_nxt->set_iter(iter_cli_nxt);
							iter_cli_nxt->set_iter(iter_sub_nxt);
						} else if (HasType(type, WA_HEAD_S)
								&& HasType(type, WA_TAIL_C)) {
							iter_sub_nxt->add_flag(type);
							iter_cli->add_flag(type);
							iter_sub_nxt->set_iter(iter_cli);
							iter_cli->set_iter(iter_sub_nxt);
						} else if (HasType(type, WA_TAIL_S)
								&& HasType(type, WA_HEAD_C)) {
							iter_sub->add_flag(type);
							iter_cli_nxt->add_flag(type);
							//	iter_sub->set_iter(iter_cli_nxt);
							//	iter_cli_nxt->set_iter(iter_sub);
						} else if (HasType(type, WA_TAIL_S)
								&& HasType(type, WA_TAIL_C)) {
							iter_sub->add_flag(type);
							iter_cli->add_flag(type);
							//	iter_sub->set_iter(iter_cli);
							//	iter_cli->set_iter(iter_sub);
						} else {
							SHOULD_NOT_REACH;
						}
					} else if (!HasType(type, WA_SUBJ)
							&& !HasType(type, WA_CLIP)) {
						// general case
						interp = CalIntersect(seg_cli, seg_sub);
						iter_cli_nxt = l_cli.insert(iter_cli_nxt,
								Node(interp, type, l_cli.end()));
						iter_sub_nxt = l_sub.insert(iter_sub_nxt,
								Node(interp, type, iter_cli_nxt));
						iter_cli_nxt->set_iter(iter_sub_nxt);
					}
				}
				// inner ++
				iter_cli = _next_ori_with_end(l_cli, iter_cli);
			}
			// outter ++
			iter_sub = _next_ori_with_end(l_sub, iter_sub);
		}
		std::cout << "  =============res ===========\n";
		List res;
		clipping_search(res);
		_show_list(res);
		if (!res.empty()) {
			Polygon poly;
			copy_list_to_polygon(res, poly);
			resp.push_back(poly);
		}
		std::cout << "  =============up res ========\n";
	}

	void _set_used(iterator& iter) {
		if (iter->has_flag(WA_IN)) {
			iter->add_flag(WA_USED_IN);
		}
		if (iter->has_flag(WA_OUT)) {
			iter->add_flag(WA_USED_OUT);
		}
	}

	void clipping_search(List& res) {
		res.clear();
		//Find first IN on subject
		iterator iter_sub = l_sub.begin();
		iterator iter_in = l_sub.end();
		for (; iter_sub != l_sub.end(); ++iter_sub) {
			if (iter_sub->has_flag(WA_IN)) {
				iter_in = iter_sub;
				break;
			}
		}
		if (iter_in == l_sub.end()) {
			return;
		}
		//
		List* plist = &l_sub;
		iterator move = iter_in;
		int count = 0;
		do {
			std::cout << " count = " << count << "\n";
			res.push_back(*move);  //record
			move->point().show();
			_set_used(move);
			ShowWAType(move->get_flag());
			iterator iter_other = iter_in->get_iter();
			int op = 0;
			if (op == 0 && plist == &l_sub
					&& (move->has_flag(WA_OUT)
							|| move->get_iter()->has_flag(WA_OUT))) {
				move = move->get_iter();
				plist = &l_cli;
				_set_used(move);
				op++;
			}
			if (op == 0 && plist == &l_cli && (move->has_flag(WA_USED_IN))) {
				move = move->get_iter();
				plist = &l_sub;
				if (move == iter_in) {
					break;
				}
				op++;
			}
			//advance
			move = _next(*plist, move);
			count++;
		} while (move != iter_in && count < 100);
		res.pop_back(); //delete the last one

	}
	void _show_list(const List& list) const {
		std::cout << "Size = " << list.size() << "\n";
		const_iterator iter = list.begin();
		int count = 0;
		for (; iter != list.end(); ++iter, ++count) {
			std::cout << count << " ";
			std::cout << iter->point().to_string();
			ShowWAType(iter->get_flag());
		}
	}
	void show(int i = 0) const {
		if (i == 0) {
			_show_list(l_cli);
		} else {
			_show_list(l_sub);
		}
	}

}
;

template<typename T>
int GnuplotActor_WAList(Gnuplot_actor& actor, std::list<_WANode_<T> >& list) {
	typedef typename std::list<_WANode_<T> >::iterator iterator;
	actor.clear();
	actor.command() = "using 1:2:3:4 title \"\"";
	actor.style() = "with vectors";
	if (list.empty()) {
		actor.data().push_back("");
		return _ERROR;
	}
	typedef typename _WANode_<T>::Point Poi;
	for (iterator iter = list.begin(); iter != list.end(); ++iter) {
		iterator iter_nxt = _GetNext(list, iter);
		const Poi& ps = iter->point();
		const Poi& pe = iter_nxt->point();
		std::stringstream sstr;
		sstr << ps.x() << " " << ps.y() << " " << pe.x() - ps.x() << " "
				<< pe.y() - ps.y();
		actor.data().push_back(sstr.str());
	}
	return _SUCCESS;
}

template<typename T>
int Intersect(const Polygon_<T>& clip, const Polygon_<T>& sub,
		std::list<Polygon_<T> >& res) {
	//typedef Polygon_<T> Polygon;
	typedef typename Polygon_<T>::Segment Segment;
	typedef typename Polygon_<T>::Point Point;
	typedef _WANode_<T> Node;
	typedef std::list<_WANode_<T> > List;
	typedef typename std::list<_WANode_<T> >::iterator iterator;
//1 Is clip inside of the subject
//  (clip is very small)
	if (IsIn(sub, clip)) {
		res.push_back(clip);
		//std::cout << " --- clip in sub\n";
		return _SUCCESS;
	}
//2 Is subject inside of the clip
//  (clip is very small)
	if (IsIn(clip, sub)) {
		res.push_back(sub);
		//std::cout << " --- sub in clip\n";
		return _SUCCESS;
	}
//3 Copy clip polygon to list
	List lc;
	//CopyPolygonToList(clip, lc);
//  Copy subject polygon to list
	List ls;
	//CopyPolygonToList(sub, ls);
//4 Copy and connect list
	List lcd;
	//CopyConnect(lc, lcd);
	List lsd;
	//CopyConnect(ls, lsd);

//5 Traveral subject
//  outter loop
//  subject
	//iterator iter_sub = ls.begin();
	//for (; iter_sub != ls.end(); ++iter_sub) {
	//iterator iter_subd = iter_sub->get_iter();
	//iterator iter_subdn = _GetNext(lsd, iter_subd);
	//Segment seg_sub = _GetSegment(ls, iter_sub);
	// inner loop
	// clip
	//iterator iter_cli = lc.begin();
	//for (; iter_cli != lc.end(); ++iter_cli) {
	//iterator iter_clid = iter_cli->get_iter();
	//iterator iter_clidn = _GetNext(lcd, iter_clid);
	//Segment seg_cli = _GetSegment(lc, iter_cli);
	//Point interp;

	//int type = IntersectWAType(seg_cli, seg_sub);
	//seg_cli.show();
	//seg_sub.show();
	//if (type != WA_ERR) { //intersect
	//ShowWAType(type);
	//general case
	//if (!HasType(type, WA_SUBJ) && !HasType(type, WA_CLIP)) {
	//	interp = CalIntersect(seg_cli, seg_sub);
	//	iter_clidn = lcd.insert(iter_clidn,
	//			Node(interp, type, lcd.end()));
	//	iter_subdn = lsd.insert(iter_subdn,
	//		Node(interp, type, iter_clidn));
	//	iter_clidn->set_iter(iter_subdn);
	//}
	//}
	//}
	//}
	// Find in point on sublist
	//List out;
	//iterator iter_subd = lsd.begin();
	//for (; iter_sub != lsd.end(); ++iter_subd) {
	//	if (HasType(iter_sub->get_flag(), WA_IN)
	//			&& !HasType(iter_sub->get_flag(), WA_IN | WA_OUT)) {
	//		//search(pin, out);
	//		Polygon polygon;
	//		CopyListToPolygon(out, polygon);
	//		res.push_back(polygon);
	//		out.clear();
	//		iter_subd = lsd.begin();
	//	}
	//}

	//show
	//std::list<Gnuplot_actor> lga;
	//Gnuplot_actor ga;
	//GnuplotActor_WAList(ga, lcd);
	//lga.push_back(ga);
	//GnuplotActor_WAList(ga, lsd);
	//lga.push_back(ga);
	//Gnuplot gp;
	//gp.set_equal_ratio();
	//
	//int wt = IntersectWAType(sc, so);
	//ShowWAType(wt);
	//gp.set_label(1, ParseWAType(wt), 0.5, 0.1);
	//gp.set_xrange(0, 5.1);
	//gp.set_yrange(0, 5.1);
	//gp.plot(lga);
}

}

#endif
