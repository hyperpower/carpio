#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "path.hpp"
#include <functional>

namespace carpio {

template<typename COO_VALUE, typename VALUE, st DIM>
struct GhostID_ {
	//typedef int (*pFun_set_bc)(Node_<COO_VALUE, VALUE, DIM>*,
	//		GhostID_<COO_VALUE, VALUE, DIM>&, utPointer);

	st root_idx;   //the root idx of the origin node
	//st idx;        //the local idx of the origin node
	Path_<DIM> path; //the path of the origin node
	st step;       //the steps of ghost node, we can choose multiple ghost Node,
				   // usually step = 0
	Direction direction; //The direction only on x, y or z
//--------------------------------------------------------------------
	//int bc_type;
	//pFun_set_bc pfun_bc;
	int shape_idx;
	int seg_idx;
};

template<typename COO_VALUE, typename VALUE, st DIM>
struct GhostID_compare_ {
	typedef GhostID_<COO_VALUE, VALUE, DIM> Gid;
	bool operator()(const Gid& lhs, const Gid& rhs) const {
		if (lhs.root_idx < rhs.root_idx) {
			return true;
		} else if (lhs.root_idx == rhs.root_idx) {
			if (lhs.path < rhs.path) {
				return true;
			} else if (lhs.path == rhs.path) {
				if (int(lhs.direction) < int(rhs.direction)) {
					return true;
				} else if (int(lhs.direction) == int(rhs.direction)) {
					return lhs.step < rhs.step;
				}
			}
		}
		return false;
	}
};

template<typename COO_VALUE, typename VALUE, int DIM>
class Ghost_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Ghost_<COO_VALUE, VALUE, DIM> Self;
	typedef const Ghost_<COO_VALUE, VALUE, DIM> const_Self;
	typedef Ghost_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>* const_pSelf;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell* pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data* pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM>* pNode;
	typedef void (*pfunction)(pNode, utPointer);
	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef void (*pfunction_conditional)(arrayList&, pNode, utPointer);
	typedef GhostID_<COO_VALUE, VALUE, DIM> GhostID;
	typedef GhostID_compare_<COO_VALUE, VALUE, DIM> GhostID_compare;

	typedef std::pair<GhostID, pNode> GhostNode;

	typedef std::map<GhostID, pNode, GhostID_compare> GhostMap;
	typedef typename GhostMap::iterator iterator;
	typedef typename GhostMap::const_iterator const_iterator;

protected:
	GhostMap _ghostmap;
	pGrid _pgrid;
	//
	/*
	 *  new ghost node
	 */
	pNode new_ghost_node(pNode pn, Direction dir) {
		int ghost_node_type = _Ghost_;
		//direction 4 5 6 7
		Cell c(*(pn->cell));
		switch (dir) {
		case _XM_:
			c.transfer(-c.get_d(_X_), 0.0, 0.0);
			break;
		case _YM_:
			c.transfer(0.0, -c.get_d(_Y_), 0.0);
			break;
		case _XP_:
			c.transfer(c.get_d(_X_), 0.0, 0.0);
			break;
		case _YP_:
			c.transfer(0.0, c.get_d(_Y_), 0.0);
			break;
		case _ZM_:
			c.transfer(0.0, 0.0, -c.get_d(_Z_));
			break;
		case _ZP_:
			c.transfer(0.0, 0.0, c.get_d(_Z_));
			break;
		default:
			ASSERT(false);
		}
		//pNode f, int nt, st level, st root_idx, st path,const Cell& c
		pNode ghostnode = new Node(  //
				nullptr,//
				ghost_node_type,  //
				pn->get_level(),  //
				pn->get_root_idx(),  //
				pn->get_idx(),       //
				pn->get_path(),      //
				c);
		ghostnode->father = pn;
		if (pn->data != nullptr) {
			ghostnode->data = new Data(*(pn->data));
		}
		return ghostnode;
	}
	/*
	 *  new ghost nodes outside of the boundary of Grid
	 */
	void new_ghost_nodes(pGrid pg) {
		for (typename Grid::iterator_face iter = pg->begin_face();
				iter != pg->end_face(); ++iter) {
			if (iter->ft() == _Boundary_) {
				pNode po = (*iter).pori();
				GhostID gid;
				gid.root_idx = po->get_root_idx(); //the root idx of the origin node
				gid.path = po->get_path(); //the local idx of the origin node
				gid.step = 0; //the steps of ghost node, we can choose multiple ghost Node,
							  // usually step =
				gid.direction = iter->dir(); //The direction only on x, y or z
				// --------------------------------
				//gid.bc_type = 0;
				//gid.pfun_bc = nullptr;
				gid.shape_idx = -1;
				gid.seg_idx = -1;
				pNode pghost = new_ghost_node(po, iter->dir());
				//change the index id ------------------------
				_ghostmap.insert(GhostNode(gid, pghost));
			}
		}
	}

public:
	Ghost_(pGrid pg) :
			_pgrid(pg) {
	}
	~Ghost_() {
		for (typename GhostMap::iterator it = _ghostmap.begin();
				it != _ghostmap.end(); ++it) {
			if (it->second != nullptr) {
				delete it->second;
				it->second = nullptr;
			}
		}
	}
	void build() {
		ASSERT(_ghostmap.empty());
		new_ghost_nodes(_pgrid);
	}
protected:
	void _connect(const GhostNode& gn){
		pNode po =  gn.second->father;
		pNode pg = gn.second;
		Direction dir_o= gn.first.direction;
		po->set_neighbor(pg,dir_o);
	}
public:
	void connect() {
		// when the ghost is built, the ghost node are already connect to original node.
		// the father of ghost node is the original one,
		// so the find_neighbor_fast() can not be used on ghost node.
		// here, we connect original node to the ghost node
		// the find_neighbor_fast() can be used on boundary node
		for (typename GhostMap::iterator it = _ghostmap.begin();
						it != _ghostmap.end(); ++it) {
			this->_connect(*it);
		}
	}

	iterator begin() {
		return _ghostmap.begin();
	}
	const_iterator begin() const {
		return _ghostmap.begin();
	}
	iterator end() {
		return _ghostmap.end();
	}
	const_iterator end() const {
		return _ghostmap.end();
	}

};

template<typename COO_VALUE, typename VALUE>
class BoundaryCondition_ {
public:
	typedef VALUE vt;
	typedef COO_VALUE cvt;
	typedef std::function<vt(cvt, cvt, cvt)> Fun;
	typedef BoundaryCondition_<cvt, vt> Self;
protected:
	// data
	// 1 Bondary conditon type
	// 2 function
	int _type;
	Fun _function;
public:
	// Constructor

	BoundaryCondition_() {
		// default boundary condition is symmetric boundary condition
		_type = 2;
		_function = [](cvt x, cvt y, cvt z) {return 0;};
	}
	BoundaryCondition_(int type, Fun fun) :
			_type(type), _function(fun) {
	}
	BoundaryCondition_(const Self& self) :
			_type(self._type), _function(self._function) {
	}
	// get
	int get_type() const {
		return _type;
	}
	vt get_val(cvt x, cvt y, cvt z) const {
		return _function(x, y, z);
	}
	// set
	void set_function(Fun fun) {
		_function = fun;
	}
	void set_default_1_bc(const vt& val) {
		_type = 1;
		_function = [&val](cvt x, cvt y, cvt z) {return val;};
	}
	void set_default_2_bc(const vt& val) {
		_type = 2;
		_function = [&val](cvt x, cvt y, cvt z) {return val;};
	}
};

template<typename COO_VALUE, typename VALUE>
class BoundaryIndex_ {
protected:
	struct BCID {
		int shape_idx;
		int seg_idx;
		st val_idx;
	};

	struct BCID_compare {
		typedef BCID BCid;
		bool operator()(const BCid& lhs, const BCid& rhs) const {
			if (lhs.shape_idx < rhs.shape_idx) {
				return true;
			} else if (lhs.shape_idx == rhs.shape_idx) {
				return lhs.seg_idx < rhs.seg_idx;
			} else if (lhs.seg_idx == rhs.seg_idx) {
				return int(lhs.val_idx) < int(rhs.val_idx);
			} else {
				return false;
			}
		}
	};
public:
	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef BoundaryCondition_<cvt, vt> BoundaryCondition;
	typedef BoundaryCondition_<cvt, vt>* pBoundaryCondition;
	//typedef BCID_compare_<cvt,vt> BCID_compare;

	typedef std::map<BCID, pBoundaryCondition, BCID_compare> BCMap;
protected:
	BCMap _BCmap;
	static const BoundaryCondition default_BC;

	typedef typename BCMap::iterator iterator;
	typedef typename BCMap::const_iterator const_iterator;

public:
	//constructor
	BoundaryIndex_() :
			_BCmap() {

	}
	//
	pBoundaryCondition find(int si, int segi, st vali) {
		BCID key;
		key.seg_idx = segi;
		key.shape_idx = si;
		key.val_idx = vali;
		iterator it = _BCmap.find(key);
		if (it != _BCmap.end()) {
			// found
			return (*it);
		} else {
			// not found
			return default_BC;
		}
	}
};

}

#endif
