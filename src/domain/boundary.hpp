#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"

namespace carpio {

template<typename COO_VALUE, typename VALUE, st DIM>
struct GhostID_ {
	typedef int (*pFun_set_bc)(Node_<COO_VALUE, VALUE, DIM>*,
			GhostID_<COO_VALUE, VALUE, DIM>&, utPointer);

	st root_idx;   //the root idx of the origin node
	st idx;        //the local idx of the origin node
	st step;       //the steps of ghost node, we can choose multiple ghost Node,
				   // usually step = 0
	Direction direction; //The direction only on x, y or z

	int bc_type;
	pFun_set_bc pfun_bc;
};

template<typename COO_VALUE, typename VALUE, st DIM>
struct GhostID_compare_ {
	typedef GhostID_<COO_VALUE, VALUE, DIM> Gid;
	bool operator()(const Gid& lhs, const Gid& rhs) const {
		if (lhs.root_idx < rhs.root_idx) {
			return true;
		} else if (lhs.root_idx == rhs.root_idx) {
			return lhs.idx < rhs.idx;
		} else if (lhs.idx == rhs.idx) {
			return int(lhs.direction) < int(rhs.direction);
		} else if (lhs.direction == rhs.direction) {
			return lhs.step < rhs.step;
		} else {
			return false;
		}
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

protected:
	typedef std::map<GhostID, pNode, GhostID_compare> GhostMap;
	GhostMap _ghostmap;
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
		pNode ghostnode = new Node(nullptr, ghost_node_type,
				pn->get_level(), pn->get_root_idx(), pn->get_idx(), c);
		ghostnode->father = pn;
		if(pn->data!=nullptr){
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
				gid.idx = po->get_idx(); //the local idx of the origin node
				gid.step = 0; //the steps of ghost node, we can choose multiple ghost Node,
							  // usually step = 0
				gid.direction = iter->dir(); //The direction only on x, y or z
				gid.bc_type = 0;
				gid.pfun_bc = nullptr;

				pNode pghost = new_ghost_node(po, iter->dir());
				//change the index id ------------------------
				//pn->data->aCenterData[Idx_IDX] = -gid.node_idx - 1; //negative
			 	_ghostmap.insert(GhostNode(gid, pghost));
			}
		}
	}

public:
	Ghost_(pGrid pg) {
		new_ghost_nodes(pg);
	}

};

template<typename VALUE>
class BoundaryCondition_ {
public:
	typedef VALUE vt;
	// Vbc
	// 1 Bondary conditon type
	// 2 value
	typedef std::pair<st, VALUE> Vbc; //variable boundary condition
	ArrayListT<Vbc> _abc;
	// ArrayListV<pFun> _bc_fun;
public:
	// Constructor
	// The number of the variable
	BoundaryCondition_(st n):_abc(n){

	}
	// Set
	/*
	 * Set boundary condtion 1
	 *    i  :  index of the variable
	 *    val:  the value
	 */
	void set_bc_1(st i, vt val){

	}
	void set_bc_2(st i, vt val){

	}
};

}

#endif
