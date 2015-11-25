#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
struct GhostID {
	typedef int (*pFun_set_bc)(Node_<COO_VALUE, VALUE, DIM>*,
				GhostID<COO_VALUE, VALUE, DIM>&, utPointer);

	st root_idx;   //the global idx of the origin node
	st idx;
	st step;       //the steps of ghost node, we can choose multiple ghost Node,
				   // usually step = 0
	Direction direction; //The direction only on x, y or z

	int bc_type;
	pFun_set_bc pfun_bc;
};

template<typename COO_VALUE, typename VALUE, int DIM>
struct GhostID_compare {
	typedef GhostID<COO_VALUE, VALUE, DIM> Gid;
	bool operator()(const Gid& lhs,
			const Gid& rhs) const {
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
class Boundary_ {
public:
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Boundary_<COO_VALUE, VALUE, DIM> Self;
	typedef const Boundary_<COO_VALUE, VALUE, DIM> const_Self;
	typedef Boundary_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef const Boundary_<COO_VALUE, VALUE, DIM>* const_pself;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell* pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data* pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM>* pNode;
	typedef void (*pfunction)(pNode, utPointer);
	typedef void (*pfunction_conditional)(arrayList&, pNode, utPointer);



};

}

#endif
