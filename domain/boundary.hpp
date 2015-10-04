#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
struct GhostID {
	typedef int (*pFun_set_bc)(Node<COO_VALUE, VALUE, DIM>*,
				GhostID<COO_VALUE, VALUE, DIM>&, utPointer);

	size_t root_idx;  //the global idx of the origin node
	size_t path;
	size_t step;   //the steps of ghost node, we can choose multiple ghost node,
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
			return lhs.path < rhs.path;
		} else if (lhs.path == rhs.path) {
			return int(lhs.direction) < int(rhs.direction);
		} else if (lhs.direction == rhs.direction) {
			return lhs.step < rhs.step;
		} else {
			return false;
		}
	}
};

template<typename COO_VALUE, typename VALUE, int DIM>
class Boundary {
public:
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const size_t NumNeighbors = NumFaces;

	typedef COO_VALUE coo_value_t;
	typedef VALUE value_t;
	typedef Boundary<COO_VALUE, VALUE, DIM> self;
	typedef const Boundary<COO_VALUE, VALUE, DIM> const_self;
	typedef Boundary<COO_VALUE, VALUE, DIM>* pself;
	typedef const Boundary<COO_VALUE, VALUE, DIM>* const_pself;
	typedef Cell<COO_VALUE, Dim> cell_t;
	typedef cell_t* pcell;
	typedef Data<VALUE, Dim> data_t;
	typedef data_t* pdata;
	typedef Node<COO_VALUE, VALUE, DIM> node;
	typedef Node<COO_VALUE, VALUE, DIM>* pnode;
	typedef void (*pfunction)(pnode, utPointer);
	typedef void (*pfunction_conditional)(arrayList&, pnode, utPointer);

};

}
}

#endif
