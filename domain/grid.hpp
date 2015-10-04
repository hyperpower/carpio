#ifndef GRID_H_
#define GRID_H_

#include "../TypeDef.h"
#include "grid_def.h"
#include "node.h"
#include "cell.h"

#include "../Algebra/Space.h"

namespace Larus {
namespace Grid {
template<typename COO_VALUE, typename VALUE, int DIM>
class Grid {
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const size_t NumNeighbors = NumFaces;

	typedef COO_VALUE coo_value_t;
	typedef VALUE value_t;
	typedef Grid<COO_VALUE, VALUE, DIM> self;
	typedef Grid<COO_VALUE, VALUE, DIM>* pself;
	typedef Cell<COO_VALUE, Dim> cell_t;
	typedef cell_t* pcell;
	typedef Data<VALUE, Dim> data_t;
	typedef data_t* pdata;
	typedef Node<COO_VALUE, VALUE, DIM> node;
	typedef Node<COO_VALUE, VALUE, DIM>* pnode;
	typedef void (*pfunction)(pnode, utPointer);
	typedef void (*pfunction_conditional)(arrayList&, pnode, utPointer);
	/*
	 *  data
	 */
	SpaceT<pnode, Dim> nodes;
	/*
	 *  constructor
	 */
	Grid() {

	}
	Grid(size_t ni, coo_value_t ox, coo_value_t dx, //
			size_t nj = 0, coo_value_t oy = 0, coo_value_t dy = 0, //
			size_t nk = 0, coo_value_t oz = 0, coo_value_t dz = 0) :
			nodes(ni, nj, nk){
		//nodes.recontruct(ni, nj, nk);
		size_t i_1d = 0;
		for (size_t i = 0; i < ni; i++) {
			for (size_t j = 0; j < ((Dim >= 2) ? nj : 1); j++) {
				for (size_t k = 0; k < ((Dim == 3) ? nk : 1); k++) {
					// new
					nodes.at_1d(i_1d) =  //
							new node(
							NULL_PTR, // father
									0, // type
									0, //level
									i_1d, //root idx
									0, //path
									ox + (i + 0.5) * dx, 0.5 * dx, //x
									oy + (j + 0.5) * dy, 0.5 * dy, //y
									oz + (k + 0.5) * dz, 0.5 * dz); //z
					i_1d ++;
				}
			}
		}
	}
protected:
	void _delete() {
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes.at_1d(i) != NULL_PTR) {
				delete nodes.at_1d(i);
			}
		}
	}
public:
	~Grid() {
		_delete();
	}
	/*
	 *  size
	 */
	inline size_t iLen() const {
		return nodes.iLen();
	}
	inline size_t jLen() const {
		return nodes.jLen();
	}
	inline size_t kLen() const {
		return (Dim < 3) ? 0 : nodes.kLen();
	}
	inline bool Empty() const {
		if (nodes.size() <= 0) {
			return true;
		} else {
			return false;
		}
	}
	size_t GetDim() const {
		return Dim;
	}
	size_t Size() const {
		return nodes.size();
	}

};
}
}
#endif
