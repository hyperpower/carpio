#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"
#include "adaptive.hpp"
#include "stencil.hpp"
#include "data.hpp"
#include "shape.hpp"
#include "boundary.hpp"

#include <cmath>

namespace carpio {
typedef Float Cvt;
typedef Float Vt;

typedef Cell_<Cvt, 2> Cell_2D;
typedef Cell_<Cvt, 3> Cell_3D;

typedef Data_<Float, 2> Data_2D;
typedef Data_<Float, 3> Data_3D;

typedef PData_<Float, Float, 2> PData_2D;
typedef PData_<Float, Float, 3> PData_3D;

typedef Node_<Float, Float, 2> Node_2D;
typedef Node_<Float, Float, 3> Node_3D;
typedef Node_2D* pNode_2D;
typedef Node_3D* pNode_3D;
typedef const Node_<Float, Float, 2>* const_pNode_2D;
typedef const Node_<Float, Float, 3>* const_pNode_3D;

typedef Grid_<Float, Float, 2> Grid_2D;
typedef Grid_<Float, Float, 3> Grid_3D;

typedef Ghost_<Float, Float, 2> Ghost_2D;
typedef Ghost_<Float, Float, 3> Ghost_3D;

typedef Stencil_<Float, Float, 2, 1> Stencil_2D1;
typedef Stencil_<Float, Float, 2, 2> Stencil_2D2;

template<typename COO_VALUE, typename VALUE, st DIM>
class Domain_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Domain_<COO_VALUE, VALUE, DIM> Self;
	typedef Domain_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Grid_<COO_VALUE, VALUE, DIM>& ref_Grid;
	typedef const Grid_<COO_VALUE, VALUE, DIM>& const_ref_Grid;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef Ghost_<COO_VALUE, VALUE, DIM> Ghost;
	typedef Ghost_<COO_VALUE, VALUE, DIM> *pGhost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM> * const_pGhost;
	typedef Ghost_<COO_VALUE, VALUE, DIM>& ref_Ghost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>& const_ref_Ghost;

	typedef Adaptive_<COO_VALUE, VALUE, DIM> Adaptive;
	typedef Adaptive_<COO_VALUE, VALUE, DIM> *pAdaptive;
	typedef BoundaryIndex_<COO_VALUE, VALUE> BoundaryIndex;
	typedef BoundaryIndex_<COO_VALUE, VALUE>* pBoundaryIndex;

	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode> *pFace;
	typedef Shape_<COO_VALUE, DIM> Shape;
	typedef Shape_<COO_VALUE, DIM>* pShape;

public:
	// data
	// 2D : the domain is bounded by a shape
	//      the shape doesn't have holes.
	pShape _pshape_bound;
	// Innner solid
	std::list<pShape> _l_pshape_inner;
	// Grid: the calculation region
	pGrid _pgrid;
	// Adaptive
	pAdaptive _padaptive;
	// BC index
	pBoundaryIndex _pbindex;
	// Ghost nodes
	pGhost _pghost;

protected:
	// Requirement
	// The inner solid should be all in the shape bound
	void _check_solids() {

	}
	void _new_grid(cvt UL) {
		// build grid ------------------
		Float max_x = _pshape_bound->max_x();
		Float max_y = _pshape_bound->max_y();
		Float min_x = _pshape_bound->min_x();
		Float min_y = _pshape_bound->min_y();
		st n_x = std::ceil((max_x - min_x) / UL);
		st n_y = std::ceil((max_y - min_y) / UL);
		_pgrid = new Grid(  //
				n_x, min_x, UL, //
				n_y, min_y, UL);
	}
	void _new_adaptive(st lmin = 1, st lmax = 1) {
		_padaptive = new Adaptive(_pgrid, lmin, lmax);
	}
	void _new_boundary_index() {
		_pbindex = new BoundaryIndex();
	}
	void _new_ghost() {
		_pghost = new Ghost(_pgrid);
	}

public:
	Domain_(pShape bound, cvt unit_length) :
			_pshape_bound(bound) {
		_new_grid(unit_length);
		_new_adaptive();
		_new_boundary_index();
		_new_ghost();
	}

	Domain_(pShape bound, cvt unit_length, st minl, st maxl) :
			_pshape_bound(bound) {
		_new_grid(unit_length);
		_new_adaptive(minl, maxl);
		_new_boundary_index();
		_new_ghost();
	}

	~Domain_() {
		if (_pghost != nullptr) {
			delete _pghost;
			_pghost = nullptr;
		}
		if (_pbindex != nullptr) {
			delete _pbindex;
			_pbindex = nullptr;
		}
		if (_padaptive != nullptr) {
			delete _padaptive;
			_padaptive = nullptr;
		}
		if (_pgrid != nullptr) {
			delete _pgrid;
			_pgrid = nullptr;
		}
	}

	void build() {
		// adapt the solid
		_padaptive->adapt_bound_solid(*_pshape_bound, 1.0);
		for (typename std::list<pShape>::iterator iter =
				_l_pshape_inner.begin(); iter != _l_pshape_inner.end();
				++iter) {
			if ((*iter) != nullptr) {
				_padaptive->adapt_inner_solid(*(*iter));
			}
		}
		// connect nodes
		_pgrid->connect_root();
		_pgrid->connect_nodes();
		// build ghost
		_pghost->build();
		_pghost->connect();
		// set shape_idx and seg_idx in ghost node
		/*
		 * Method: calculates the ghost node belongs to which boundary segment.
		 * find the segment in shape which intersect with the line between ghost
		 * node and the boundary node
		 * choose nearest segment
		 * The line is axi align
		 * function : vt Intersect(Segment seg, Axis axis, vt ori)
		 */
		_pghost->set_boundary_index(0,*_pshape_bound);

	}


public:

	pGrid p_grid() {
		return _pgrid;
	}
	const_pGrid p_grid() const {
		return _pgrid;
	}
	ref_Grid grid() {
		return *_pgrid;
	}
	const_ref_Grid grid() const {
		return *_pgrid;
	}
	pGhost p_ghost() {
		return _pghost;
	}
	const_pGhost p_ghost() const {
		return _pghost;
	}
	ref_Ghost ghost() {
		return *_pghost;
	}
	const_ref_Ghost ghost() const {
		return *_pghost;
	}

};

}

#endif
