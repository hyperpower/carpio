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
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef Ghost_<COO_VALUE, VALUE, DIM> Ghost;
	typedef Ghost_<COO_VALUE, VALUE, DIM> *pGhost;
	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode> *pFace;
	typedef Shape_<COO_VALUE, DIM> Shape;
	typedef Shape_<COO_VALUE, DIM>* pShape;

public:
	// data
	// 2D : the domain is bounded by a shape
	//      the shape doesn't have holes.
	pShape pshape_bound;
	// Innner solid
	std::list<pShape> l_pshape_inner;
	// Grid: the calculation region
	pGrid pgrid;
	// Ghost nodes
	pGhost pghost;

protected:
	// Requirement
	// The inner solid should be all in the shape bound
	void _check_solids(){

	}
public:
	Domain_(pShape bound, std::list<pShape> l_ph, pGrid pg){

	}



};

}

#endif
