#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_


#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"
#include "adaptive.hpp"
#include "stencil.hpp"
#include "data.hpp"

namespace carpio {

typedef Cell_<Float, 2> Cell_2D;
typedef Cell_<Float, 3> Cell_3D;

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



}

#endif
