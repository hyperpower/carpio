#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_


#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"
#include "adaptive.hpp"


namespace carpio {

typedef Cell_<Float, 2> Cell_2D;
typedef Cell_<Float, 3> Cell_3D;

typedef Node_<Float, Float, 2> Node_2D;
typedef Node_<Float, Float, 3> Node_3D;

typedef Grid_<Float, Float, 2> Grid_2D;
typedef Grid_<Float, Float, 3> Grid_3D;


}

#endif
