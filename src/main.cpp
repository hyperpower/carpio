#include <iostream>
#include "algebra/polynomial.hpp"
#include "domain/node.hpp"
//include "domain/io_grid.hpp"
//#include "_test/_test_solver.hpp"
//#include "_test/_test_grid.hpp"
//#include "_test/_test_ns.hpp"
//#include "_test/_test_equation.hpp"
//#include "_test/test_advection.hpp"
//#include "_test/_test_interpolate.hpp"
//#include "_test/_test_poisson.hpp"
//#include "_test/_test_polynomial.hpp"
#include "_test/_test_expression.hpp"
#include <gtest/gtest.h>



using namespace std;
using namespace carpio;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int res = RUN_ALL_TESTS();
    return res;
}
