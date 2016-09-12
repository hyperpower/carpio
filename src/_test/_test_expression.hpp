#ifndef _TEST_EXPRESSION_H_
#define _TEST_EXPRESSION_H_

#include "../calculation/expression.hpp"
#include "gtest/gtest.h"

namespace carpio {

TEST(Expression, case1) {
	Expression2_<Float, Float, 3> exp;
	exp.insert(1,nullptr);
	exp.insert(1,nullptr);
	exp.show();

}



}

#endif
