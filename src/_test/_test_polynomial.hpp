#ifndef _TEST_POLYNOMIAL_H_
#define _TEST_POLYNOMIAL_H_

#include "../algebra/algebra.hpp"
#include "gtest/gtest.h"

namespace carpio {
TEST(test_polynomial,simple1) {
	typedef Polynomial_<Float, std::string, int> Poly;
	Poly poly;
	Poly::Term t1(1.0, "b", 0);
	Poly::Term t2(4.0, "a", 0);
	Poly::Term t5(4.0, "c", 0);
	Poly::Term t3(5.0, "a", 2);
	Poly::Term t4(6.0, "a", 3);
	poly.insert(t1);
	poly.insert(t2);
	poly.insert(t3);
	poly.insert(t4);
	poly.insert(t5);
	//poly.concise();
	EXPECT_EQ(3, poly.size());
	poly.show();
}

TEST(test_polynomial,const1) {
	typedef Polynomial_<Float, std::string, int> Poly;
	Poly poly;
	Poly::Term t1(3.0, "a", 1);
	Poly::Term t2(4.0, "a", 1);
	Poly::Term t3(5.0, "a", 0);
	Poly::Term t4(6.0, "a", 0);
	poly.insert(t1);
	poly.insert(t2);
	poly.insert(t3);
	poly.insert(t4);
	//EXPECT_EQ(3, poly.size());
	poly.show();
}

TEST(test_polynomial,plus_ploy) {
	typedef Polynomial_<Float, std::string, int> Poly;
	Poly poly;
	Poly::Term t1(1.0, "a", 1);
	Poly::Term t2(2.0, "b", 1);
	Poly::Term t3(3.0, "d", 0);
	poly.insert(t1);
	poly.insert(t2);
	poly.insert(t3);
	Poly poly2;
	Poly::Term t12(1.0, "a", 1);
	Poly::Term t22(2.0, "b", 1);
	Poly::Term t32(2.0, "b", 0);
	poly2.insert(t12);
	poly2.insert(t22);
	poly2.insert(t32);

	//EXPECT_EQ(3, poly.size());
	poly.show();
	poly2.show();
	poly.plus(poly2);
	poly.show();
	std::cout<<"concise ------ \n"<<std::endl;
	poly.concise();
	poly.show();
}

TEST(test_polynomial,plus_ploy2) {
	typedef Polynomial_<Float, int, int> Poly;
	Poly poly;
	Poly::Term t1(1.0, 10, 1);
	Poly::Term t2(2.0, 11, 0);
	poly.insert(t1);
	poly.insert(t2);
	Poly poly2;
	Poly::Term t12(1.0, 13, 1);
	Poly::Term t22(2.0, 14, 0);
	poly2.insert(t12);
	poly2.insert(t22);

	//EXPECT_EQ(3, poly.size());
	poly.show();
	poly2.show();
	poly.plus(poly2);
	poly.show();
	std::cout<<"concise ------ \n"<<std::endl;
	poly.concise();
	poly.show();
}
}
#endif
