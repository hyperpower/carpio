#ifndef _TEST_POLYNOMIAL_H_
#define _TEST_POLYNOMIAL_H_

#include "../algebra/algebra.hpp"
#include "../algebra/vector_list.hpp"
#include "../utility/Clock.h"
#include "../utility/random.h"
//#include "../calculation/expression.hpp"
#include "gtest/gtest.h"

namespace carpio {
TEST(DISABLED_test_polynomial,simple1_DISABLED) {
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

TEST(DISABLED_test_polynomial,const1_DISABLED) {
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

TEST(DISABLED_test_polynomial,plus_ploy_DISABLED) {
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
	std::cout << "concise ------ \n" << std::endl;
	poly.concise();
	poly.show();
}

TEST(DISABLED_test_polynomial,plus_ploy2) {
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
	std::cout << "concise ------ \n" << std::endl;
	poly.concise();
	poly.show();
}

TEST(DISABLED_test_polynomial,vector_list_base) {
	Clock t;
	VectorList_<int> vl(10);
	std::cout << " size :" << vl.size() << std::endl;
	std::cout << " cap  :" << vl.capacity() << std::endl;
	//vl[1] = 5;
	vl.push_back(3);
	for (VectorList_<int>::iterator iter = vl.begin(); iter != vl.end();
			++iter) {
		std::cout << *iter << std::endl;
	}

}

Polynomial_<Float, int, int> get_a_ploy() {
	typedef Polynomial_<Float, int, int> Poly;
	Poly res;
	Poly::Term t1(1.0, 10, 1);
	Poly::Term t2(2.0, 11, 0);
	res.insert(t1);
	res.insert(t2);

	return res;
}

TEST(test_polynomial,move_contructor) {
	typedef Polynomial_<Float, int, int> Poly;
	Poly p = get_a_ploy();
	p.show();
}

void poly_add() {
	Polynomial2_<std::string, Float> poly("z");
	poly.insert(1.0, "a");
	poly.insert(2.0, "b");
	poly.insert(3.0, "a");
	poly.show();

	Polynomial2_<std::string, Float> poly2("z");
	poly2.insert(1.0, "a");
	poly2.insert(2.0, "b");
	poly2.insert(3.0, "c");
	poly2.show();

	poly.plus(poly2);
	poly.times(0);
	poly.concise();
	poly.show();
}


TEST(test_polynomial,random_number) {
	Random::randomSeed();
	for (int i = 0; i < 10; i++) {
		std::cout << Random::nextFloat() << "\n";
		std::cout << Random::nextString(50) << "\n";
	}

	//Polynomial2_<std::string, Float> poly("z");
	//Expression2_<Float, Float, 3> exp;
	poly_add();

}

}
#endif
