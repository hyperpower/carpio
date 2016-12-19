#ifndef __TEST_TS_H_
#define __TEST_TS_H_

#include "../io/gnuplot.h"
#include "../ts/ts_tri_moller.h"

#include "gtest/gtest.h"
#include <math.h>



namespace carpio {

TEST(TS, test1) {
	std::cout<<"TS Test \n";
}

TEST(TS, test_oritation){

		Float v0[3];
		Float v1[3];
		Float v2[3];

		Float d[3];

		std::cout << " test oritation--------------------------- \n";
		v0[0] = 0;
		v0[1] = 0;
		v0[2] = 0;
		v1[0] = 1;
		v1[1] = 0;
		v1[2] = 0;
		v2[0] = 0;
		v2[1] = 1;
		v2[2] = 0;

		d[0] = 0;
		d[1] = 0;
		d[2] = 1;

		Float* a = v2;
		Float* b = v0;
		Float* c = v1;
		std::cout << SIGN3(ORIENT3D(d, b, c, a)) << std::endl;

}

}

#endif
