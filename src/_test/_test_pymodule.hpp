#ifndef _TEST_PYMODULE_HPP_
#define _TEST_PYMODULE_HPP_

#include "../run/pymodule.h"

using namespace std;
namespace carpio {

void test_py_1() {
	cout << " ----------------- \n";
	PyParse pyp("/home/zhou/workspace/carpio/Debug/_para.py");
	cout << " ----------------- \n";
}

}

#endif

