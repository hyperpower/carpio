#ifndef _TEST_PATH_HPP_
#define _TEST_PATH_HPP_

#include "../domain/path.hpp"
#include "../carpio_define.hpp"

namespace carpio {
void test_path_1(){
	Path_<3> path(3);
	path.set();
	std::cout<<"size  :  "<<path.size()<<std::endl;

}

}

#endif
