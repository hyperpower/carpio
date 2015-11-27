#ifndef IO_DEFINE_H_
#define IO_DEFINE_H_

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "../carpio_define.hpp"
#include <string>
#include <fstream>

namespace carpio {
/*
 *  out put a and be to string,
 *  class V1 and v2 must overload operator<<
 *  sep is the separator
 *  For example:
 *  a = 1.3
 *  b = 1.4
 *  sep = " "
 *  return  1.3 1.4
 *
 */
template<class V1, class V2>
std::string ToString(V1 a, V2 b, const std::string sep) {
	std::ostringstream sst;
	sst << a << sep << b;
	return sst.str();
}

}

#endif /* IO_H_ */
