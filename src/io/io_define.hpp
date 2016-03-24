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

template<class V1, class V2, class V3>
std::string ToString(V1 a, V2 b, V3 c, const std::string sep) {
	std::ostringstream sst;
	sst << a << sep << b;
	return sst.str();
}

template<class V1, class V2, class V3, class V4, class V5, class V6, class V7>
std::string ToString(V1 a, V2 b, V3 c, V4 d, V5 e, V6 f, V7 g,
		const std::string sep) {
	std::ostringstream sst;
	sst << a << sep; //1
	sst << b << sep; //2
	sst << c << sep; //3
	sst << d << sep; //4
	sst << e << sep; //5
	sst << f << sep; //6
	sst << g;        //7
	return sst.str();
}

}

#endif /* IO_H_ */
