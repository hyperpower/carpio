
#ifndef GEOMETRY_DEFINE_H_
#define GEOMETRY_DEFINE_H_

#include "../carpio_define.hpp"

#include <math.h>

namespace carpio {

template <class TYPE>
inline TYPE Distance(
		const TYPE &x1,
		const TYPE &y1,
		const TYPE &x2,
		const TYPE &y2){
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}


}
#endif
