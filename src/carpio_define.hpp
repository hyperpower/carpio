#ifndef _DEFINE_HPP_
#define _DEFINE_HPP_

#include <array>
#include <deque>
#include <forward_list>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <assert.h>


namespace carpio{

#define ASSERT(expr) assert(expr)
#define ASSERT_MSG(expr, msg) assert((expr)&&(msg))
#define CAST(type, p)          ((type)p)
#define CAST_REF(type, p)      (*((type)p))
#define _IF_TRUE_RETRUN(expr)  if(expr){return;};
#define _IF_FALSE_RETRUN(expr)  if(false==(expr)){return;};
// value type
typedef std::size_t st; //size type
typedef double Float;
typedef void* utPointer;
typedef const void* const_utPointer;


// std  container
template<typename T>
using Vector = std::vector<T>;

//return code
#define _SUCCESS   0
#define _ERROR     1
#define _WARNING   2


}


#endif
