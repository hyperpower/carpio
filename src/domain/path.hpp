#ifndef PATH_HPP_
#define PATH_HPP_

#include "../carpio_define.hpp"
#include "domain_define.hpp"
#include "../utility/bit_array.h"
#include <functional>

#include <math.h>
#include <inttypes.h>

namespace carpio {
template<st DIM>
class Path_ {
public:
	static const st Dim = DIM;
protected:
	// protected data
	BIT_ARRAY* _bar;

	// protected function
	void _create(st nbits) {
		_bar = bit_array_create(bit_index_t(nbits));
	}

public:
	// Constructor - create a new bit array of length nbits
	Path_() {
		_create(0);
	}
	Path_(st n) {
		_create(n);
	}
	// Destructor - free the memory used for a bit array
	~Path_() {
		bit_array_free(_bar);
	}
	// Size
	st size() const {
		return st(bit_array_length(_bar));
	}
	// get dim
	st get_dim() const {
		return Dim;
	}
	// is Empty
	bool empty() const {
		if (_bar == nullptr || this->size() == 0) {
			return true;
		} else {
			return false;
		}
	}
	// set if [idx] == 0 then set to 1,  if [idx] == 1 then no change
	void set() {
		bit_array_set_all(_bar);
	}
	void set(st idx) {
		bit_array_set_bit(_bar, bit_index_t(idx));
	}
	void clear() {
		bit_array_clear_all(_bar);
	}
	void clear(st idx) {
		bit_array_clear_bit(_bar, bit_index_t(idx));
	}
	void toggle(){
		bit_array_toggle_all(_bar);
	}
	void toggle(st idx){
		bit_array_toggle_bit(_bar, bit_index_t(idx));
	}


};

}

#endif
