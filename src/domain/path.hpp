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
	typedef bit_index_t bi;
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
		bit_array_clear_all(_bar);
	}
	Path_(st n) {
		_create(n);
		bit_array_clear_all(_bar);
	}
	Path_(const Path_& p) {
		_create(0);
		bit_array_copy_all(_bar, p._bar);
	}
	Path_& operator=(const Path_& p) {
		if (this == &p) {
			return *this;
		}
		bit_array_copy_all(_bar, p._bar);
		return *this;
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
	char get(st idx) const {
		return bit_array_get_bit(_bar, bit_index_t(idx));
	}
	// is Empty
	bool empty() const {
		if (_bar == nullptr || this->size() == 0) {
			return true;
		} else {
			return false;
		}
	}
	bool is_valid() const {
		// 1 if this is empty -> not valid
		if (empty()) {
			return false;
		}
		// 2 we assume that the root node path equals to 0
		if (this->size() == 1) {
			if (this->get(0) == 0) {
				return true;
			} else {
				return false;
			}
		}
		// 3 the size of the normal node path should be divided by DIM
		return 0 == (this->size() % Dim);
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
	void toggle() {
		bit_array_toggle_all(_bar);
	}
	void toggle(st idx) {
		bit_array_toggle_bit(_bar, bit_index_t(idx));
	}
	void resize(st new_size) {
		char ret = bit_array_resize(_bar, bi(new_size));
		ASSERT_MSG(ret == 1, " >! resize failure.");
	}
	// copy
	// this function will resize this
	void copy(const Path_& src, st bidx, st len) {
		bit_array_copy(_bar, 0, src._bar, bi(bidx), bi(len));
	}
	void append(const Path_& app) {
		st os = this->size();
		st ns = os + app.size();
		this->resize(ns);
		bit_array_copy(_bar, bi(os), app._bar, 0, app.size());
	}
	void show(int config = 0) const {
		std::cout << "(";
		if (config != 0) {
			for (st i = 0; i < this->size(); i++) {
				if (i % Dim == 0) {
					std::cout << " ";
				}
				char c = this->get(i);
				if (c == 0) {
					std::cout << "0";
				} else if (c == 1) {
					std::cout << "1";
				} else {
					std::cout << "x";
				}

			}
			std::cout << ")";
		} else {
			std::cout << "Dim = " << Dim << ", size = " << this->size()
					<< ", (";
			for (st i = 0; i < this->size(); i++) {
				if (i % Dim == 0) {
					std::cout << " ";
				}
				char c = this->get(i);
				if (c == 0) {
					std::cout << "0";
				} else if (c == 1) {
					std::cout << "1";
				} else {
					std::cout << "x";
				}
			}
			std::cout << ")\n";
		}

	}

};

}

#endif
