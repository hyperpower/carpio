/*
 * Expression.h
 *
 *  Created on: Feb 2, 2015
 *      Author: zhou
 */

#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "../carpio_define.hpp"
#include <utility>

#include <stdint.h>

template<class T>
struct IsZero_: std::unary_function<T, bool> {
	bool operator()(T number) const{
		return (number == 0);
	}
};

namespace carpio {
template<class COE, class TERM, class EXP,                                    //
		class ISZERO_COE = IsZero_<COE>,                                     //
		class ISZERO_EXP = IsZero_<EXP>,                                     //
		class COMPARE_TERM = std::less<TERM>,                                 //
		class COMPARE_EXP = std::less<EXP> >                                 //
class Polynomial {
protected:
	typedef std::pair<TERM, EXP> Key;
	typedef COE Value;
	typedef ISZERO_COE is_zero_coe;
	typedef ISZERO_EXP is_zero_exp;
	ISZERO_COE _is_zero_coe;
	ISZERO_EXP _is_zero_exp;

	struct _compare {
		COMPARE_TERM _compare_term;
		COMPARE_EXP _compare_exp;

		ISZERO_EXP __is_zero_exp;
		bool operator()(const Key& lhs, const Key& rhs) const {
			if(__is_zero_exp(lhs.second) && __is_zero_exp(rhs.second)){
				return false;
			}
			if (_compare_term(lhs.first, rhs.first)) {
				return true;
			} else if (!_compare_term(rhs.first, lhs.first)) {
				return _compare_exp(lhs.second, rhs.second);
			}
			return false;
		}
	};


	typedef std::map<Key, Value, _compare> Map;
	Map _map;
public:

	typedef typename Map::const_reference const_reference;
	typedef typename Map::pointer pointer;
	typedef typename Map::const_pointer const_pointer;
	typedef typename Map::iterator iterator;
	typedef typename Map::const_iterator const_iterator;

	typedef typename Map::difference_type difference_type;
	typedef typename Map::size_type size_type;

	class Term: public std::pair<const Key, Value> {
	public:
		Term() {
		}
		Term(const Key& key, const Value& v) :
				std::pair<const Key, Value>(key, v) {
		}
		Term(const std::pair<const Key, Value>& p) :
				std::pair<const Key, Value>(p.first, p.second) {
		}
		Term(const COE& coe, const TERM& term, const EXP& exp) :
				std::pair<const Key, Value>(std::pair<TERM, EXP>(term, exp),
						coe) {
		}

		const COE& coe() const {
			return this->second;
		}
		COE& coe() {
			return this->second;
		}
		const TERM& term() const {
			return this->first.first;
		}
		TERM& term() {
			return this->first.first;
		}
		const TERM& exp() const {
			return this->first.second;
		}
		TERM& exp() {
			return this->first.second;
		}

	};

public:
	Polynomial() :
			_map() {
	}
	Polynomial(const Polynomial& rhs) :
			_map(rhs._map) {
	}
	Polynomial operator=(const Polynomial& rhs) {
		this->_map = rhs._map;
		return *this;
	}
	/*
	 *  iterator
	 */
	iterator begin() {
		return _map.begin();
	}
	const_iterator begin() const {
		return _map.begin();
	}
	iterator end() {
		return _map.end();
	}
	const_iterator end() const {
		return _map.end();
	}
	/*
	 * capacity
	 */
	bool empty() const {
		return _map.empty();
	}
	size_type size() const {
		return _map.size();
	}
	Value& operator[](const Key& k) {
		return _map[k];
	}
	const Value& operator[](const Key& k) const {
		return _map[k];
	}
protected:
	typedef std::pair<iterator, bool> _ret;
public:

	_ret insert(const std::pair<const Key, Value> &x) {
		iterator it = _map.find(x.first);
		if (it != _map.end()){
			it->second = it->second + x.second;
		}
		return _map.insert(x);
	}

	void plus(const Polynomial &poly) {
		for (Polynomial::const_iterator iter = poly.begin(); iter != poly.end();
				++iter) {
			_ret ret = this->insert((*iter));
			if (ret.second == false) {
				(ret.first)->second += iter->second;
			}
		}
	}

	void minus(const Polynomial &poly) {
		for (Polynomial::const_iterator iter = poly.begin(); iter != poly.end();
				++iter) {
			Key k = iter->first;
			Value v = iter->second;
			Term term(k, -v);       //minus
			_ret ret = this->insert(term);
			if (ret.second == false) {
				(ret.first)->second = (ret.first)->second - v;
			}
		}
	}
	void times(const COE& rhs) {   //overload operator*
		for (Polynomial::const_iterator iter = this->begin();
				iter != this->end(); ++iter) {
			iter->second = iter->second * rhs;
		}
	}
	void divide(const COE& rhs) {  //overload operator/
		for (Polynomial::const_iterator iter = this->begin();
				iter != this->end(); ++iter) {
			iter->second = iter->second / rhs;
		}
	}

	void trim_zero() {
		iterator it_b = _map.begin();
		iterator it_e = _map.end();
		iterator it_t;

		while (it_b != it_e) {
			if (_is_zero_coe(it_b->second)) {  // Criteria checking here
				it_t = it_b;            // Keep a reference to the iter
				++it_b;                 // Advance in the map
				_map.erase(it_t);       // Erase it !!!
			} else {
				++it_b;                 // Just move on ...
			}
		}
	}

	void merge_const() {
		iterator it_b = _map.begin();
		iterator it_e = _map.end();
		iterator it_t;
		size_type flag = 0;
		std::pair<Key, Value> ct;
		while (it_b != it_e) {
			if (_is_zero_exp(it_b->first.second)) {  // Criteria checking here
				if (flag == 0) {
					ct.first.first = it_b->first.first;
					ct.first.second = it_b->first.second;
					ct.second = it_b->second;
				} else {
					ct.second = ct.second + it_b->second;
				}
				flag++;
				it_t = it_b;            // Keep a reference to the iter
				++it_b;                 // Advance in the map
				_map.erase(it_t);       // Erase it !!!
			} else {
				++it_b;                 // Just move on ...
			}
		}
		if (flag > 0) {
			_map.insert(ct);
		}
	}

	void concise(){
		trim_zero();
		merge_const();
	}

	void show() const;
};
//Polynomial operator-(const Polynomial &, const Polynomial &);
//Polynomial operator+(const Polynomial &, const Polynomial &);
//Polynomial operator*(const Float&, const Polynomial&);

}

#endif /* ALGEBRA_EXPRESSION_H_ */
