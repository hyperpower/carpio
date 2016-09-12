/************************
 //  \file   MatrixT.h
 //  \brief
 // 
 //  \author zhou
 //  \date   23 janv. 2014 
 ***********************/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <iterator>

#include "array_list.hpp"

namespace carpio {

template<typename T>
class MatrixT {
public:
	// type definitions===================
	typedef T value_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef st size_type;
	typedef st difference_type;
protected:
	st m_iLen;
	st m_jLen;
	ArrayListT<T> *m_mp;
public:
	//constructor==========================
	MatrixT();
	MatrixT(const MatrixT<T>& a);
	MatrixT(size_type iLen, size_type jLen);
	MatrixT(size_type iLen, size_type jLen, T **value);
	void reconstruct(size_type iLen, size_type jLen);
	//=============================================
	MatrixT<T>& operator=(const MatrixT<T> &a);
	//=============================================
	~MatrixT();
	//Capacity=====================================
	size_type size() const;
	size_type size_i() const;
	size_type size_j() const;
	bool empty() const;
	//Element access===============================
	ArrayListT<T>& operator[](size_type index);
	const ArrayListT<T>& operator[](size_type index) const;
	reference operator()(size_type i, size_type j);
	const_reference operator()(size_type i, size_type j) const;
	reference at(size_type i, size_type j);
	const_reference at(size_type i, size_type j) const;
	T get(size_type i, size_type j);
	T* getpValue(size_type i, size_type j);
	void set(size_type i, size_type j, const T& value);
	void set_row(size_type i, const T& value);
	void set_col(size_type i, const T& value);
	void assign(const T& value);
	//=============================================

	void swap(size_type i1, size_type j1, size_type i2, size_type j2);

	size_type count_equal(T value); //overload ==

	inline bool check_idx_i(size_type);
	inline bool check_idx_j(size_type);

	//Low efficient function for small matrix
	void append_row(const ArrayListT<T> &);
	void appendCol(const ArrayListT<T> &);
	void insert_row_back(size_type, const ArrayListT<T>&);
	void insert_col_back(size_type, const ArrayListT<T>&);
	void delete_row(size_type);
	void delete_col(size_type);

protected:
	//class iterator: public std::iterator<std::input_iterator_tag, KeyType>
	//<class Category, class T, class Distance = ptrdiff_t,
	//       class Pointer = T*, class Reference = T&>
	template<class _T, class _Ptr, class _Ref, class _pMatrix>
	class _iterator: public std::iterator<std::bidirectional_iterator_tag, _T,
			ptrdiff_t, _Ptr, _Ref> {
	private:
		typedef _T value_type;
		typedef _Ptr pointer;
		typedef _Ref reference;

		_pMatrix pm;
		_Ptr ptr;
		size_type i;
		size_type j;
	public:
		_iterator() :
				pm(NULL), ptr(nullptr), i(0), j(0) {
			/* empty */
		}

		_iterator(_pMatrix pm, _Ptr p, size_type i, size_type j) {
			this->pm = pm;
			this->ptr = p;
			this->i = i;
			this->j = j;
		}

		_iterator(const _iterator& it) {
			pm = it.pm;
			ptr = it.ptr;
			i = it.i;
			j = it.j;
		}

		_iterator& operator ++() {
			i++;
			if (i >= pm->size_i()) {
				i = 0;
				j++;
				if (j >= pm->size_j()) {
					ptr = nullptr;
					return *this;
				}
			}
			ptr = &pm->operator()(i, j);
			return *this;
		}

		_iterator& operator --() {
			i--;
			if (i < 0) {
				i = pm->size_i() - 1;
				j--;
				if (j < 0) {
					i = 0;
					j = pm->size_j();
					ptr = nullptr;
					return (*this);
				}
			}
			ptr = pm->operator()(i, j);
			return *this;
		}

		_iterator operator --(int) {
			_iterator copy(*this);
			operator--();
			return copy;
		}

		_iterator operator ++(int) {
			_iterator copy(*this);
			operator++();
			return copy;
		}

		bool operator ==(const _iterator& rhs) {
			return pm == rhs.pm && ptr == rhs.ptr && i == rhs.i && j == rhs.j;
		}

		bool operator !=(const _iterator& rhs) {
			return !(*this == rhs);
		}

		reference operator *() const {
			return (*ptr);
		}

		pointer operator ->() const {
			return ptr;
		}
	};

public:
	typedef _iterator<T, T*, T&, MatrixT<T>*> iterator;
	typedef _iterator<T, const T*, const T&, const MatrixT<T>* > const_iterator;

	iterator begin() {
		return iterator(this, &(this->operator()(0, 0)), 0, 0);
	}

	const_iterator begin() const {
		return const_iterator(this, &(this->operator()(0, 0)), 0, 0);
	}

	iterator end() {
		return iterator(this, nullptr, 0, this->m_jLen);
	}
	const_iterator end() const {
		return const_iterator(this, nullptr, 0, this->m_jLen);
	}

};

template<typename T>
inline bool MatrixT<T>::check_idx_i(st i) {
	return (i >= 0 && i < m_iLen) ? true : false;
}

template<typename T>
inline bool MatrixT<T>::check_idx_j(st j) {
	return (j >= 0 && j < m_jLen) ? true : false;
}

template<typename T>
MatrixT<T>::MatrixT() {
	m_iLen = 0;
	m_jLen = 0;
	m_mp = nullptr;
}

template<typename T>
MatrixT<T>::MatrixT(const MatrixT<T>& a) {
	m_iLen = a.size_i();
	m_jLen = a.size_j();
	m_mp = new ArrayListT<T> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j] = a[i][j];
		}
	}
}

template<typename T>
MatrixT<T>::MatrixT(st iLen, st jLen) {
	m_iLen = iLen;
	m_jLen = jLen;
	//m_total = m_iLen * m_jLen;
	m_mp = new ArrayListT<T> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
}

template<typename T>
MatrixT<T>::MatrixT(st iLen, st jLen, T **value) {
	m_iLen = iLen;
	m_jLen = jLen;
	//m_total = m_iLen * m_jLen;
	m_mp = new ArrayListT<T> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j] = value[i][j];
		}
	}
}
template<typename T>
void MatrixT<T>::reconstruct(st iLen, st jLen) {
	m_iLen = iLen;
	m_jLen = jLen;
	//m_total = m_iLen * m_jLen;
	if (m_mp != nullptr) {
		delete[] m_mp;
	}
	m_mp = new ArrayListT<T> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i].reconstruct(m_jLen);
	}
}
template<typename T>
MatrixT<T>::~MatrixT() {
	delete[] m_mp;
}
template<typename T>
const ArrayListT<T>& MatrixT<T>::operator[](st index) const {
	ASSERT(index >= 0 && index < m_iLen);
	return m_mp[index];
}
template<typename T>
ArrayListT<T>& MatrixT<T>::operator[](st index) {
	ASSERT(index >= 0 && index < m_iLen);
	return m_mp[index];
}
template<typename T>
bool MatrixT<T>::empty() const {
	return m_mp == nullptr;
}

template<typename T>
MatrixT<T>& MatrixT<T>::operator=(const MatrixT<T> &a) {
	if (this == &a) {
		return *this;
	}
	if (m_iLen == a.size_i() && m_jLen == a.size_j()) {
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				this->m_mp[i][j] = a[i][j];
			}
		}
	} else {
		delete[] this->m_mp;
		m_iLen = a.size_i();
		m_jLen = a.size_j();
		//m_total = m_iLen * m_jLen;
		m_mp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				this->m_mp[i][j] = a[i][j];
			}
		}
	}
	return *this;
}

template<typename T>
st MatrixT<T>::size() const {
	return m_iLen * m_jLen;
}
template<typename T>
st MatrixT<T>::size_i() const {
	return m_iLen;
}
template<typename T>
st MatrixT<T>::size_j() const {
	return m_jLen;
}
template<typename T>
T MatrixT<T>::get(st i, st j) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}
template<typename T>
typename MatrixT<T>::reference MatrixT<T>::operator()(size_type i,
		size_type j) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}
template<typename T>
typename MatrixT<T>::const_reference MatrixT<T>::operator()(size_type i,
		size_type j) const {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}

template<typename T>
T& MatrixT<T>::at(st i, st j) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}

template<typename T>
const T& MatrixT<T>::at(st i, st j) const {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	return m_mp[i][j];
}

template<typename T>
T* MatrixT<T>::getpValue(st i, st j) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	return &m_mp[i][j];
}

template<typename T>
void MatrixT<T>::set(st i, st j, const T& value) {
	ASSERT(i >= 0 && i < m_iLen);
	ASSERT(j >= 0 && j < m_jLen);
	m_mp[i][j] = value;
}

template<typename T>
void MatrixT<T>::set_row(st i, const T& value) {
	ASSERT(i >= 0 && i < m_iLen);
	for (size_type j = 0; j < m_jLen; j++) {
		m_mp[i][j] = value;
	}
}

template<typename T>
void MatrixT<T>::set_col(st j, const T& value) {
	ASSERT(j >= 0 && j < m_jLen);
	for (size_type i = 0; i < m_iLen; i++) {
		m_mp[i][j] = value;
	}
}

template<typename T>
void MatrixT<T>::assign(const T& value) {
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[i][j] = value;
		}
	}
}

template<typename T>
void MatrixT<T>::swap(st i1, st j1, st i2, st j2) {
	ASSERT(i1 >= 0 && i1 < m_iLen && i2 >= 0 && i2 < m_iLen);
	ASSERT(j1 >= 0 && j1 < m_jLen && j2 >= 0 && j2 < m_jLen);
	if (i1 == i2 && j1 == j2) {
		return;
	}
	T tmp;
	tmp = m_mp[i1][j1];
	m_mp[i1][j1] = m_mp[i2][j2];
	m_mp[i2][j2] = tmp;
}

template<typename T>
st MatrixT<T>::count_equal(T value) {
	size_type num = 0;
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			if (m_mp[i][j] == value) {
				num++;
			}
		}
	}
	return num;
}
template<typename T>
void MatrixT<T>::append_row(const ArrayListT<T> &a) {
	ASSERT(a.size()==m_jLen||m_mp==NULL);
	if (m_mp == NULL) {
		m_iLen = 1;
		m_jLen = a.size();
		m_mp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i].reconstruct(m_jLen);
		}
		for (size_type j = 0; j < m_jLen; j++) {
			m_mp[0][j] = a[j];
		}
	} else {
		m_iLen += 1;
		ArrayListT<T> *tmp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen - 1; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type j = 0; j < m_jLen; j++) {
			tmp[m_iLen - 1][j] = a[j];
		}
		ArrayListT<T> *tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}
template<typename T>
void MatrixT<T>::appendCol(const ArrayListT<T> &a) {
	ASSERT(a.size()==m_iLen||m_mp==NULL);
	if (m_mp == NULL) {
		m_iLen = a.size();
		m_jLen = 1;
		m_mp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i].reconstruct(1);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			m_mp[i][0] = a[i];
		}
	} else {
		m_jLen += 1;
		ArrayListT<T> *tmp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < m_jLen - 1; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i][m_jLen - 1] = a[i];
		}
		ArrayListT<T> *tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}
template<typename T>
void MatrixT<T>::insert_row_back(MatrixT<T>::size_type ii,
		const ArrayListT<T>& a) {
	ASSERT(ii >= 0 && ii < m_iLen);
	ASSERT(a.Len() == m_jLen);
	m_iLen += 1;
	ArrayListT<T> *tmp = new ArrayListT<T> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		tmp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i <= ii; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			tmp[i][j] = m_mp[i][j];
		}
	}
	for (size_type j = 0; j < m_jLen; j++) {
		tmp[ii + 1][j] = a[j];
	}
	for (size_type i = ii + 2; i < m_iLen; i++) {
		for (size_type j = 0; j < m_jLen; j++) {
			tmp[i][j] = m_mp[i - 1][j];
		}
	}
	ArrayListT<T> *tmp2 = m_mp;
	m_mp = tmp;
	tmp = tmp2;
	delete[] tmp;
}
template<typename T>
void MatrixT<T>::insert_col_back(MatrixT<T>::size_type jj,
		const ArrayListT<T>& a) {
	ASSERT(jj >= 0 && jj < m_jLen);
	ASSERT(a.Len() == m_iLen);
	m_jLen += 1;
	ArrayListT<T> *tmp = new ArrayListT<T> [m_iLen];
	for (size_type i = 0; i < m_iLen; i++) {
		tmp[i].reconstruct(m_jLen);
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = 0; j <= jj; j++) {
			tmp[i][j] = m_mp[i][j];
		}
	}
	for (size_type i = 0; i < m_iLen; i++) {
		tmp[i][jj + 1] = a[i];
	}
	for (size_type i = 0; i < m_iLen; i++) {
		for (size_type j = jj + 2; j < m_jLen; j++) {
			tmp[i][j] = m_mp[i][j - 1];
		}
	}
	ArrayListT<T> *tmp2 = m_mp;
	m_mp = tmp;
	tmp = tmp2;
	delete[] tmp;
}
template<typename T>
void MatrixT<T>::delete_row(MatrixT<T>::size_type ii) {
	ASSERT(ii >= 0 && ii < m_iLen);
	m_iLen -= 1;
	if (m_iLen == 0) {
		m_iLen = 0;
		m_jLen = 0;
		delete[] m_mp;
		m_mp = NULL;
	} else {
		ArrayListT<T> *tmp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < ii; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type i = ii + 1; i < m_iLen + 1; i++) {
			for (size_type j = 0; j < m_jLen; j++) {
				tmp[i - 1][j] = m_mp[i][j];
			}
		}
		ArrayListT<T> *tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}
template<typename T>
void MatrixT<T>::delete_col(MatrixT<T>::size_type jj) {
	ASSERT(jj >= 0 && jj < m_iLen);
	m_jLen -= 1;
	if (m_jLen == 0) {
		m_iLen = 0;
		m_jLen = 0;
		delete[] m_mp;
		m_mp = NULL;
	} else {
		ArrayListT<T> *tmp = new ArrayListT<T> [m_iLen];
		for (size_type i = 0; i < m_iLen; i++) {
			tmp[i].reconstruct(m_jLen);
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = 0; j < jj; j++) {
				tmp[i][j] = m_mp[i][j];
			}
		}
		for (size_type i = 0; i < m_iLen; i++) {
			for (size_type j = jj + 1; j < m_jLen + 1; j++) {
				tmp[i][j - 1] = m_mp[i][j];
			}
		}
		ArrayListT<T> *tmp2 = m_mp;
		m_mp = tmp;
		tmp = tmp2;
		delete[] tmp;
	}
}

//end of class MatrixT===========================
//===============================================
//===============================================

template<typename T>
class MatrixV: public MatrixT<T> {
public:
	// type definitions===================
	typedef T value_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef st size_type;
	typedef st difference_type;
	//constructor==========================
	MatrixV();
	MatrixV(size_type iLen, size_type jLen);
	MatrixV(size_type iLen, size_type jLen, T **value);
	//void reconstruct(size_type size_i, size_type size_j);
	//~MatrixV();
	//=============================================
	MatrixV<T> operator+(const MatrixV<T> &a);
	MatrixV<T> operator-(const MatrixV<T> &a);
	MatrixV<T> operator*(const MatrixV<T> &a);
	//show ========================================
	void show() const;

};

template<typename T>
MatrixV<T>::MatrixV() :
		MatrixT<T>() {
}

template<typename T>
MatrixV<T>::MatrixV(st iLen, st jLen) :
		MatrixT<T>(iLen, jLen) {
	this->assign(0);
}

template<typename T>
MatrixV<T>::MatrixV(st iLen, st jLen, T **value) :
		MatrixT<T>(iLen, jLen) {
	for (size_type i = 0; i < this->m_iLen; i++) {
		for (size_type j = 0; j < this->m_jLen; j++) {
			this->m_mp[i][j] = value[i][j];
		}
	}
}

template<typename T>
MatrixV<T> MatrixV<T>::operator+(const MatrixV<T> &a) {
	ASSERT(a.iLen() == this->iLen());
	ASSERT(a.jLen() == this->jLen());
	MatrixT<T> sum(this->m_iLen, this->m_jLen);
	for (size_type i = 0; i < this->m_iLen; i++) {
		for (size_type j = 0; j < this->m_jLen; j++) {
			sum[i][j] = this->m_mp[i][j] + a[i][j];
		}
	}
	return sum;
}
template<typename T>
MatrixV<T> MatrixV<T>::operator-(const MatrixV<T> &a) {
	ASSERT(a.iLen() == this->iLen());
	ASSERT(a.jLen() == this->jLen());
	MatrixT<T> sum(this->m_iLen, this->m_jLen);
	for (size_type i = 0; i < this->m_iLen; i++) {
		for (size_type j = 0; j < this->m_jLen; j++) {
			sum[i][j] = this->m_mp[i][j] - a[i][j];
		}
	}
	return sum;
}
template<typename T>
MatrixV<T> MatrixV<T>::operator*(const MatrixV<T> &a) {
	ASSERT(a.iLen() == this->jLen());
	size_type nrow = this->m_iLen;
	size_type ncol = a.jLen();
	MatrixT<T> res(nrow, ncol);
	for (size_type i = 0; i < nrow; i++) {
		for (size_type j = 0; j < ncol; j++) {
			res[i][j] = 0;
			for (size_type k = 0; k < this->m_jLen; k++) {
				res[i][j] += this->m_mp[i][k] * a[k][j];
			}
		}
	}
	return res;
}
template<typename T>
void MatrixV<T>::show() const {
	std::cout << "> Matrix " << this->m_iLen << " x " << this->m_jLen << "\n";
	std::cout << "> ";
	for (int i = 0; i < this->m_iLen; i++) {
		for (int j = 0; j < this->m_jLen; j++) {
			std::cout << std::scientific << this->m_mp[i][j] << "  ";
		}
		std::cout << std::endl;
		std::cout << "> ";
	}
	std::cout << "< ----------\n";
}

//===============================================
//Function outside of the class==================
//===============================================
template<typename T>
ArrayListV<T> operator*(const MatrixV<T> &m, const ArrayListV<T> &a) {
	ASSERT(m.jLen() == a.size());
	arrayList res(a.size());
	for (int i = 0; i < m.iLen(); i++) {
		for (int j = 0; j < m.jLen(); j++) {
			res[i] += m[i][j] * a[j];
		}
	}
	return res;
}

typedef MatrixV<Float> Matrix;

template<typename T, st DIM1, st DIM2>
class MatrixS {
public:
	// type definitions===================
	typedef MatrixS<T, DIM1, DIM2> self;
	typedef T value_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef st size_type;
	typedef st difference_type;

	T elems[DIM1 > 0 ? DIM1 : 1][DIM2 > 0 ? DIM2 : 1];

public:
	//constructor==========================
	MatrixS() {
	}
	MatrixS(const T& a) {
		this->assign(a);
	}
	MatrixS(const MatrixS<T, DIM1, DIM2>& a) {
		typedef size_type st;
		for (st i = 0; i < DIM1; i++) {
			for (st j = 0; j < DIM2; j++) {
				elems[i][j] = a[i][j];
			}
		}
	}
	//=============================================
	self& operator=(const self &a) {
		if (this == &a) {
			return *this;
		} else {
			typedef size_type st;
			for (st i = 0; i < DIM1; i++) {
				for (st j = 0; j < DIM2; j++) {
					this->elems[i][j] = a[i][j];
				}
			}
		}
		return *this;
	}
	//=============================================
	~MatrixS() {
	}
	//Capacity=====================================
	static inline size_type size() {
		return DIM1 * DIM2;
	}
	static inline size_type iLen() {
		return DIM1;
	}
	static inline size_type jLen() {
		return DIM2;
	}
	//Element access===============================
	reference operator()(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	const_reference operator()(size_type i, size_type j) const {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	reference at(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	const_reference at(size_type i, size_type j) const {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	T get(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	T* getpValue(size_type i, size_type j) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		return elems[i][j];
	}
	void set(size_type i, size_type j, const T& value) {
		ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
		elems[i][j] = value;
	}
	void set_row(size_type i, const T& value) {
		ASSERT_MSG(i < DIM1, "out of range");
		typedef size_type st;
		for (st j = 0; j < DIM2; j++) {
			this->elems[i][j] = value;
		}
	}
	void set_col(size_type j, const T& value) {
		ASSERT_MSG(j < DIM2, "out of range");
		typedef size_type st;
		for (st i = 0; i < DIM2; i++) {
			this->elems[i][j] = value;
		}
	}
	void assign(const T& value) {
		typedef size_type st;
		for (st i = 0; i < DIM1; i++) {
			for (st j = 0; j < DIM2; j++) {
				this->elems[i][j] = value;
			}
		}
	}
	inline void ones() {
		assign(1);
	}
	inline void zeros() {
		assign(0);
	}
	//
	void show() const {
		std::cout << "> MatrixS " << DIM1 << " x " << DIM2 << "\n";
		std::cout << "> ";
		for (int i = 0; i < DIM1; i++) {
			for (int j = 0; j < DIM2; j++) {
				std::cout << std::scientific << this->elems[i][j] << "  ";
			}
			std::cout << std::endl;
			std::cout << "> ";
		}
		std::cout << "< ----------\n";
	}
};

typedef MatrixS<Float, 2, 2> Matrix2x2;
typedef MatrixS<Float, 3, 3> Matrix3x3;
typedef MatrixS<Float, 4, 4> Matrix4x4;

template<typename T, st DIM1, st DIM2>
static inline void zeros(MatrixS<Float, DIM1, DIM2>& m) {
	m.assign(0);
}

template<typename T, st DIM1, st DIM2>
static inline void ones(MatrixS<Float, DIM1, DIM2>& m) {
	m.assign(1);
}

/*
 * \brief calculate the determinant of a 2x2 matrix.
 *
 *        Adapted from:
 *        Matrix Inversion
 *        by Richard Carling
 *        from "Graphics Gems", Academic Press, 1990
 */
template<typename T>
static inline T det2x2(const T& a, const T& b, const T& c, const T& d) {
	return a * d - b * c;
}

template<typename T>
inline T determinant(const MatrixS<T, 2, 2>& m) {
	T a, b, c, d;
	a = m(0, 0);
	b = m(0, 1);
	c = m(1, 0);
	d = m(1, 1);
	return det2x2(a, b, c, d);
}

/*
 * \brief calculate the determinant of a 3x3 matrix
 *        in the form
 *
 *        | a1,  b1,  c1 |
 *        | a2,  b2,  c2 |
 *        | a3,  b3,  c3 |
 *
 *        Adapted from:
 *        Matrix Inversion
 *        by Richard Carling
 *        from "Graphics Gems", Academic Press, 1990
 */
template<typename T>
static inline T det3x3(const T& a1, const T& a2, const T& a3, const T& b1,
		const T& b2, const T& b3, const T& c1, const T& c2, const T& c3) {
	T ans3;
	ans3 = a1 * det2x2(b2, b3, c2, c3) - b1 * det2x2(a2, a3, c2, c3)
			+ c1 * det2x2(a2, a3, b2, b3);
	return ans3;
}

template<typename T>
inline T determinant(const MatrixS<T, 3, 3>& m) {
	const T& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2);
	const T& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2);
	const T& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2);
	return det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

/**
 * \brief matrix_determinant:
 * \param m Matrix 4 x 4 .
 *
 * \returns: the value of det(\p m).
 */
template<typename T>
inline T determinant(const MatrixS<T, 4, 4>& m) {
	T ans4;
	const T& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2), d1 = m(0, 3);
	const T& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2), d2 = m(1, 3);
	const T& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2), d3 = m(2, 3);
	const T& a4 = m(3, 0), b4 = m(3, 1), c4 = m(3, 2), d4 = m(3, 3);

	ans4 = a1 * det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
			- b1 * det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4)
			+ c1 * det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4)
			- d1 * det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

	return ans4;
}

template<typename T>
inline void adjoint(const MatrixS<T, 4, 4>& m, MatrixS<T, 4, 4>& ma) {
	const T& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2), d1 = m(0, 3);
	const T& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2), d2 = m(1, 3);
	const T& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2), d3 = m(2, 3);
	const T& a4 = m(3, 0), b4 = m(3, 1), c4 = m(3, 2), d4 = m(3, 3);
	/* row column labeling reversed since we transpose rows & columns */
	// the matrix without transpose is called cofactor matrix
	// the matrix with    transpose is called adjoint  matrix
	ma(0, 0) = det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
	ma(1, 0) = -det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
	ma(2, 0) = det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
	ma(3, 0) = -det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

	ma(0, 1) = -det3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
	ma(1, 1) = det3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
	ma(2, 1) = -det3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
	ma(3, 1) = det3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);

	ma(0, 2) = det3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
	ma(1, 2) = -det3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
	ma(2, 2) = det3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
	ma(3, 2) = -det3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);

	ma(0, 3) = -det3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
	ma(1, 3) = det3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
	ma(2, 3) = -det3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
	ma(3, 3) = det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

template<typename T>
inline void adjoint(const MatrixS<T, 3, 3>& m, MatrixS<T, 3, 3>& ma) {
	/* row column labeling reversed since we transpose rows & columns */
	// the matrix without transpose is called cofactor matrix
	// the matrix with    transpose is called adjoint  matrix
	ma(0, 0) = (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1));
	ma(0, 1) = (m(2, 1) * m(0, 2) - m(0, 1) * m(2, 2));
	ma(0, 2) = (m(0, 1) * m(1, 2) - m(1, 1) * m(0, 2));
	ma(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2));
	ma(1, 1) = (m(0, 0) * m(2, 2) - m(2, 0) * m(0, 2));
	ma(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2));
	ma(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
	ma(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1));
	ma(2, 2) = (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0));
}

template<typename T>
inline int inverse(const MatrixS<T, 4, 4>& m, MatrixS<T, 4, 4>& m_inv) {
	typedef MatrixS<T, 4, 4> matrix;
	T det = determinant(m);
	if (det == 0.) {
		return -1;
	}
	adjoint(m, m_inv);
	for (typename matrix::size_type i = 0; i < 4; ++i) {
		for (typename matrix::size_type j = 0; j < 4; ++j) {
			m_inv(i, j) /= det;
		}
	}
	return 1;
}

template<typename T>
inline int inverse(const MatrixS<T, 3, 3>& m, MatrixS<T, 3, 3>& m_inv) {
	typedef MatrixS<T, 3, 3> matrix;
	T det = determinant(m);
	if (det == 0.) {
		return -1;
	}
	adjoint(m, m_inv);
	for (typename matrix::size_type i = 0; i < m.iLen(); ++i) {
		for (typename matrix::size_type j = 0; j < m.jLen(); ++j) {
			m_inv(i, j) /= det;
		}
	}
	return 1;
}

template<typename T>
void transpose(const MatrixS<T, 4, 4> m, MatrixS<T, 4, 4> mt) {
	mt(0, 0) = m(0, 0);
	mt(1, 0) = m(0, 1);
	mt(2, 0) = m(0, 2);
	mt(3, 0) = m(0, 3);
	mt(0, 1) = m(1, 0);
	mt(1, 1) = m(1, 1);
	mt(2, 1) = m(1, 2);
	mt(3, 1) = m(1, 3);
	mt(0, 2) = m(2, 0);
	mt(1, 2) = m(2, 1);
	mt(2, 2) = m(2, 2);
	mt(3, 2) = m(2, 3);
	mt(0, 3) = m(3, 0);
	mt(1, 3) = m(3, 1);
	mt(2, 3) = m(3, 2);
	mt(3, 3) = m(3, 3);
}

template<typename T>
void transpose(const MatrixS<T, 3, 3> m, MatrixS<T, 3, 3> mt) {
	mt(0, 0) = m(0, 0);
	mt(1, 0) = m(0, 1);
	mt(2, 0) = m(0, 2);
	mt(0, 1) = m(1, 0);
	mt(1, 1) = m(1, 1);
	mt(2, 1) = m(1, 2);
	mt(0, 2) = m(2, 0);
	mt(1, 2) = m(2, 1);
	mt(2, 2) = m(2, 2);
}

template<typename T>
void transpose(const MatrixS<T, 2, 2> m, MatrixS<T, 2, 2> mt) {
	mt(0, 0) = m(0, 0);
	mt(1, 0) = m(0, 1);
	mt(0, 1) = m(1, 0);
	mt(1, 1) = m(1, 1);
}

}	//namespace=====================================

#endif /* MATRIXT_H_ */
