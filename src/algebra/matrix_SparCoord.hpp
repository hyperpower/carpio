/************************
 //  \file   MatrixSCO.h
 //  \brief
 // 
 //  \author zhou
 //  \date   16 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOORD_H_
#define MATRIXSPARCOORD_H_

#include "../carpio_define.hpp"
#include "array_list.hpp"
#include "matrix.hpp"

namespace carpio {

template<class VALUE> class MatrixSCC_;
template<class VALUE> class MatrixSCR_;

template<class VALUE>
class MatrixSCO_ {
public:
	typedef VALUE vt;
private:
	ArrayListV<vt> val_;    // data values (nz_ elements)
	ArrayListV<st> rowind_;    // row_ind (nz_ elements)
	ArrayListV<st> colind_;    // col_ind (nz_ elements)

	st nz_;                   // number of nonzeros
	st dim_[2];               // number of rows, cols
public:

	MatrixSCO_(void) :
			val_(0), rowind_(0), colind_(0), nz_(0) {
		dim_[0] = 0;
		dim_[1] = 0;
	}

	MatrixSCO_(const Matrix &M) {
		dim_[0] = M.size_i();
		dim_[1] = M.size_j();
		// count non-zero
		st countnz = 0;
		for (st i = 0; i < M.size_i(); i++) {
			for (st j = 0; j < M.size_j(); j++) {
				if (M[i][j] != 0) {
					countnz++;
				}
			}
		}
		nz_ = countnz;
		val_.reconstruct(nz_);
		rowind_.reconstruct(nz_);
		colind_.reconstruct(nz_);
		st idx = 0;
		for (st i = 0; i < M.size_i(); i++) {
			for (st j = 0; j < M.size_j(); j++) {
				if (M[i][j] != 0) {
					val_[idx] = M[i][j];
					rowind_[idx] = i;
					colind_[idx] = j;
					idx++;
				}
			}
		}
	}

	MatrixSCO_(const MatrixSCO_ &S) :
			val_(S.val_), rowind_(S.rowind_), colind_(S.colind_), nz_(S.nz_) {
		dim_[0] = S.dim_[0];
		dim_[1] = S.dim_[1];
	}

	MatrixSCO_(st M, st N, st nz, vt *val, st *r, st *c) :
			val_(val, nz), rowind_((*r), nz), colind_(*c, nz), nz_(nz) {
		dim_[0] = M;
		dim_[1] = N;
	}

	MatrixSCO_(const MatrixSCR_<vt> &R) :
			val_(R.NumNonzeros()), rowind_(R.NumNonzeros()), colind_(
					R.NumNonzeros()), nz_(R.NumNonzeros()) {
		dim_[0] = R.getiLen();
		dim_[1] = R.getjLen();
		int count = 0;
		int i, j;
//  Loop through rows...
		for (i = 1; i <= dim_[0]; i++) {
			for (j = count; j < R.row_ptr(i); j++) {
				val_[count] = R.val(count);
				colind_[count] = R.col_ind(count);
				rowind_[count] = i - 1;
				count++;
			}
		}
	}

	MatrixSCO_(const MatrixSCC_<vt> &C) :
			val_(C.NumNonzeros()), rowind_(C.NumNonzeros()), colind_(
					C.NumNonzeros()), nz_(C.NumNonzeros()) {
		dim_[0] = C.getiLen();
		dim_[1] = C.getjLen();

		int count = 0;
		int i, j;
//  Loop through columns...
		for (j = 1; j <= dim_[1]; j++) {
			for (i = count; i < C.col_ptr(j); i++) {
				val_[count] = C.val(count);
				rowind_[count] = C.row_ind(count);
				colind_[count] = j - 1;
				count++;
			}
		}
	}

	MatrixSCO_<vt>& operator=(const MatrixSCO_<vt> &C) {
		dim_[0] = C.dim_[0];
		dim_[1] = C.dim_[1];
		nz_ = C.nz_;
		val_ = C.val_;
		rowind_ = C.rowind_;
		colind_ = C.colind_;
		return *this;
	}

	MatrixSCO_<vt>& newsize(st M, st N, st nz) {
		dim_[0] = M;
		dim_[1] = N;
		nz_ = nz;
		val_.reconstruct(nz);
		rowind_.reconstruct(nz);
		colind_.reconstruct(nz);
		return *this;
	}

	MatrixSCO_<vt>& resize(st M, st N, st nz) {
		dim_[0] = M;
		dim_[1] = N;
		nz_ = nz;
		val_.resize(nz);
		rowind_.resize(nz);
		colind_.resize(nz);
		return *this;
	}

//slow---

	vt operator()(st i, st j) const {
		assert(i >= 0 && i < dim_[0]);
		assert(j >= 0 && j < dim_[1]);
		for (st t = 0; t < nz_; t++) {
			if (rowind_(t) == i && colind_(t) == j) {
				return val_(t);
			}
		}
		return 0.0;
	}

	ArrayListV<vt> operator*(const ArrayListV<vt> &x) const {
		st M = dim_[0];
		st N = dim_[1];
//  Check for compatible dimensions:
		assert(x.size() == N);
		ArrayListV<vt> res(M);
		for (st i = 0; i < nz_; i++) {
			res[rowind_[i]] += x[colind_[i]] * val_[i];
		}
		return res;
	}

	ArrayListV<vt> transMult(const ArrayListV<vt> &x) const {
		st tM = dim_[1];
		st tN = dim_[0];
		assert(!(x.Len() == tN));
		ArrayListV<vt> res(tM);
		for (st i = 0; i < nz_; i++) {
			res[colind_[i]] += x[rowind_[i]] * val_[i];
		}
		return res;
	}

	vt max() const {
		return val_.findMax();
	}

	vt min() const {
		return val_.findMin();
	}

	void show(st a) const {
		//if a==0 show matrix in Coordinate formate
		if (a == 0) {
			std::cout << "RowIdx " << "ColIdx " << "Value " << std::endl;
			for (st i = 0; i < nz_; i++) {
				std::cout << std::scientific << rowind_[i] << "  ";
				std::cout << std::scientific << colind_[i] << "  ";
				std::cout << std::scientific << val_[i] << "\n";
			}
		} else {
			for (st i = 0; i < this->iLen(); i++) {
				for (st j = 0; j < this->jLen(); j++) {
					bool isnz = 0;
					for (st t = 0; t < nz_; t++) {
						if (rowind_(t) == i && colind_(t) == j) {
							std::cout << std::scientific << val_(t) << "  ";
							isnz = 1;
						}
					}
					if (isnz == 0) {
						std::cout << std::scientific << 0.0 << "  ";
					}
				}
				std::cout << std::endl;
			}
		}

	}

	vt& val(st i) {
		return val_(i);
	}

	st& row_ind(st i) {
		return rowind_(i);
	}

	st& col_ind(st i) {
		return colind_(i);
	}

	const vt& val(st i) const {
		return val_(i);
	}

	const st& row_ind(st i) const {
		return rowind_(i);
	}

	const st& col_ind(st i) const {
		return colind_(i);
	}

	st size() const {
		return dim_[0] * dim_[1];
	}

	st iLen() const {
		return dim_[0];
	}

	st jLen() const {
		return dim_[1];
	}

	st NumNonzeros() const {
		return nz_;
	}

};

}
// This is the end of namespace

#endif /* MATRIXSPARCOORD_H_ */
