/************************
 //  \file   MatrixSCR.h
 //  \brief
 // 
 //  \author zhou
 //  \date   25 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOMPROW_H_
#define MATRIXSPARCOMPROW_H_

#include "../carpio_define.hpp"
#include "array_list.hpp"
#include "matrix.hpp"

namespace carpio {

template<class VALUE> class MatrixSCO_;
template<class VALUE> class MatrixSCC_;

template<class VALUE>
class MatrixSCR_ {
public:
	typedef VALUE vt;
private:
	ArrayListV<vt> val_;       // data values (nz_ elements)
	ArrayListV<st> rowptr_;    // row_ptr     (dim_[0]+1 elements)
	ArrayListV<st> colind_;    // col_ind     (nz_ elements)

	st nz_;                   // number of nonzeros
	st dim_[2];               // number of rows, cols

public:
	MatrixSCR_(void) :
			val_(0), rowptr_(0), colind_(0), nz_(0) {
		dim_[0] = 0;
		dim_[1] = 0;
	}
	//MatrixSCC(const Matrix &M);
	MatrixSCR_(const MatrixSCR_<vt> &S) :
			val_(S.val_), rowptr_(S.rowptr_), colind_(S.colind_), nz_(S.nz_) {
		dim_[0] = S.dim_[0];
		dim_[1] = S.dim_[1];
	}
	MatrixSCR_(const MatrixSCC_<vt> &C) :
			val_(C.NumNonzeros()),  //
			rowptr_(C.iLen() + 1),  //
			colind_(C.NumNonzeros()),  //
			nz_(C.NumNonzeros())  //
	{
		dim_[0] = C.getiLen();
		dim_[1] = C.getjLen();

		st i, j;

		ArrayListV < st > tally(C.iLen() + 1, 0);
		//      First pass through nonzeros.  Tally entries in each row.
		//      And calculate rowptr array.
		for (i = 0; i < nz_; i++) {
			tally[C.row_ind(i)]++;
		}
		rowptr_[0] = 0;
		for (j = 0; j < dim_[0]; j++)
			rowptr_(j + 1) = rowptr_(j) + tally(j);
		//      Make copy of rowptr for use in second pass.
		tally = rowptr_;
		//      Second pass through nonzeros.   Fill in index and value entries.
		st count = 0;
		for (i = 1; i <= dim_[1]; i++) {
			for (j = count; j < C.col_ptr(i); j++) {
				val_[tally(C.row_ind(j))] = C.val(j);
				colind_[tally(C.row_ind(j))] = i - 1;
				tally[C.row_ind(count)]++;
				count++;
			}
		}

	}

	MatrixSCR_(const MatrixSCO_<vt> &CO) :
			val_(CO.NumNonzeros()), rowptr_(CO.iLen() + 1), colind_(
					CO.NumNonzeros()), nz_(CO.NumNonzeros()) {
		dim_[0] = CO.iLen();
		dim_[1] = CO.jLen();

		st i;
		ArrayListV < st > tally(CO.iLen() + 1, 0);
		//      First pass through nonzeros.  Tally entries in each row.
		//      And calculate rowptr array.
		for (i = 0; i < nz_; i++) {
			tally[CO.row_ind(i)]++;
		}
		rowptr_(0) = 0;
		for (i = 0; i < dim_[0]; i++) {
			rowptr_(i + 1) = rowptr_(i) + tally(i);
		}
		//      Make copy of rowptr for use in second pass.
		tally = rowptr_;
		//      Second pass through nonzeros.   Fill in index and value entries.
		for (i = 0; i < nz_; i++) {
			val_[tally(CO.row_ind(i))] = CO.val(i);
			colind_[tally(CO.row_ind(i))] = CO.col_ind(i);
			tally[CO.row_ind(i)]++;
		}

	}

	MatrixSCR_(st M, st N, st nz, vt *val, st *r, st *c) :
			val_(val, nz), rowptr_(*r, M + 1), colind_(*c, nz), nz_(nz) {
		dim_[0] = M;
		dim_[1] = N;
	}

	MatrixSCR_(st M, st N, st nz, const ArrayListV<vt> &val,
			const ArrayListV<st> &r, const ArrayListV<st> &c) :
			val_(val), rowptr_(r), colind_(c), nz_(nz) {
		dim_[0] = M;
		dim_[1] = N;
	}

	vt& val(st i) {
		return val_(i);
	}

	st& col_ind(st i) {
		return colind_(i);
	}

	st& row_ptr(st i) {
		return rowptr_(i);
	}

	const vt& val(st i) const {
		return val_(i);
	}

	const st& col_ind(st i) const {
		return colind_(i);
	}

	const st& row_ptr(st i) const {
		return rowptr_(i);
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

	vt max() const {
		return val_.findMax();
	}

	vt min() const {
		return val_.findMin();
	}

	MatrixSCR_<vt>& operator=(const MatrixSCR_<vt> &R) {
		dim_[0] = R.dim_[0];
		dim_[1] = R.dim_[1];
		nz_ = R.nz_;
		val_ = R.val_;
		rowptr_ = R.rowptr_;
		colind_ = R.colind_;
		return *this;
	}

	MatrixSCR_<vt>& newsize(st M, st N, st nz) {
		dim_[0] = M;
		dim_[1] = N;
		nz_ = nz;
		val_.reconstruct(nz);
		rowptr_.reconstruct(M + 1);
		colind_.reconstruct(nz);
		return *this;
	}

	vt operator()(st i, st j) const {
		ASSERT(i >= 0 && i < dim_[0]);
		ASSERT(j >= 0 && j < dim_[1]);
		for (st t = rowptr_(i); t < rowptr_(i + 1); t++) {
			if (colind_(t) == j)
				return val_(t);
		}
		return 0.0;
	}

	ArrayListV<vt> operator*(const ArrayListV<vt> &x) const {
		st M = dim_[0];
		st N = dim_[1];
		//  Check for compatible dimensions:
		ASSERT(x.size() == N);

		ArrayListV<vt> res(M);
		for (st i = 0; i < M; ++i) {
			for (st j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
				res[i] += x[colind_[j]] * val_[j];
			}
		}
		return res;
	}

	ArrayListV<vt> transMult(const ArrayListV<vt> &x) const {
		st Mt = dim_[1];
		st Nt = dim_[0];
		//  Check for compatible dimensions:
		ASSERT(x.size() == Nt);
		ArrayListV<vt> res(Mt);
		for (st i = 0; i < Mt; ++i) {
			for (st j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
				res[i] += x[colind_[j]] * val_[j];
			}
		}
		return res;
	}

	void trans() {
		MatrixSCC_<vt> C(dim_[1], dim_[0], nz_, val_, colind_, rowptr_);
		dim_[0] = C.getiLen();
		dim_[1] = C.getjLen();

		st i, j;

		ArrayListV < st > tally(C.getiLen() + 1, 0);
		//      First pass through nonzeros.  Tally entries in each row.
		//      And calculate rowptr array.
		for (i = 0; i < nz_; i++) {
			tally[C.row_ind(i)]++;
		}
		rowptr_[0] = 0;
		for (j = 0; j < dim_[0]; j++)
			rowptr_(j + 1) = rowptr_(j) + tally(j);
		//      Make copy of rowptr for use in second pass.
		tally = rowptr_;
		//      Second pass through nonzeros.   Fill in index and vt entries.
		st count = 0;
		for (i = 1; i <= dim_[1]; i++) {
			for (j = count; j < C.col_ptr(i); j++) {
				val_[tally(C.row_ind(j))] = C.val(j);
				colind_[tally(C.row_ind(j))] = i - 1;
				tally[C.row_ind(count)]++;
				count++;
			}
		}
	}

	MatrixSCR_<vt> getTrans() const {
		MatrixSCR_ res(dim_[0], dim_[1], nz_, val_, rowptr_, colind_);
		res.trans();
		return res;
	}

	void show(st a) const {
		if (a == 0) {
			std::cout << "RowPtr " << "ColInd " << "vt " << std::endl;
			for (st i = 0; i < dim_[0]; i++) {
				for (st ii = rowptr_[i]; ii < rowptr_[i + 1]; ii++) {
					std::cout << std::scientific << i << "  ";
					std::cout << std::scientific << colind_[ii] << "  ";
					std::cout << std::scientific << val_[ii] << "\n";
				}
			}
		} else {
			for (st i = 0; i < dim_[0]; i++) {
				for (st j = 0; j < dim_[1]; j++) {
					bool flag = 0;
					for (st ii = rowptr_[i]; ii < rowptr_[i + 1]; ii++) {
						if (colind_[ii] == j) {
							std::cout << std::scientific << val_[ii] << "  ";
							flag = 1;
						}
					}
					if (flag == 0) {
						std::cout << std::scientific << 0.0 << "  ";
					}
				}
				std::cout << std::endl;
			}
		}
	}
};

}
// This is the end of namespace

#endif /* MATRIXSPARCOMPROW_H_ */
