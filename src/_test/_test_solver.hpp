// *
// * test_solver.h
// *
// *  Created on: May 3, 2015
// *      Author: zhou
//

#ifndef _TEST_SOLVER_H_
#define _TEST_SOLVER_H_

#include "../io/mmio.h"
#include "../algebra/matrix.hpp"
#include "../algebra/solver_matrix.hpp"
#include "gtest/gtest.h"

namespace carpio {

void test_gauss_e() {
	Matrix A(3, 3);
	A[0][0] = 2;
	A[0][1] = 1;
	A[0][2] = -1; //
	A[1][0] = -3;
	A[1][1] = -1;
	A[1][2] = 2; //
	A[2][0] = -2;
	A[2][1] = 1;
	A[2][2] = 2; //

	arrayList b(3);
	arrayList x(3);
	b[0] = 8;
	b[1] = -11;
	b[2] = -3;

	A.show();
	std::cout << " =====" << std::endl;
	b.show();

	solver_gaussian_elimination(A, b);
	A.show();
	b.show();

}

TEST(test_solver,Gauss) {
	Matrix A(3, 3);
	A[0][0] = 2;
	A[0][1] = 1;
	A[0][2] = -1; //
	A[1][0] = -3;
	A[1][1] = -1;
	A[1][2] = 2; //
	A[2][0] = -2;
	A[2][1] = 1;
	A[2][2] = 2; //
	arrayList b(3);
	arrayList x(3);
	b[0] = 8;
	b[1] = -11;
	b[2] = -3;
	ASSERT_EQ(-11, b[1]);
	//ASSERT_EQ(-3, b[2]);
	solver_gaussian_elimination(A, b);
	ASSERT_EQ(3, b.size());
	ASSERT_DOUBLE_EQ(3.0, b[1]);
	ASSERT_DOUBLE_EQ(-1.0, b[2]);
}

TEST(test_solver,read_matrix) {
	string workdir = "./src/_test/input_files/";
	MatrixSCO_<Float> mf;
	string fn_matrix = "685_bus";
	mm_read_mtx_sparse(workdir + fn_matrix + ".mtx", mf);
	ASSERT_EQ(3249, mf.NumNonzeros()) ;
	cout << "matrix information oo==============\n";
	cout << "i = " << mf.iLen() << "   j = " << mf.jLen() << endl;

	MatrixSCR_<Float> mfr(mf);
	cout << "matrix information  ==============\n";
	cout << "i = " << mfr.iLen() << "   j = " << mfr.jLen() << endl;

	cout << "matrix information  ==============\n";
	cout << "  " << mf(0, 0) << "    " << mfr(0, 0) << endl;
	ArrayListV<Float> b(mfr.iLen());
	b.assign(2.5);
	ArrayListV<Float> x(mfr.iLen());
	x.assign(1);

	//arrayListV<Float> ax = mf*x;
	//ax.show();
	ASSERT_DOUBLE_EQ(mf(20, 22), mfr(20, 22));

	//set up ========
	int max_iter = 1000;
	Float tol = 1e-6;
	std::list<Float> lr;  //list residual
	//solver =======================
	int res1 = IC_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	cout << "return code" <<res1<<endl;
	cout << "max iter = " << max_iter << endl;
	cout << "tol      = " << tol << endl;
	//gnuplot_show_ylog(lr);
	cout << "solver jacobi " << endl;
	x.assign(1);
	lr.clear();
	max_iter = 1000;
	tol = 1e-6;
	int res2 = Dia_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	cout << "return code" <<res2<<endl;
	cout << "max iter = " << max_iter << endl;
	cout << "tol      = " << tol << endl;
	//gnuplot_show_ylog(lr);
	//
	//gnuplot_show(mfr);
	//int i=0;
	//for(ListT<Float>::iterator it=lr.begin(); it!=lr.end(); it++){
	//	cout<<i<< "  "<<(*it)/nrm2(b)<<endl;
	//	i++;
	//}

	//cout << max_iter << "  " << tol << endl;
	//to row comp
}

} //end namespace

#endif /* TEST_SOLVER_H_ */
