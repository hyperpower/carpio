// *
// * test_solver.h
// *
// *  Created on: May 3, 2015
// *      Author: zhou
//

#ifndef _TEST_SOLVER_H_
#define _TEST_SOLVER_H_

#include "../io/mmio.h"
#include "../utility/format.h"
#include "../algebra/algebra.hpp"
#include "gtest/gtest.h"
#include "../io/gnuplot.h"
#include "../io/io_gnuplot_domain.h"

namespace carpio {
using namespace std;

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
	string fn_matrix = "pde225";
	mm_read_mtx_sparse(workdir + fn_matrix + ".mtx", mf);
	//mm_read_mtx_sparse(".mtx", mf);
	//ASSERT_EQ(3249, mf.NumNonzeros());
	cout << "matrix information oo==============\n";
	fmt::print(" i = {:>5} j = {:>5}\n", mf.iLen(), mf.jLen());

	MatrixSCR_<Float> mfr(mf);
	//mf.show(0);
	// read b
	ArrayListV<Float> b = mfr.sum_row();
	//string fn_array = "e05r0100_rhs1";
	//mm_read_array(workdir + fn_array + ".mtx", b);
	//b.show();
	ArrayListV<Float> x(mfr.iLen());
	x.assign(0);

	//arrayListV<Float> ax = mf*x;
	//ax.show();
	ASSERT_DOUBLE_EQ(mf(20, 22), mfr(20, 22));

	//set up ========
	int max_iter = 1000;
	Float tol    = 1e-6;
	std::list<Float> lr;  //list residual
	//solver =======================
	//int res1 = IC_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	//EXPECT_EQ(0, res1) << " >! Solver return " << res1;
	//cout << "max iter = " << max_iter << endl;
	//cout << "tol      = " << tol << endl;
	//cout << "x 10       " << x[10] << endl;

	//gnuplot_show_ylog(lr);
	cout << "solver jacobi " << endl;
	x.assign(0);
	lr.clear();
	max_iter = 5000;
	tol = 1e-12;
	int res2 = BiCGSTAB(mfr, x, b, max_iter, tol, lr);

	fmt::print("{}", "BiCGSTAB -----------------\n");
	fmt::print("return code : {:>10d}\n", res2);
	fmt::print("max iter    : {:>10d}\n", max_iter);
	fmt::print("max iter    : {:>10d}\n", lr.size());
	fmt::print("tol         : {:>10f}\n", tol);

	//res2 = BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	//cout << "BiCG stab -----------------\n";
	//cout << "return code" << res2 << endl;
	//cout << "max iter = " << max_iter << endl;
	//cout << "tol      = " << tol << endl;
	//cout << "x 10       " << x[10] << endl;
	gnuplot_show_ylog(lr);
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
#ifdef VIENNACL_WITH_OPENCL

TEST(test_solver,cl_matrix) {

	string workdir = "./src/_test/input_files/";
	MatrixSCO_<Float> mf;
	string fn_matrix = "steam3";
	mm_read_mtx_sparse(workdir + fn_matrix + ".mtx", mf);
	//ASSERT_EQ(3249, mf.NumNonzeros());
	cout << "matrix information oo==============\n";
	cout << "i = " << mf.iLen() << "   j = " << mf.jLen() << endl;

	MatrixSCR_<Float> mfr(mf);
	cout << "matrix information  ==============\n";
	cout << "i = " << mfr.iLen() << "   j = " << mfr.jLen() << endl;

	cout << "matrix information  ==============\n";
	cout << "  " << mf(0, 0) << "    " << mfr(0, 0) << endl;

	ArrayListV<Float> b(mfr.iLen());
	b.assign(1);
	ArrayListV<Float> x(mfr.iLen());
	x.assign(1);

	Float tol = 1e-6;
	int mi = 1000;
	Solver_<Float>::Solve(mfr, x, b, tol, mi);
	std::cout << "x 10    " << x[10] << "\n";

	//viennacl::vector<Float> vclb(mfr.iLen());
	//viennacl::vector<Float> vclx(mfr.iLen());
	//Solver_<Float>::MatSCR_vcl mat_vcl(mfr.size1(), mfr.size2(),
	//		mfr.NumNonzeros());
	//copy(b.begin(), b.end(), vclb.begin());
	//copy(x.begin(), x.end(), vclx.begin());
	//Solver_<Float>::Copy(mfr, mat_vcl, mfr.NumNonzeros());

	//copy(vclx.begin(), vclx.end(), x.begin());

}
#endif

} //end namespace

#endif /* TEST_SOLVER_H_ */
