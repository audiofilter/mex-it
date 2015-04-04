#include "Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrixd;

#include <iostream>

void mex_function(const Matrixd &x, const Matrixd &y, Matrixd &out1) {
	out1 = x * y;

#ifdef DEBUG
	std::cout << "In mex_function x: = " << x.rows() << "*" << x.cols() << "\n" << x << "\n";
	std::cout << "In mex_function y: = " << y.rows() << "*" << y.cols() << "\n" << y << "\n";

	std::cout << "In mex_function out1: = " << out1.rows() << "*" << out1.cols() << "\n" << out1 << "\n";
#endif
}
#include "mex_wrap.cxx"

#include <iostream>

using namespace mex_binding;
using namespace std;

#define M 4
#define N 3

int main() {

	int nrhs = 2;
	int nlhs = 1;

	/* (4,3) matrix * (3,4) matrix => (4,4) matrix */

	double A[N*M] = {1,1,1,1,2,2,2,2,3,3,3,3};
	double B[N*M] = {-1,1,-1,1,-2,2,-2,2,-3,3,-3,3};
	double mult[M*M] = {-2,-2,-2,-2,3,3,3,3,-7,-7,-7,-7,6,6,6,6};
	double mex_mult[M*M];

	mxArray *in1 = mxCreateDoubleMatrix(M,N, mxREAL);
	mxArray *in2 = mxCreateDoubleMatrix(N,M, mxREAL);

	memcpy(mxGetPr(in1),A, N*M*sizeof(double));
	memcpy(mxGetPr(in2),B, N*M*sizeof(double));

	mxArray *out1 = mxCreateDoubleMatrix(M,M, mxREAL);

	const mxArray *prhs[nrhs];
	mxArray *plhs[nlhs];
	
	prhs[0] = in1;
	prhs[1] = in2;
	plhs[0] = out1;

	call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);

	memcpy(mex_mult,mxGetPr(plhs[0]),M*M*sizeof(double));

#ifdef DEBUG
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> m(M,M);
	memcpy(m.data(),mxGetPr(plhs[0]),M*M*sizeof(double));
	cout << __FILE__ << " => mex_function: \nMatrix out = " << m << "\n";
#endif

	// Go through all elements and check
	for (int i=0;i<M*M;i++) {
		assert(mult[i] == mex_mult[i]);
	}




}
