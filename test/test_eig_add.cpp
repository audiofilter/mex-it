#include <iostream>

#include "Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrixd;

#include <iostream>

static void mex_function(const Matrixd &x, const Matrixd &y, Matrixd &out1, double &out2) {

#ifdef DEBUG
	std::cout << "In mex_function x: = " << x.rows() << "*" << x.cols() << "\n" << x << "\n";
	std::cout << "In mex_function y: = " << y.rows() << "*" << y.cols() << "\n" << y << "\n";
#endif
	out1 = x + y;
	out2 = (x(0,0)-y(0,0));

#ifdef DEBUG
	std::cout << "In mex_function out1: = " << out1.rows() << "*" << out1.cols() << "\n" << out1 << "\n";
	std::cout << "In mex_function out2: = \n" << out2 << "\n";
#endif
}
#include "mex-it.h"

using namespace mex_binding;
using namespace std;

#define M 4
#define N 3

int main() {
	int nrhs = 2;
	int nlhs = 2;

	double A[N*M];
	double B[N*M];
	double sum[N*M];
	double mex_sum[N*M];

	int c=0;
	for (int i=0;i<N;i++) {
		for (int j=0;j<M;j++) {
			A[c] = i+1;
			B[c] = (j%2==0) ? -i-1 : i+1;
			sum[c] = A[c] + B[c];
			c++;
		}
	}

	mxArray *in1 = mxCreateDoubleMatrix(M,N, mxREAL);
	mxArray *in2 = mxCreateDoubleMatrix(M,N, mxREAL);

	memcpy(mxGetPr(in1),A, N*M*sizeof(double));
	memcpy(mxGetPr(in2),B, N*M*sizeof(double));

	mxArray *out1 = mxCreateDoubleMatrix(M,N, mxREAL);
	mxArray *out2 = mxCreateDoubleScalar(-1.0);

	const mxArray *prhs[nrhs];
	mxArray *plhs[nlhs];
	
	prhs[0] = in1;
	prhs[1] = in2;
	plhs[0] = out1;
	plhs[1] = out2;

	// Inputs and outputs are 'Mex' matrices here, but Eigen:Matrices in function
	call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);

	double result = mxGetScalar(plhs[1]);

	memcpy(mex_sum,mxGetPr(plhs[0]),N*M*sizeof(double));

	// Go through all element and check
	for (int i=0;i<M*N;i++) {
		assert(sum[i] == mex_sum[i]);
	}

	assert(result == 2);

}
