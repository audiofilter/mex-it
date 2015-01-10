#include "Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrixd;

#include <iostream>

void mex_function(const Matrixd &x, const Matrixd &y, Matrixd &out1) {

	std::cout << "In mex_function x: = " << x.rows() << "*" << x.cols() << "\n" << x << "\n";
	std::cout << "In mex_function y: = " << y.rows() << "*" << y.cols() << "\n" << y << "\n";

	out1 = x * y;

	std::cout << "In mex_function out1: = " << out1.rows() << "*" << out1.cols() << "\n" << out1 << "\n";
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

	double A[N*M];
	double B[N*M];

	int c=0;
	for (int i=0;i<N;i++) {
		for (int j=0;j<M;j++) {
			A[c] = i+1;
			B[c] = (j%2==0) ? -i-1 : i+1;
			c++;
		}
	}

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

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> m(M,M);
	memcpy(m.data(),mxGetPr(plhs[0]),M*M*sizeof(double));

	cout << __FILE__ << " => mex_function: \nMatrix out = " << m << "\n";

}
