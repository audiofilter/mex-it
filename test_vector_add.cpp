#include <vector>
typedef std::vector<double> Vectord;

#include <iostream>

void mex_function(const Vectord &x, const Vectord &y, Vectord &out1, double &out2) {
	out1.resize(x.size());
	for (int i=0;i<x.size();i++) {
		out1[i] = x[i] + y[i];
	}
	out2 = x[0]-y[0];
#ifdef DEBUG
	std::cout << "In mex_function x[" << x.size() << "] = ";
	for (int i=0;i<x.size()-1;i++) std::cout << x[i] << ",";
	std::cout << x[x.size()-1] << "\n";

	std::cout << "In mex_function y[" << y.size() << "] = ";
	for (int i=0;i<y.size()-1;i++) std::cout << y[i] << ",";
	std::cout << y[y.size()-1] << "\n";

	std::cout << "In mex_function out1[" << out1.size() << "] = ";
	for (int i=0;i<out1.size()-1;i++) std::cout << out1[i] << ",";
	std::cout << out1[out1.size()-1] << "\n";

	std::cout << "In mex_function out2: = " << out2 << "\n";
#endif
}

#include "mex_wrap.cxx"

using namespace mex_binding;
using namespace std;

#define M 4

int main() {

	int nrhs = 2;
	int nlhs = 2;

	double A[] = {1,2,3,4};
	double B[] = {0,1,2,3};
	double sum[] = {1,3,5,7};
	double mex_sum[] = {0,0,0,0};

	mxArray *in1 = mxCreateDoubleMatrix(M,1, mxREAL);
	mxArray *in2 = mxCreateDoubleMatrix(M,1, mxREAL);

	memcpy(mxGetPr(in1),A, M*sizeof(double));
	memcpy(mxGetPr(in2),B, M*sizeof(double));

	mxArray *out1 = mxCreateDoubleMatrix(M,1, mxREAL);
	mxArray *out2 = mxCreateDoubleScalar(-1.0); // will be overwritten!

	const mxArray *prhs[nrhs];
	mxArray *plhs[nlhs];
	
	prhs[0] = in1;
	prhs[1] = in2;
	plhs[0] = out1;
	plhs[1] = out2;

	call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);

	double result = mxGetScalar(plhs[1]);

	std::vector<double> m(M);
	memcpy(m.data(),mxGetPr(plhs[0]),M*sizeof(double));

	memcpy(mex_sum,mxGetPr(plhs[0]),M*sizeof(double));

	for (int i=0;i<M;i++) {
		assert(mex_sum[i] == sum[i]);
	}
	assert(result == 1);


#ifdef DEBUG
	cout << __FILE__ << " => mex_function: \n\t\tMatrix out[" << M << "] = " << "\n\t\t";
	for (int i=0;i<M-1;i++) std::cout << m[i] << ",";
	std::cout << m[M-1] << "\n";

	cout << __FILE__ << " => mex_function: \n\t\tout2 is: " << result << endl;
#endif
	
}
