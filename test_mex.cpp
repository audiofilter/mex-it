#include "mex_function.h"
#include "mex_stuff.h"

using namespace mex_binding;
using namespace std;

#include "mex_wrap.cxx"

#include <iostream>
int main() {

	int nrhs = 3;
	int nlhs = 1;

	double x=3.0;
	double y=4.0;
	double z=5.0;

	mxArray *in1 = mxCreateDoubleScalar(x);
	mxArray *in2 = mxCreateDoubleScalar(y);
	mxArray *in3 = mxCreateDoubleScalar(z);

	mxArray *out1 = mxCreateDoubleScalar(-1.0);

	const mxArray *prhs[nrhs];
	mxArray *plhs[nlhs];
	
	prhs[0] = in1;
	prhs[1] = in2;
	prhs[2] = in3;
	plhs[0] = out1;

	call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);

	double result = mxGetScalar(plhs[0]);

	cout << __FILE__ << " => mex_function: result of (" << x << " + " << y << ") * " << z << " is: " << result << endl;

}
