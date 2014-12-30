#include "mex_stuff.h"
#include "mex.h"
#include <iostream>

using namespace mex_binding;
using namespace std;

void mex_function(const double &x, const double &y, const double& z, double& result) {
	result = int((x + y)*z);
	cout << __FILE__ << " => mex_function: result of (" << x << " + " << y << ") * " << z << " is: " << result << endl;
}

#include "mex_wrap.cxx"

int main() {

	int nrhs = 3;
	int nlhs = 1;

	mxArray *in1 = mxCreateDoubleScalar(3.0);
	mxArray *in2 = mxCreateDoubleScalar(4.0);
	mxArray *in3 = mxCreateDoubleScalar(5.0);

	mxArray *out1 = mxCreateDoubleScalar(-1.0);

	const mxArray *prhs[nrhs];
	mxArray *plhs[nlhs];
	
	prhs[0] = in1;
	prhs[1] = in2;
	prhs[2] = in3;
	plhs[0] = out1;

	call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);

}
