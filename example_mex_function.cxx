#include "mex_stuff.h"

using namespace mex_binding;
using namespace std;

void mex_function(const double &x, const double &y, double& result) {
	result = x + y;
	///	cout << "result of " << x << " + " << y << " is: " << result << endl;
}

#include "mex_wrap.cxx"
