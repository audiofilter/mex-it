#include "mex_function.h"

void mex_function(const double &x, const double &y, const double &z, double& result) {
	result = (x + y)*z;
}

