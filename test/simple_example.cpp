// We just need the mex_function header file
// This could be put into mex-it.h but this way allows naming of the mex function in matlab to this file name (without the .cxx)

void mex_function(const double &x, const double &y, const double &z, int& result) {
	result = static_cast<int>((x + y)*z);
}

#include "mex-it.h"
