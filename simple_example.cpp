// We just need the mex_function header file
// This could be put into mex_wrap.cxx but this way allows naming of the mex function in matlab to this file name (without the .cxx)
#include "mex_stuff.h"

void mex_function(const double &x, const double &y, const double &z, int& result) {
	result = (x + y)*z;
}

using namespace mex_binding;
using namespace std;

#include "mex_wrap.cxx"
