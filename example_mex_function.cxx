// We just need the mex_function header file
// This could be put into mex_wrap.cxx but this way allows naming of the mex function in matlab to this file name (without the .cxx)
#include "mex_function.h"
#include "mex_stuff.h"

using namespace mex_binding;
using namespace std;

#include "mex_wrap.cxx"
