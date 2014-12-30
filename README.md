mex-it
======

C++ 11 code to automatically create Matlab mex function based on generic C++ functions


The methodology is based on that used in the dlib library 
* supply the mex_function and let the compiler do the rest.


However, the code was completely reworked to use as many classes from std C++11 as possible and to make use of variadic templates to greatly reduce the amount of code needed.
* Currently code to callback matlab is not part of this library.

This supports basic C++ data types. Soon support for Eigen Matrices may be added.


####How do the compiler know the inputs and outputs?
* Non-const references mean outputs
* Everything else are inputs

### Build status 

* Either Cmake or setup your own mex script
* Mac OS X with Clang 6
* Linux GCC 4.9.2
Currently doesn't work on Travis due to limitations


### Example
Create a c++ implementation file like this

```c++
void mex_function(const double &x, const double &y, const double &z, double& result) {
	result = (x + y)*z;
}
#include "mex_stuff.h"
using namespace mex_binding;
using namespace std;
#include "mex_wrap.cxx"
```