mex-it
======

C++ 11 code to automatically create Matlab mex function based on generic C++ functions


The methodology is based on that used in the dlib library 
* supply the mex_function and let the compiler do the rest.


However, the code was completely reworked to use as many classes from std C++11 as possible and to make use of variadic templates to greatly reduce the amount of code needed.
* Currently code to callback matlab is not part of this library.

This supports basic C++ data types. 
Has limited support for Eigen Matrices (Matrix type of double + ColMajor)


####How do the compiler know the inputs and outputs?
* Non-const references mean outputs
* Everything else are inputs

### Build status 

* Either Cmake or setup your own mex script
* Mac OS X with Clang 6
* Linux GCC 4.9.2
* Visual Studio 2013 

### Example
Create a c++ implementation file like this

```c++
void mex_function(const double &x, const double &y, const double &z, double& result) {
	result = (x + y)*z;
}
#include "mex_wrap.cxx"
```

### Then build with CMake or using mex within Matlab

For single file example, you can do this in matlab (for recent GCC/Clang)

mex CXXFLAGS="\$CXXFLAGS -std=c++11" simple_example.cpp

### Requirements
Eigen needed for Eigen examples
*	Assumes Eigen is in either /usr/local/include/eigen3 or /usr/include/eigen3, please edit FindEigen.cmake for other paths
