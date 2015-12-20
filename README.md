mex-it
======

C++ 11 code to automatically create Matlab mex function based on generic C++ functions.

###Only a single header file "mex-it.h" is required to be included as in example below


The methodology is based on that used in the dlib library 
* supply the mex_function and let the compiler do the rest.


However, the code was completely reworked to use as many classes from std C++11 as possible and to make use of variadic templates to greatly reduce the amount of code needed.
* Currently code to callback matlab is not part of this library.

This supports basic C++ data types. 
Has limited support for Eigen Matrices (Matrix type of double + ColMajor), plus Eigen Vectors (see examples)


####How does the compiler know the inputs and outputs?
* Non-const references mean outputs
* Everything else are inputs

### Build status 

* Either Cmake or setup your own mex script
* Mac OS X with Clang 6
- Linux : Travis  [![Build Status](https://travis-ci.org/audiofilter/mex-it.png)](https://travis-ci.org/audiofilter/mex-it)
* Visual Studio 2013 (32-bit for Octave use): AppVeyor [![Build status](https://ci.appveyor.com/api/projects/status/4n3dshqn0oh24g0y?svg=true)](https://ci.appveyor.com/project/audiofilter/mex-it)

### Example
Create a c++ implementation file like this

```c++
void mex_function(const double &x, const double &y, const double &z, double& result) {
	result = (x + y)*z;
}
#include "mex-it.h"
```

### Then build with CMake or using mex within Matlab

For single file example, you can do this in matlab (for recent GCC/Clang)

```sh
mex CXXFLAGS="\$CXXFLAGS -std=c++11" simple_example.cpp
```
### Requirements
Eigen needed for Eigen examples
*	Assumes Eigen is in either /usr/local/include/eigen3 or /usr/include/eigen3, please edit FindEigen.cmake for other paths

### Other examples

C++ files starting with test_ are examples that can't generate mex files but do test the interfaces in pure C++ using mex header
files

* test_eig_add.cpp     -> Test of adding 2 Eigen Matrices of type 'double'
* test_eig_mult.cpp    -> Test of multipling 2 Eigen Matrices of type 'double'
* test_vector_add.cpp  -> Test of adding std::vector<double>
* test_mex.cpp -> Uses include of 'mex_function.h' to test the function in 'mex_function.cpp'

Mex example

* eigen_add.cpp -> For creating mex function adds 2 eigen matrices
* simple_example.cpp -> A very basic example that returns (x+y)*z for inputs x,y, and z
* vector_example.cpp -> An example that returns (x+y)*z for inputs x,y, and z
* eigen_vector_example.cpp  -> An example of a vector add for Eigen (double) vectors
* eigen_example.cpp  -> Same example as in test_eig_add.cpp but here can generate a mex function


Matlab/Octave Scripts to test
* test_eigen_vector_example.m    -> Test of adding 2 Eigen Vectors of type 'double' - also checks if wrong type is used
* test_eigen_add.m               -> Test of adding 2 Eigen Matrices of type 'double' - also checks if wrong type is used
* test_vector_example.m          -> Test of adding 2 std::vector<double> - also checks if wrong type is used
