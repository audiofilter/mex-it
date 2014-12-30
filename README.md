mex-it
======

C++ 11 code to automatically create Matlab mex function based on generic C++ functions


Currently this is not yet functional

The methodology is based on that used in the dlib library -> supply the mex_function and let the compiler do the rest.
However, the code was completely reworked to use as many classes from std C++11 as possible and to make use of variadic templates to greatly reduce the amount of code needed.
Currently code to callback matlab is not part of this library.
This supports basic C++ data types. Soon support for Eigen Matrices will be added.

How do the compiler know the inputs and outputs?
		Non-const references mean outputs
		Everything else are inputs

