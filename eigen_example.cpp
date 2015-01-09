#include "mex_stuff.h"
#include "Eigen/Core"
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrixf;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrixd;

//  This file defines a function callable from MATLAB once you mex it.

#include <iostream>

void mex_function(const Matrixd &x, const Matrixd &y, Matrixd &out1, double &out2, double some_number) {
	out1 = x + y;
	out2 = (x(0,0)-y(0,0));
	
	// we can also use cout to print things as usual:
	std::cout << "some_number: " << some_number << std::endl;
}
#include "mex_wrap.cxx"
