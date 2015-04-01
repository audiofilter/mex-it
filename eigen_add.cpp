#include "Eigen/Core"
#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor> Matrixd;

//  This file defines a function callable from MATLAB once you mex it.
void mex_function(const Matrixd &x, const Matrixd &y, Matrixd &out1) {
	out1 = x + y;
}
#include "mex_wrap.cxx"
