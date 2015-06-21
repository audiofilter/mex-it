#include "Eigen/Core"
#include <iostream>

typedef Eigen::VectorXd Vec_d;

void mex_function(const Vec_d &x, const Vec_d &y, Vec_d &out1) {

	out1 = x + y;
}
#include "mex_wrap.cxx"
