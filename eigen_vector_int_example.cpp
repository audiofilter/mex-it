#include "Eigen/Core"
#include <iostream>

typedef Eigen::VectorXd Vec_d;
typedef Eigen::Matrix<int,Eigen::Dynamic,1> Vec_i;

void mex_function(const Vec_i& x, const Vec_d& y, Vec_d& sum) {
    if (x.size() < y.size()) {
        std::cout << "Need x size to be equal or larger to y size\n";
        return;
    }
    sum = y;
    for (size_t i=0;i<y.size();i++) {
      sum[i] = sum[i] + (double)x[i];
      std::cout << x[i] << " + " << y[i] << " = " << sum[i] << "\n";
		}
}

#include "mex-it.h"
