#include <iostream>
#include <vector>

void mex_function(const std::vector<int>& x, const std::vector<double>& y, std::vector<double>& sum) {
    if (x.size() < y.size()) {
        std::cout << "Need x size to be equal or larger to y size\n";
        return;
    }
    for (size_t i=0;i<y.size();i++) {
			sum.push_back( x[i] + y[i] );
      std::cout << x[i] << " + " << y[i] << " = " << sum[i] << "\n";
		}
}

#include "mex-it.h"
