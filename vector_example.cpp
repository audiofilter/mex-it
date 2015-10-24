#include <iostream>
#include <vector>

void mex_function(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& sum) {
    if (x.size() < y.size()) {
        std::cout << "Need x size to be equal or larger to y size\n";
        return;
    }
    for (size_t i=0;i<y.size();i++) {
			sum.push_back( x[i] + y[i] );
		}
}

#include "mex-it.h"
