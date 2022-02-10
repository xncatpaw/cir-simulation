/**
 * @file test_rng.cpp
 * @author bx.zh 
 * @brief Used to test the random generator class.
 * 
 * @copyright Copyleft (c) all wrongs reserved. 
 */

#include <iostream>
#include <vector>
#include "../include/simu_type.hpp"
#include "../include/random.hpp"

int main()
{
    auto gauss_gen = cir::Gaussian();
    std::vector<int> shape({2, 3});
    std::cout << "shape is " << shape[0] << ", " << shape[1] << std::endl;
    double *rv_buf = new double[6]();
    gauss_gen.gen(rv_buf, shape, cir::PRECISE);
    for(auto i=0; i<6; ++i)
        std::cout << rv_buf[i] << " ";
    std::cout << std::endl;

    return 0;
}