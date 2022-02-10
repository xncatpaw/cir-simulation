/**
 * @file test_rng.cpp
 * @author bx.zh 
 * @brief Used to test the random generator class.
 * 
 * @copyright Copyleft (c) all wrongs reserved. 
 */

#include <iostream>
#include <vector>
#include <random>
#include "../include/simu_type.hpp"
#include "../include/random.hpp"

int main(int argc, char* argv[])
{
    auto gauss_gen = cir::Gaussian();
    // std::vector<int> shape({2, 3});
    size_t size = 6;
    auto _type = cir::PRECISE;
    if(argc >= 2)
    { 
        std::string str_size(argv[1]);
        size = size_t(std::stoi(str_size)); 
    }
    if(argc >= 3)
    {
        std::string str_type(argv[2]);
        _type = static_cast<cir::SimuType>(std::stoi(str_type));
    }
    std::cout << "size : " << size << ", ";
    std::cout << "type : " << _type << std::endl;
    double *rv_buf = new double[size]();
    gauss_gen.gen(rv_buf, size, _type);
    for(size_t i=0; i<size; ++i)
        std::cout << rv_buf[i] << " ";
    std::cout << std::endl;


    // std::cout << "Test the uniforme generator" << std::endl;
    // for(size_t i=0; i<size; ++i)
    //     std::cout << gauss_gen.uni_distr(gauss_gen.generator) << " ";
    // std::cout << std::endl;

    return 0;
}