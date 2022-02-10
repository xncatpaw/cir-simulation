/**
 * @file random.hpp
 * @author bx.zh 
 * @brief defines the random number geenrators.
 * @date 2022-02-10
 * 
 * @copyright Copyleft (c) all wrongs reserved.
 */
#pragma once
#include <random>
#include <vector>
#include "simu_type.hpp"

namespace cir
{
    /**
     * @brief The Gaussian distribution generator class.
     */
    class Gaussian
    {
        private:
            std::default_random_engine generator;
            std::normal_distribution<double> gauss_distr;
        
        public:
            Gaussian();
            void operator()(double *p_out, std::vector<int> shape, SimuType type);
            void gen(double *p_out, std::vector<int> shape, SimuType type);
    };
}