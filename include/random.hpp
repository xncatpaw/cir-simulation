/**
 * @file random.hpp
 * @author bx.zh 
 * @brief defines the random number generators.
 * @date 2022-02-10
 * 
 * @copyright Copyleft (c) all wrongs reserved.
 */
#pragma once
#include <random>
// #include <vector>
#include "simu_type.hpp"

namespace cir
{
    using param_t = std::poisson_distribution<int>::param_type;

    /**
     * @brief The Gaussian distribution generator class.
     */
    class Gaussian
    {
    private:
        static double _SQRT_3;
        static double _SQRT_6;
        static double _SQRT_3p6;
        static double _SQRT_3m6;
        static double _PROB_2nd;
        static double _PROB_3rd;
        std::default_random_engine generator;
        std::normal_distribution<double> gauss_distr;
        std::uniform_real_distribution<double> uni_distr;
    
    public:
        Gaussian();
        Gaussian(unsigned seed);
        Gaussian(const Gaussian &); // Copy constructor.

        void operator()(double *p_out, const size_t& length, SimuType type=PRECISE);
        double operator()();
        void gen(double *p_out, const size_t& length, SimuType type);
        double gen();
        double gen_unif();
    };


    /**
     * @brief Class used to generate the non-central chi squared distribution.
     */
    class NonCentralChi2
    {
    private:
        double _nu;
        std::default_random_engine generator;
        std::normal_distribution<double> _gauss;
        std::chi_squared_distribution<double> _chi_2;
        std::poisson_distribution<int> _poisson;
    
    public:
        NonCentralChi2(double nu);
        double gen(double lambda);
        double operator()(double lambda);
    };
}