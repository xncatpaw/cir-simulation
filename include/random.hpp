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
    };


    /**
     * @brief Class used to generate the non-central chi squared distribution.
     */
    // class NonCentralChi_2
    // {
    //     private:
    //         uint _p, _q;
    //         double _lam;
        
    //     public:
    //         NonCentralChi_2(uint p, uint q, double lam);
    // };
}