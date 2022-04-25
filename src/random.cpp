/**
 * @file random.cpp
 * @author bx.zh 
 * @brief defines the random number generators.
 * @date 2022-02-10
 * 
 * @copyright Copyleft (c) all wrongs reserved.
 */
#include <vector>
#include <random>
#include <assert.h>
#include <chrono>
#include <cmath>
#include <numeric>
#include "../include/simu_type.hpp"
#include "../include/random.hpp"

namespace cir
{
    // Static constants.
    double Gaussian::_SQRT_3 = std::sqrt(double(3.0));
    double Gaussian::_SQRT_6 = std::sqrt(double(6.0));
    double Gaussian::_SQRT_3p6 = std::sqrt(double(3.0) + Gaussian::_SQRT_6);
    double Gaussian::_SQRT_3m6 = std::sqrt(double(3.0) - Gaussian::_SQRT_6);
    double Gaussian::_PROB_2nd = double(1.0) / double(6.0);
    double Gaussian::_PROB_3rd = (Gaussian::_SQRT_6 - 2.0) / (4.0 * Gaussian::_SQRT_6);


    Gaussian::Gaussian() : generator(std::chrono::steady_clock::now().time_since_epoch().count()),
                            gauss_distr(), uni_distr()
    { }


    Gaussian::Gaussian(unsigned seed) : generator(seed), gauss_distr(), uni_distr()
    { }


    Gaussian::Gaussian(const Gaussian &other) : generator(other.generator), 
                                                gauss_distr(other.gauss_distr), 
                                                uni_distr(other.uni_distr)
    { }


    void Gaussian::operator()(double *p_out, const size_t& length, SimuType type)
    {
        gen(p_out, length, type);
    }


    void Gaussian::gen(double *p_out, const size_t& length, SimuType type)
    {
        // auto ord = shape.size();
        assert(length > 0);
        
        // Choose the distribution generator.
        if (type==PRECISE) // Use the exact standard normal distribution.
        {
            // Generate.
            for(size_t i=0; i<length; ++i)
            { *(p_out+i) = gauss_distr(generator); }
        }
        else if (type==SECOND_ORDER) // Use the \pm \sqrt(3) and 0.
        { 
            // The 2nd order discretized gaussian.
            // P(X= {0, -\sqrt(3), +\sqrt(3)}) = {2/3, 1/6, 1/6}. 
            for(size_t i=0; i<length; ++i) 
            {
                double new_val = 0.0;
                double _tmp = uni_distr(generator);
                if (_tmp <= _PROB_2nd)
                { new_val = _SQRT_3; }
                else if (_tmp >= 1-_PROB_2nd)
                { new_val = -1.0 * _SQRT_3; }
                
                *(p_out+i) = new_val;
            }
        }
        else
        { }

        return;
    }


    double Gaussian::operator()()
    {
        return gen();
    }


    double Gaussian::gen()
    {
        return gauss_distr(generator);
    }


    double Gaussian::gen_unif()
    {
        return uni_distr(generator);
    }


    NonCentralChi2::NonCentralChi2(double nu) : _nu(nu), 
            generator(std::chrono::steady_clock::now().time_since_epoch().count()),
            _gauss(), _chi_2(nu), _poisson()
    { }


    double NonCentralChi2::gen(double lambda)
    {
        _poisson.param(param_t(lambda/2));
        int N = _poisson(generator);
        double Y = 0.0;
        for(auto i=0; i<N; ++i)
        {
            auto G = _gauss(generator);
            Y += G*G;
            G = _gauss(generator);
            Y += G*G;
        }
        double Z = _chi_2(generator);
        return Y + Z;
    }


    double NonCentralChi2::operator()(double lambda)
    {
        return gen(lambda);
    }
} 
