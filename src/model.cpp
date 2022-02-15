/**
 * @file model.hpp
 * @author bx.zh 
 * @brief Defined the CIR simulation model here.
 * @date 2022-02-15
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 */

#include <random>
#include <assert.h>
#include <cmath>
#include "../include/simu_type.hpp"
#include "../include/random.hpp"
#include "../include/model.hpp"

namespace cir
{
    CIR::CIR(double k, double a, double sigma) : _k(k), _a(a), _sigma(sigma),
            _sigma_sqr(sigma*sigma)
    {
        _nu = 4 * _a / _sigma_sqr;
    }


    double CIR::psi(double k, double t)
    {
        double res = 0;
        if(k==0)
        { res = t; }
        else
        { res = (1 - std::exp(-k*t)) / k; }
        return res;
    }


    void CIR::psi(double *p_k, double *p_t, double *p_out)
    {
        double res = 0;
        if((*p_k)==0)
        { res = *(p_t); }
        else
        { res = (1 - std::exp(-(*p_k) * (*p_t))) / (*p_k); }
        *p_out = res;
    }


}
