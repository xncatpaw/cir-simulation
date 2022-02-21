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


    void CIR::_step_imp_3(double* p_X_crt, double* p_X_nxt, double dW, double h)
    {
        auto _1kh = _k * h + 1.0; // 1 + kh.
        auto term_sigdW = _sigma * dW;
        auto tmp = term_sigdW * term_sigdW; // \sigma^2 * dW^2.
        tmp += 4.0 * _1kh * (*p_X_crt + (_a - _sigma_sqr/2.0) * h);
        tmp = std::sqrt(tmp) + term_sigdW;
        tmp = tmp / (2.0 * _1kh);
        *p_X_nxt = (tmp * tmp);
    }


    void CIR::_step_imp_4(double* p_X_crt, double* p_X_nxt, double dW, double h)
    {
        auto _1kh = _k * h /2.0 + 1.0; // 1 + kh/2. 
        auto term_sigdW = _sigma/2.0 * dW;
        auto term_X_crt_sqrt = std::sqrt(*p_X_crt);

        auto tmp = term_sigdW + term_X_crt_sqrt;
        tmp = tmp * tmp + 4.0 * _1kh * (_a - _sigma_sqr/4.0) / 2.0 * h;
        tmp = std::sqrt(tmp) + term_sigdW + term_X_crt_sqrt;
        tmp = tmp / (_1kh * 2.0);
        *p_X_nxt = (tmp * tmp);
    }
}
