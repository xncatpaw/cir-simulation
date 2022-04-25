/**
 * @file hes.cpp
 * @author bx.zh 
 * @brief Defined the Heston simulation model here.
 * @date 2022-03-04
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */
#include <assert.h>
#include <cmath>
#include "../include/simu_type.hpp"
#include "../include/random.hpp"
#include "../include/cir.hpp"
#include "../include/hes.hpp"
namespace cir
{
    Heston::Heston(double mu, double k, double a, double sigma, double rho) :
        _vol(k, a, sigma), _mu(mu), _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
    { }


    Heston::Heston(double mu, double k, double a, double sigma, double rho, unsigned seed) :
        _vol(k, a, sigma), _gauss(seed), _mu(mu), _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
    { }


    Heston::Heston(double mu, double rho, const CIR & other) : 
         _vol(other), _gauss(), _mu(mu), _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
    { }


    Heston::~Heston() { }


    double Heston::_step_eul(double S_cur, double V_cur, double dWS, double h)
    {
        auto idx = (_mu - V_cur/2)*h + std::sqrt(V_cur)*dWS;
        auto S_nxt = S_cur * std::exp(idx);
        return S_nxt;
    }


    double Heston::_step_cmb(double S_cur, double V_cur, double V_nxt, double dWS, double h,
                            double gamma_1, double gamma_2)
    {
        auto V_cmb = gamma_1*V_cur + gamma_2*V_nxt;
        auto idx = (_mu - V_cmb/2)*h + std::sqrt(V_cmb)*dWS;
        auto S_nxt = S_cur * std::exp(idx);
        return S_nxt;
    }


    void Heston::_step(HESStepFnT hes_step, 
                    double* p_S_cur, double* p_S_nxt, double* p_V_cur, double* p_V_nxt,
                    double dWS,  double h, double* p_other_arg)
    {
        *p_S_nxt = hes_step(*p_S_cur, *p_V_cur, *p_V_nxt, dWS, h, p_other_arg);
    }


    HESStepFnT Heston::_pick_func(Scheme scheme)
    {
        HESStepFnT res;
        switch (scheme)
        {
        case EUL:
            res = [this](double S_cur, double V_cur, double V_nxt,\
                        double dWS, double h, double* p_other_arg) \
                    {return this->_step_eul(S_cur, V_cur, dWS, h);};
            break;
        case CMB:
            res = [this](double S_cur, double V_cur, double V_nxt,\
                        double dWS, double h, double* p_other_arg) \
                    {return this->_step_cmb(S_cur, V_cur, V_nxt, dWS, h, 
                                            p_other_arg[0], p_other_arg[1]);};
            break;
        }
        return res;
    }


    bool Heston::check_cond(double h, double* p_S0, double* p_other_arg, Scheme scheme)
    {
        bool res = true;
        switch (scheme)
        {
        case EUL:
            break;
        case CMB:
            break;
        case IMP_3:
            res = false;
            break;
        case IMP_4:
            res = false;
            break;
        case EXP:
            res = false;
            break;
        case TG:
            res = false;
            break;
        case QE:
            res = false;
            break;
        case EXT:
            res = false;
            break;
        }
        return res;
    }


    void Heston::gen(double *p_S_out, double *p_V_out, double h, size_t n, size_t num,
            double *p_S0, double *p_V0, Scheme hes_scheme, Scheme cir_scheme,
            double lambda, double gamma_1, double gamma_2, bool trace)
    {
        auto std_ = std::sqrt(h);
        double dW = 0.0, dB = 0.0, dWS = 0.0, dWV = 0.0;
        double gamma[] = {gamma_1, gamma_2};
        assert(check_cond(h, p_S0, gamma, hes_scheme));
        assert(_vol.check_cond(h, p_V0, &lambda, cir_scheme));

        auto step_func_hes = _pick_func(hes_scheme);
        auto step_func_cir = _vol._pick_func(cir_scheme);
        for(size_t i=0; i<num; ++i)
        {
            auto p_S_tmp = trace ? p_S_out+i*(n+1) : p_S_out+i;
            auto p_V_tmp = trace ? p_V_out+i*(n+1) : p_V_out+i;
            p_S_tmp[0] = p_S0[i];
            p_V_tmp[0] = p_V0[i];
            for(size_t j=0; j<n; ++j)
            {
                dW = _gauss.gen() * std_;
                dB = _gauss.gen() * std_;
                dWS = dW;
                dWV = _rho * dW + _rhocmp * dB;
                if (trace)
                {
                    _vol._step(step_func_cir, p_V_tmp+j, p_V_tmp+j+1, dWV, h, &lambda);
                    _step(step_func_hes, p_S_tmp+j, p_S_tmp+j+1, p_V_tmp+j, p_V_tmp+j+1, dWS, h, gamma);
                }
                else
                {
                    double V_crt = *p_V_tmp;
                    _vol._step(step_func_cir, p_V_tmp, p_V_tmp, dWV, h, &lambda);
                    _step(step_func_hes, p_S_tmp, p_S_tmp, &V_crt, p_V_tmp, dWS, h, gamma);
                }
            }
        }
    }


    void Heston::operator()(double *p_S_out, double *p_V_out, double h, size_t n, size_t num,
            double *p_S0, double *p_V0, Scheme hes_scheme, Scheme cir_scheme,
            double lambda, double gamma_1, double gamma_2, bool trace)
    {
        gen(p_S_out, p_V_out, h, n, num, p_S0, p_V0, hes_scheme, cir_scheme, 
            lambda, gamma_1, gamma_2, trace);
    }
}