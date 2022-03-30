/**
 * @file hes.cpp
 * @author bx.zh 
 * @brief Defined the Heston simulation model here.
 * @date 2022-03-04
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */

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
         _vol(other), _mu(mu), _gauss(), _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
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


    HESStepFnT Heston::_pick_func(HESScheme scheme)
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


    bool Heston::check_cond(double h, double* p_S0, double* p_other_arg, HESScheme scheme)
    {
        bool res = true;
        switch (scheme)
        {
        case EUL:
            break;
        case CMB:
            break;
        }
        return res;
    }
}