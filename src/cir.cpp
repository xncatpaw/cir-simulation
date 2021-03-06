/**
 * @file cir.cpp
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
#include "../include/cir.hpp"

namespace cir
{
    CIR::CIR(double k, double a, double sigma) : _k(k), _a(a), _sigma(sigma),
            _sigma_sqr(sigma*sigma), _gauss(), _chi_2(4*_a/_sigma_sqr)
    {
        _theta = _k!=0 ? _a/_k : 0;
        _nu = 4 * _a / _sigma_sqr;
    }


    CIR::CIR(double k, double a, double sigma, unsigned seed) : _k(k), _a(a), _sigma(sigma),
            _sigma_sqr(sigma*sigma), _gauss(seed), _chi_2(4*_a/_sigma_sqr)
    {
        _theta = _k!=0 ? _a/_k : 0;
        _nu = 4 * _a / _sigma_sqr;
    }


    CIR::CIR(const CIR& other) : _k(other._k), _a(other._a), _theta(other._theta),
             _sigma(other._sigma), _sigma_sqr(other._sigma_sqr), 
             _nu(other._nu),  _gauss(other._gauss), _chi_2(other._nu) 
    {}


    CIR::~CIR() { }


    double CIR::ita(double h) const
    {
        auto tmp = std::exp(-1.0 * _k * h);
        auto res = 4 * _k * tmp / (_sigma_sqr * (1-tmp));
        return res;
    }


    double CIR::_step_eul(double X_crt, double dW, double h)
    {
        auto X_nxt = (_a - _k*X_crt) * h + _sigma*std::sqrt(X_crt)*dW;
        X_nxt = X_nxt > 0.0 ? X_nxt : 0.0;
        return X_nxt;
    }


    double CIR::_step_imp_3(double X_crt, double dW, double h)
    {
        auto _1kh = _k * h + 1.0; // 1 + kh.
        auto term_sigdW = _sigma * dW;
        auto tmp = term_sigdW * term_sigdW; // \sigma^2 * dW^2.
        tmp += 4.0 * _1kh * (X_crt + (_a - _sigma_sqr/2.0) * h);
        tmp = std::sqrt(tmp) + term_sigdW;
        tmp = tmp / (2.0 * _1kh);
        auto X_nxt = (tmp * tmp);
        return X_nxt;
    }


    double CIR::_step_imp_4(double X_crt, double dW, double h)
    {
        auto _1kh = _k * h /2.0 + 1.0; // 1 + kh/2. 
        auto term_sigdW = _sigma/2.0 * dW;
        auto term_X_crt_sqrt = std::sqrt(X_crt);

        auto tmp = term_sigdW + term_X_crt_sqrt;
        tmp = tmp * tmp + 4.0 * _1kh * (_a - _sigma_sqr/4.0) / 2.0 * h;
        tmp = std::sqrt(tmp) + term_sigdW + term_X_crt_sqrt;
        tmp = tmp / (_1kh * 2.0);
        auto X_nxt = (tmp * tmp);
        return X_nxt;
    }


    double CIR::_step_tg(double X_crt, double dW, double h)
    {
        auto tmp_exp = std::exp(-h*_k);
        auto m_ = _theta + (X_crt - _theta) * tmp_exp;
        auto s2_ = X_crt*_sigma_sqr*tmp_exp*(1-tmp_exp)/_k \
            + _theta*_sigma_sqr*(1-tmp_exp)*(1-tmp_exp)/(2*_k);
        
        /* TODO */
        return 0;
    }


    double CIR::_step_qe(double X_crt, double dW, double h, double lambda)
    {
        auto tmp_exp = std::exp(-h*_k);
        auto m_ = _theta + (X_crt - _theta) * tmp_exp;
        auto s2_ = X_crt*_sigma_sqr*tmp_exp*(1-tmp_exp)/_k \
            + _theta*_sigma_sqr*(1-tmp_exp)*(1-tmp_exp)/(2*_k);
        auto psi_ = s2_ / (m_*m_);
        double res = 0.0;
        if (psi_ <= lambda)
        {
            auto inv_psi_ = m_*m_ / s2_;
            auto b2_ = 2.0*inv_psi_ - 1.0 + std::sqrt(2*inv_psi_)*std::sqrt(2*inv_psi_-1);
            auto a_ = m_ / b2_;
            auto tmp_ = std::sqrt(b2_) + dW/std::sqrt(h);
            res = a_ * tmp_ * tmp_;
        }
        else
        {
            auto p_ = (psi_-1) / (psi_+1);
            auto beta_ = (1-p_) / m_;
            auto u = _gauss.gen_unif();
            res = u<=p_ ? 0 : std::log((1-p_)/(1-u)) / beta_;
        }
        return res;
    }


    double CIR::_step_exp(double X_crt, double dW, double h, double lambda)
    {
        auto _1kh = 1.0 - _k*h/2.0;
        auto tmp = _1kh * std::sqrt(X_crt) ;
        tmp += (_sigma * dW) / (2.0 * _1kh);
        tmp = tmp * tmp;
        tmp = tmp + (_a - _sigma_sqr/4.0) * h;
        tmp = tmp + lambda * (dW*dW - h);

        auto X_nxt = tmp;
        return X_nxt;
    }


    double CIR::_step_ext(double X_crt, double dW, double h)
    {
        auto _itah = ita(h);
        auto lambda = X_crt * _itah;
        auto X_nxt = _chi_2(lambda) * std::exp(-1.0*_k*h) / _itah;
        return X_nxt;
    }


    void CIR::_step(CIRStepFnT func, double* p_X_crt, double* p_X_nxt, double dW, double h, double* p_other_arg)
    {
        *p_X_nxt = func(*p_X_crt, dW, h, p_other_arg);
    }


    void CIR::gen(double* p_out, double h, size_t n, size_t num, double* p_x0, Scheme scheme, double lambda, bool trace)
    {
        auto std_ = std::sqrt(h); // The standared deviation of dW.
        double dW = 0.0;
        assert (check_cond(h, p_x0, &lambda, scheme));
        auto step_func = _pick_func(scheme);
        for(size_t i=0; i<num; ++i)
        {
            auto p_tmp = trace ? p_out+i*(n+1) : p_out+i; // The shifted pointer to current sample.
                                                            // If trace, shall shift i*(n+1), if not, shall shift i.
            p_tmp[0] = p_x0[i]; // The first value.
            for(size_t j=0; j<n; ++j)
            {
                dW = _gauss.gen() * std_;
                if(trace)
                { _step(step_func, p_tmp+j, p_tmp+j+1, dW, h, &lambda); }
                else
                { _step(step_func, p_tmp, p_tmp, dW, h, &lambda); }
            }
        }
    }


    void CIR::operator()(double* p_out, double T, size_t n, size_t num, double* p_x0, Scheme scheme, double lambda, bool trace)
    {
        assert(T>=0 && n>0 && num>0 );
        double h = T / (double (n));
        gen(p_out, h, n, num, p_x0, scheme, lambda, trace);
    }


    CIRStepFnT CIR::_pick_func(Scheme scheme)
    {
        CIRStepFnT res;
        switch (scheme)
        {
        case EUL:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_eul(X_crt, dW, h);};
            break;
        case IMP_3:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_imp_3(X_crt, dW, h);};
            break;
        case IMP_4:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_imp_4(X_crt, dW, h);};
            break;
        case TG:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_tg(X_crt, dW, h);};
            break;
        case QE:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_qe(X_crt, dW, h, *p_other_arg);};
            break;
        case EXP:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_exp(X_crt, dW, h, *p_other_arg);};
            break;
        case EXT:
            res = [this](double X_crt, double dW, double h, double* p_other_arg) \
                    {return this->_step_ext(X_crt, dW, h);};
        default:
            break;
        }
        return res;
    }


    bool CIR::check_cond(double h, double* p_x0, double* p_other_arg, Scheme scheme)
    {
        bool res = true;
        double lambda = p_other_arg ? * p_other_arg : 0;
        switch (scheme)
        {
        case IMP_3:
            res = (_nu>0.5) && (_k>=0 || h<=(-1.0/_k));
            break;
        case IMP_4:
            res = (_nu>1) && (_k>=0 || h<=(-2.0/_k));
            break;
        case TG:
            break;
        case QE:
            // double lambda = *p_other_arg;
            res = (lambda>=1 && lambda<=2);
            break;
        case EXP:
            // double lambda = *p_other_arg;
            res = (lambda>=0 && lambda<=(_a-_sigma_sqr/4.0));
            break;
        case CMB:
            res = false;
            break;
        case EUL:
            res = true;
        case EXT:
            res = true;
        default:
            break;
        }
        return res;
    }
    // double CIR::psi(double k, double t)
    // {
    //     double res = 0;
    //     if(k==0)
    //     { res = t; }
    //     else
    //     { res = (1 - std::exp(-k*t)) / k; }
    //     return res;
    // }


    // void CIR::psi(double *p_k, double *p_t, double *p_out)
    // {
    //     double res = 0;
    //     if((*p_k)==0)
    //     { res = *(p_t); }
    //     else
    //     { res = (1 - std::exp(-(*p_k) * (*p_t))) / (*p_k); }
    //     *p_out = res;
    // }

}


