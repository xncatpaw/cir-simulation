/**
 * @file cir.hpp
 * @author bx.zh 
 * @brief Defined the CIR simulation model here.
 * @date 2022-02-15
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */

#pragma once
// #include <memory>
#include <functional>
#include <random>
#include "simu_type.hpp"
#include "random.hpp"

namespace cir
{
    // typedef void (CIR::*StepFn)(double* p_X_crt, double* p_X_nxt, double dW, double h, double* p_other_arg);
    // typedef std::function<void (CIR, double*, double*, double, double, double*)> StepFnT;
    typedef std::function<double (double X_crt, double dW, double h, double*)> CIRStepFnT;
    class CIR
    {
    private:
        // CIR eq has two forms:
        //  dX_t = (a - kX_t)dt + \sigma \sqrt(X_t) dW_t,
        //  dX_t = k(\theta - X_t)dt + \sigma\sqrt(X_t) dW_t.
        double _k, _a, _theta, _sigma; 
        double _sigma_sqr, _nu;
        Gaussian _gauss;
        NonCentralChi2 _chi_2;
        // std::chi_squared_distribution<double> _chi_2;
        // std::poisson_distribution<double> _poisson;

    private:
        double ita(double h) const;

    public:
        /**
         * @brief Construct a new CIR object
         *  dX_t = (a-kX_t)dt + sigma sqrt(X_t) dW_t.
         */
        CIR(double k, double a, double sigma);
        CIR(double k, double a, double sigma, unsigned seed);
        CIR(const CIR &);
        ~CIR();
        // double psi(double k, double t);
        // void psi(double * p_k, double *p_t, double *p_out);

        /**
         * @brief Used to generate the process, with or without trace.
         * @param p_out Pointer to the output buffer. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param T The time range.
         * @param n The number of steps.
         * @param num The number of samples to simulate.
         * @param p_x0 The pointer to initial values. It shall of size n. 
         * @param scheme Simulation scheme to use, shall be IMP_3, IMP_4 or EXP.
         * @param lambda The lambda param for scheme QE or EXP.
         *          For EXP, 0 <= lambda <= a-\sigma^2/4, 
         *          For QE, 1 <= lambda <= 2.
         * @param trace Whether keep or not the trace of process.
         */
        void operator()(double* p_out, double T, size_t n, size_t num, double* p_x0, Scheme scheme, 
                        double lambda=0.0, bool trace=false);
        // void operator()(double* p_out, double h, size_t n, size_t num, double x0, SimuScheme scheme);

        /**
         * @brief Used to generate the process, with or without trace.
         * @param p_out Pointer to the output buffer. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param h The one step size.
         * @param n The number of steps.
         * @param num The number of samples to simulate.
         * @param p_x0 The pointer to initial values. It shall of size num. 
         * @param scheme Simulation scheme to use, shall be EUL, IMP_3, IMP_4, EXP, TG, QE, or EXT.
         * @param lambda The lambda param for scheme EXP.
         * @param trace Whether keep or not the trace of process.
         */
        void gen(double* p_out, double h, size_t n, size_t num, double* p_x0, Scheme scheme, double lambda=0.0, bool trace=false);

        /**
         * @brief Used to execute one step simulation of Euler method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_eul(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of implicite-3 method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_imp_3(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of implicite-4 method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_imp_4(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of truncated Gaussian method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_tg(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of Quadratic Exponential method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_qe(double X_crt, double dW, double h, double lambda);

        /**
         * @brief Used to execute one step simulation of explicit method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @param lambda The lambda param for EXP or QE scheme.
         *          For EXP, 0 <= lambda <= a-\sigma^2/4, 
         *          For QE, 1 <= lambda <= 2.
         * @return The value of X_{i+1}.
         */
        double _step_exp(double X_crt, double dW, double h, double lambda);

        /**
         * @brief Used to execute one step simulation of exact distribution method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_ext(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation with the specified method.
         * @param p_Xcrt Pointer to the value of X_{i}. 
         * @param p_Xnxt Pointer to the value of X_{i+1}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @param p_other_arg The pointer to extra pararm(s), e.g. lambda for scheme EXP. 
         */
        void _step(CIRStepFnT func, double* p_X_crt, double* p_X_nxt, double dW, double h, double* p_other_arg=nullptr);

        /**
         * @brief Used to pick the selected step function.
         * 
         * @param scheme The scheme to simulate, could be IMP_3, IMP_4, TG or EXP, defined in Scheme.
         * @return CIRStepFnT The function used to execute.
         */
        CIRStepFnT _pick_func(Scheme scheme);

        /**
         * @brief Used to check whether the condition hypothesis is satisfied.
         * 
         * @param h 
         * @param scheme 
         * @return true
         * @return false 
         */
        bool check_cond(double h, double* p_x0, double* p_other_arg, Scheme scheme);
};
}