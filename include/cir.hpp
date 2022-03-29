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
#include "simu_type.hpp"
#include "random.hpp"

namespace cir
{
    // typedef void (CIR::*StepFn)(double* p_X_crt, double* p_X_nxt, double dW, double h, double* p_other_arg);
    // typedef std::function<void (CIR, double*, double*, double, double, double*)> StepFnT;
    typedef std::function<double (double X_crt, double dW, double h, double*)> StepFnT;
    class CIR
    {
    private:
        // CIR eq has two forms:
        //  dX_t = (a - kX_t)dt + \sigma \sqrt(X_t) dW_t,
        //  dX_t = k(\theta - X_t)dt + \sigma\sqrt(X_t) dW_t.
        double _k, _a, _theta, _sigma; 
        double _nu, _sigma_sqr;
        Gaussian _gauss;

    public:

    public:
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
         * @param lambda The lambda param for scheme EXP.
         * @param trace Whether keep or not the trace of process.
         */
        void operator()(double* p_out, double T, size_t n, size_t num, double* p_x0, CIRScheme scheme, 
                        double lambda=0.0, bool trace=false);
        // void operator()(double* p_out, double h, size_t n, size_t num, double x0, SimuScheme scheme);

        /**
         * @brief Used to generate the process, with or without trace.
         * @param p_out Pointer to the output buffer. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param h The one step size.
         * @param n The number of steps.
         * @param num The number of samples to simulate.
         * @param p_x0 The pointer to initial values. It shall of size n. 
         * @param scheme Simulation scheme to use, shall be IMP_3, IMP_4 or EXP.
         * @param lambda The lambda param for scheme EXP.
         * @param trace Whether keep or not the trace of process.
         */
        void gen(double* p_out, double h, size_t n, size_t num, double* p_x0, CIRScheme scheme, double lambda=0.0, bool trace=false);

        /**
         * @brief Used to execute one step simulation of implicite-3 method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_imp_3(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of implicite-3 method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_imp_4(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of implicite-3 method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_tg(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of implicite-3 method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @return The value of X_{i+1}.
         */
        double _step_qe(double X_crt, double dW, double h);

        /**
         * @brief Used to execute one step simulation of explicit method.
         * @param Xcrt The value of X_{i}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @param lambda The lambda param in explicit method. 0 <= lambda <= a - sigma^2/4. 
         * @return The value of X_{i+1}.
         */
        double _step_exp(double X_crt, double dW, double h, double lambda);

        /**
         * @brief Used to execute one step simulation with the specified method.
         * @param p_Xcrt Pointer to the value of X_{i}. 
         * @param p_Xnxt Pointer to the value of X_{i+1}. 
         * @param dW The value of W_{i+1} - W_{i}.
         * @param h The value of one step size.
         * @param lambda The lambda param in explicit method. 0 <= lambda <= a - sigma^2/4. 
         */
        void _step(StepFnT func, double* p_X_crt, double* p_X_nxt, double dW, double h, double lambda=0.0);

        /**
         * @brief Used to pick the selected step function.
         * 
         * @param scheme The scheme to simulate, could be IMP_3, IMP_4, TG or EXP, defined in CIRScheme.
         * @return StepFnT The function used to execute.
         */
        StepFnT _pick_func(CIRScheme scheme);

        /**
         * @brief Used to check whether the condition hypothesis is satisfied.
         * 
         * @param h 
         * @param scheme 
         * @return true
         * @return false 
         */
        bool check_cond(double h, double* p_x0, double lambda, CIRScheme scheme);
};
}