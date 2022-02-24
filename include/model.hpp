/**
 * @file model.hpp
 * @author bx.zh 
 * @brief Defined the CIR simulation model here.
 * @date 2022-02-15
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */

#pragma once
#include <memory>
#include "simu_type.hpp"
#include "random.hpp"

namespace cir
{
    class CIR
    {
        private:
            double _k, _a, _sigma;
            double _nu, _sigma_sqr;
            Gaussian _gauss;

        // protected:
        public:
            /**
             * @brief Used to execute one step simulation of implicite-3 method.
             * @param p_Xcrt Pointer to the value of X_{i}. 
             * @param p_Xnxt Pointer to the value of X_{i+1}. 
             * @param dW The value of W_{i+1} - W_{i}.
             * @param h The value of one step size.
             */
            void _step_imp_3(double* p_X_crt, double* p_X_nxt, double dW, double h);
            /**
             * @brief Used to execute one step simulation of implicite-4 method.
             * @param p_Xcrt Pointer to the value of X_{i}. 
             * @param p_Xnxt Pointer to the value of X_{i+1}. 
             * @param dW The value of W_{i+1} - W_{i}.
             * @param h The value of one step size.
             */
            void _step_imp_4(double* p_X_crt, double* p_X_nxt, double dW, double h);
            void _step_exp(double* p_X_crt, double* p_X_nxt, double dW, double h, double lambda);

        public:
            CIR(double k, double a, double sigma);
            CIR(double k, double a, double sigma, unsigned seed);
            double psi(double k, double t);
            void psi(double * p_k, double *p_t, double *p_out);

            void operator()(double* p_out, double T, size_t n, size_t num, double x0, SimuScheme scheme);
            // void operator()(double* p_out, double h, size_t n, size_t num, double x0, SimuScheme scheme);

            /**
             * @brief Used to generate the process, with or without trace.
             * 
             * @param p_out Pointer to the output buffer. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
             * @param h The one step size.
             * @param n The number of steps.
             * @param num The number of samples to simulate.
             * @param p_x0 The pointer to initial values. It shall of size n. 
             * @param scheme Simulation scheme to use, shall be IMP_3, IMP_4 or EXP.
             * @param lambda The lambda param for scheme EXP.
             * @param trace Whether keep or not the trace of process.
             */
            void gen(double* p_out, double h, size_t n, size_t num, double* p_x0, SimuScheme scheme, double lambda=0.0, bool trace=false);
    };
}