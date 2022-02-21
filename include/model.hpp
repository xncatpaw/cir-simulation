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

namespace cir
{
    class CIR
    {
        private:
            double _k, _a, _sigma;
            double _nu, _sigma_sqr;

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
            // void _step_exp(std::shared_ptr<double>, size_t);

        public:
            CIR(double k, double a, double sigma);
            double psi(double k, double t);
            void psi(double * p_k, double *p_t, double *p_out);

            void operator()(double* p_out, double T, size_t n, size_t num, double x0, SimuType type);
            void gen(double* p_out, double T, size_t n, size_t num, double x0, SimuType type);
    };
}