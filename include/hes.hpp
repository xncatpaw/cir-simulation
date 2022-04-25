/**
 * @file hes.hpp
 * @author bx.zh 
 * @brief Defined the Heston simulation model here.
 * @date 2022-03-04
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */

#pragma once
#include "simu_type.hpp"
#include "random.hpp"
#include "cir.hpp"

namespace cir
{
    typedef std::function<double (double S_cur, double V_cur, double V_nxt, \
                        double dWS, double h, double* p_other_arg)> HESStepFnT;
    class Heston
    {
    private:
        CIR _vol;
        Gaussian _gauss;
        double _mu, _rho, _rhocmp;

    public:
        /**
         * @brief Construct a new Heston object
         *  dS_t = S_t(mu dt + sqrt(V_t)dW^S_t),
         *  dV_t = (a-kV_t)dt  = sigma sqrt(V_t)dW^V_t.
         */
        Heston(double mu, double k, double a, double sigma, double rho=0.0);
        Heston(double mu, double k, double a, double sigma, double rho, unsigned seed);
        Heston(double mu, double rho, const CIR &);
        ~Heston();

        /**
         * @brief Used to execute one step simulation of Euler method.
         * 
         * @param S_cur The value of S_{i}.
         * @param V_cur The value of V_{i}.
         * @param dWS The value of W^S_{i+1} - W^S_{i}.
         * @param h The value of one step size.
         * @return double The value of S_{i+1}.
         */
        double _step_eul(double S_cur, double V_cur, double dWS, double h);

        /**
         * @brief Used to execute one step simulation of Combinaison method.
         * 
         * @param S_cur The value of S_{i}.
         * @param V_cur The value of V_{i}.
         * @param V_nxt The value of V_{i+1}.
         * @param dWS The value of W^S_{i+1} - W^S_{i}.
         * @param h The value of one step size.
         * @param gamma_1 The coef of V_{i}.
         * @param gamma_2 The coef of V_{i+1}.
         * @return double The value of S_{i+1}.
         */
        double _step_cmb(double S_cur, double V_cur, double V_nxt, double dWS, double h, 
                        double gamma_1, double gamma_2);

        /**
         * @brief Used to pick the selected step function.
         * 
         * @param scheme The scheme to simulate, could be IMP_3, IMP_4, TG or EXP, defined in Scheme.
         * @return StepFnT The function used to execute.
         */

        /**
         * @brief Used to execute one step simulation with indicated step func.
         * 
         * @param hes_step The Heston model step func.
         * @param p_S_cur The pointer to value of S_{i}.
         * @param p_S_nxt The pointer to value of S_{i+1}.
         * @param p_V_cur The pointer to value of V_{i}.
         * @param p_V_nxt The pointer to value of V_{i+1}.
         * @param dWS The value of W^S_{i+1} - W^S_{i}.
         * @param h The value of one step size.
         * @param p_other_arg The pointer to other arg(s).
         */
        void _step(HESStepFnT hes_step, 
                    double* p_S_cur, double* p_S_nxt, double* p_V_cur, double* p_V_nxt,
                    double dWS,  double h, double* p_other_arg);

        HESStepFnT _pick_func(Scheme scheme);

        /**
         * @brief Used to check whether the condition hypothesis is satisified.
         * 
         * @param h 
         * @param p_S0 
         * @param p_other_arg 
         * @param scheme 
         * @return true 
         * @return false 
         */
        bool check_cond(double h, double* p_S0, double* p_other_arg, Scheme scheme);
    };
    
}