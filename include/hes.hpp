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
         * @brief Used to generate the process, with or without trace.
         * @param p_S_out Pointer to the output buffer of S value. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param p_V_out Pointer to the output buffer of V value. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param h The one step size.
         * @param n The number of steps.
         * @param num The number of samples to simulate.
         * @param p_S0 The pointer to initial values of S. It shall of size num. 
         * @param p_V0 The pointer to initial values of V. It shall of size num. 
         * @param hes_scheme Simulation scheme to use, shall be EUL or CMB.
         * @param cir_scheme Simulation scheme to use, shall be EUL, IMP_3, IMP_4, EXP, TG, QE, or EXT.
         * @param lambda The lambda param for scheme QE or EXP.
         * @param trace Whether keep or not the trace of process.
         */
        void gen(double *p_S_out, double *p_V_out, double h, size_t n, size_t num,
                double *p_S0, double *p_V0, Scheme hes_scheme, Scheme cir_scheme,
                double lambda=0.0, double gamma_1=0.5, double gamma_2=0.5, bool trace=false);


        /**
         * @brief Used to generate the process, with or without trace.
         * @param p_S_out Pointer to the output buffer of S value. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param p_V_out Pointer to the output buffer of V value. If keep trace, it shall have size at least num * (n+1). If not, at least num. 
         * @param h The one step size.
         * @param n The number of steps.
         * @param num The number of samples to simulate.
         * @param p_S0 The pointer to initial values of S. It shall of size num. 
         * @param p_V0 The pointer to initial values of V. It shall of size num. 
         * @param hes_scheme Simulation scheme to use, shall be EUL or CMB.
         * @param cir_scheme Simulation scheme to use, shall be EUL, IMP_3, IMP_4, EXP, TG, QE, or EXT.
         * @param lambda The lambda param for scheme QE or EXP.
         * @param trace Whether keep or not the trace of process.
         */
        void operator()(double *p_S_out, double *p_V_out, double h, size_t n, size_t num,
                double *p_S0, double *p_V0, Scheme hes_scheme, Scheme cir_scheme,
                double lambda=0.0, double gamma_1=0.5, double gamma_2=0.5, bool trace=false);


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
                    double dWS,  double h, double* p_other_arg=nullptr);

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