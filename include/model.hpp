/**
 * @file model.hpp
 * @author bx.zh 
 * @brief Defined the CIR simulation model here.
 * @date 2022-02-15
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */

#pragma onece
#include "simu_type.hpp"

namespace cir
{
    class CIR
    {
        private:
            double _k, _a, _sigma;
            double _nu, _sigma_sqr;

        public:
            CIR(double k, double a, double sigma);
            double psi(double k, double t);
            void psi(double * p_k, double *p_t, double *p_out);

            void operator()(double* p_out, double T, size_t n, size_t num, double x0, SimuType type);
            void gen(double* p_out, double T, size_t n, size_t num, double x0, SimuType type);

    };
}