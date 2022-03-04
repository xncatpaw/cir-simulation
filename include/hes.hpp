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
    class Heston
    {
    private:
        CIR _vol;
        Gaussian _gauss_X;
        double _mu;

    public:
        Heston(double mu, double k, double a, double sigma);
        Heston(double mu, double k, double a, double sigma, unsigned seed_X, unsigned seed_V);
        Heston(double mu, const CIR &);
        ~Heston();

    };
    
}