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
        Gaussian _gauss;
        double _mu, _rho, _rhocmp;




    public:
        Heston(double mu, double k, double a, double sigma, double rho=0.0);
        Heston(double mu, double k, double a, double sigma, double rho, unsigned seed);
        Heston(double mu, double rho, const CIR &);
        ~Heston();



    };
    
}