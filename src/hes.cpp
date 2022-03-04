/**
 * @file hes.cpp
 * @author bx.zh 
 * @brief Defined the Heston simulation model here.
 * @date 2022-03-04
 * 
 * @copyright Copyleft (c) 2022: all wrongs reserved.
 * 
 */

#include "../include/simu_type.hpp"
#include "../include/random.hpp"
#include "../include/cir.hpp"
#include "../include/hes.hpp"

namespace cir
{
    Heston::Heston(double mu, double k, double a, double sigma, double rho) :
            _vol(k, a, sigma), _mu(mu), 
            _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
    { }

    
    Heston::Heston(double mu, double k, double a, double sigma, double rho, unsigned seed) :
            _vol(k, a, sigma), _gauss(seed), _mu(mu), 
            _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
    { }


    Heston::Heston(double mu, double rho, const CIR & other) : _vol(other), _mu(mu), _gauss(),
            _rho(rho), _rhocmp(std::sqrt(1-rho*rho))
    { }

    Heston::~Heston() { }
}