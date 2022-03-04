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
    Heston::Heston(double mu, double k, double a, double sigma) : 
            _vol(k, a, sigma), _gauss_X(), _mu(mu)
    { }

    
    Heston::Heston(double mu, double k, double a, double sigma, unsigned seed_X, unsigned seed_V) :
            _vol(k, a, sigma, seed_V), _gauss_X(seed_X), _mu(mu)
    { }


    Heston::Heston(double mu, const CIR & other) : _vol(other), _mu(mu)
    { }


    Heston::~Heston() { }
}