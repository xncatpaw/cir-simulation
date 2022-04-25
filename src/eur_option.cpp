#include "../include/eur_option.hpp"
#include <complex>
#include <math.h>
#include <boost/math/quadrature/trapezoidal.hpp>


using boost::math::quadrature::trapezoidal;

namespace cir
{
    EurOption::EurOption(double S0, double K, double T, double t, double r, double lbd, double v0, double mu, double k, double a, double sigma, double rho) :
        _S0(S0), _K(K), _T(T), _t(t), _r(r), _lbd(lbd), _v0(v0), _mu(mu), _k(k), _a(a), _sigma(sigma), _rho(rho), i(0,1)
    {}

    // EurOption::EurOption(double S0, double K, double T, double t, double r, double mu, double k, double a, double sigma, double rho, unsigned seed) :
    //     _underlying(mu, k , a, sigma, rho, seed), _S0(S0), _K(K), _T(T), _t(t), _r(r), i(std::complex<double>(0,1))
    // {}
    
    // EurOption::EurOption(double S0, double K, double T, double t, double r, double mu, double rho, const CIR & other) :
    //     _underlying(mu, rho, other), _S0(S0), _K(K), _T(T), _t(t), _r(r), i(std::complex<double>(0,1))
    // {}
    
    EurOption::~EurOption() {}

    double EurOption::closed_price(double upper_bound)
    {
        return _S0 * probability(1, upper_bound) - _K * dicount_bound_price() * probability(2, upper_bound);
    }

    double EurOption::probability(int j, double upper_bound)
    {
        double error;
        double L1;
        size_t max_refinements = 12;
        // unsigned int max_depth = 0;
        double tolerance = 1e-9;
        auto lambda_integrand = [j, this](double phi) {return integrand(j, std::log(_S0), _v0, phi);};
        double integral = trapezoidal(lambda_integrand, 1e-10, upper_bound, tolerance, max_refinements, &error, &L1);
        return .5 + integral / M_PI;
    }

    double EurOption::integrand(int j, double x, double v, double phi)
    {
        std:: complex<double> aux = std::exp(-i * phi * std::log(_K)) * f(j, x, v, phi) / i / phi;
        return aux.real();
    }

    std::complex<double> EurOption::f(int j, double x, double v, double phi)
    {
        return std::exp(C(j, phi) + D(j, phi) * v + i * phi * x);
    }

    std::complex<double> EurOption::C(int j, double phi)
    {
        double b = _k + _lbd + (j - 2) * _rho * _sigma;
        double u = -j + 1.5;
        double tau = _T - _t;
        std::complex<double> aux1 = _rho * _sigma * phi * i;
        std::complex<double> d = std::sqrt( (aux1 - b) * (aux1 - b) - _sigma * _sigma * phi * (2. * u * i - phi) );
        std::complex<double> g = (b - aux1 + d) / (b - aux1 - d);
        return _r * phi * i * tau + _a / _sigma / _sigma * ( (b - aux1 + d) * tau - 2. * std::log( (1. - g * std::exp(d * tau)) / (1. - g)) );
    }

    std::complex<double> EurOption::D(int j, double phi)
    {
        double b = _k + _lbd + (j - 2) * _rho * _sigma;
        double u = -j + 1.5;
        double tau = _T - _t;
        std::complex<double> aux1 = _rho * _sigma * phi * i;
        std::complex<double> d = std::sqrt( (aux1 - b) * (aux1 - b) - _sigma * _sigma * phi * (2. * u * i - phi) );
        std::complex<double> g = (b - aux1 + d) / (b - aux1 - d);
        return (b - aux1 + d) / _sigma / _sigma * (1. - std::exp(d * tau)) / (1. - g * std::exp(d * tau)); 
    }

    double EurOption::dicount_bound_price()
    {
        double tau = _T - _t;
        return std::exp(-_r * tau);
    }
}
