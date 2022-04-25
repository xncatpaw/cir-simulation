#pragma once
#include <complex>

namespace cir
{
    class EurOption
    {
        private:
        double _S0, _K, _T, _t, _r;
        double _lbd;
        double _v0, _mu, _k, _a, _sigma, _rho;
        double _hi_bound;
        const std::complex<double> i;

        public:
        EurOption(double S0, double K, double T, double t, double r, double lbd, double v0, double mu, double k, double a, double sigma, double rho=0.0);
        // EurOption::EurOption(double S0, double K, double T, double t, double r, double mu, double k, double a, double sigma, double rho, unsigned seed);
        // EurOption::EurOption(double S0, double K, double T, double t, double r, double mu, double rho, const CIR &);
        ~EurOption();
        double closed_price(double upper_bound);
        double probability(int j, double upper_bound);
        double integrand(int j, double x, double v, double phi);
        std::complex<double> f(int j, double x, double v, double phi);
        std::complex<double> C(int j, double phi);
        std::complex<double> D(int j, double phi);
        double dicount_bound_price();
    };
}