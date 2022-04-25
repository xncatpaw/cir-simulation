#include <iostream>
#include <chrono>
#include "../include/eur_option.hpp"

int main()
{
    // static const int _N=100;
    // static const int _NUM=10;

    double S0 = 100.0, K = 100.0, T = 1.0, t = 0.0, r = 0.05;
    double lbd = 5.0;
    double v0 = 1.0, mu = 0.0, k = 1.0, a = 1.0, sigma = 1.0, rho = 0.3;

    cir::EurOption eur_option(S0, K, T, t, r, lbd, v0, mu, k, a, sigma, rho=rho);

    std::cout << eur_option.dicount_bound_price() << std::endl;
    std::cout << eur_option.D(1, 1) << std::endl;
    std::cout << eur_option.C(1, 1) << std::endl;
    std::cout << eur_option.f(1, 4.6, 1, 1) << std::endl;
    for (double i = 1e-10; i < 10; ++i)
    {
        std::cout << eur_option.integrand(1, 4.6, 1, i) << " ";
    }
    std::cout << std::endl;
    std::cout << eur_option.closed_price(10) << std::endl;

    return 0;
}