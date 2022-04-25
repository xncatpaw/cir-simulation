#include <iostream>
#include <random>
#include <chrono>
#include "../include/random.hpp"
#include "../include/cir.hpp"
#include "../include/hes.hpp"
// #include <boost/math/distributions/non_central_chi_squared.hpp>


int main()
{
    static const int _N=100;
    static const int _NUM=10;

    double sigma=1.0, a=1.0, k=1.0, v0=1.0;
    double mu=0.0, rho=0.3;
    double T = 1.0;

    cir::CIR cir_model(k, a, sigma);
    cir::Heston hes_model(mu, rho, cir_model);

    double *p_S = new double[(_N+1)*_NUM];
    double *p_V = new double[(_N+1)*_NUM];
    double *p_S0 = new double[_NUM];
    double *p_V0 = new double[_NUM];

    for(auto j=0; j<_NUM; ++j)
    {
        p_S0[j] = 100;
        p_V0[j] = v0;
    }

    hes_model.gen(p_S, p_V, T/_N, _N, _NUM, p_S0, p_V0, cir::CMB, cir::EXP, 
                    0.0, 0.5, 0.5, false);
    // cir_model.gen(p_V, T/_N, _N, _NUM, p_V0, cir::EXP, 0, true);

    for(auto i=0; i<=_N; ++i)
    {
        for(auto j=0; j<_NUM; ++j)
        {
            std::cout << "(" << p_S[i*_NUM+j] << ", " << p_V[i*_NUM+j] << ") ";
        }
        std::cout << std::endl;
    }


    delete[] p_S;
    delete[] p_V;
    delete[] p_S0;
    delete[] p_V0;
    return 0;
}