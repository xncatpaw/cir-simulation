#include <iostream>
#include <random>
#include <chrono>
#include "../include/random.hpp"
// #include <boost/math/distributions/non_central_chi_squared.hpp>


int main()
{
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());
    std::chi_squared_distribution<double> chi2(0.3);
    cir::NonCentralChi2 nc_chi2(0.3);

    for(auto i=0; i<10; ++i)
    {
        std::cout << nc_chi2(0.5) << std::endl;
    }
    return 0;
}