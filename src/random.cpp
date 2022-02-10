#include <vector>
#include <random>
#include <assert.h>
#include <chrono>
#include "../include/simu_type.hpp"
#include "../include/random.hpp"

namespace cir
{
    Gaussian::Gaussian() : generator(std::chrono::steady_clock::now().time_since_epoch().count()),
                            gauss_distr()
    { }

    void Gaussian::gen(double *p_out, std::vector<int> shape, SimuType type)
    {
        // auto ord = shape.size();
        auto tot_num = 1; // Total number of r.v. to generate.
        for(auto num : shape)
        {
            assert (num > 0);
            tot_num *= num;
        }
        
        // Choose the distribution generator.
        if (type==PRECISE)
        {
            // Generate.
            for(auto i=0; i<tot_num; ++i)
            { *(p_out+i) = gauss_distr(generator); }
        }
        else if (type==SECOND_ORDER)
        { /* code */ }
        else
        { }
        

        return;
    }
} 
