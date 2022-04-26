#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <tuple>
#include "../include/random.hpp"
#include "../include/cir.hpp"
#include "../include/hes.hpp"

typedef std::tuple<cir::Scheme, double, std::string> ParamType;

int main(int argc, char* argv[])
{
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_stop = std::chrono::high_resolution_clock::now();

    int _N = 100;
    int _NUM = 10000;

    // double sigma = 1.0, v0 = 0.04, k = 0.5, a = 0.02;
    double sigma=0.2, v0=0.04, k=0.5, a=0.03, T=10.0, s0=100, rho=-0.9, u=0;
    // double sigma=0.2, v0=0.04, k=0.3, a=0.03, T=15.0, s0=100, rho=-0.5, u=0;
    // double sigma=0.2, v0=0.09, k=1, a=0.09, T=5.0, s0=100, rho=-0.5, u=0;
    double h = T/_N; 
    double lambda_exp = a - (sigma*sigma/4);

    std::ofstream f_out_V("../csv/hes_V.csv");
    std::ofstream f_out_S("../csv/hes_S.csv");

    double *p_V = new double[_NUM];
    double *p_S = new double[_NUM];
    double *p_V0 = new double[_NUM];
    double *p_S0 = new double[_NUM];
    for(auto j=0; j<_NUM; ++j)
    {
        p_V0[j] = v0;
        p_S0[j] = s0;
    }

    cir::CIR cir_model(k, a, sigma);
    cir::Heston hes_model(u, rho, cir_model);

    ParamType arr_cir_parm[] = \
    {
        std::make_tuple(cir::EUL, 0.0, "EUL"),
        std::make_tuple(cir::IMP_3, 0.0, "IMP_3"),
        std::make_tuple(cir::IMP_4, 0.0, "IMP_4"),
        std::make_tuple(cir::EXP, 0.0, "EXP_0"),
        std::make_tuple(cir::EXP, lambda_exp, "EXP_1"),
        std::make_tuple(cir::QE, 1.0, "QE_1"),
        std::make_tuple(cir::QE, 1.5, "QE_1_5"),
        std::make_tuple(cir::QE, 2, "QE_2"),
        std::make_tuple(cir::EXT, 0, "EXT"),
    }; 

    ParamType arr_hes_parm[] = \
    {
        std::make_tuple(cir::EUL, 0.0, "EUL"),
        std::make_tuple(cir::CMB, 0.0, "CMB-0"),
        std::make_tuple(cir::CMB, 0.5, "CMB-0_5"),
        std::make_tuple(cir::CMB, 1, "CMB-1")
    };
    

    f_out_V << "scheme, time(us), ";     
    f_out_S << "scheme, time(us), ";     
    for(auto j=0; j<_NUM; ++j)
    {
        f_out_V << j << ", "; 
        f_out_S << j << ", "; 
    }
    f_out_V << std::endl;
    f_out_S << std::endl;


    for(auto hes_parm : arr_hes_parm)
    {
        auto hes_scheme = std::get<0>(hes_parm);
        auto gamma_1 = std::get<1>(hes_parm);
        auto gamma_2 = 1-gamma_1;
        auto hes_name = std::get<2>(hes_parm);

        for(auto cir_parm : arr_cir_parm)
        {
            auto cir_scheme = std::get<0>(cir_parm);
            auto lambda = std::get<1>(cir_parm);
            auto cir_name = std::get<2>(cir_parm);
            // if(scheme==cir::QE or scheme==cir::EXP or scheme==cir::TG)
            // { continue; }
            auto t_start = std::chrono::high_resolution_clock::now();
            hes_model.gen(p_S, p_V, h, _N, _NUM, p_S0, p_V0, 
                        hes_scheme, cir_scheme, lambda, gamma_1, gamma_2);
            auto t_stop = std::chrono::high_resolution_clock::now();
            auto dur_us = (t_stop - t_start).count();

            auto name = hes_name + "-" + cir_name;
            std::cout << name <<" " << dur_us << std::endl;
            f_out_V << name << ", " << dur_us << ", " ;
            f_out_S << name << ", " << dur_us << ", " ;
            for(auto j=0; j<_NUM; ++j)
            { 
                f_out_V << p_V[j] << ", "; 
                f_out_S << p_S[j] << ", "; 
            }
            f_out_V << std::endl;
            f_out_S << std::endl;
        }
    }


    f_out_V.close();
    f_out_S.close();
    delete[] p_V;
    delete[] p_V0;
    delete[] p_S;
    delete[] p_S0;
    return 0;
}
