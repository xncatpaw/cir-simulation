#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <tuple>
#include "../include/random.hpp"
#include "../include/cir.hpp"


typedef std::tuple<cir::Scheme, double, std::string> ParamType;

int main(int argc, char* argv[])
{
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_stop = std::chrono::high_resolution_clock::now();

    int _N = 100;
    int _NUM = 10000;

    // double sigma = 1.0, v0 = 0.04, k = 0.5, a = 0.02;
    double sigma=1.0, v0=1.0, k=1.0, a=1.0;
    double T = 1.0, h = T/_N;

    std::string out_path("../csv/cir_10000.csv");
    std::ofstream f_out(out_path);

    double *p_V = new double[_NUM];
    double *p_V0 = new double[_NUM];
    for(auto j=0; j<_NUM; ++j)
    {
        p_V0[j] = v0;
    }

    cir::CIR cir_model(k, a, sigma);

    ParamType arr_parm[] = \
    {
        std::make_tuple(cir::EUL, 0.0, "EUL"),
        std::make_tuple(cir::IMP_3, 0.0, "IMP_3"),
        std::make_tuple(cir::IMP_4, 0.0, "IMP_4"),
        std::make_tuple(cir::EXP, 0.0, "EXP_0"),
        std::make_tuple(cir::EXP, 0.5, "EXP_5"),
        std::make_tuple(cir::QE, 1.0, "QE_1"),
        std::make_tuple(cir::QE, 1.5, "QE_1_5"),
        std::make_tuple(cir::QE, 2, "QE_2"),
        std::make_tuple(cir::EXT, 0, "EXT"),
    }; 
    

    f_out << "scheme, time(us), ";     
    for(auto j=0; j<_NUM; ++j)
    { f_out << j << ", "; }
    f_out << std::endl;

    for(auto parm : arr_parm)
    {
        auto scheme = std::get<0>(parm);
        auto lambda = std::get<1>(parm);
        auto name = std::get<2>(parm);
        // if(scheme==cir::QE or scheme==cir::EXP or scheme==cir::TG)
        // { continue; }
        auto t_start = std::chrono::high_resolution_clock::now();
        cir_model.gen(p_V, h, _N, _NUM, p_V0, scheme, lambda);
        auto t_stop = std::chrono::high_resolution_clock::now();
        auto dur_us = (t_stop - t_start).count();

        std::cout << name <<" " << dur_us << std::endl;
        f_out << name << ", " << dur_us << ", " ;
        for(auto j=0; j<_NUM; ++j)
        { f_out << p_V[j] << ", "; }
        f_out << std::endl;
    }

    t_start = std::chrono::high_resolution_clock::now();
    cir_model.gen(p_V, T, 1, _NUM, p_V0, cir::EXT);
    t_stop = std::chrono::high_resolution_clock::now();
    auto dur_us = (t_stop - t_start).count();

    std::cout << "EXT-1" <<" " << dur_us << std::endl;
    f_out << "EXT-1" << ", " << dur_us << ", " ;
    for(auto j=0; j<_NUM; ++j)
    { f_out << p_V[j] << ", "; }
    f_out << std::endl;

    f_out.close();
    delete[] p_V;
    delete[] p_V0;
    return 0;
}