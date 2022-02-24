/**
 * @file test_step.cpp
 * @author  bx.zh
 * @brief Used to test the functions in class CIR.
 * 
 * @copyright Copyleft (c) all wrongs reserved. 
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <chrono>
#include "../include/model.hpp"


int main(int argc, char* argv[])
{
    // std::shared_ptr<double> p_X(new double[101]);
    double *p_X_3 = new double[101];
    double *p_X_4 = new double[101];
    double *p_X_e0 = new double[101];
    double *p_X_e1 = new double[101];

    double T = 1.0;
    size_t n = 10;
    double sigma=1.0, a=1.0, k=1.0, x0=1.0, h = T / n, lambda=0.5;
    double h_sqrt = std::sqrt(h);

    cir::CIR model_test(k, a, sigma);
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());
    std::normal_distribution<double> gauss_distr(0.0, h_sqrt);

    std::cout << "Test with trace" << std::endl;
    p_X_3[0] = x0;
    p_X_4[0] = x0;
    p_X_e0[0] = x0;
    p_X_e1[0] = x0;
    for(size_t i=0; i<n; ++i)
    {
        double dW = gauss_distr(generator);
        model_test._step_imp_3(p_X_3+i, p_X_3+i+1, dW, h);
        model_test._step_imp_4(p_X_4+i, p_X_4+i+1, dW, h);
        model_test._step_exp(p_X_e0+i, p_X_e0+i+1, dW, h, 0);
        model_test._step_exp(p_X_e1+i, p_X_e1+i+1, dW, h, lambda);
    }
    for(size_t i=0; i<n; ++i)
    {
        std::cout << std::fixed;
        std::cout << std::setprecision(5) << p_X_3[i] << "  " << p_X_4[i] << "  " \
                  << p_X_e0[i] << " " << p_X_e1[i] << std::endl;
    }


    std::cout << std::endl;
    std::cout << "Test without trace" << std::endl;
    p_X_3[0] = x0;
    p_X_4[0] = x0;
    for(size_t i=0; i<n; ++i)
    {
        std::cout << std::setprecision(5) << p_X_3[0] << "  " << p_X_4[0] << "  " \
                  << p_X_e0[0] << " " << p_X_e1[0] << std::endl;
        double dW = gauss_distr(generator);
        model_test._step_imp_3(p_X_3, p_X_3, dW, h);
        model_test._step_imp_4(p_X_4, p_X_4, dW, h);
    }
}