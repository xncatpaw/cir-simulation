/**
 * @file test.cpp
 * @brief Used to test whether the eigen dependency is satisfied.
 *        If this file cannot be complied normally, please install
 *        the eigen library.
 * @date 2022-02-07
 * 
 */
#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
 
int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
}