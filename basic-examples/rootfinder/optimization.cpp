#include <iomanip>
#include <iostream>
#include <clstatphys/rootfinder/bisection_method.hpp>
#include <clstatphys/rootfinder/false_position_method.hpp>
#include <clstatphys/rootfinder/newton_method.hpp>

double f(double x) { return 3.293*x*x-5.33*x+1; }   
double df(double x) { return 6.586*x-5.33; }   

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << std::setprecision(10);
  std::cout << "[test for optimzation.hpp]\n";
  std::cout << "[bisection method]\n";
  optimization::BisectionMethod optimizer_bi;
  int iteration;
  iteration = optimizer_bi.find_zero(f,0.0,1.0);
  if(iteration < 0){
    std::cout << "Initial enclosure failuser " << std::endl;
  } else {
    std::cout << iteration << " " << optimizer_bi.zero() << std::endl;
  }

  std::cout << "[false position method]\n";
  optimization::FalsePositionMethod optimizer_fp;
  iteration = optimizer_fp.find_zero(f,0.0,1.0);
  if(iteration < 0){
    std::cout << "Initial enclosure failuser " << std::endl;
  } else {
    std::cout << iteration << " " << optimizer_fp.zero() << std::endl;
  }

  std::cout << "[Newton method]\n";
  optimization::NewtonMethod optimizer_newton;
  iteration = optimizer_newton.find_zero(f,df,0.0);
  if(iteration < 0){
    std::cout << "Initial enclosure failuser " << std::endl;
  } else {
    std::cout << iteration << " " << optimizer_newton.zero() << std::endl;
  }


#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
