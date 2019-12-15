#ifndef EXACTSOLUTION_PARTITION_FUNCTION_CLASSICAL_XY_CHAIN_HPP
#define EXACTSOLUTION_PARTITION_FUNCTION_CLASSICAL_XY_CHAIN_HPP

#include <cmath>
#include <string>
#include <boost/math/special_functions/bessel.hpp>
#include <functional>


namespace exactsolution{

class PartitionFunctionClassicalXYChain{
public :
  static std::string name() { return "Exact partition function , classical XY, open boundary"; }
  PartitionFunctionClassicalXYChain(double J,unsigned int num) : J_(J), num_(num) {}
    
  double operator()(double T) const {
    double value = 0.0;
    value = 2.0 * M_PI * boost::math::cyl_bessel_i(0,J_/T);
    value = 2.0 * M_PI * std::pow(value , num_-1);
    value += std::pow( std::sqrt(2.0 * M_PI * T), num_);
    return value;
  }

  double logZ(double T) const {
    double value = std::log(2.0 * M_PI);
    value += (double)(num_ - 1) * std::log(2.0 * M_PI * boost::math::cyl_bessel_i(0,J_/T));
    value += (double)(num_) * 0.5 * std::log(2.0 * M_PI * T );
    return value;
  }

private :
  double J_;
  unsigned int num_;
};  

} //end namespce

#endif //EXACTSOLUTION_PARTITION_FUNCTION_CLASSICAL_XY_CHAIN_HPP
