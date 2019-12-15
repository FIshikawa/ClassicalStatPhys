#ifndef EXACTSOLUTION_PARTITION_FUNCTION_HARMONIC_CHAIN_HPP
#define EXACTSOLUTION_PARTITION_FUNCTION_HARMONIC_CHAIN_HPP

#include <cmath>
#include <string>

namespace exactsolution{

class PartitionFunctionHarmonicChain{
public :
  static std::string name() { return "Exact partition function , harmonic chain, open boundary"; }
  PartitionFunctionHarmonicChain(double J,unsigned int num) : J_(J), num_(num) {}
    
  double operator()(double T) const {
    double value = 0.0;
    value = std::sqrt(2.0 * M_PI * T / J_);
    value = std::pow(value , num_-1);
    value += std::pow( std::sqrt(2.0 * M_PI * T), num_);
    return value;
  }
  
  double logZ(double T) const {
    double value = 0.0;
    value = (double)(num_ - 1) * 0.5 * std::log(2.0 * M_PI * T / J_);
    value += (double)(num_) * 0.5 * std::log(2.0 * M_PI * T );
    return value;
  }

private :
  double J_;
  unsigned int num_;
}; 

} //end namespace

#endif //EXACTSOLUTION_PARTITION_FUNCTION_HARMONIC_CHAIN_HPP
