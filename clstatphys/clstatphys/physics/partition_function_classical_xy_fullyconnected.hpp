#ifndef EXACTSOLUTION_PARTITION_FUNCTION_CLASSICAL_XY_FULLYCONNECTED_HPP
#define EXACTSOLUTION_PARTITION_FUNCTION_CLASSICAL_XY_FULLYCONNECTED_HPP

#include <cmath>
#include <string>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>
#include <clstatphys/rootfinder/false_position_method.hpp>
#include <clstatphys/rootfinder/newton_method.hpp>


namespace exactsolution{

class PartitionFunctionClassicalXYFullyConnected{
public :
  static std::string name() { return "Exact partition function, classical XY, Fully Connected"; }
  PartitionFunctionClassicalXYFullyConnected(double J,unsigned int num) : J_(J), num_(num) {}
    
  double operator()(double T){
    double value = 0.0;
    value = logZ(T);
    value = std::exp(value);
    return value;
  }

  double logZ(double T) {
    double value = 0.0;
    double M_temp = M(T);
    value -= 0.5 * std::log(2.0 * M_PI * T) - 0.5 * J_ / T - 0.5 * M_temp * M_temp / T * J_  + F(M_temp / T * J_);
    value *= num_;
    return value;
  }

  double T(double U){
    double T0 = 0.0;
    double T1 = 0.75;
    if(U >= 0.25 + 0.5 * J_){
      return 2.0 * U - J_;
    } else{
      //False Position Method
      optimization::FalsePositionMethod optimizer;
      double J = J_;
      auto diff_U = [&J,&U,this](double x) { return U - ( 0.5 * x + 0.5 * J * ( 1 - this->M(x) * this->M(x)));};
      optimizer.find_zero(diff_U, T0, T1);
      return optimizer.zero();
    }
  }

  double M(double T) {
    double M_tmp = 0.0;
    if(T >= J_ * 0.5) return 0.0;
    else if( T == 0.0) return 1.0;
    else{
      M_tmp = std::sqrt(1 - 2.0 * T / J_);
      optimization::NewtonMethod optimizer;
      double J = J_;
      auto diff_M = [&J,&T,this](double x) { return x - this->dF(x * J / T);};
      auto d_diff_M = [&J,&T,this](double x) { return 1 - J / T * this->ddF(x * J / T);};
      optimizer.find_zero(diff_M, d_diff_M, M_tmp);
      return optimizer.zero();
    }
  }

  double F(double value){
    return std::log(2.0 * M_PI * boost::math::cyl_bessel_i(0,value));
  }

  double dF(double value){
    double value_output = 0.0;
    value_output = boost::math::cyl_bessel_i(1.0,value) / boost::math::cyl_bessel_i(0.0,value);
    return value_output;
  }

  double ddF(double value){
    double value_output = 0.0;
    value_output = 0.5 * ( boost::math::cyl_bessel_i(0.0,value) + boost::math::cyl_bessel_i(2.0,value) ) * boost::math::cyl_bessel_i(0.0,value); 
    value_output -=  boost::math::cyl_bessel_i(1.0,value) * boost::math::cyl_bessel_i(1.0,value);
    value_output /= boost::math::cyl_bessel_i(0.0,value) * boost::math::cyl_bessel_i(0.0,value);
    return value_output;
  }

private :
  double J_;
  unsigned int num_;
  mutable double U_;
  mutable double T_;
}; //end harmonic_chain 

}
#endif //EXACTSOLUTION_PARTITION_FUNCTION_CLASSICAL_XY_FULLYCONNECTED_HPP
