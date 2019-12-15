#ifndef INTEGRATOR_NERI_YOSHIDA_HPP
#define INTEGRATOR_NERI_YOSHIDA_HPP

// ref .American Journal of Physics 73, 938 (2005) 
//      and original paper Phys. Lett. A 150, 262–268 (1990)

#include <string>
#include <vector>
#include <cmath>

namespace integrator {

class Yoshida4th{
public:
  static std::string name() { return "4th order Yoshida method (Neri method)"; }
  Yoshida4th(int dim) : dim_(dim), k1_(dim), k2_(dim), k3_(dim), k4_(dim){}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const int  dim2 = dim_ / 2;
    const double h2 = h * h;
    double b1 = 1.0/(4.0 - 2.0*std::pow(2.0,1.0/3.0));  
    double b2 = (1.0 - std::pow(2.0,1.0/3.0))/(4.0 - 2.0*std::pow(2.0,1.0/3.0));  
    double c1 = 1.0 / (2.0 - std::pow(2.0,1.0/3.0));
    double c2 = 1.0 / (1.0 - std::pow(2.0,2.0/3.0));
    f(t, y, k1_);
    for (int i = 0; i < dim2; ++i){
      y[i] = y[i] +  h * c1 * k1_[i] + h2 * c1 * b1 * k1_[i+dim2];
      y[i + dim2] = y[i + dim2] + h * b1 * k1_[i + dim2];
    }
    f(t, y, k2_);
    for (int i = 0; i < dim2; ++i){
      y[i] = y[i] +  h * c2 * k2_[i] + h2 * c2 * b2 * k2_[i+dim2];
      y[i + dim2] = y[i + dim2] + h * b2 * k2_[i + dim2];
    }
    f(t, y, k3_);
    for (int i = 0; i < dim2; ++i){
      y[i] = y[i] +  h * c1 * k3_[i] + h2 * c1 * b2 * k3_[i+dim2];
      y[i + dim2] = y[i + dim2] + h * b2 * k3_[i + dim2];
    }
    f(t, y, k4_);
    for (int i = dim2; i < dim_; ++i) y[i] = y[i] + h * b1 * k4_[i];
  }
private:
  int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> k3_;
  mutable std::vector<double> k4_;
};

} // namespace integrator

#endif // INTEGRATOR_NERI_YOSHIDA_HPP
