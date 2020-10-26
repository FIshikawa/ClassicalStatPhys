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
  static std::string name() { return "4th order Yoshida method (Neri method)";}
  Yoshida4th(int dim) : dim_(dim), k1_(dim), k2_(dim), k3_(dim){}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const int  dim2 = dim_ / 2;
    const double h2 = h * h;
    double c1 = 1.0/(4.0 - std::pow(2.0,4.0/3.0));  
    double c2 = (1.0 - std::pow(2.0,1.0/3.0))/(4.0 - std::pow(2.0,4.0/3.0));
    double b1 = 1.0 / (2.0 - std::pow(2.0,1.0/3.0));
    double b2 = 1.0 / (1.0 - std::pow(2.0,2.0/3.0));
    for (int i = 0; i < dim2; ++i) y[i] = y[i] +  h * c1 * y[i + dim2];
    f(t, y, k1_);
    for (int i = 0; i < dim2; ++i){
      y[i + dim2] = y[i + dim2] + h * b1 * k1_[i + dim2];
      y[i] = y[i] +  h * c2 * y[i + dim2];
    }
    f(t, y, k2_);
    for (int i = 0; i < dim2; ++i){
      y[i + dim2] = y[i + dim2] + h * b2 * k2_[i + dim2];
      y[i] = y[i] +  h * c2 * y[i + dim2];
    }
    f(t, y, k3_);
    for (int i = 0; i < dim2; ++i){
      y[i + dim2] = y[i + dim2] + h * b1 * k3_[i + dim2];
      y[i] = y[i] +  h * c1 * y[i + dim2];
    }
  }
private:
  int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> k3_;
};

} // namespace integrator

#endif // INTEGRATOR_NERI_YOSHIDA_HPP
