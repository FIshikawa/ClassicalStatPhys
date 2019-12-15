#ifndef INTEGRATOR_CAUCHY_EULER_HPP
#define INTEGRATOR_CAUCHY_EULER_HPP

#include <cmath>
#include <string>
#include <vector>
#include <random>

namespace integrator{

class CauchyEuler{
public:
  static std::string name() { return "Cauchy Rnadom Disspated Dynamics, Euler-Maruyama Scheme (Euler)"; }
  CauchyEuler(double T, double gamma, unsigned int dim) : dim_(dim), gamma_(gamma), T_(T), k_(dim){}

  template <class Rand, class F>
  void step(double t, double h, std::vector<double>& y, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0)
    std::cauchy_distribution<> cauchy_dist(0.0,1.0);
    unsigned int num_ = dim_ / 2;
    double *x = &y[0];
    double *v = &y[num_];
    double *fx = &k_[0];
    double *fv = &k_[num_]; 
    double W = gamma_ * T_;
    f(t,y,k_);
    for(int i = 0; i < num_; ++i){
      v[i] = v[i] + (-1.0 * gamma_ * v[i] +  fv[i] + W * cauchy_dist(mt))*h;
      x[i] = x[i] + fx[i]*h; 
    }
  }

protected:
  double gamma_;
  double T_;
  unsigned int dim_;
  mutable std::vector<double> k_;
}; 

} //end namespace

#endif //INTEGRATOR_CAUCHY_EULER_HPP
