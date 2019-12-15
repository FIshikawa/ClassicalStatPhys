#ifndef INTEGRATOR_DISSIPATED_EULER_HPP
#define INTEGRATOR_DISSIPATED_EULER_HPP

#include <cmath>
#include <string>
#include <vector>
#include <random>

namespace integrator{

class DissipatedEuler{
public:
  static std::string name() { return "Disspated Dynamics, Euler-Maruyama Scheme (Euler)"; }
  DissipatedEuler(double T, double gamma, unsigned int dim) : dim_(dim), gamma_(gamma), T_(T), k_(dim){}

  template <class Rand, class F>
  void step(double t, double h, std::vector<double>& y, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0)
    std::normal_distribution<> normal_dist(0.0,1.0);
    unsigned int num_ = dim_ / 2;
    double *x = &y[0];
    double *v = &y[num_];
    double *fx = &k_[0];
    double *fv = &k_[num_]; 
    double W = std::sqrt( 2 * gamma_ * T_/ h); 
    f(t,y,k_);
    for(int i = 0; i < num_; ++i){
      v[i] = v[i] + (-1.0 * gamma_ * v[i] +  fv[i] + W * normal_dist(mt))*h;
      x[i] = x[i] + fx[i]*h; 
    }
  }

protected:
  double gamma_;
  double T_;
  unsigned int dim_;
  mutable std::vector<double> k_;
}; 

}//end namespace

#endif //INTEGRATOR_DISSIPATED_EULER_HPP
