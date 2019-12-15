#ifndef INTEGRATOR_DISSIPATED_RUNGEKUTTA_HPP
#define INTEGRATOR_DISSIPATED_RUNGEKUTTA_HPP

#include <cmath>
#include <string>
#include <vector>
#include <random>

namespace integrator{

class DissipatedRungeKutta{
public:
  static std::string name() { return "Disspated Dynamics, Heun Scheme (Runge-Kutta)"; }
  DissipatedRungeKutta(unsigned int dim, double temperture, double gamma) : 
    dim_(dim), gamma_(gamma), temperture_(temperture),  k1_(dim), k2_(dim),yt_(dim) {}

  template <class Rand, class F>
  void step(double t, double h, std::vector<double>& y, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0)
    std::normal_distribution<> normal_dist(0.0,1.0);
    unsigned int num_ = dim_ / 2;
    double *x = &y[0];
    double *v = &y[num_];
    double *x_t = &yt_[0];
    double *v_t = &yt_[num_];
    double *fx = &k1_[0];
    double *fv = &k1_[num_]; 
    double *fx_t = &k2_[0];
    double *fv_t = &k2_[num_]; 
    double D = std::sqrt( 2 * gamma_ * temperture_/ h); 

    f(t,y,k1_);

    for(int i = 0; i < num_; ++i){
      v_t[i] = v[i] 
               + h * ( -1.0 * gamma_ * v[i] +  fv[i] + D * normal_dist(mt) );
      x_t[i] = x[i] + fx[i] * h; 
    }

    f(t,yt_,k2_);

    for(int i = 0; i < num_; ++i){
      v[i] = v[i] 
             + h * ( 0.5 * (-1.0 * gamma_ * v[i] + -1.0 * gamma_ * v_t[i] + fv_t[i] + fv[i]) + D * normal_dist(mt) );
      x[i] = x[i] + 0.5 * ( fx_t[i] + fx[i] ) * h; 
    }

  }

protected:
  double gamma_;
  double temperture_;
  unsigned int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> yt_;
}; 

} //end namespace

#endif //INTEGRATOR_DISSIPATED_RUNGEKUTTA_HPP
