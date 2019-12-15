#ifndef ENSEMBLE_STEADY_CAUCHY_HPP
#define ENSEMBLE_STEADY_CAUCHY_HPP

#include <cmath>
#include <string>
#include <random>

namespace ensemble{

class SteadyCauchy {
public:
  static std::string name() { return "Steady Cauchy moment distribution "; }
  SteadyCauchy(double J, double T, int num) : J_(J),T_(T),num_(num){}

 template <class Rand>
 void equilibrate_velocity(std::vector<double>& z, double T,  Rand & mt) const {
   std::cauchy_distribution<> moment_dist(0.0,1.0);
   double *v = &z[num_]; 
   double v_total = 0.0;
   for(int i = 0; i < num_; ++i){
     v[i] = moment_dist(mt) * T;
     v_total += v[i];
   }
   v_total = v_total/num_;
   for(int i = 0; i < num_; ++i) v[i] -= v_total;
 }
  template <class Rand, class F>
  void montecalro(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::cauchy_distribution<> cauchy_dist(0.0,1.0);
    double *x = &z[0];
    double *v = &z[num_];
    std::vector<double> z_t(2*num_,0.0);
    double *x_t = &z_t[0];
    double *v_t = &z_t[num_];
    std::vector<double> force(2*num_);
    std::vector<double> force_t(2*num_);
    double *fx = &force[0];
    double *fv = &force[num_]; 
    double *fx_t = &force_t[0];
    double *fv_t = &force_t[num_]; 
    double gamma = 1.5;
    double dt = 0.001;
    double W = gamma * T_; 
    for(int step = 0; step < 100; ++step){
      f(0,z,force);
      for(int i = 0; i < num_; ++i){
        v_t[i] = v[i] + (-1.0 * gamma * v[i] +  fv[i] + W * cauchy_dist(mt) )*dt;
        x_t[i] = x[i] + fx[i]*dt; 
      }
      f(0,z_t,force_t);
      for(int i = 0; i < num_; ++i){
        v[i] = v[i] + ( 0.5 * (-1.0 * gamma * v[i] + -1.0 * gamma * v_t[i] + fv_t[i] + fv[i]) + W * cauchy_dist(mt) )*dt;
        x[i] = x[i] + 0.5*(fx_t[i]+fx[i])*dt; 
      }
    }
    counter += num_;
  }

private:
  double J_,T_;
  int num_;
}; 

}

#endif //ENSEMBLE_STEADY_CAUCHY_HPP
