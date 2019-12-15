#ifndef ENSEMBLE_LANGEVIN_THERMAL_BATH_HPP
#define ENSEMBLE_LANGEVIN_THERMAL_BATH_HPP

#include <cmath>
#include <string>
#include <vector>
#include <random>

namespace ensemble{

class LangevinThermalBath {
public:
  static std::string name() { return "Langevin thermal bath"; }
  LangevinThermalBath(double J, double T, int num) : T_(T),num_(num) {}

  template <class Rand, class F>
  void montecalro(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::normal_distribution<> normal_dist(0.0,1.0);
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
    double D = std::sqrt( 2 * gamma * T_/ dt); 
    for(int i = 0; i < num_; ++i) v[i] = std::sqrt(T_) * normal_dist(mt);
    for(int step = 0; step < 1000; ++step){
      f(0,z,force);
      for(int i = 0; i < num_; ++i){
        v_t[i] = v[i] + (-1.0 * gamma * v[i] +  fv[i] + D * normal_dist(mt))*dt;
        x_t[i] = x[i] + fx[i]*dt; 
      }
      f(0,z_t,force_t);
      for(int i = 0; i < num_; ++i){
        v[i] = v[i] + ( 0.5 * (-1.0 * gamma * v[i] + -1.0 * gamma * v_t[i] + fv_t[i] + fv[i]) + D * normal_dist(mt))*dt;
        x[i] = x[i] + 0.5*(fx_t[i]+fx[i])*dt; 
      }
    }
    counter += num_;
  }

  template <class Rand>
  void equilibrate_velocity(std::vector<double>& z, Rand & mt, double T = T_) const {
    double beta = std::sqrt(T);
    std::normal_distribution<> normal_dist(0.0,1.0);
    double *v = &z[num_]; 
    double v_total = 0.0;
    for(int i = 0; i < num_; ++i){
      v[i] = beta * normal_dist(mt);
      v_total += v[i];
    }
    v_total = v_total/ num_;
    for(int i = 0; i < num_; ++i) v[i] -= v_total;
  }

private:
  double T_;
  int num_;

}; 

} //end namespace

#endif //ENSEMBLE_LANGEVIN_THERMAL_BATH_HPP
