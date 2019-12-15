#ifndef ENSEMBLE_HARMONIC_CHAIN_PERIODIC_BOUNDARY_HPP
#define ENSEMBLE_HARMONIC_CHAIN_PERIODIC_BOUNDARY_HPP

#include <cmath>
#include <string>
#include <random>

namespace ensemble{

class HarmonicChainPeriodicBoundary {
public:
  static std::string name() { return "Harmonic Chain Open Boundary Equilbrium"; }
  HarmonicChainPeriodicBoundary(double J, double T, int num) : J_(J),T_(T),num_(num){}

  template <class Rand>
  void set_initial_state(std::vector<double>& z, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::normal_distribution<> xdis(0.0,1.0);
    double sqT = std::sqrt(T_); //sqrt of inverse temperture
    double sqJ = std::sqrt(J_); //sqrt of inverse temperture
    double Q = 0.0;
    std::vector<double> q(num_);
    double *x = &z[0];
    for(int j = 0; j < num_ ; ++j) x[j] = 0.0;//initialize
    for(int j = 0; j < num_ ; ++j) q[j] = 1.0 * sqT / sqJ * xdis(mt); //set x
    for(int j = 0; j < num_ ; ++j) Q += q[j]; //set x
    Q = Q / num_;
    for(int j = 0; j < num_ ; ++j) q[j] -= Q; //set x
    for(int i = 1; i < num_ ; ++i) x[i] = x[i-1] + q[i-1];
  }

  template <class Rand>
  void equilibrate_velocity(std::vector<double>& z, Rand & mt, double T = -1.0) const {
    if(T < 0.0) T = T_;
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
  double J_, T_;
  int num_;
};

} //end namespace

#endif //ENSEMBLE_HARMONIC_CHAIN_PERIODIC_BOUNDARY_HPP
