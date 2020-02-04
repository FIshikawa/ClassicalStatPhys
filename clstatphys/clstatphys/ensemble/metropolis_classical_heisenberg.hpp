#ifndef ENSEMBLE_METROPOLIS_CLASSICAL_HEISENBERG_HPP
#define ENSEMBLE_METROPOLIS_CLASSICAL_HEISENBERG_HPP

#include <cmath>
#include <string>
#include <random>
#include <utility>

namespace ensemble{

class MetropolisClassicalHeisenberg {
public:
  static std::string name() { return "classical Heisenberg model , Metropolis"; }
  MetropolisClassicalHeisenberg(double J, double T, unsigned int num) : num_(num),T_(T) {} 

  template <class Rand, class F>
  void montecarlo(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::uniform_int_distribution<> dist(0,num_-1);
    std::uniform_real_distribution<> Spin_Rand_x(0,2*M_PI);
    std::uniform_real_distribution<> Spin_Rand_y(-1.0,1.0);
    std::uniform_real_distribution<> realdist(0,1.0);
    double *x = &z[0];
    double *y = &z[num_];
    int i = dist(mt);
    double dP = 0.0;
    double P_accept;
    double PastSpin_x = x[i];
    double PastSpin_y = y[i];
    double NewSpin_x =  Spin_Rand_x(mt);
    double NewSpin_y =  std::acos(Spin_Rand_y(mt));
    x[i] = NewSpin_x;
    y[i] = NewSpin_y;
    double new_dE = f.target_potential_energy(i,z,0.0);
    x[i] = PastSpin_x;
    y[i] = PastSpin_y;
    double past_dE = f.target_potential_energy(i,z,0.0);
    
    dP = std::exp(-1.0/T_ * new_dE)/std::exp(-1.0/T_ * past_dE);
    P_accept = realdist(mt);
   
    if(dP > P_accept){
      counter += 1;
      x[i] = NewSpin_x;
      y[i] = NewSpin_y;
    }
  }

  template <class Rand>
  void equilibrate_velocity(std::vector<double>& z, Rand & mt, double T = T_) const {
    double beta = std::sqrt(T);
    std::normal_distribution<> normal_dist(0.0,1.0);
    double *v = &z[2*num_]; 
    double *w = &z[3*num_]; 
    double v_total = 0.0;
    double w_total = 0.0;
    for(int i = 0; i < num_; ++i){
      v[i] = beta * normal_dist(mt);
      w[i] = beta * normal_dist(mt);
      v_total += v[i];
      w_total += w[i];
    }
    v_total = v_total/ num_;
    w_total = w_total/ num_;
    for(int i = 0; i < num_; ++i) v[i] -= v_total;
    for(int i = 0; i < num_; ++i) w[i] -= w_total;
  }


private:
  int num_;
  double T_;

}; 

} //end namespace

#endif //ENSEMBLE_METROPOLIS_CLASSICAL_HEISENBERG_HPP
