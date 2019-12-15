#ifndef ENSEMBLE_METROPOLIS_CLASSICAL_XY_HPP
#define ENSEMBLE_METROPOLIS_CLASSICAL_XY_HPP

#include <cmath>
#include <string>
#include <random>

namespace ensemble{

class MetropolisClassicalXY {
public:
  static std::string name() { return "classical XY model , Metropolis"; }
  MetropolisClassicalXY(double J, double T, int num) : num_(num),T_(T) {} 

  template <class Rand, class F>
  void montecalro(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::uniform_int_distribution<> dist(0,num_-1);
    std::uniform_real_distribution<> Spin_Rand(0,2*M_PI);
    std::uniform_real_distribution<> realdist(0,1.0);
    int i = dist(mt);
    double dP = 0.0;
    double P_accept;
    double PastSpin = z[i];
    double NewSpin =  Spin_Rand(mt);
    z[i] = NewSpin;
    double new_dE = f.target_potential_energy(i,z,0.0);
    z[i] = PastSpin;
    double past_dE = f.target_potential_energy(i,z,0.0);
    
    dP = std::exp(-1.0/T_ * new_dE)/std::exp(-1.0/T_ * past_dE);
    P_accept = realdist(mt);
   
    if(dP > P_accept){
      counter += 1;
      z[i] = NewSpin;
    }
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
  int num_;
  double T_;
}; 

} //end namespace

#endif //ENSEMBLE_METROPOLIS_CLASSICAL_XY_HPP

