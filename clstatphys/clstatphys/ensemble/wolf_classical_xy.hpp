#ifndef ENSEMBLE_WOLF_CLASSICAL_XY_HPP
#define ENSEMBLE_WOLF_CLASSICAL_XY_HPP

#include <cmath>
#include <string>
#include <random>
#include <stack>

namespace ensemble{

class WolfClassicalXY {
public:
  static std::string name() { return "classical XY model equilibrium , Single cluster update, Wolf"; }
  WolfClassicalXY(double J, double T, int num) : J_(J), T_(T), num_(num) {}

  template <class Rand, class F>
  void montecalro(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::uniform_int_distribution<> dist(0,num_-1);
    std::uniform_real_distribution<> Spin_Rand(0,2*M_PI);
    std::uniform_real_distribution<> realdist(0,1.0);
    std::stack<int> spin_checklist;
    double *spin = &z[0];
    int check_tag = dist(mt);
    double phi = Spin_Rand(mt); //project angle
    int spin_original = sign(std::cos(spin[check_tag] - phi));
    spin[check_tag] = -spin[check_tag] + 2.*phi + M_PI;
    spin_checklist.push(check_tag);
    int cluster_size = 0;
    while(!spin_checklist.empty()){
      ++cluster_size;
      int center_tag = spin_checklist.top();
      double spin_center = std::cos(spin[center_tag] - phi);
      spin_checklist.pop();
      int Nd_t = f.Nd();
      for(int k = 0; k < Nd_t; ++k){
        int nearest_tag = f.table(center_tag,k);
        double spin_nearest = std::cos(spin[nearest_tag] - phi);
        double J_int = J_ * std::fabs(spin_nearest * spin_center);
        spin_nearest = sign(spin_nearest);
        double P_accept = 1 - std::exp(-2.*J_int/T_);
        if(spin_nearest == spin_original && realdist(mt) < P_accept){
          spin_checklist.push(nearest_tag);
          spin[nearest_tag] = -spin[nearest_tag] + 2.*phi + M_PI;
        }
      }
    } 
    counter += cluster_size;
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
  int sign(double hoge)const{ return (hoge > 0) - (hoge < 0);} 
  int num_;
  double J_, T_;
}; //equil end

} //end namespace

#endif //ENSEMBLE_WOLF_CLASSICAL_XY_HPP
