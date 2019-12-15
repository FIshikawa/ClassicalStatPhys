#ifndef ENSEMBLE_WATER_BAG_HPP
#define ENSEMBLE_WATER_BAG_HPP

#include <cmath>
#include <string>
#include <random>

namespace ensemble{

class WaterBag {
public:
  static std::string name() { return "WaterBag Distribution"; }
  WaterBag(double T, int num) : T_(T), num_(num){}

  template <class Rand>
  void set_initial_state(std::vector<double>& z, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::uniform_real_distribution<> p_Rand(0.0,1.0);
    double K = 0.0;
    double P = 0.0;
    double *x = &z[0];
    double *p = &z[num_];
    for(int i = 0; i < num_ ; ++i) x[i] = 0.0;
    for(int i = 0; i < num_ ; ++i){
      p[i] = p_Rand(mt);
      P += p[i];
    }
    for(int i = 0; i < num_ ; ++i){
      p[i] -= P / (double)num_;
      K +=  p[i] * p[i];
    }
    K = 0.5 * K / (double)num_;
    if(K != 0.0) for(int i = 0; i < num_ ; ++i) p[i] = p[i] / std::sqrt(2.0 * K / T_);
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
  double T_;
  int num_;
};

}//end namespace

#endif //ENSEMBLE_WATER_BAG_HPP
