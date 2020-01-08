#ifndef ENSEMBLE_POINT_HEATING_HPP
#define ENSEMBLE_POINT_HEATING_HPP

#include <cmath>
#include <string>
#include <random>

namespace ensemble{

class PointHeating{
public:
  static std::string name() { return "Point Heating Distribution"; }
  PointHeating(double T1, double T2, int num1, int num_all) : 
    T1_(T1), num1_(num1), T2_(T2), num_all_(num_all){}

  template <class Rand>
  void set_initial_state(std::vector<double>& z, Rand & mt) const {
    equilibrate_velocity(z,mt);
    equilibrate_position(z,mt);
  }

  template <class Rand>
  void equilibrate_velocity(std::vector<double>& z, Rand & mt) const {
    std::normal_distribution<> v_rand_1(0.0,std::sqrt(T1_));
    std::normal_distribution<> v_rand_2(0.0,std::sqrt(T2_));
    double *v = &z[num_all_]; 
    double v_total = 0.0;
    for(int i = 0; i < num_all_; ++i) v[i] = v_rand_2(mt);
    for(int i = 0; i < num1_; ++i) v[i] = v_rand_1(mt);
    for(int i = 0; i < num_all_; ++i) v_total += v[i];
    v_total = v_total/ num_all_;
    for(int i = 0; i < num_all_; ++i) v[i] -= v_total;
  }

  template <class Rand>
  void equilibrate_position(std::vector<double>& z, Rand & mt) const {
    std::normal_distribution<> x_rand_1(0.0,std::sqrt(T1_));
    std::normal_distribution<> x_rand_2(0.0,std::sqrt(T2_));
    double *x = &z[0];
    for(int i = 0; i < num_all_; ++i) x[i] = x_rand_2(mt);
    for(int i = 0; i < num1_; ++i) x[i] = x_rand_1(mt);
  }


private:
  double T1_, T2_;
  int num_all_, num1_;
};

}//end namespace

#endif //ENSEMBLE_WATER_BAG_HPP
