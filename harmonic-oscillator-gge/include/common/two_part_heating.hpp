#ifndef ENSEMBLE_POINT_HEATING_HPP
#define ENSEMBLE_POINT_HEATING_HPP

#include <cmath>
#include <string>
#include <random>

namespace ensemble{

class TwoPartHeating{
public:
  static std::string name() { return "Two Part Heating Distribution"; }
  TwoPartHeating(int left1, int right1, int left2, int right2,
      double T1, double T2, int num_all) : 
    left1_(left1), right1_(right1), left2_(left2), right2_(right2),
    T1_(T1), T2_(T2), num_all_(num_all) {}

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
    for(int i = left1_; i < right1_; ++i) v[i] = v_rand_1(mt);
    for(int i = left2_; i < right2_; ++i) v[i] = v_rand_1(mt);

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
    for(int i = left1_; i < right1_; ++i) x[i] = x_rand_1(mt);
    for(int i = left2_; i < right2_; ++i) x[i] = x_rand_1(mt);
  }


private:
  double T1_, T2_;
  int num_all_, left1_, right1_, left2_, right2_;
};

}//end namespace

#endif //ENSEMBLE_WATER_BAG_HPP
