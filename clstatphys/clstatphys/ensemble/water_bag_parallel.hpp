#ifndef ENSEMBLE_WATER_BAG_PARALLEL_HPP
#define ENSEMBLE_WATER_BAG_PARALLEL_HPP

#include <omp.h>
#include <cmath>
#include <string>
#include <random>
#include <memory>

namespace ensemble{

class WaterBagParallel{
public:
  static std::string name() { return "WaterBag Distribution : Parallelized"; }
  WaterBagParallel(double T, int num) : T_(T),num_(num){}
  template <class Rand>
  void set_initial_state(std::vector<double>& z, std::vector<std::shared_ptr<Rand> > & mts) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    double K = 0.0;
    double P = 0.0;
    double *x = &z[0];
    double *p = &z[num_];

    #pragma omp parallel for 
    for(int i = 0; i < num_ ; ++i){
      x[i] = 0.0;
    }

    #pragma omp parallel
    {
    int tid = omp_get_thread_num();
    #pragma omp for
    for(int i = 0; i < num_ ; ++i){
      auto  mt = *mts[tid];
      std::uniform_real_distribution<> p_Rand(0.0,1.0);
      p[i] = p_Rand(mt);
    }
    }

    #pragma omp parallel for reduction(+:P)
    for(int i = 0; i < num_ ; ++i){
      P += p[i];
    }
    #pragma omp parallel for reduction(+:K)
    for(int i = 0; i < num_ ; ++i){
      p[i] -= P / (double)num_;
      K +=  p[i] * p[i];
    }

    K = 0.5 * K / (double)num_;
    #pragma omp parallel for 
    for(int i = 0; i < num_ ; ++i){
      p[i] = p[i] / std::sqrt(2.0 * K / T_);
    }
  }

private:
  int num_;
  double T_;
}; 

} //end namespace

#endif //ENSEMBLE_WATER_BAG_PARALLEL_HPP
