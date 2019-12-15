#ifndef ENSEMBLE_NORMALMODE_ENSEMBLE_PERIODIC_BOUNDARY_HPP
#define ENSEMBLE_NORMALMODE_ENSEMBLE_PERIODIC_BOUNDARY_HPP

#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <clstatphys/physics/periodic_boundary_normalmode_fftw.hpp>

namespace ensemble{

class NormalModeEnsemblePeriodicBoundaryFFTW{
public:
  static std::string name() { return "Nomal Mode Ensembe under Periodic Boundary"; }
  NormalModeEnsemblePeriodicBoundaryFFTW(int num, int k_initial, int N_k, double E) : 
    num_(num), z_k_(2*num), k_initial_(k_initial), N_k_(N_k/2), Ek_(E/N_k) {}
  template <class Rand>
  void set_initial_state(std::vector<double>& z, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::uniform_real_distribution<> uniform_rand(0.0,2.0*M_PI);
    double *p_k = &z_k_[num_];
    double *x_k = &z_k_[0];
    double *p = &z[num_];
    double *x = &z[0];
    for(int i = 0; i < k_initial_; ++i){
      x_k[i] = 0.0;
      p_k[i] = 0.0;
    }
    for(int i = k_initial_ ; i < k_initial_ + N_k_; ++i){ 
      double phi = uniform_rand(mt);
      p_k[i] = std::cos(phi) * std::sqrt(2.0 * Ek_);
      double omega = 2.0 * std::sin(i * M_PI /num_);
      x_k[i] = -1.0 * std::sin(phi) * std::sqrt(2.0 * Ek_) / omega ;
    }
    for(int i = k_initial_ + N_k_; i < num_ - N_k_ - (k_initial_ - 1); ++i){
      p_k[i] = 0.0;
      x_k[i] = 0.0;
    }
    for(int i = num_  - N_k_ - (k_initial_ - 1); i < num_ - (k_initial_ - 1); ++i){
      double phi = uniform_rand(mt);
      p_k[i] = std::cos(phi) * std::sqrt(2.0 * Ek_);
      double omega = 2.0 * std::sin(i * M_PI /num_);
      x_k[i] = -1.0 * std::sin(phi) * std::sqrt(2.0 * Ek_) / omega ;
    }
    for(int i = num_ - (k_initial_ - 1) ; i < num_ ; ++i){
      p_k[i] = 0.0;
      x_k[i] = 0.0;
    }
    PeriodicBoundaryNormalModeFFTW(z_k_, z);
  }
private:
  int num_;
  int k_initial_;
  int N_k_;
  double Ek_;
  mutable std::vector<double> z_k_;
};

} // namespace initializer

#endif // ENSEMBLE_NORMAL_MODE_ENSEMBLE_HPP
