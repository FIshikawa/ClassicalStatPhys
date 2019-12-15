#ifndef ENSEMBLE_LANGEVIN_THERMOSTAT_HPP
#define ENSEMBLE_LANGEVIN_THERMOSTAT_HPP

#include <cmath>
#include <string>
#include <vector>
#include <random>

#include <integrator/dissipated_runge_kutta.hpp>

namespace ensemble{

class LangevinThermostat {
public:
  static std::string name() { return "Langevin thermostat"; }
  LangevinThermostat(int num_particles, double temperture, double dt, double gamma,  double relax_time) : 
    num_particles_(num_particles),dt_(dt),total_step_(static_cast<int>(relax_time/dt)),
    integrator_(2*num_particles, temperture, gamma) {}

  template <class Rand, class F>
  void sample(std::vector<double>& z, F const& f, Rand & mt) const {
    for(int step = 0; step < total_step_;) montecarlo(z, step, f, mt);
  }

  template <class Rand, class F>
  void montecarlo(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    integrator_.step(0.0, dt_, z, f, mt);
    counter += 1;
  }

  template <class Rand>
  void equilibrate_velocity(std::vector<double>& z, Rand & mt, double temperture = -1.0) const {
    if(temperture < 0.0) temperture = temperture_;
    double beta = std::sqrt(temperture);
    std::normal_distribution<> normal_dist(0.0,1.0);
    double *v = &z[num_particles_]; 
    for(int i = 0; i < num_particles_; ++i) v[i] = beta * normal_dist(mt);
  }


private:
  int num_particles_, total_step_;
  double temperture_, dt_;
  integrator::DissipatedRungeKutta integrator_;

}; 

} //end namespace

#endif //ENSEMBLE_LANGEVIN_THERMOSTAT_HPP
