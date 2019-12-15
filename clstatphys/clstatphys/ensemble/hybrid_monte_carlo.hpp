#ifndef ENSEMBLE_HYBRID_MONTECARLO_HPP
#define ENSEMBLE_HYBRID_MONTECARLO_HPP

#include <cmath>
#include <string>
#include <random>

#include <integrator/yoshida_4th.hpp>

namespace ensemble{

class HybridMonteCarlo{
public:
  static std::string name() { return "Hybrid Monte Carlo"; }
  HybridMonteCarlo(int num_particles, double temperture, double dt, double relax_time, int total_accept) : 
    num_particles_(num_particles), temperture_(temperture), dt_(dt), total_step_(static_cast<int>(relax_time/dt)), 
    total_accept_(total_accept), integrator_(2*num_particles) {} 

  template <class Rand, class F>
  void sample(std::vector<double>& z, F const& f, Rand & mt) const {
    int counter = 0;
    while(counter < total_accept_) montecarlo(z, counter, f, mt);
  }

  template <class Rand, class F>
  void montecarlo(std::vector<double>& z, int& counter, F const& f, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    equilibrate_velocity(z, mt);
    std::vector<double> z_t = z;
    for(int step = 0; step < total_step_; ++step) integrator_.step(0.0, dt_, z_t, f);

    double new_dE = f.energy(0.0,z_t);
    double past_dE = f.energy(0.0,z);

    double dP = std::exp(-1.0 / temperture_ * (new_dE - past_dE));

    std::uniform_real_distribution<> realdist(0,1.0);
    double P_accept = realdist(mt);
   
    if(dP > P_accept){
      z = z_t;
      counter += 1;
    }
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
  int num_particles_, total_accept_, total_step_;
  double temperture_, dt_;
  integrator::Yoshida4th integrator_;
}; 

} //end namespace

#endif //ENSEMBLE_METROPOLIS_CLASSICAL_XY_HPP

