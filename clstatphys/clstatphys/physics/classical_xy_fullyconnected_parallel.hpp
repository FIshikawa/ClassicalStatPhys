#ifndef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_PARALLEL_HPP
#define HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_PARALLEL_HPP

#include <string>
#include <cmath>
#include <vector>
#include <omp.h>

namespace hamiltonian{

class ClassicalXYFullyConnectedParallel { 
public:
  // z[0...n-1]: position, z[n...2n-1]: momentum
  static std::string name() { return "Classical XYspin Fully Connected Dynamics Parallelized : J = J_input / N"; }
  ClassicalXYFullyConnectedParallel(int num_particles, double J, int N_adj) : n_(num_particles), J_(J / (double)num_particles), Nd_(N_adj){s_total_x_ = 0.0; s_total_y_ = 0.0;}

  double energy(double t, std::vector<double> const& z) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double ene = 0;

    #pragma omp parallel for reduction(+:ene)
    for (int l = 0; l < n_; ++l){
      ene += 0.5 * v[l]*v[l];
      ene += 0.5 * J_ * ( n_ - std::cos(x[l]) * s_total_x_ - std::sin(x[l]) * s_total_y_ );
    }
    return ene;
  }

  double potential_energy(double  t, std::vector<double> const& z) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double ene = 0;

    #pragma omp parallel for reduction(+:ene)
    for (int l = 0; l < n_ ; ++l){
      ene += 0.5 * J_ * ( n_ - std::cos(x[l]) * s_total_x_ - std::sin(x[l]) * s_total_y_ );
    }
    return ene;
  }

  double kinetic_energy(double /* t */, std::vector<double> const& z) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double ene = 0;

    #pragma omp parallel for reduction(+:ene)
    for (int l = 0; l < n_ ; ++l){
      ene += 0.5 * v[l]*v[l];
    }
    return ene;
  }

  virtual double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[n_];
    double ene = 0;
    ene += 0.5 * v[l]*v[l];
    return ene;
  }

  virtual double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    double ene = 0;
    ene =  1.0 * J_ * ( n_ - std::cos(x[l]) * s_total_x_ - std::sin(x[l]) * s_total_y_ );
    return ene;
  }

  // "force" calculation
  virtual void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double *fx = &force[0];
    double *fv = &force[n_];

    #pragma omp parallel for 
    for(int i = 0; i < n_ ; ++i){ 
      fx[i] = v[i];
    }

    double s_total_x = 0.0;
    double s_total_y = 0.0;

    #pragma omp parallel for reduction(+:s_total_x, s_total_y)
    for(int l = 0; l < n_; ++l){
        s_total_x += std::cos(x[l]) ;
        s_total_y += std::sin(x[l]) ;
    }

    #pragma omp parallel for 
    for(int i = 0; i < n_; ++i){
      fv[i] = - 1.0 * J_  *  (  std::sin(x[i]) * s_total_x - std::cos(x[i]) * s_total_y );
    }
  }

  void set_meanfield(double s_total_x, double s_total_y) {
    s_total_x_ = s_total_x ;
    s_total_y_ = s_total_y ;
  }

  double check_meanfield(int i){
    if(i == 0) return s_total_x_;
    else return s_total_y_;
  }

  int Nd()const { return Nd_ ;}

protected:
  int n_;
  double J_;
  int Nd_;
  double s_total_y_;
  double s_total_x_;
};

} //end namespace

#endif //HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_PARALLEL_HPP
