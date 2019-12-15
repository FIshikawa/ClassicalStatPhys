#ifndef NORMALMODE_ENERGY_PERIODIC_FFTW_HPP
#define NORMALMODE_ENERGY_PERIODIC_FFTW_HPP

#include <cmath>
#include <vector>
#include <physics/periodic_boundary_normalmode_fftw.hpp>

void NormalModeEnergyPeriodicFFTW(std::vector<double>& z, std::vector<double>& Ek){
  int num_particles = z.size() / 2;
  const double *x = &z[0];
  const double *p = &z[num_particles];

  std::vector<double> z_k(2 * num_particles);

  PeriodicBoundaryNormalModeFFTW(z, z_k);

  const double *x_k = &z_k[0];
  const double *p_k = &z_k[num_particles];

  double omega_square;
  for(int i = 0 ; i < num_particles ; ++i){
    omega_square = 2.0 * ( 1.0 - std::cos(2.0 * M_PI * i / num_particles));
    Ek[i] = 0.5 * ( p_k[i] * p_k[i] + omega_square * x_k[i] * x_k[i] );
  }

}

#endif //NORMALMODE_ENERGY_HPP
