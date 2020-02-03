#ifndef NORMALMODE_ENERGY_FFTW_HPP
#define NORMALMODE_ENERGY_FFTW_HPP

#include <cmath>
#include <vector>
#include <clstatphys/physics/fixed_boundary_normalmode_fftw.hpp>

void NormalModeEnergyFFTW(std::vector<double> const & z, std::vector<double>& Ek){
  int num_particles = z.size() / 2;
  const double *x = &z[0];
  const double *p = &z[num_particles];

  std::vector<double> z_k(2 * num_particles);

  FixedBoundaryNormalModeFFTW(z, z_k);

  const double *x_k = &z_k[0];
  const double *p_k = &z_k[num_particles];

  double omega_square;
  for(int i = 0 ; i < num_particles ; ++i){
    omega_square = 2.0 * ( 1.0 - std::cos(M_PI * (i+1) / (num_particles+1)));
    Ek[i] = 0.5 * ( p_k[i] * p_k[i] + omega_square * x_k[i] * x_k[i] );
  }

}

#endif //NORMALMODE_ENERGY_HPP
