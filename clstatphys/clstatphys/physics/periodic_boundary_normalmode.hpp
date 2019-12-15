#ifndef PERIODIC_BOUNDARY_NORMALMODE_HPP
#define PERIODIC_BOUNDARY_NORMALMODE_HPP

#include <cmath>
#include <vector>

void PeriodicBoundaryNormalMode(const std::vector<double>& z, std::vector<double>& z_k){
  int num_particles = z.size() / 2;
  const double *x = &z[0];
  const double *p = &z[num_particles];
  double *x_k = &z_k[0];
  double *p_k = &z_k[num_particles];
  double normalize_constant = std::sqrt(num_particles);

  for(int i = 0; i < num_particles; ++i){
    p_k[i] = 0;
    x_k[i] = 0;
    for(int j = 0 ; j < num_particles; ++j){
      p_k[i] += p[j] * (std::cos(2.0*M_PI*i*j/num_particles) + std::sin(2.0*M_PI*i*j/num_particles));
      x_k[i] += x[j] * (std::cos(2.0*M_PI*i*j/num_particles) + std::sin(2.0*M_PI*i*j/num_particles));
    }
    p_k[i] /= normalize_constant;
    x_k[i] /= normalize_constant;
  }

}

#endif //FIXED_BOUNDARY_NORMALMODE_HPP
