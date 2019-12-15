#ifndef FIXED_BOUNDARY_NORMALMODE_HPP
#define FIXED_BOUNDARY_NORMALMODE_HPP

#include <vector>

void FixedBoundaryNormalMode(const std::vector<double>& z, std::vector<double>& z_k){
  int num_particles = z.size() / 2;
  const double *x = &z[0];
  const double *p = &z[num_particles];
  double *x_k = &z_k[0];
  double *p_k = &z_k[num_particles];
  double normalize_constant = std::sqrt(2.0/(num_particles+1));

  for(int i = 0; i < num_particles; ++i){
    p_k[i] = 0;
    x_k[i] = 0;
    for(int j = 0 ; j < num_particles; ++j){
      p_k[i] += normalize_constant * p[j] * std::sin(M_PI*(i+1)*(j+1)/(num_particles+1));
      x_k[i] += normalize_constant * x[j] * std::sin(M_PI*(i+1)*(j+1)/(num_particles+1));
    }
  }

}

#endif //FIXED_BOUNDARY_NORMALMODE_HPP
