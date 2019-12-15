#ifndef FIXED_BOUNDARY_NORMALMODE_FFTW_HPP
#define FIXED_BOUNDARY_NORMALMODE_FFTW_HPP

#include <vector>
#include <fftw3.h>

void FixedBoundaryNormalModeFFTW(const std::vector<double>& z, std::vector<double>& z_k){
  int num_particles = z.size() / 2;
  const double *x = &z[0];
  const double *p = &z[num_particles];
  double *x_k = &z_k[0];
  double *p_k = &z_k[num_particles];

  double *p_in, *p_out, *x_in, *x_out;
  fftw_plan plan_momentum, plan_displace;


  x_in  = (double *) fftw_malloc(sizeof(double) * num_particles);
  x_out = (double *) fftw_malloc(sizeof(double) * num_particles);
  p_in  = (double *) fftw_malloc(sizeof(double) * num_particles);
  p_out = (double *) fftw_malloc(sizeof(double) * num_particles);
  plan_momentum = fftw_plan_r2r_1d(num_particles, p_in, p_out, FFTW_RODFT00, FFTW_ESTIMATE);
  plan_displace = fftw_plan_r2r_1d(num_particles, x_in, x_out, FFTW_RODFT00, FFTW_ESTIMATE);


  for(int i = 0; i < num_particles; ++i){
    p_in[i] = p[i];
    x_in[i] = x[i];
  }

  fftw_execute(plan_momentum);
  fftw_execute(plan_displace);

  double normalize_constant = std::sqrt(2*(num_particles+1));
  for(int i = 0; i < num_particles; ++i){
    p_k[i] = p_out[i] / normalize_constant; 
    x_k[i] = x_out[i] / normalize_constant; 
  }

  fftw_destroy_plan(plan_momentum);
  fftw_destroy_plan(plan_displace);
  fftw_free(p_in);
  fftw_free(p_out);
  fftw_free(x_in);
  fftw_free(x_out);
}

#endif //FIXED_BOUNDARY_NORMALMODE_FFTW_HPP
