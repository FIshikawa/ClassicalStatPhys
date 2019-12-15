#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include <clstatphys/physics/periodic_boundary_normalmode.hpp>

namespace {
class PeriodicBoundaryNormalModeFFTW: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    precision = std::numeric_limits<double>::epsilon();
    num_particles = 1e+2; 
    normalize_constant = std::sqrt(1.0 / num_particles);
    z.resize(2 * num_particles);
    z_k.resize(2 * num_particles);
    z_inverse.resize(2 * num_particles);
    double *x = &z[0];
    double *p = &z[num_particles];
    for(int i = 0 ; i < num_particles; ++i){
      x[i] = normalize_constant * std::cos( 6.0 * M_PI * i / num_particles);
      p[i] = normalize_constant * std::cos( 4.0 * M_PI * i / num_particles); 
    }
  }
  double precision;
  int num_particles;
  double normalize_constant;
  std::vector<double> z;
  std::vector<double> z_k;
  std::vector<double> z_inverse;
};
}

TEST_F(PeriodicBoundaryNormalModeFFTW, InverseTest){
  PeriodicBoundaryNormalMode(z, z_k);
  PeriodicBoundaryNormalMode(z_k, z_inverse);
  double *x = &z[0];
  double *p = &z[num_particles];
  double *x_inverse = &z_inverse[0];
  double *p_inverse = &z_inverse[num_particles];
  for(int i = 0 ; i < num_particles; ++i){
    EXPECT_NEAR(x[i], x_inverse[i], precision * num_particles);
    EXPECT_NEAR(p[i], p_inverse[i], precision * num_particles);
  }
}

TEST_F(PeriodicBoundaryNormalModeFFTW, DeltaFunctionDisplaceTest){
  PeriodicBoundaryNormalMode(z, z_k);
  double *x_k = &z_k[0];
  for(int i = 0 ; i < num_particles; ++i){
    if(i == 3) EXPECT_NEAR(x_k[i], 0.5, precision * 100);
    else if(i == num_particles - 3) EXPECT_NEAR(x_k[i], 0.5, precision * 100);
    else EXPECT_NEAR(x_k[i], 0.0, precision * num_particles);
  }
}

TEST_F(PeriodicBoundaryNormalModeFFTW, DeltaFunctionMomentSinTest){
  PeriodicBoundaryNormalMode(z, z_k);
  double *p_k = &z_k[num_particles];
  for(int i = 0 ; i < num_particles; ++i){
    if(i == 2) EXPECT_NEAR(p_k[i], 0.5, precision * 100);
    else if(i == num_particles - 2) EXPECT_NEAR(p_k[i], 0.5, precision * 100);
    else EXPECT_NEAR(p_k[i], 0.0, precision * num_particles);
  }
}
