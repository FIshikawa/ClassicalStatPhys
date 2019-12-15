#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include <physics/periodic_boundary_normalmode_fftw.hpp>

namespace {
class PeriodicBoundaryNormalModeFFTWTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    precision = std::numeric_limits<double>::epsilon();
    num_particles = 1e+2; 
    normalize_constant = std::sqrt(num_particles);
    z.resize(2 * num_particles);
    z_k.resize(2 * num_particles);
    z_inverse.resize(2 * num_particles);
    double *x = &z[0];
    double *p = &z[num_particles];
    for(int i = 0 ; i < num_particles; ++i){
      x[i] = std::cos( 6.0 * M_PI * i / num_particles) / normalize_constant ;
      p[i] = std::cos( 4.0 * M_PI * i / num_particles) / normalize_constant ; 
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

TEST_F(PeriodicBoundaryNormalModeFFTWTest, InverseTest){
  PeriodicBoundaryNormalModeFFTW(z, z_k);
  PeriodicBoundaryNormalModeFFTW(z_k, z_inverse);
  double *x = &z[0];
  double *p = &z[num_particles];
  double *x_inverse = &z_inverse[0];
  double *p_inverse = &z_inverse[num_particles];
  for(int i = 0 ; i < num_particles; ++i){
    EXPECT_NEAR(x[i], x_inverse[i], precision * num_particles);
    EXPECT_NEAR(p[i], p_inverse[i], precision * num_particles);
  }
}

TEST_F(PeriodicBoundaryNormalModeFFTWTest, DeltaFunctionDisplaceTest){
  PeriodicBoundaryNormalModeFFTW(z, z_k);
  double *x_k = &z_k[0];
  for(int i = 0 ; i < num_particles; ++i){
    if(i == 3) EXPECT_DOUBLE_EQ(x_k[i], 0.5);
    else if(i == num_particles - 3) EXPECT_DOUBLE_EQ(x_k[i], 0.5);
    else EXPECT_NEAR(x_k[i], 0.0, precision * num_particles);
  }
}

TEST_F(PeriodicBoundaryNormalModeFFTWTest, DeltaFunctionMomentSinTest){
  PeriodicBoundaryNormalModeFFTW(z, z_k);
  double *p_k = &z_k[num_particles];
  for(int i = 0 ; i < num_particles; ++i){
    if(i == 2) EXPECT_DOUBLE_EQ(p_k[i], 0.5);
    else if(i == num_particles - 2) EXPECT_DOUBLE_EQ(p_k[i], 0.5);
    else EXPECT_NEAR(p_k[i], 0.0, precision * num_particles);
  }
}
