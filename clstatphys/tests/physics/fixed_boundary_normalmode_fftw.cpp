#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include <clstatphys/physics/fixed_boundary_normalmode_fftw.hpp>

namespace {
class FixedBoundaryNormalModeFFTWTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    precision = std::numeric_limits<double>::epsilon();
    num_particles = 1e+2; 
    normalize_constant = std::sqrt(2.0 / (num_particles+1));
    z.resize(2 * num_particles);
    z_k.resize(2 * num_particles);
    z_inverse.resize(2 * num_particles);
    double *x = &z[0];
    double *p = &z[num_particles];
    for(int i = 0 ; i < num_particles; ++i){
      x[i] = normalize_constant * std::sin( 6.0 * M_PI * (i + 1) / (num_particles + 1));
      p[i] = normalize_constant * std::sin( 3.0 * M_PI * (i + 1) / (num_particles + 1)); 
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

TEST_F(FixedBoundaryNormalModeFFTWTest, InverseTest){
  FixedBoundaryNormalModeFFTW(z, z_k);
  FixedBoundaryNormalModeFFTW(z_k, z_inverse);
  double *x = &z[0];
  double *p = &z[num_particles];
  double *x_inverse = &z_inverse[0];
  double *p_inverse = &z_inverse[num_particles];
  for(int i = 0 ; i < num_particles; ++i){
    EXPECT_NEAR(x[i], x_inverse[i], precision * 10);
    EXPECT_NEAR(p[i], p_inverse[i], precision * 10);
  }
}

TEST_F(FixedBoundaryNormalModeFFTWTest, DeltaFunctionDisplaceTest){
  FixedBoundaryNormalModeFFTW(z, z_k);
  double *x_k = &z_k[0];
  for(int i = 0 ; i < num_particles; ++i){
    if(i == 5) EXPECT_DOUBLE_EQ(x_k[i], 1.0);
    else EXPECT_NEAR(x_k[i], 0.0, precision * 10);
  }
}

TEST_F(FixedBoundaryNormalModeFFTWTest, DeltaFunctionMomentSinTest){
  FixedBoundaryNormalModeFFTW(z, z_k);
  double *p_k = &z_k[num_particles];
  for(int i = 0 ; i < num_particles; ++i){
    if(i == 2) EXPECT_DOUBLE_EQ(p_k[i], 1.0);
    else EXPECT_NEAR(p_k[i], 0.0, precision * 10);
  }
}
