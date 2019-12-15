#include <gtest/gtest.h>
#include <main/quasi_integrable_thermalization_divergence_functions.hpp>

namespace {
class QuasiIntegrableThermalizationDivergenceFunctionsTest : public ::testing::TestWithParam<double>{
  protected:  
  virtual void SetUp(){
    k_initial = 1;
    num_particles = 1024;
    N_normalmode = num_particles / 16;
    harmonic_spectrum.resize(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) harmonic_spectrum[i] = 0.0;
    toda_spectrum.resize(num_particles);
  }
  int N_normalmode,num_particles;
  int k_initial;
  double E_initial;
  std::vector<double> harmonic_spectrum,toda_spectrum;
};

TEST_P(QuasiIntegrableThermalizationDivergenceFunctionsTest, HarmonicTest) {
  SetHarmonicDist(harmonic_spectrum, N_normalmode, k_initial);
  double normalize_const = 0;
  for(int i = 0; i < num_particles; ++i){
    normalize_const += harmonic_spectrum[i];
    if(i < k_initial)
      EXPECT_DOUBLE_EQ(harmonic_spectrum[i],0.0);
    else if(k_initial <= i && i < k_initial + N_normalmode/2) 
      EXPECT_DOUBLE_EQ(harmonic_spectrum[i],1.0/N_normalmode);
    else if(k_initial + N_normalmode/2 <= i && i < num_particles - N_normalmode/2 - (k_initial-1))
      EXPECT_DOUBLE_EQ(harmonic_spectrum[i],0.0);
    else if(num_particles - N_normalmode/2 - (k_initial-1) <= i && i < num_particles - (k_initial-1))
      EXPECT_DOUBLE_EQ(harmonic_spectrum[i],1.0/N_normalmode);
    else
      EXPECT_DOUBLE_EQ(harmonic_spectrum[i],0.0);
  }
  EXPECT_NEAR(normalize_const,1.0,1E-8);
}

TEST_P(QuasiIntegrableThermalizationDivergenceFunctionsTest, TodaTest) {
  E_initial = GetParam() * 1024;
  int counter = SetTodaDist(toda_spectrum, E_initial);
  ASSERT_EQ(counter,1024);
  double normalize_const = 0;
  for(int i = 0; i < num_particles; ++i){
    normalize_const += toda_spectrum[i];
  }
  EXPECT_NEAR(normalize_const,1.0,1E-8);
}

INSTANTIATE_TEST_SUITE_P(InstantiationName,
                            QuasiIntegrableThermalizationDivergenceFunctionsTest,
                            ::testing::Values(0.1,0.01,0.001,0.0025,0.005,0.05,0.075));

}//end namespace
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
