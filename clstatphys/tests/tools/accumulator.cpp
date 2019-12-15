#include <gtest/gtest.h>
#include <tools/accumulator.hpp>

namespace {
class AccumulatorTest : public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set constants 
    n_bin = 100;
    N_loop = 1E+02;
    exact_mean = 0.0;
    exact_variance = 1.0;
    // create data set
    std::size_t seed = 1234;
    std::mt19937 mt(seed);
    std::normal_distribution<> dist(exact_mean, exact_variance);
    data.resize(N_loop);
    for(int i = 0; i < N_loop; ++i) data[i] = dist(mt); 
  }
  int  n_bin;
  int  N_loop;
  std::vector<double> data;
  double exact_mean;
  double exact_variance;
};

TEST_F(AccumulatorTest, FundamentalTest) {
  tools::Accumulator accumulator;
  tools::Accumulator accumulator_t;
  for(int i = 0; i < N_loop; ++i) accumulator << data[i]; 
  for(int i = 0; i < N_loop; ++i) accumulator_t << data[i]; 
  EXPECT_DOUBLE_EQ(accumulator.count(), N_loop);
  double sum1 = accumulator.sum1();
  accumulator = accumulator + accumulator_t;
  EXPECT_DOUBLE_EQ(accumulator.count(), 2*N_loop);
  EXPECT_DOUBLE_EQ(accumulator.sum1(), sum1+accumulator_t.sum1());
}

TEST_F(AccumulatorTest, ExactTest) {
  tools::Accumulator accumulator;
  for(int i = 0; i < N_loop; ++i) accumulator << data[i]; 
  double variance = accumulator.variance();
  double mean = accumulator.mean();
  double n = accumulator.count();
  EXPECT_NEAR(mean, exact_mean, 1.0 / std::sqrt(N_loop));
  EXPECT_NEAR(variance, exact_variance, 1.0 / std::sqrt(N_loop));
  EXPECT_NEAR(mean*n, accumulator.sum1(), n / std::sqrt(N_loop));
  EXPECT_NEAR(variance*n, accumulator.sum2(), n / std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.kurtosis(),3.0,1.0/std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.kurtosis_excess(),0.0,1.0/std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.skewness(),0.0,1.0/std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.average(), exact_mean, 1.0 / std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.error(), exact_variance/(N_loop-1), 1.0 / std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.central_moment1(),1.0,std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.central_moment2(),1.0,std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.central_moment3(),1.0,std::sqrt(N_loop));
  EXPECT_NEAR(accumulator.central_moment4(),1.0,std::sqrt(N_loop));
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
