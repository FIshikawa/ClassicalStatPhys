#include <gtest/gtest.h>
#include <clstatphys/tools/histogram.hpp>

namespace {
class HistogramTest : public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set constants 
    n_bin = 100;
    N_loop = 1E+02;
    exact_mean = 0.0;
    exact_variance = 1.0;
    data.resize(N_loop);
    hist.resize(n_bin);
    bins.resize(n_bin);
    histogram.initialize(n_bin);
    // check initialize
    ASSERT_EQ(n_bin, histogram.size());
    ASSERT_DOUBLE_EQ(0, histogram.min());
    ASSERT_DOUBLE_EQ(0, histogram.max());
    // create data set
    std::size_t seed = 1234;
    std::mt19937 mt(seed);
    std::normal_distribution<> dist(exact_mean, exact_variance);
    for(int i = 0; i < N_loop; ++i) data[i] = dist(mt); 
    // create histogram
    histogram(data);
    ASSERT_EQ(N_loop, histogram.total_count());
    // output hist and bins
    histogram.output(hist,bins);
  }
  tools::Histogram histogram;
  int  n_bin;
  int  N_loop;
  std::vector<double> data,bins,hist;
  double exact_mean;
  double exact_variance;
};

TEST_F(HistogramTest, FundamentalTest) {
  EXPECT_DOUBLE_EQ(bins[1], bins[0] + histogram.interval());
  EXPECT_DOUBLE_EQ(bins[0], histogram.min());
  EXPECT_DOUBLE_EQ(bins[bins.size()-1], histogram.max());
}

TEST_F(HistogramTest, HistTest) {
  double normalize_const_output = histogram.normalize_const(); 
  double normalize_const_from_hist = 0.0;
  double variance = 0.0;
  double mean = 0.0;
  double bin_width = histogram.interval();
//double bin_width = bins[1] - bins[0];
  for(int i = 0 ; i < n_bin ; ++i){
    normalize_const_from_hist += hist[i] * bin_width;
    variance += hist[i]/normalize_const_output * bins[i] * bins[i] * bin_width;
    mean += hist[i]/normalize_const_output * bins[i] * bin_width;
  }
  variance -= mean*mean;
  EXPECT_DOUBLE_EQ(normalize_const_from_hist, normalize_const_output);
  EXPECT_NEAR(mean, histogram.mean(),bin_width);
  EXPECT_NEAR(variance, histogram.variance(),bin_width);

}

TEST_F(HistogramTest, ExactTest) {
  double N_loop_double = static_cast<double>(N_loop);
  double expect_mean_error = std::sqrt(histogram.unbiassed_variance() / N_loop_double);
  EXPECT_FALSE((exact_mean - histogram.mean()) < expect_mean_error * 0.1);
  EXPECT_FALSE((exact_variance - histogram.unbiassed_variance()) < 1/std::sqrt(N_loop_double) * 0.1 );
  EXPECT_NEAR(exact_mean, histogram.mean(), expect_mean_error);
  EXPECT_NEAR(exact_variance, histogram.unbiassed_variance(), 1/std::sqrt(N_loop_double)); 
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
