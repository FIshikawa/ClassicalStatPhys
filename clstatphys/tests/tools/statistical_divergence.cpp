#include <gtest/gtest.h>
#include <random>
#include <limits>
#include <tools/statistical_divergence.hpp>

namespace {
class StatisticalDivergenceTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set constants 
    N_loop = 40;
    p.resize(N_loop);
    q.resize(N_loop);
    // create data set
    double sum_value_p = 0;
    double sum_value_q = 0;
    for(int i = 0; i < N_loop; ++i){
      if( i < 3*N_loop/4) p[i] = 1.0/30; 
      if( N_loop/4 <= i) q[i] = 1.0/30; 
      sum_value_p += p[i];
      sum_value_q += q[i];
    }
    ASSERT_DOUBLE_EQ(sum_value_p,1.0);
    ASSERT_DOUBLE_EQ(sum_value_q,1.0);
  }
  int  N_loop;
  std::vector<double> p,q;
};


TEST_F(StatisticalDivergenceTest, ExactTest) {
  tools::StatisticalDivergence f_distance;
  f_distance.calc(p,q);
  double prec = std::sqrt(2.0 * std::numeric_limits<double>::epsilon());
  double js_divergence = f_distance.jensen_shannon_divergence();
  double h_distance = f_distance.hellinger_distance();
  double l1_norm = f_distance.l1_norm();
  double l2_norm = f_distance.l2_norm();
  std::cout << "js  : " << js_divergence << std::endl
            << "h   : " << h_distance << std::endl
            << "l1  : " << l1_norm<< std::endl
            << "l2  : " << l2_norm<< std::endl;
  EXPECT_DOUBLE_EQ(l1_norm, 20.0/30);
  EXPECT_DOUBLE_EQ(l2_norm, std::sqrt(20)/30);
  EXPECT_DOUBLE_EQ(h_distance,std::sqrt(1.0-20.0/30));
  EXPECT_DOUBLE_EQ(js_divergence, 1.0/3.0*std::log(2));
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
