#include <gtest/gtest.h>
#include <clstatphys/tools/false_position_method.hpp>

double f(double x) { return 3*x*x-5*x+1; }

TEST(FlasePositionTest, FlasePosition1) {
  optimization::FalsePositionMethod optimizer;
  int iteration;
  iteration = optimizer.find_zero(f, 0, 1);
  EXPECT_TRUE(iteration > 0);
  EXPECT_DOUBLE_EQ((5 - std::sqrt(13)) / 6, optimizer.zero());
}

TEST(FlasePositionTest, FlasePosition2) {
  optimization::FalsePositionMethod optimizer;
  int iteration;
  iteration = optimizer.find_zero(f, 1, 10);
  EXPECT_TRUE(iteration > 0);
  EXPECT_DOUBLE_EQ((5 + std::sqrt(13)) / 6, optimizer.zero());
}

TEST(FlasePositionTest, FlasePosition3) {
  optimization::FalsePositionMethod optimizer;
  int iteration;
  iteration = optimizer.find_zero(f, 10, 1);
  EXPECT_TRUE(iteration > 0);
  EXPECT_DOUBLE_EQ((5 + std::sqrt(13)) / 6, optimizer.zero());
}

TEST(FlasePositionTest, FlasePosition4) {
  optimization::FalsePositionMethod optimizer;
  int iteration;
  iteration = optimizer.find_zero(f, 3, 4);
  EXPECT_FALSE(iteration > 0);
}
