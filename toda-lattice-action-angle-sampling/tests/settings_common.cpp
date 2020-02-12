#include <gtest/gtest.h>
#include <common/settings_common.hpp>

TEST(SettingsTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10", "10", "100","10"};
  int argc = 10;
  int input_counter = 0;
  SettingsCommon settings_common(argc, argv, input_counter);
  ASSERT_EQ(10, settings_common.Ns);
  ASSERT_EQ(10, settings_common.Ns_observe);
  ASSERT_EQ(10, settings_common.N_mc);
  ASSERT_EQ(100, settings_common.N_loop);
  ASSERT_EQ(10, settings_common.num_iterations);
  settings_common.declare(std::cout);
}

