#include <gtest/gtest.h>
#include <common/settings_common.hpp>

TEST(SettingsTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10","100","10.0","100","100","linear","10"};
  int argc = 11;
  int input_counter = 0;
  SettingsCommon settings_common(argc, argv, input_counter);
  ASSERT_EQ(10, settings_common.Ns);
  ASSERT_EQ(10, settings_common.Ns_observe);
  ASSERT_EQ(100, settings_common.N_time);
  ASSERT_EQ(10.0, settings_common.t);
  ASSERT_EQ(100, settings_common.N_loop);
  ASSERT_EQ(100, settings_common.N_time_measure);
  ASSERT_EQ("linear", settings_common.plot_scale);
  ASSERT_EQ(10, settings_common.num_iterations);
  settings_common.declare(std::cout);
}

