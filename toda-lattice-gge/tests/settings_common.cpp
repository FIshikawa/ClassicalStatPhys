#include <gtest/gtest.h>
#include <common/settings_common.hpp>

//TEST(SettingsTest, ConstructTest){
//  SettingCommon setting_common();
//}

TEST(SettingsTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10","100","10.0","1.0","1","4","100","100","100","linear"};
  int argc = 10;
  int input_counter = 0;
  SettingsCommon setting_common(argc, argv, input_counter);
  ASSERT_EQ(10, setting_common.Ns);
  ASSERT_EQ(10, setting_common.Ns_observe);
  ASSERT_EQ(100, setting_common.N_time);
  ASSERT_EQ(10.0, setting_common.t);
  ASSERT_EQ(1.0, setting_common.E_initial);
  ASSERT_EQ(1, setting_common.k_initial);
  ASSERT_EQ(4, setting_common.N_normalmode);
  ASSERT_EQ("linear", setting_common.plot_scale);
  setting_common.declare(std::cout);
}

