#include <gtest/gtest.h>
#include <specific/harmonic_oscillator_periodic_boundary_thermalization_from_point_heating.hpp>

TEST(SettingsTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10","100","10.0","1.0","100","100","linear",
                    "2.0", "2.0", "5"};
  int argc = 13;
  int input_counter = 0;
  Settings settings(argc, argv, input_counter);
  settings.declare(std::cout);
  ASSERT_EQ(10, settings.Ns);
  ASSERT_EQ(10, settings.Ns_observe);
  ASSERT_EQ("linear", settings.plot_scale);
  ASSERT_EQ(2.0, settings.J);
  ASSERT_EQ(2.0, settings.heating_temperture);
  ASSERT_EQ(5, settings.num_heating_particles);

  auto hamiltonian = settings.hamiltonian();
  ASSERT_EQ(2, hamiltonian.Nd());
  auto ensembler = settings.ensembler();
}

