#include <gtest/gtest.h>
#include <specific/toda_lattice_periodic_boundary_thermalization_from_normalmode_ensemble.hpp>

TEST(SettingsTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10","100","10.0","1.0","1","4","100","100","100","linear","1.0","1.0"};
  int argc = 10;
  int input_counter = 0;
  Settings settings(argc, argv, input_counter);
  settings.declare(std::cout);
  ASSERT_EQ(10, settings.Ns);
  ASSERT_EQ("linear", settings.plot_scale);
  ASSERT_EQ(4, settings.N_normalmode);
  ASSERT_EQ(1.0, settings.J);
  ASSERT_EQ(1.0, settings.alpha);

  auto hamiltonian = settings.hamiltonian();
  ASSERT_EQ(2, hamiltonian.Nd());
  auto ensembler = settings.ensembler();
}

