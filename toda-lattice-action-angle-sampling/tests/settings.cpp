#include <gtest/gtest.h>
#include <specific/toda_lattice_periodic_boundary_hybrid_monte_carlo.hpp>

TEST(SettingsTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10", "10", "100","10",
                    "1.0", "1", "4",
                    "10", "10", "0.1", "1.0",
                    "2.0", "2.0"};
  int argc = 32;
  int input_counter = 0;
  Settings settings(argc, argv, input_counter);
  settings.declare(std::cout);

  ASSERT_EQ(1.0, settings.E_initial);
  ASSERT_EQ(1, settings.k_initial);
  ASSERT_EQ(4, settings.N_normalmode);

  ASSERT_EQ(10, settings.total_accept);
  ASSERT_EQ(10, settings.relax_time);
  ASSERT_EQ(0.1, settings.dt_relax);
  ASSERT_EQ(1.0, settings.temperture);

  ASSERT_EQ(2.0, settings.J);
  ASSERT_EQ(2.0, settings.alpha);

  auto hamiltonian = settings.hamiltonian();
  ASSERT_EQ(2, hamiltonian.Nd());
  auto ensembler = settings.ensembler();
  auto monte_carlo_sampler = settings.monte_carlo_sampler();
}

