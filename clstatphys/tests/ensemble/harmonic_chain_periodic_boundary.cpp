#include <gtest/gtest.h>
#include <ensemble/harmonic_chain_periodic_boundary.hpp>
#include <physics/harmonic_oscillator.hpp>
#include <lattice/chain.hpp>
#include <limits>

namespace {
class HarmonicChainPeriodicBoundaryEnsemblePeriodicTest : public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set lattice
    int Ns = 1e+2;
    lattice::Chain lattice(Ns);
    num_particles = lattice.set_num_particles(Ns); 
    ASSERT_EQ(num_particles, Ns);
    N_adj = lattice.number_adjacent();
    ASSERT_EQ(N_adj,2);
    pair_table = std::vector<std::vector<int> >(num_particles,std::vector<int>(N_adj));
    lattice.create_table(pair_table);
    ASSERT_EQ(pair_table[0][0],num_particles-1);
    ASSERT_EQ(pair_table[0][1],1);
    ASSERT_EQ(pair_table[num_particles-1][0],num_particles-2);
    ASSERT_EQ(pair_table[num_particles-1][1],0);
    // set constants
    J = 1.0;
    temperture = 1.0;
    z.resize(2 * num_particles);
  }
  std::vector<double> z;
  int num_particles, N_adj;
  double temperture, J;
  std::vector<std::vector<int> > pair_table;
};

TEST_F(HarmonicChainPeriodicBoundaryEnsemblePeriodicTest, FundamentalTest) {
  hamiltonian::HarmonicOscillator hamiltonian(num_particles,J,pair_table,N_adj);
  // set random numbers
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  // create data set
  ensemble::HarmonicChainPeriodicBoundary initializer(J, temperture, num_particles);

  double internal_enegy = 0.0;
  double potential_energy = 0.0;
  double kinetic_energy = 0.0;
  int num_iteration = 1e+4;

  for(int i = 0; i < num_iteration; ++i){
    initializer.set_initial_state(z, mt);
    initializer.equilibrate_velocity(z,mt);

    internal_enegy += hamiltonian.energy(0.0,z)/num_particles/num_iteration;
    potential_energy += hamiltonian.potential_energy(0.0,z)/num_particles/num_iteration;
    kinetic_energy += hamiltonian.kinetic_energy(0.0,z)/num_particles/num_iteration;
  }

  double sum_x = 0.0;
  for(int i = 0 ; i < num_iteration; ++i) sum_x += z[i];
  std::cout << "sum : " << sum_x << std::endl;

  std::cout << "internal_enegy : " << internal_enegy << std::endl;
  std::cout << "kinetic_energy: " << kinetic_energy << std::endl;
  std::cout << "potential_energy: " << potential_energy << std::endl;
  double expected_internal = 0.5 * 2.0 * temperture;
  std::cout << "expected_internal : " << expected_internal << std::endl;

  double precision = 1.0/std::sqrt(1.0*num_iteration);

  EXPECT_NEAR(internal_enegy/expected_internal,1.0, 1.0/num_particles + precision);
}

}
