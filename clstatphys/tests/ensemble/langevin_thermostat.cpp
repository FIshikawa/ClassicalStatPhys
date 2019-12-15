#include <gtest/gtest.h>
#include <ensemble/langevin_thermostat.hpp>
#include <physics/harmonic_oscillator_fixed_end.hpp>
#include <lattice/chain.hpp>
#include <limits>

namespace {
class LangevinThermostatTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set by lattice::Chain
    int Ns = 2;
    lattice::Chain lattice(Ns);
    num_particles = lattice.set_num_particles(Ns); 
    ASSERT_EQ(num_particles, Ns);
    N_adj = lattice.number_adjacent();
    ASSERT_EQ(N_adj,2);
    pair_table = std::vector<std::vector<int> >(num_particles,std::vector<int>(N_adj));
    lattice.create_table(pair_table);
    // set interaction constants 
    J = 1.0;
    temperture = 1.0;
    dt = 0.01;
    relax_time = 1.0;
    num_iteration = 1e+4;
    gamma = 2.0;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      z[i+num_particles] = 1.0;
      if(i % 2 == 0) z[i] = 2.0;
      else z[i] = 0;
    }
  }
  std::vector<double> z;
  int num_particles, N_adj, total_accept, num_iteration;
  double temperture, J, relax_time, dt, gamma;
  std::vector<std::vector<int> > pair_table;
};

TEST_F(LangevinThermostatTest, BasicTest) {
  hamiltonian::HarmonicOscillatorFixedEnd hamiltonian(num_particles,J,pair_table,N_adj);
  // set random numbers
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  // create data set
  ensemble::LangevinThermostat sampler(num_particles, temperture, dt, gamma, relax_time);

  double internal_enegy = 0.0;
  double potential_energy = 0.0;
  double kinetic_energy = 0.0;

  for(int mc_step = 0; mc_step < num_iteration; ++mc_step){
    sampler.sample(z, hamiltonian, mt);

    internal_enegy += hamiltonian.energy(0.0,z)/num_particles/num_iteration;
    potential_energy += hamiltonian.potential_energy(0.0,z)/num_particles/num_iteration;
    kinetic_energy += hamiltonian.kinetic_energy(0.0,z)/num_particles/num_iteration;
  }

  double expected_energy = temperture;

  double error = std::sqrt(1.0/num_iteration) * 5;

  std::cout << "expected_error : " << error << std::endl;
  std::cout << "internal_enegy : " << internal_enegy << std::endl;
  std::cout << "kinetic_energy: " << kinetic_energy << std::endl;
  std::cout << "potential_energy: " << potential_energy << std::endl;

  double expected_internal = 0.5 * 2.0 * temperture;
  std::cout << "expected_internal : " << expected_internal << std::endl;

  EXPECT_NEAR(potential_energy/(expected_internal/2.0),1.0,error);
  EXPECT_NEAR(kinetic_energy/(expected_internal/2.0),1.0,error);
  EXPECT_NEAR(internal_enegy/expected_internal,1.0,error);
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
