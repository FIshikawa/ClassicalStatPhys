#include <gtest/gtest.h>
#include <physics/harmonic_oscillator_fixed_end.hpp>
#include <lattice/chain.hpp>
#include <integrator/dissipated_runge_kutta.hpp>

namespace {
class DissipatedRungeKuttaTest: public ::testing::Test {
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
    gamma = 1.0;
    dt = 0.01;
    N_time = static_cast<int>(1.0/dt);
    num_iteration = 1e+4;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      z[i+num_particles] = 1.0;
      if(i % 2 == 0) z[i] = 2.0;
      else z[i] = 0;
    }
  }

  int  num_particles, N_adj, N_time, num_iteration;
  double J, temperture, gamma, dt;
  std::vector<std::vector<int> > pair_table;
  std::vector<double> z;
};

TEST_F(DissipatedRungeKuttaTest, BasicTest) {
  // set integrator 
  integrator::DissipatedRungeKutta integrator(2*num_particles, temperture, gamma);

  hamiltonian::HarmonicOscillatorFixedEnd hamiltonian(num_particles,J,pair_table,N_adj);

   // set random numbers
  std::size_t seed = 1234;
  std::mt19937 mt(seed);

  double internal_enegy = 0.0;
  double potential_energy = 0.0;
  double kinetic_energy = 0.0;
  for(int step = 0; step < 10*N_time; ++step) integrator.step(0.0,dt,z,hamiltonian,mt);

  for(int mc_step = 0; mc_step < num_iteration; ++mc_step){
    for(int step = 0; step < N_time; ++step) integrator.step(0.0,dt,z,hamiltonian,mt);

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
