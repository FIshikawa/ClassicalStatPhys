#include <gtest/gtest.h>
#include <physics/sinh_lattice_fixed_end_parallel.hpp>
#include <lattice/chain_open_boundary.hpp>
#include <integrator/runge_kutta_4th.hpp>

namespace {
class SinhLatticeFixedEndParallelTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set by lattice::Chain
    int Ns = 5;
    lattice::ChainOpenBoundary lattice(Ns);
    num_particles = lattice.set_num_particles(Ns); 
    ASSERT_EQ(num_particles, Ns);
    N_adj = lattice.number_adjacent();
    ASSERT_EQ(N_adj,2);
    pair_table = std::vector<std::vector<int> >(num_particles,std::vector<int>(N_adj));
    lattice.create_table(pair_table);
    ASSERT_EQ(pair_table[0][0],0);
    ASSERT_EQ(pair_table[0][1],1);
    ASSERT_EQ(pair_table[num_particles-1][0],num_particles-2);
    ASSERT_EQ(pair_table[num_particles-1][1],num_particles-1);
    // set interaction constants 
    J = 1.0;
    alpha = 1.0;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      z[i+num_particles] = 1.0;
      if(i % 2 == 0) z[i] = 1.0;
      else z[i] = 0;
    }
  }

  int  num_particles, N_adj;
  double J,alpha;
  std::vector<std::vector<int> > pair_table;
  std::vector<double> z;
};

TEST_F(SinhLatticeFixedEndParallelTest, ForceTest) {
  hamiltonian::SinhLatticeFixedEndParallel hamiltonian(num_particles,J,alpha,pair_table,N_adj);
  // set integrator 
  integrator::RungeKutta4th integrator(2*num_particles);
  double dt = 1.0e-2;
  double pt = 0.0;
  double N_time = 1e+4;
  double error_rungekutta = dt * dt * dt * dt * dt * 10; //10 is assumed constant of error 
  double energy_init = hamiltonian.energy(pt,z);
  double potential_init = hamiltonian.potential_energy(pt,z);
  double kinetic_init = hamiltonian.kinetic_energy(pt,z);
  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z,hamiltonian);
    pt += dt;
  }
  EXPECT_FALSE(energy_init == hamiltonian.energy(pt,z));
  EXPECT_NEAR(energy_init/hamiltonian.energy(pt,z),1.0,error_rungekutta*N_time);
  EXPECT_FALSE(potential_init == hamiltonian.potential_energy(pt,z)); 
  EXPECT_FALSE(kinetic_init == hamiltonian.kinetic_energy(pt,z)); 
  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z,hamiltonian);
    pt += dt;
  }
  EXPECT_FALSE(energy_init == hamiltonian.energy(pt,z));
  EXPECT_NEAR(energy_init/hamiltonian.energy(pt,z),1.0,error_rungekutta*N_time);
  EXPECT_FALSE(potential_init == hamiltonian.potential_energy(pt,z)); 
  EXPECT_FALSE(kinetic_init == hamiltonian.kinetic_energy(pt,z)); 
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
