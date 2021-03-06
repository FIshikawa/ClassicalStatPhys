#include <gtest/gtest.h>
#include <limits>
#include <clstatphys/ensemble/point_heating.hpp>
#include <clstatphys/physics/harmonic_oscillator_fixed_end.hpp>
#include <clstatphys/lattice/chain.hpp>

namespace {
class PointHeatingTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set by lattice::Chain
    int Ns = 4;
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

TEST_F(PointHeatingTest, BasicTest) {
  hamiltonian::HarmonicOscillatorFixedEnd hamiltonian(num_particles,J,pair_table,N_adj);
  // set random numbers
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  // create data set
  ensemble::PointHeating sampler(temperture*0.5, temperture, num_particles/2, num_particles);

  std::vector<double> variance_position(num_particles);
  double total_velocity = 0.0;

  for(int mc_step = 0; mc_step < num_iteration; ++mc_step){
    sampler.set_initial_state(z, mt);
    for(int i = 0; i < num_particles; ++i){
      variance_position[i] += z[i]*z[i]/num_iteration;
      total_velocity += z[i+num_particles]/num_iteration;
    }
  }

  EXPECT_NEAR(total_velocity,0.0,1e-16);

  double error = 2.0/std::sqrt(num_iteration);
  std::cout << "error : " << error << std::endl;
  std::cout << "mean poition" << std::endl;
  for(int i = 0; i < num_particles/2; ++i){
    EXPECT_NEAR(variance_position[i],temperture*0.5,error);
    std::cout << i << " : " <<  variance_position[i] << std::endl;
  }
  for(int i = num_particles/2; i < num_particles; ++i){
    EXPECT_NEAR(variance_position[i],temperture,error);
    std::cout << i << " : " <<  variance_position[i] << std::endl;
  }
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
