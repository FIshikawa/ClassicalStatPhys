#include <gtest/gtest.h>
#include <physics/toda_lax_form.hpp>
#include <physics/todalattice.hpp>
#include <lattice/chain.hpp>
#include <integrator/runge_kutta_4th.hpp>

namespace {
class TodaLaxFormTest: public ::testing::Test {
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
    ASSERT_EQ(pair_table[0][0],num_particles - 1);
    ASSERT_EQ(pair_table[0][1],1);
    ASSERT_EQ(pair_table[num_particles-1][0],num_particles-2);
    ASSERT_EQ(pair_table[num_particles-1][1],0);
    // set interaction constants 
    J = 0.25;
    alpha = 2.0;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      if(i % 2 == 0){
        z[i] = 1.0;
        z[i+num_particles] = 1.0;
      }
      else{
        z[i] = 0;
        z[i+num_particles] = 1.0;
      }
    }
    N_conservations = 20;
  }
  int N_conservations;
  int  num_particles, N_adj;
  double J,alpha;
  std::vector<std::vector<int> > pair_table;
  std::vector<double> z;
};

TEST_F(TodaLaxFormTest, ConservationTest) {
  hamiltonian::TodaLattice hamiltonian(num_particles,J,alpha,pair_table,N_adj);
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,pair_table,N_adj);
  // set integrator 
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  std::cout << " L : " << std::endl 
            << L << std::endl;

  integrator::RungeKutta4th integrator(2*num_particles);
  std::vector<double> conservations_init(N_conservations,0.0);
  std::vector<double> conservations(N_conservations,0.0);
  std::vector<double> eigenvalues(num_particles,0.0);
  std::vector<double> z_temp(z);
  double energy_total = hamiltonian.energy(0.0,z_temp);
  double dt = 1.0e-2;
  double pt = 0.0;
  double N_time = 1e+5;
  double error_rungekutta = dt * dt * dt * dt * dt * energy_total * energy_total; //10 is assumed constant of error 
  std::cout << "initial energy : " << energy_total << std::endl;
  toda_lax_form.conservations_with_eigenvalues(z_temp, conservations_init, eigenvalues);
  for(int i = 0; i < conservations.size(); ++i){
    std::cout <<"i : " << i << ", conservations_init : " << conservations_init[i] << std::endl;
  }
  for(int i = 0; i < num_particles; ++i){
    std::cout <<"i : " << i << ", eigenvalues_init : " << eigenvalues[i] << std::endl;
  }
  double momentum_total = num_particles * 1.0;
  EXPECT_DOUBLE_EQ(energy_total+num_particles*J, conservations_init[1]*2.0/(alpha*alpha));
  EXPECT_DOUBLE_EQ(momentum_total, conservations_init[0]*2.0/alpha);

  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z_temp,hamiltonian);
    pt += dt;
  }
  toda_lax_form.conservations_with_eigenvalues(z_temp, conservations, eigenvalues);
  EXPECT_NEAR(2.0*conservations[0]/momentum_total/alpha,1.0,error_rungekutta*N_time);
  EXPECT_NEAR(2.0*conservations[1]/(energy_total + num_particles*J)/(alpha*alpha),1.0,error_rungekutta*N_time);
  L = toda_lax_form.L_matrix(z_temp);
  std::cout << " L : " << std::endl 
            << L << std::endl;
  for(int i = 0; i < conservations.size(); ++i){
    EXPECT_NEAR(conservations[i]/conservations_init[i],1.0,error_rungekutta*N_time);
  }

  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z_temp,hamiltonian);
    pt += dt;
  }
  toda_lax_form.conservations_with_eigenvalues(z_temp, conservations, eigenvalues);
  EXPECT_NEAR(2.0*conservations[0]/momentum_total/alpha,1.0,error_rungekutta*N_time);
  EXPECT_NEAR(2.0*conservations[1]/(energy_total + num_particles*J)/(alpha*alpha),1.0,error_rungekutta*N_time);
  L = toda_lax_form.L_matrix(z_temp);
  std::cout << " L : " << std::endl 
            << L << std::endl;
  for(int i = 0; i < conservations.size(); ++i){
    EXPECT_NEAR(conservations[i]/conservations_init[i],1.0,2.0*error_rungekutta*N_time);
  }
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
