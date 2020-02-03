#include <gtest/gtest.h>
#include <clstatphys/lattice/chain.hpp>
#include <clstatphys/integrator/runge_kutta_4th.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <clstatphys/physics/toda_lax_form.hpp>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/physics/toda_action_variables.hpp>

namespace {
class TodaActionVariablds: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set by lattice::Chain
    int Ns = 32;
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

    // set interaction constants 
    J = 1.0;
    alpha = 1.0;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      if(i % 2 == 0){
        z[i] = -1;
        z[i+num_particles] = -1;
      }
      else{
        z[i] = 1;
        z[i+num_particles] = 1;
      }
    }
  }
  int  num_particles, N_adj;
  double J,alpha;
  std::vector<std::vector<int> > pair_table;
  std::vector<double> z;
};

TEST_F(TodaActionVariablds, IntegrationTest) {
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  std::cout << " L : " << std::endl 
            << L << std::endl;

  integrable::TodaDiscriminant discriminant(num_particles,L, "periodic");

  std::vector<double> roots(2*num_particles);
  discriminant.total_roots(roots);
  std::cout << "roots : " << std::endl;
  for(int i = 0 ; i < 2*num_particles; ++i) std::cout << roots[i] << " ";
  std::cout << std::endl;
  std::cout << "check value at roots: " << std::endl;
  for(int i = 0 ; i < 2*num_particles; ++i) std::cout << discriminant(roots[i]) << " ";
  std::cout << std::endl;

  std::vector<double> action_variables(num_particles-1);
  int num_iterations = 16;
  TodaActionVariables(action_variables, discriminant, num_iterations);

  for(int i = 0;i < num_particles-1;++i) std::cout << action_variables[i] << " ";
  std::cout << std::endl;

}

TEST_F(TodaActionVariablds, TimeDevelopTest) {
  hamiltonian::TodaLattice hamiltonian(num_particles,J,alpha,pair_table,N_adj);
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");
  integrator::RungeKutta4th integrator(2*num_particles);

  int num_iterations = 16;
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  integrable::TodaDiscriminant discriminant(num_particles,L, "periodic");
  std::vector<double> action_variables_init(num_particles-1);
  TodaActionVariables(action_variables_init, discriminant, num_iterations);

  std::cout << "action variables init: " << std::endl << "   ";
  for(int i = 0; i < num_particles-1;++i) std::cout << action_variables_init[i] << " ";
  std::cout << std::endl;

  double dt = 1.0e-2;
  double pt = 0.0;
  double N_time = 1e+5;
  double error_rungekutta = dt * dt * dt * dt * dt * 20; //10 is assumed constant of error 


  std::vector<double> z_temp(z);
  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z_temp,hamiltonian);
    pt += dt;
  }

  L = toda_lax_form.L_matrix(z_temp);
  discriminant = integrable::TodaDiscriminant(num_particles, L, "periodic");
  std::vector<double> action_variables_mid(num_particles-1);
  TodaActionVariables(action_variables_mid, discriminant, num_iterations);

  std::cout << "action variables mid: " << std::endl << "   ";
  for(int i = 0; i < num_particles-1;++i) std::cout << action_variables_mid[i] << " ";
  std::cout << std::endl;

  for(int i = 0; i < num_particles-1; ++i){
    EXPECT_NEAR(action_variables_init[i],action_variables_mid[i],error_rungekutta*N_time);
  }
  std::cout << std::endl;

  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z_temp,hamiltonian);
    pt += dt;
  }

  L = toda_lax_form.L_matrix(z_temp);
  discriminant = integrable::TodaDiscriminant(num_particles, L, "periodic");
  std::vector<double> action_variables_fin(num_particles-1);
  TodaActionVariables(action_variables_fin, discriminant, num_iterations);

  std::cout << "action variables fin: " << std::endl << "   ";
  for(int i = 0; i < num_particles-1;++i) std::cout << action_variables_fin[i] << " ";
  std::cout << std::endl;

  for(int i = 0; i < num_particles-1; ++i){
    EXPECT_NEAR(action_variables_init[i],action_variables_fin[i],error_rungekutta*N_time);
  }
  std::cout << std::endl;


}



}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
