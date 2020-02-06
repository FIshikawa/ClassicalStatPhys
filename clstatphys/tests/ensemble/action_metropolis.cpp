#include <gtest/gtest.h>
#include <clstatphys/physics/toda_lax_form.hpp>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/physics/toda_action_variables.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <clstatphys/ensemble/action_metropolis.hpp>

namespace {
class MetropolisActionTodaTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    num_particles = 12;
    // set interaction constants 
    num_iteration = 16;
    J = 1.0;
    alpha = 1.0;
    dx = 0.1;
    dp = 0.1;
    betas = std::vector<double>(num_particles-1, 1.0);
    num_mc_step = 1e+4;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      z[i+num_particles] = 1.0;
      if(i % 2 == 0) z[i] = 2.0;
      else z[i] = 0;
    }
  }
  std::vector<double> z, betas;
  int num_particles, num_iteration, num_mc_step;
  double dp, dx;
  double J,alpha;
};

TEST_F(MetropolisActionTodaTest, BasicTest) {
  // set random numbers
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  // create data set
  ensemble::MetropolisActionToda sampler(
                                         num_particles, 
                                         num_iteration, 
                                         betas,
                                         J,
                                         alpha,
                                         dp,
                                         dx
                                         );
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");


  std::vector<double> actions_init(num_particles-1);
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  integrable::TodaDiscriminant discriminant(num_particles,L, "periodic");
  TodaActionVariables(actions_init, discriminant, num_iteration);
  std::cout << "actions init :" << std::endl;
  for(int i = 0; i < num_particles-1;++i) 
    std::cout << "   i : " << i << " , J : " << actions_init[i] << std::endl;

  int counter_init = 0;
  for(int mc_step = 0; mc_step < 100; ++mc_step) sampler.montecarlo(z, counter_init, mt);
  std::cout << "initialized number : " << counter_init << std::endl;

  std::vector<double> actions(num_particles-1);
  std::vector<double> actions_temp(num_particles-1);

  int counter_measured = 0;
  for(int mc_step = 0; mc_step < num_mc_step; ++mc_step){
    sampler.montecarlo(z, counter_measured, mt);

    L = toda_lax_form.L_matrix(z);
    discriminant = integrable::TodaDiscriminant(num_particles,L, "periodic");
    TodaActionVariables(actions_temp, discriminant, num_iteration);

    for(int i = 0; i < num_particles-1;++i) actions[i] += actions_temp[i] / num_mc_step;
  }
  std::cout << "measured number : " << counter_measured<< std::endl;

  double error = std::sqrt(1.0/num_mc_step) * 5;

  std::cout << "actions end:" << std::endl;
  for(int i = 0; i < num_particles-1;++i) 
    std::cout << "   i : " << i << " , J : " << actions[i] << " , beta : " << betas[i] << std::endl;
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
