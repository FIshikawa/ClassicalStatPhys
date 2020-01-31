#include <gtest/gtest.h>
#include <clstatphys/physics/toda_lax_form.hpp>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/physics/toda_action_variables.hpp>

namespace {
class TodaActionVariablds: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set by lattice::Chain
    num_particles = 5;
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
  int  num_particles;
  double J,alpha;
;
  std::vector<double> z;
};

TEST_F(TodaActionVariablds, IntegrationTest) {
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");
  // set integrator 
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


}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
