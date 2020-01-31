#include <gtest/gtest.h>
#include <clstatphys/physics/toda_lax_form.hpp>
#include <clstatphys/physics/toda_discriminant.hpp>

namespace {
class TodaDiscriminant: public ::testing::Test {
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
      z[i+num_particles] = 1.0;
      if(i % 2 == 0) z[i] = 1.0;
      else z[i] = -1;
    }
  }
  int  num_particles;
  double J,alpha;
;
  std::vector<double> z;
};

TEST_F(TodaDiscriminant, ZeroRootsPlusTest) {
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");
  // set integrator 
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  std::cout << " L : " << std::endl 
            << L << std::endl;

  integrable::TodaDiscriminant discriminant(num_particles,L, "periodic");

  std::vector<double> w(num_particles);
  rokko::dlvector w_temp(num_particles);
  int info = rokko::lapack::syev('V', 'U', L, w_temp);
  std::cout << " eigvalues : " << std::endl
            << w_temp << std::endl;
  for(int i = 0; i < num_particles; ++i) w[i] = w_temp[i];

  for(int i = 0; i < num_particles; ++i) EXPECT_NEAR(2.0,discriminant(w[i]),1e-8);

}

TEST_F(TodaDiscriminant, ZeroRootsMinusTest) {
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");
  // set integrator 
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  L(0,num_particles-1) *= -1;
  L(num_particles-1,0) *= -1;
  std::cout << " L : " << std::endl 
            << L << std::endl;

  integrable::TodaDiscriminant discriminant(num_particles,L, "periodic");

  std::vector<double> w(num_particles);
  rokko::dlvector w_temp(num_particles);
  int info = rokko::lapack::syev('V', 'U', L, w_temp);
  std::cout << " eigvalues : " << std::endl
            << w_temp << std::endl;

  for(int i = 0; i < num_particles; ++i) w[i] = w_temp[i];

  for(int i = 0; i < num_particles; ++i) EXPECT_NEAR(2.0,discriminant(w[i]),1e-8);
}

TEST_F(TodaDiscriminant, ZeroRootsTotal) {
  integrable::TodaLaxForm toda_lax_form(num_particles,J,alpha,"periodic");
  rokko::dlmatrix L = toda_lax_form.L_matrix(z);
  // set integrator 
  integrable::TodaDiscriminant discriminant(num_particles, L, "periodic");

  std::vector<double> w(2*num_particles);
  discriminant.total_roots(w);
  std::cout << " roots values : " << std::endl;
  for(int i = 0; i < 2*num_particles; ++i){
    std::cout << w[i] << " ";
    if(discriminant(w[i]) > 0) EXPECT_NEAR(2.0,discriminant(w[i]),1e-12);
    else EXPECT_NEAR(-2.0,discriminant(w[i]),1e-12);
  }
  std::cout << std::endl;
}


}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
