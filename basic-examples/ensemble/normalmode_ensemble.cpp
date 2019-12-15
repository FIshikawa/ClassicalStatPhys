#include<iostream>
#include<vector>
#include<random>

#include<clstatphys/ensemble/normalmode_ensemble.hpp>
#include<clstatphys/lattice/chain.hpp>
#include<clstatphys/physics/fpu_fixed_end.hpp>
#include<clstatphys/physics/normalmode_energy.hpp>

int main() {
  std::cout << "[test :" << ensemble::NormalModeEnsemble::name()  << "]\n";
  //class declare 
  int Ns = 100;
  lattice::Chain lattice(Ns);
  int num_particles = lattice.set_num_particles(Ns); 
  int N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
  lattice.create_table(pair_table);
  double J = 1.0;
  double alpha = 1.0;
  double beta = 1.0;
  hamiltonian::FPUFixedEnd hamiltonian(num_particles,J,alpha,beta,pair_table,N_adj);
  //class declare end

  //target declare
  double epsilon = 1e-2;
  int k_initial = 0;
  int N_k = num_particles/10;
  double E = epsilon * num_particles;
  ensemble::NormalModeEnsemble initializer(num_particles,k_initial,N_k,E);
  //end target declare

  //parameter set 
  std::vector<double> z(2*num_particles);
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  //parmeter set end

  //initialize 
  initializer.set_initial_state(z,mt);
  std::vector<double> Ek(num_particles);
  NormalModeEnergy(z,Ek);
  double E_total = 0.0;
  //initialize end

  //output 
  std::cout << " #output k Energy/epsiclon " << std::endl;
  for(int i = 1 ; i < num_particles ; ++i){
    std::cout << i << " " << Ek[i]/epsilon << std::endl;
    E_total += Ek[i];
  }
  std::cout << "# E_0 : " << E << ", E_total : " << E_total << ", E_ham : " << hamiltonian.energy(0.0,z) << std::endl;
  //end output 

  return 0;
} 
