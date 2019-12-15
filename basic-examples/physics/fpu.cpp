#include<iostream>
#include<vector>
#include<random>

#include<clstatphys/integrator/runge_kutta_4th.hpp>
#include<clstatphys/physics/fpu_fixed_end.hpp>
#include<clstatphys/lattice/chain.hpp>

int main() {
  //class declare 
  int Ns = 10;
  lattice::Chain lattice(Ns);
  int num_particles = lattice.set_num_particles(Ns); 
  integrator::RungeKutta4th integrator(2*num_particles);
  //integrator::velver integrator(2*num_particles);
  int N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
  lattice.create_table(pair_table);
  //class declare end

  //target declare
  std::cout << "[test :" << hamiltonian::FPUFixedEnd::name()  << "]\n";
  double J = 1.0;
  double alpha = 1.0;
  double beta = 1.0;
  hamiltonian::FPUFixedEnd hamiltonian(num_particles,J,alpha,beta,pair_table,N_adj);
  //end target declare

  //parameter set 
  double dt = 0.001;
  double t = 1.0;
  int N_time = (int) (t / dt);
  std::vector<double> z(2*num_particles);
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  std::uniform_real_distribution<> uniform_rand(0.0,1.0);
  for(int i = 0; i < 2*num_particles; ++i) z[i] = uniform_rand(mt);
  z[0] = 0.0;
  //parmeter set end

  //time develop
  double pt = 0.0;
  double E = hamiltonian.energy(pt,z);
  double K = hamiltonian.kinetic_energy(pt,z);
  double P = hamiltonian.potential_energy(pt,z);
  std::cout << "pt : " << pt << ", q[0] : " << z[0] << ", q[1] : " << z[1] << ", E : " << E << ", K : " << K << ". P : " << P << std::endl; 
  for(int step = 0; step < N_time ; ++step){
    pt += dt;
    integrator.step(pt, dt, z, hamiltonian);
  }//end time step
  E = hamiltonian.energy(pt,z);
  K = hamiltonian.kinetic_energy(pt,z);
  P = hamiltonian.potential_energy(pt,z);
  std::cout << "pt : " << pt << ", q[0] : " << z[0] << ", q[1] : " << z[1] << ", E : " << E << ", K : " << K << ". P : " << P << std::endl; 
  //time develop end
    
  return 0;
} 
