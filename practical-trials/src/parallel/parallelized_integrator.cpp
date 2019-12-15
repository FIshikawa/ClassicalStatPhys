#include<omp.h>
#include<vector>
#include<random>
#include<fstream>
#include<iostream>
#include<string>
#include<chrono>

#include "Lattice.hpp"
#include "hamiltonian.hpp"
#include "integrator.hpp"

int main(int argc, char **argv){
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[test parallelaizeation]\n";
  int N_parallel_max = omp_get_num_procs();
  std::cout << "[limit threads : " << N_parallel_max << "]\n" ;

  //set output 
  std::cout<< "open ouput_file" << std::endl;
  std::ofstream ofs("result.dat",std::ios::out);
  std::ofstream ofsr("result_TimeAtNsMax.dat",std::ios::out);
  ofs.close();
  ofs.open("result.dat",std::ios::app);
  ofs << std::scientific;
  ofsr << std::scientific;
  ofs << "#output Ns ";
  ofsr << "#output Np time" << std::endl;

  for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){
    ofs << "time_Np_" << std::to_string(N_parallel) << " ";
  }
  ofs << std::endl;
  //set output end 

  std::cout<< "start Ns iteration" << std::endl;

  int N_iterate = 1000;
  int Ns_init = 1e4;
  int Ns_max = 1e6;
  int dNs = (Ns_max - Ns_init)/N_iterate;

  for(int i = 0; i < N_iterate; ++i){
    int Ns = Ns_init + dNs * i;

    //std::cout<< "create lattice" << std::endl;
    Lattice::FullyConnected lattice(Ns);
    int N_adj = lattice.number_adjacent();
    //std::cout << "N_adj : " << N_adj << std::endl;
    int num_particles = lattice.set_num_particles(Ns);
    //std::cout << "num_particles : " << num_particles << std::endl;
    std::vector<double> z(2*num_particles);
    std::size_t seed = 1234;
    std::mt19937 mt(seed);
    std::uniform_real_distribution<> uniform_rand(0,2*M_PI);

    omp_set_num_threads(N_parallel_max);
    #pragma omp parallel 
    {
   // #pragma omp single
   // {
   // std::cout << "[initialize z : now threads : " << omp_get_num_threads() << "]\n";
   // }
    #pragma omp for
    for(int i = 0 ; i < 2*num_particles; ++i){
      z[i] = uniform_rand(mt);
    }
    }

    //target class def
    //std::cout<< "generate instance by the classes" << std::endl;
    double J = 1.0;
    double T = 1.0;
    hamiltonian::classicalXY_FullyConnected_Parallel ham(num_particles,J,N_adj);
    integrator::velver_parallel integrator(2*num_particles);
    ham.set_meanfield((double)num_particles, 0.0); // only Fully-connected
    //target class def end

    ofs << Ns;
    for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){
      omp_set_num_threads(N_parallel);
      double return_value = 0.0;
      
     // #pragma omp parallel
     // {
     // #pragma omp single
     // { 
     //   std::cout << "[in parallel : now threads : " << omp_get_num_threads() << "]\n";
     // }
     // }

      double z_temp = z[0];
      std::chrono::system_clock::time_point  start, end;
      start = std::chrono::system_clock::now();
      //target behavior
      integrator.step(0,0.01,z,ham);
      //target hehavior end
      end = std::chrono::system_clock::now();
      double result_value = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); 
      return_value = z_temp - z[0];

     // std::cout << " return value : " << return_value << std::endl;
      ofs << " " << result_value; 

      if(i == N_iterate - 1) ofsr << N_parallel << " " << result_value << std::endl;
    }
    ofs << std::endl;
  }
  std::cout << "[finish test]\n"; 
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
