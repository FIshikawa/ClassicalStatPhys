#include<omp.h>
#include<vector>
#include<random>
#include<iostream>
#include<fstream>
#include<string>
#include<memory>
#include<chrono>

#include "initializer.hpp"
#include "Lattice.hpp"

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

    std::cout<< "create lattice" << std::endl;
    Lattice::FullyConnected lattice(Ns);
    int N_adj = lattice.number_adjacent();
    std::cout << "N_adj : " << N_adj << std::endl;
    int num_particles = lattice.set_num_particles(Ns);
    std::cout << "num_particles : " << num_particles << std::endl;
    std::vector<double> z(2*num_particles);


    //set random value for  parallelization
    unsigned long seed = 1234;
    std::minstd_rand seed_gen(seed);
    std::vector<std::shared_ptr<std::mt19937> > mts(N_parallel_max);
    omp_set_num_threads(N_parallel_max);
    #pragma omp parallel 
    {
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    for (int p = 0; p < num_threads; ++p) {
      if (p == tid)
        mts[p].reset(new std::mt19937(seed_gen()));
        #pragma omp barrier
    }
    }
    //set random value for  parallelization end 

    auto summation = [&z](){ 
      double total = 0;
      #pragma omp parallel for reduction(+:total)
      for(int i= 0 ; i < z.size(); ++i){
        total += z[i] * z[i];
      }
      return total;
    };

    //target class def
    std::cout<< "generate instance by the classes" << std::endl;
    double J = 1.0;
    double T = 1.0;
    initializer::WaterBag_Prallelized initializer(J,T,num_particles);
    int counter = 0;
    //target class def end

    ofs << Ns;
    for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){
      omp_set_num_threads(N_parallel);
      double return_value = 0.0;
      
      #pragma omp parallel
      {
      #pragma omp single
      { 
        std::cout << "[in parallel : now threads : " << omp_get_num_threads() << "]\n";
      }
      }

      std::chrono::system_clock::time_point  start, end;
      start = std::chrono::system_clock::now();
      //target behavior 
      initializer.set_initial_state(z,mts);
      //target behavior end 
      end = std::chrono::system_clock::now();
      double result_value = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); 
      return_value = summation();

      std::cout << " return value : " << return_value << std::endl;
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
