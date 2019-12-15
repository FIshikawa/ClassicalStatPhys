#include<omp.h>
#include<vector>
#include<random>
#include<fstream>
#include<string>
#include<chrono>
#include<memory>
#include<iostream>
#include<iomanip>

int main(int argc, char **argv){

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::string condition_dat = "condition.dat";
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
  int Ns_init = 1e5;
  int Ns_max =  1e7;
  int dNs = (Ns_max - Ns_init)/N_iterate;

  for(int i = 0; i < N_iterate ; ++i){
    int Ns = Ns_init + dNs * i;
    int num_particles = Ns;
    std::vector<double> z(2*num_particles);
    unsigned long seed = 1234;
    std::minstd_rand seed_gen(seed);
    std::vector<std::shared_ptr<std::mt19937> > mts(N_parallel_max);

    omp_set_num_threads(N_parallel_max);
    #pragma omp parallel 
    {
    #pragma omp single
    {
    std::cout << "[initialize z : now threads : " << omp_get_num_threads() << "]\n";
    }
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    for (int p = 0; p < num_threads; ++p) {
      if (p == tid)
        mts[p].reset(new std::mt19937(seed_gen()));
        #pragma omp barrier
    }
    }
 
    #pragma omp parallel 
    {
    int tid = omp_get_thread_num();
    std::uniform_real_distribution<> uniform_rand(0.0,1.0);
    std::mt19937& mt = *mts[tid];
    #pragma omp for
    for(int i = 0 ; i < 2*num_particles; ++i){
      z[i] = uniform_rand(mt);
    }
    }

    //target class def
    auto do_experiment = [&z,&num_particles](){
      double total_z = 0.0;
      #pragma omp parallel for reduction(+:total_z)
      for(int i = 0 ; i < 2*num_particles; ++i){
        total_z += z[i];
      }
      return total_z;
    };
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
      return_value = do_experiment();
      //target behavior end
      end = std::chrono::system_clock::now();
      double result_value = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); 

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
