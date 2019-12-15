#include<vector>
#include<unordered_map>
#include<random>
#include<omp.h>
#include<iostream>
#include<memory>

#include"accumulator.hpp"

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[test alea ]\n";

  //set parameter
  int N_time = 2;
  int N_parallel = omp_get_num_procs();
  int N_iterate = 1e6;
  //end 

  //set results and error
  std::unordered_map<std::string, std::vector<double> > results{};
  std::vector<std::string> physical_value_list = {"Spin", "Energy"};
  for(auto key : physical_value_list){
    results[key].resize(N_time);
    results["error_" + key].resize(N_time);
    for(int i = 0 ; i < N_time ; ++i){
      results[key][i] = 0.0;
      results["error_" + key][i] = 0.0;
    }
  }
  //end set 

  //set random value for  parallelization
  unsigned long seed = 1234;
  std::minstd_rand seed_gen(seed);
  std::vector<std::shared_ptr<std::mt19937> > mts(N_parallel);
  omp_set_num_threads(N_parallel);
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

  #pragma omp parallel
  {
  //set physical value measure
  std::unordered_map<std::string, std::vector<stat::accumulator> > physical_values;
  for(auto key : physical_value_list){
    //physical_values[key].resize(N_time);
    physical_values[key] = std::vector<stat::accumulator>(N_time);
  }
  //end set physical value measure 

  //set random
  int tid = omp_get_thread_num();
  std::uniform_real_distribution<> uniform_rand(0.0,1.0);
  std::mt19937& mt = *mts[tid];
  //end set

  //substitute values
  for(int i = 0; i < N_iterate; ++i){
    for(int step = 0 ; step < N_time; ++step){
      physical_values["Spin"][step] << uniform_rand(mt);
      physical_values["Energy"][step] << uniform_rand(mt);
    }
  }
  //end substitute

  //output resluts
  int num_threads = omp_get_num_threads();
  for (int p = 0; p < num_threads; ++p) {
      if (p == tid){
        std::cout << "parallel : " << tid << std::endl; 
        for(auto key  : physical_value_list){
          std::cout << " key : " << key << std::endl;
          for(int i = 0 ; i < N_time; ++i){
            std::cout << "  time : " << i << ", value : " << physical_values[key][i].mean() << " +- " << physical_values[key][i].error() << std::endl; 
            results[key][i] += physical_values[key][i].mean() / N_parallel; 
            results["error_" + key][i] += physical_values[key][i].error() / N_parallel;
          }
        }
      }
     #pragma omp barrier
  }
  //end 
  }//end parallelized

  std::cout <<"total results" << std::endl;
  for(auto element : results){
    std::cout << " result : key : " << element.first << std::endl;
    for(int i = 0 ; i < N_time ; ++i) std::cout << "  time : " << i << ", value : " <<  element.second[i] << std::endl; 
  }


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
