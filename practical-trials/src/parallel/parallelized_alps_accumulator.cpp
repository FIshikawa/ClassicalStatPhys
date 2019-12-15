#include<vector>
#include<unordered_map>
#include<random>
#include<alps/accumulators.hpp>
#include<omp.h>
#include<iostream>
#include<memory>

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
  std::unordered_map<std::string, double> results{};
  std::vector<std::string> physical_value_list = {"Spin", "Energy"};
  for(auto key : physical_value_list){
    results[key] = 0.0;
    results["error_" + key] = 0.0;
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

  std::vector<alps::accumulators::accumulator_set> physical_values(N_parallel);  
  #pragma omp parallel
  {
  int tid = omp_get_thread_num();
  //set physical value measure
  int num_threads = omp_get_num_threads();
  for (int p = 0; p < num_threads; ++p) {
    if (p == tid){
      for(auto key : physical_value_list){
        physical_values[tid] << alps::accumulators::LogBinningAccumulator<double>(key);
      }
    #pragma omp barrier
    }
  }
  //end set physical value measure 

  //set random
  std::uniform_real_distribution<> uniform_rand(0.0,1.0);
  std::mt19937& mt = *mts[tid];
  //end set

  //substitute values
  for(int i = 0; i < N_iterate; ++i){
      physical_values[tid]["Spin"] << uniform_rand(mt);
      physical_values[tid]["Energy"] << uniform_rand(mt);
  }
  //end substitute


  //output resluts
  for (int p = 0; p < num_threads; ++p) {
      if (p == tid){
        alps::accumulators::result_set physical_values_finalize(physical_values[tid]);
        std::cout << "parallel : " << tid << std::endl; 
        for(auto key  : physical_value_list){
          std::cout << " key : " << key << std::endl;
            std::cout <<  "  value : " << physical_values_finalize[key].mean<double>() << " +- " << physical_values_finalize[key].error<double>() << std::endl; 
            results[key] += physical_values_finalize[key].mean<double>() / N_parallel; 
            results["error_" + key] += physical_values_finalize[key].error<double>() / N_parallel;
        }
      }
     #pragma omp barrier
  }
  //end 
  }//end parallelized

  std::cout <<"total results" << std::endl;
  for(auto element : results){
    std::cout << " result : key : " << element.first << std::endl;
    std::cout << "  value : " <<  element.second << std::endl; 
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
