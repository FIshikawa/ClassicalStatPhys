#include<vector>
#include<alps/accumulators.hpp>
#include<iostream>
#include<random>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[test alps accumulator]\n";

  //set random
  unsigned long seed = 1234;
  std::mt19937 mt(seed);
  std::uniform_real_distribution<> uniform_rand(0.0,1.0);
  //end random

  //set accumulator_set 
  alps::accumulators::accumulator_set  accum_set;
  accum_set  << alps::accumulators::LogBinningAccumulator<double>("square");
  //end accumu set 

  //iteration
  for(int i = 0 ; i < 1e6; ++i){
    //accum_set["square"] << 1.0;
    accum_set["square"] << uniform_rand(mt) * uniform_rand(mt);
  }

  alps::accumulators::result_set results(accum_set);
  std::cout << "reslut : " << results["square"] << std::endl;
  std::cout << "reslut : " << results["square"].mean<double>() << " +- " << results["square"].error<double>() << std::endl;


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
