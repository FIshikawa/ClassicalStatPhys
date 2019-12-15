#include<vector>
#include<alps/alea.hpp>
#include<unordered_map>
#include<random>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[test alea ]\n";
  unsigned long seed = 1234;
  std::minstd_rand seed_gen(seed);
  std::mt19937 mt(seed_gen());
  std::uniform_real_distribution<> uniform_rand(0,1.0);

  std::unordered_map<std::string, std::vector<alps::alea::cov_acc<double> > > physical_values{};  
  std::vector<std::string> physical_value_list = {"Spin", "Energy"};
  std::cout << "[check physical_value dict ]\n"; 
  for(auto key : physical_value_list){
    physical_values[key].resize(2);
    for(int i = 0; i < 2 ; ++i){
      physical_values[key][i].set_size(1);
    }
  }

  for(int i = 0; i < 1e6; ++i){
    physical_values["Spin"][0] << uniform_rand(mt);
    physical_values["Energy"][0] << 1.0;
  }
  std::unordered_map<std::string, std::vector<alps::alea::cov_result<double> > > physical_values_result{}; 
  for(auto key  : physical_value_list){
    physical_values_result[key].resize(2);
    for(int i = 0; i < 2 ; ++i){
      physical_values_result[key][i]= physical_values[key][i].finalize();
    }
    std::cout << "[" << key << "] "
              << "value : " << physical_values_result[key][0].mean() << " +- " << physical_values_result[key][0].var()
              << " cov : " << physical_values_result[key][0].cov() << std::endl;
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
