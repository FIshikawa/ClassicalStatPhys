#include <iostream>
#include <clstatphys/physics/partition_function_classical_xy_fullyconnected.hpp>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  double J = 1.0;
  int num_particles = 10;
  double E,T;
  exactsolution::PartitionFunctionClassicalXYFullyConnected calc_partition(J, num_particles);
  std::cout << "[test for exactcalc.hpp]\n";
  std::cout << "J = " << J << std::endl
            << "num_paritcles = " << num_particles << std::endl;
  for(int i = 0 ; i < 200; ++i){
    T = 0.01 * i + 0.01;
    E = 0.5 * T;
    //std::cout << " T : " << T 
    //          << " exact M = " << calc_partition.M(T)
    std::cout << " E : " << E
              << " exact T = " << calc_partition.T(E)
              << std::endl;
  }
  std::cout << "[end test]\n"; 
  



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
