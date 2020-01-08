#include <gtest/gtest.h>
#include <specific/harmonic_oscillator_periodic_boundary_thermalization_from_point_heating.hpp>
#include <common/sampling_process.hpp>
#include <common/physical_quantities.hpp>

TEST(SamplingProcessTest, MethodsTest){
  char *argv[32] = {"./test", "10", "4","100","10.0","1.0","1","4","10","10","10","linear","4.0","2.0"};
  int argc = 14;
  int input_counter = 0;  
  Settings settings(argc, argv, input_counter);
  settings.declare(std::cout);
  PhysicalQuantities physical_quantities(settings);
  SamplingProcess(settings, physical_quantities);
}

