#include <gtest/gtest.h>
#include <specific/toda_lattice_periodic_boundary_thermalization_from_normalmode_ensemble.hpp>
#include <common/sampling_process.hpp>
#include <common/physical_quantities.hpp>
#include <common/finalize.hpp>

TEST(FinalizeTest, MethodsTest){
  char *argv[32] = {"./test", "10", "10","100","10.0","100","100","linear","10",
                    "1.0","1","4","2.0","2.0"};
  int argc = 15;
  int input_counter = 0;  
  Settings settings(argc, argv, input_counter);
  settings.declare(std::cout);
  PhysicalQuantities physical_quantities(settings);
  SamplingProcess(settings, physical_quantities);
  Finalize(settings, physical_quantities);
}

