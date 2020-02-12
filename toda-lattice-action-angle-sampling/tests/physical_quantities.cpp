#include <gtest/gtest.h>
#include <clstatphys/tools/data_recorder.hpp>
#include <specific/toda_lattice_periodic_boundary_hybrid_monte_carlo.hpp>
#include <common/physical_quantities.hpp>
#include <common/sampling_process.hpp>


TEST(PhysicalQuantitesTest, BasicTest){
  char *argv[32] = {"./test", "10", "10", "10", "100", "10",
                    "1.0", "1", "4",
                    "10", "10", "0.1", "1.0",
                    "2.0", "2.0"};
  int argc = 32;
  int input_counter = 0;
  Settings settings(argc, argv, input_counter);
  settings.declare(std::cout);
  PhysicalQuantities physical_quantities(settings);
  std::cout << "complete : construction " << std::endl;

  int num_particles = settings.num_particles;
  std::vector<double> z(2*num_particles,1.0);
  int measure_step = 0;
  physical_quantities.measure(z,measure_step);
  std::cout << "complete : measure step " << std::endl;

  PhysicalQuantities physical_quantities_temp(settings);
  SamplingProcess(settings, physical_quantities_temp);
  std::cout << "complete : sampling process" << std::endl;
  tools::DataRecorder dataput(settings.condition_dat);  
  physical_quantities_temp.output(settings, dataput);
}

