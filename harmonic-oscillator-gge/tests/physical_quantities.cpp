#include <gtest/gtest.h>
#include <clstatphys/tools/data_recorder.hpp>
#include <specific/harmonic_oscillator_periodic_boundary_thermalization_from_point_heating.hpp>
#include <common/physical_quantities.hpp>
#include <common/sampling_process.hpp>


TEST(PhysicalQuantitesTest, BasicTest){
  char *argv[32] = {"./test", "10", "10","100","10.0","1.0","100","100","linear",
                    "2.0","2.0","5"};
  int argc = 14;
  int input_counter = 0;
  Settings settings(argc, argv, input_counter);
  PhysicalQuantities physical_quantities(settings);

  int num_particles = settings.num_particles;
  std::vector<double> z(2*num_particles,1.0);
  int measure_step = 0;
  physical_quantities.measure(z,measure_step);

  PhysicalQuantities physical_quantities_temp(settings);
  SamplingProcess(settings, physical_quantities_temp);
  tools::DataRecorder dataput(settings.condition_dat);  
  physical_quantities_temp.output(settings, dataput);
}

