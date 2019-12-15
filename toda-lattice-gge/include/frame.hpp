#include <clstatphys/tools/data_recorder.hpp>

using RandomGenerator = std::mt19937_64;

int main(int args, char **argv){
  // set parameters
  Parameters parameters;
  parameters.set(args, argv);

  // set dataput
  tools::DataRecorder dataput(parameters.common.condition_dat); 
  dataput.time_tag() << "start time " << std::endl;

  parameters.declare(dataput);

  // set modules 
  Modules modules(parameters);

  // set physical quantities
  PhysicalQuantities physical_quantities(parameters);

  // randam generetor set 
  std::size_t seed = 1234; 

  // seed set
  std::minstd_rand seed_gen(seed);
  RandomGenerator mt(seed_gen());

  // sampling proecss
  SamplingProcess(parameters, modules, physical_quantities, mt, dataput);

  // set result data
  ResultData result_data(parameters, physical_quantities, dataput);

  // finalize 
  Finalize(parameters, physical_quantities, result_data, dataput);
  
  // output results
  OutputResults(parameters, result_data, dataput);

  dataput.time_tag() << " end time " << std::endl; 

  return 0;
}
