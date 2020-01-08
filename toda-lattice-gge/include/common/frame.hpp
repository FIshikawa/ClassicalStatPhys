#include <mpi.h>
#include <clstatphys/tools/data_recorder.hpp>

int main(int argc, char **argv){
  int input_counter = 1;
  // set parameters
  Settings settings(argc, argv, input_counter);
  int process_id = settings.process_id;
  tools::DataRecorder dataput(settings.condition_dat);
  dataput << "[process id : " << process_id << "] works " << std::endl;

  if(process_id == 0) dataput.time_tag() << " start time " << std::endl;

  // declare parameters
  settings.declare(dataput);

  // set physical quantities
  if(process_id == 0) dataput.time_tag() << " define physical quantities : start" << std::endl;
  PhysicalQuantities physical_quantities(settings);
  if(process_id == 0) dataput.time_tag() << " define physical quantities : end" << std::endl;

  // sampling proecss
  if(process_id == 0) dataput.time_tag() << " sampling : start" << std::endl;
  SamplingProcess(settings, physical_quantities);
  if(process_id == 0) dataput.time_tag() << " sampling : finish" << std::endl;

  // finalize 
  if(process_id == 0) dataput.time_tag() << " finalize : start" << std::endl;
  Finalize(settings, physical_quantities);
  if(process_id == 0) dataput.time_tag() << " finalize : end" << std::endl;
  
  // output results
  if(process_id == 0) dataput.time_tag() << " output : start" << std::endl;
  if(process_id == 0) physical_quantities.output(settings,dataput);
  if(process_id == 0) dataput.time_tag() << " output : end" << std::endl;

  if(process_id == 0) dataput.time_tag() << " end time " << std::endl; 

  return 0;
}
