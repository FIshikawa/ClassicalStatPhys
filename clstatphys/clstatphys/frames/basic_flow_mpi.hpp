#include <mpi.h>
#include <clstatphys/tools/data_recorder.hpp>

using RandomGenerator = std::mt19937_64;

struct SettingMPI{
  int N_mpi_parallel;
  int process_id;
  int name_length;
  int num_threads;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  inline void set(int args, char **argv){
    MPI_Comm_size(MPI_COMM_WORLD, &N_mpi_parallel);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    MPI_Get_processor_name(processor_name, &name_length);
  }

  inline void declare(tools::DataRecorder & dataput){
    dataput << "Number of MPI parallelization : N_mpi_parallel =" << N_mpi_parallel << std::endl;
  }
};

struct Settings{
  SettingCommon common;
  SettingDepend depend;
  SettingMPI    mpi;

  inline void set(int args, char **argv){
    int input_counter = 0;
    common.mpi(args, argv);
    common.set(args, argv, input_counter);
    depend.set(args, argv, common, input_counter);
  }

  inline void declare(tools::DataRecorder & dataput){
    common.declare(dataput);
    mpi.declare(dataput);
    depend.declare(dataput);
  }
};

struct Results{
  ResultsCommon results_common;
  ResultsDepend results_depend;

  inline void set(){
    results_common.set();
    results_depend.set();
  }

  template <struct Settings>
  inline void output(Settings & settings, tools::DataRecorder & dataput){
    if(settings.mpi.process_id == 0){
      results_common.output(dataput);
      results_depend.output(dataput);
    }
  }
};

int main(int args, char **argv){
  int mpi_error = MPI_Init(&argc, &argv);
  // set parameters
  Settings settings;
  settings.set(args, argv);

  // set dataput
  tools::DataRecorder dataput(settings.common.condition_dat); 
  if(settings.mpi.process_id == 0)
    dataput.time_tag() << "start time " << std::endl;

   // declare configurations
  settings.declare(dataput);

  // set physical quantities
  PhysicalQuantities physical_quantities(settings);

  // randam generetor set 
  std::size_t seed = 1234; 

  // seed set
  std::minstd_rand seed_gen(seed);
  RandomGenerator mt(seed_gen());

  // sampling proecss
  SamplingProcess(settings, physical_quantities, mt, dataput);

  // set result data
  Results results;
  results.set(settings, physical_quantities, dataput);

  // finalize 
  Finalize(settings, physical_quantities, result_data, dataput);
  
  // output results
  result_data.output(settings, dataput);

  if(parameters.mpi.process_id == 0)
    dataput.time_tag() << " end time " << std::endl; 

  mpi_error = MPI_Finalize();
  return 0;
}
