#ifndef SETTINGS_COMMON_HPP
#define SETTINGS_COMMON_HPP

#include <mpi.h>
#include <boost/lexical_cast.hpp>
#include <clstatphys/physics/toda_lax_form.hpp>
#include <clstatphys/lattice/chain.hpp>

using Lattice = lattice::Chain;

struct SettingsCommon{
  int num_iterations = 32; 
  int Ns = 10; // lenght of chain
  int N_loop = 1; //loop number
  int Ns_observe = Ns; //observed particle 
  int N_mc = 10; //numer of time step
  int order = 0;
  int num_particles = 10; //nummer of particles
  int N_mpi_parallel = 2; // number of mpi parallel
  int N_each = 10;// number of each parallel
  int N_adj = 2;// number of adjacent spins
  int N_total_data = 10; // number of data
  int process_id = 1;
  int name_length;
  int num_threads;
  int mpi_error = 1;
  int mpi_error_end = 1;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  std::string condition_dat{"condi.dat"};
  std::string result_directory{"./"};


  SettingsCommon(int argc, char **argv, int & input_counter){
    if(mpi_error > 0){
      mpi_error = MPI_Init(&argc, &argv);
      MPI_Comm_size(MPI_COMM_WORLD, &N_mpi_parallel);
      MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
      MPI_Get_processor_name(processor_name, &name_length);
    }
    set(argc, argv, input_counter);
  }
  SettingsCommon() = default;

  inline void set(int argc, char **argv, int & input_counter){
    if (argc > input_counter) result_directory = boost::lexical_cast<std::string>(argv[input_counter]); ++input_counter;
    if (argc > input_counter) Ns =               boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) Ns_observe =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (Ns < Ns_observe){
      std::cerr << "Ns_observe should be lower than Ns" << std::endl;
      std::exit(1);
    }
    if (argc > input_counter) N_mc =             boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) num_iterations =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;

    num_particles = Lattice::set_num_particles(Ns); //nummer of particles
    Lattice lattice_t(Ns);
    N_adj = lattice_t.number_adjacent();

    N_each = N_loop / N_mpi_parallel; // number of each parallel
    condition_dat = result_directory + condition_dat;

    N_total_data = N_mc;
  }

  template <class Dataput>
  inline void declare(Dataput & dataput){
    dataput << "<<Common Settings>>" << std::endl 
            << "  Sampling GGE of Toda lattice by Monte Carlo" << std::endl
            <<  "  " << Lattice::name() << std::endl

            << "  Number of patricles : num_particles = " << num_particles << std::endl
            << "  Length of whole system : Ns = " << Ns << std::endl
            << "  Length of observe system : Ns_observe = " << Ns_observe << std::endl
            << "  Number of result data (time step): N_total_data = " << N_total_data << std::endl
            << "  Number of monte carlo steps : N_mc = " << N_mc << std::endl
            << "  Number of MPI parallelization : N_mpi_parallel =" << N_mpi_parallel << std::endl
            << "  Number of loop : N_loop =" << N_loop << std::endl
            << "  Each thread loop : N_each = " << N_each << std::endl 
            << "  Number of iteration for action variables : num_iterations = " << num_iterations << std::endl
            << "  Result directory : " << result_directory << std::endl;
  }

  integrable::TodaLaxForm toda_lax_form(){
    Lattice lattice_t(Ns);
    std::vector<std::vector<int> > pair_table(
                                              num_particles,
                                              std::vector<int>(N_adj)
                                              );
    lattice_t.create_table(pair_table);
    return integrable::TodaLaxForm(num_particles,1.0,1.0,"periodic");
  }

  inline void finalize(){mpi_error = MPI_Finalize();};

}; // end struct define

#endif
