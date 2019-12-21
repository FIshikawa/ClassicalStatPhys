#include <boost/lexical_cast.hpp>
#include <clstatphys/integrator/yoshida_4th.hpp>
#include <clstatphys/tools/time_measure.hpp>
#include <clstatphys/physics/toda_lax_form.hpp>
#include <clstatphys/lattice/chain.hpp>

using Lattice = lattice::Chain;
using Integrator = integrator::Yoshida4th;

struct SettingsCommon{
  int Ns = 10; // lenght of chain
  int N_loop = 1; //loop number
  int Ns_observe = Ns; //observed particle 
  int N_time = 1000; //numer of time step
  int N_normalmode = 5; //numer of time step
  int N_time_measure = 10; //measuring of per time or linear
  int k_initial = 0; //start wave vector filled by initialization 
  int n_bin = 1000;
  int order = 0;
  int num_particles = 10; //nummer of particles
  int N_mpi_parallel = 2; // number of mpi parallel
  int N_each = 10;// number of each parallel
  int N_adj = 2;// number of adjacent spins
  int N_total_data = 100; // number of data
  double dt = 0.01;//time step length
  double E_initial = 1.0; // initial energy
  double t = 10; // Whole time
  std::string plot_scale{"linear"};
  std::string condition_dat{"condi.dat"};
  std::string result_directory{"./"};
  tools::TimeMeasure time_measure_;
  

  SettingsCommon(int argc, char **argv, int & input_counter){
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
    if (argc > input_counter) N_time =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) t =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) E_initial =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) k_initial =        boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) N_normalmode =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (Ns < N_normalmode + k_initial){
      std::cerr << "k_initial + N_noramalmode should be lower than Ns" << std::endl;
      std::exit(1);
    } 
    if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) N_time_measure =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) plot_scale =       boost::lexical_cast<std::string>(argv[input_counter]);++input_counter;

    num_particles = Lattice::set_num_particles(Ns); //nummer of particles
    Lattice lattice_t(Ns);
    N_adj = lattice_t.number_adjacent();

    N_each = N_loop / N_mpi_parallel; // number of each parallel
    dt = t / N_time; //time step length
    condition_dat = result_directory + condition_dat;

    time_measure_ = tools::TimeMeasure(t, dt, N_time_measure, plot_scale);
    N_total_data = time_measure_.number_of_total_data();

  }

  template <class Dataput>
  inline void declare(Dataput & dataput){
    dataput << "Declare Configulation " << std::endl 
            << "Explore Toda lattice Generalized Gibbs Model" << std::endl
            <<  Integrator::name() << std::endl

            << "Number of patricles : num_particles = " << num_particles << std::endl
            << "Length of whole system : Ns = " << Ns << std::endl
            << "Length of observe system : Ns_observe = " << Ns_observe << std::endl
            << "Number of time for measure: N_time_measure = " << N_time_measure << std::endl
            << "Number of result data (time step): N_total_data = " << N_total_data << std::endl
            << "Order of result data : order = " << time_measure_.order() << std::endl
            << "Developing ime : t = " << t << std::endl
            << "Number of time steps : N_time = " << N_time << std::endl
            << "Time interbal : dt = " << dt << std::endl
            << "Energy initial : E_initial =" << E_initial << std::endl
            << "Number of non-zero energy normal modes : N_normalmode =" << N_normalmode << std::endl
            << "Start wave vector filled by initialization : k_initial =" << k_initial << std::endl
            << "Number of MPI parallelization : N_mpi_parallel =" << N_mpi_parallel << std::endl
            << "Number of loop : N_loop =" << N_loop << std::endl
            << "Each thread loop : N_each = " << N_each << std::endl 
            << "Plot scale : plot_scale = " << plot_scale<< std::endl
            << "Result directory : " << result_directory << std::endl;
  }

  Integrator integrator(){return Integrator(2*num_particles);}

  tools::TimeMeasure time_measure(){
    return tools::TimeMeasure(t, dt, N_time_measure, plot_scale);
  }

  integrable::TodaLaxForm toda_lax_form(){
    Lattice lattice_t(Ns);
    std::vector<std::vector<int> > pair_table(
                                              num_particles,
                                              std::vector<int>(N_adj)
                                              );
    lattice_t.create_table(pair_table);
    return integrable::TodaLaxForm(num_particles,1.0,1.0,pair_table,N_adj);
  }
};
