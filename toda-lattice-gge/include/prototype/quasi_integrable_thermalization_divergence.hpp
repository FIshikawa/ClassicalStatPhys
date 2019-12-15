#include <omp.h>
#include <mpi.h>
#include <cstdlib>
#include <memory>
#include <main/quasi_integrable_thermalization_divergence_functions.hpp>
#include <boost/lexical_cast.hpp>
#include <tools/data_recorder.hpp>
#include <tools/auto_correlation_function.hpp>
#include <tools/histogram.hpp>
#include <tools/accumulator.hpp>
#include <tools/statistical_divergence.hpp>
#include <integrator/yoshida_4th.hpp>

using Integrator = integrator::Yoshida4th;
// NormalModeEnergy, Ensembler, Hamiltonian, Lattice are defined in each src header 


int main(int argc, char **argv){
  //MPI configuration 
  int mpi_error;
  int N_mpi_parallel;
  int process_id;
  int name_length;
  int num_threads;
  mpi_error = MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &N_mpi_parallel);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processor_name, &name_length);

  int Ns = 1024; // lenght of chain
  int N_loop = 1; //loop number
  int Ns_observe = Ns; //observed particle 
  int N_time = 10; //numer of time step
  int N_time_resolve = 1; //numer of time step
  int N_normalmode = Ns/16; //numer of time step
  int k_initial = 1; //start wave vector filled by initialization 
  int n_bin = 1000;
  double E_initial = 0.1*Ns; // initial energy
  double t = 10; // Whole time
  std::string condition_dat = "condi.dat";
  std::string result_directory("./");

  // set common parameters
  int input_counter = 1;
  if (argc > input_counter) result_directory = boost::lexical_cast<std::string>(argv[input_counter]); ++input_counter;
  if (argc > input_counter) Ns =               boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) Ns_observe =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if(Ns < Ns_observe){
    std::cerr << "Ns_observe should be lower than Ns" << std::endl;
    std::exit(1);
  }
  if (argc > input_counter) N_time =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) t =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) E_initial =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) k_initial =        boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_normalmode =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if(Ns < N_normalmode + k_initial){
    std::cerr << "k_initial + N_normalmode should be lower than Ns" << std::endl;
    std::exit(1);
  } 
  if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time_resolve =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;

  // set private parameters
  Parameters parameters;
  ParameterSet(argc, argv, parameters, input_counter);

  int num_particles = Lattice::set_num_particles(Ns); //nummer of particles
  int N_time_measure = N_time / N_time_resolve; //number of time span of histogram
  int N_each = N_loop / N_mpi_parallel; // number of each parallel
  double dt = t / N_time; //time step length
  condition_dat = result_directory + condition_dat;

  tools::DataRecorder dataput(condition_dat); 
  if(process_id == 0){
    dataput.time_tag() << "start time " << std::endl;

    // declare common settings
    dataput << "Declare Configulation " << std::endl 
            << "Explore Non-thermal Fixed Point : (Quasi) Integrable Models" << std::endl
            <<  Hamiltonian::name() << std::endl
            <<  Ensembler::name() << std::endl
            <<  Integrator::name() << std::endl
            <<  Lattice::name() << std::endl;

    // declare private settings 
    ParameterDeclare(parameters,dataput);
    dataput << "Number of patricles : num_particles = " << num_particles << std::endl
            << "Length of whole system : Ns = " << Ns << std::endl
            << "Length of observe system : Ns_observe = " << Ns_observe << std::endl
            << "Number of time resolve for measure: N_time_resolve = " << N_time_resolve << std::endl
            << "Number of time for measure: N_time_measure = " << N_time_measure << std::endl
            << "Developing ime : t = " << t << std::endl
            << "Number of time steps : N_time = " << N_time << std::endl
            << "Time interbal : dt = " << dt << std::endl
            << "Energy initial : E_initial =" << E_initial << std::endl
            << "Number of non-zero energy normal modes : N_normalmode =" << N_normalmode << std::endl
            << "Start wave vector filled by initialization : k_initial =" << k_initial << std::endl
            << "Number of MPI parallelization : N_mpi_parallel =" << N_mpi_parallel << std::endl
            << "Number of loop : N_loop =" << N_loop << std::endl
            << "Each thread loop : N_each = " << N_each << std::endl 
            << "Number of time resolve for hist : N_time_resolve = " << N_time_resolve << std::endl
            << "Result directory : " << result_directory << std::endl;
  }

  dataput << "[process : " << process_id << "] processor name : " << processor_name << std::endl;

  //set results and name vectors
  std::unordered_map<std::string, std::vector<double> > results;  
  std::vector<std::string> physical_value_list = {
                                                  "Velocity",
                                                  "Velocity_total",
                                                  "Displacement",
                                                  "Displacement_total",
                                                  "Energy_total",
                                                  "Kinetic_Energy",
                                                  "Potential_Energy",
                                                  "SpectralEntropy",
                                                  "EffectiveFractionOfModes",
                                                  "JSDivergence_from_Toda",
                                                  "JSDivergence_from_Harmonic",
                                                  "JSDivergence_from_Thermal",
                                                  "L2_from_Toda",
                                                  "L2_from_Harmonic",
                                                  "L2_from_Thermal",
                                                  "HD_from_Toda",
                                                  "HD_from_Harmonic",
                                                  "HD_from_Thermal"
                                                  };
  for(const auto& key : physical_value_list){
    results[key].resize(N_time_measure);
    results["error_" + key].resize(N_time_measure);
    for(int i = 0 ; i < N_time_measure ; ++i){
      results[key][i] = 0.0;
      results["error_"+key][i] = 0.0;
    }
  }

  //set correlatin function 
  std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d{};  
  results_2d["Normalmode_Energy"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));
  results_2d["error_Normalmode_Energy"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));
  results_2d["Normalmode_Dist"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));
  results_2d["error_Normalmode_Dist"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));

  if(process_id == 0) dataput.time_tag() << "start sampling process : time development " << std::endl;

  //lattice set
  Lattice lattice(Ns);
  int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj)); //interaction table 
  lattice.create_table(pair_table);
  std::vector<double> z(2*num_particles);
  
  //randam generetor set 
  std::size_t seed = 1234 * (1 + process_id);

  //ensember set 
  using RandomGenerator = std::mt19937_64;
  std::minstd_rand seed_gen(seed);
  RandomGenerator mt(seed_gen());
  Ensembler ensembler(num_particles,k_initial,N_normalmode,E_initial);

  //hamiltonian set 
  Hamiltonian hamiltonian = HamiltonianSet(parameters, num_particles, pair_table, N_adj);
  //integratro set 
  Integrator integrator(2*num_particles);

  //measure physical values set
  std::unordered_map<std::string, std::vector<tools::Accumulator> > physical_values;
  for(const auto& key : physical_value_list){
    physical_values[key].resize(N_time_measure);
    for(int i = 0; i < N_time_measure; ++i) physical_values[key][i].reset(); 
  }

  //set normalmode energy
  std::vector<std::vector<tools::Accumulator> > normalmode_energy(N_time_measure,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_time_measure ; ++step){
    normalmode_energy[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) normalmode_energy[step][i].reset();
  }
  std::vector<std::vector<tools::Accumulator> > normalmode_dist(N_time_measure,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_time_measure ; ++step){
    normalmode_dist[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) normalmode_dist[step][i].reset();
  }

  int observe = (num_particles-1) / 2;
  int measure_step = 0;
  double *x = &z[0];
  double *v = &z[num_particles];
  double pt = 0.0;
  double v_total,x_total,sum_spectral,average_spectral,spectral_entropy;
  double js_divergece_from_toda, js_divergece_from_harmonic, js_divergece_from_thermal,
         h_distance_from_toda, h_distance_from_harmonic, h_distance_from_thermal,   
         l2_norm_from_toda, l2_norm_from_harmonic, l2_norm_from_thermal;

  std::vector<double> normalmode_energy_temp(num_particles);
  std::vector<double> normalmode_dist_temp(num_particles);
  tools::StatisticalDivergence f_divergence;

  // set reference specturm distributions 
  std::vector<double> thermal_spectrum(num_particles,1.0/(num_particles-1));
  thermal_spectrum[0] = 0.0;
  std::vector<double> harmonic_spectrum(num_particles,0.0);
  std::vector<double> toda_spectrum(num_particles,0.0);
  SetHarmonicDist(harmonic_spectrum, N_normalmode, k_initial); 
  SetTodaDist(toda_spectrum, E_initial);

  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){

    // set initial state
    ensembler.set_initial_state(z,mt);

    //time develop
    for(int step = 0; step < N_time; ++step){  
      //set total values 
      if((step % N_time_resolve) == 0){
        measure_step = step / N_time_resolve;

        v_total = 0.0;
        x_total = 0.0;
        for(int i = 0; i < num_particles; ++i){
          v_total += v[i] / num_particles;
          x_total += x[i] / num_particles;
        } 
        //set normalmode energy
        NormalModeEnergy(z,normalmode_energy_temp);

        // calc spectral entropy
        sum_spectral = 0;
        spectral_entropy= 0;
        average_spectral= 0;
        for(int i = 0; i < num_particles; ++i) sum_spectral += normalmode_energy_temp[i];
        for(int i = 0; i < num_particles; ++i){
          average_spectral = normalmode_energy_temp[i] / sum_spectral;
          normalmode_dist_temp[i] = average_spectral;
          if(average_spectral > 0) spectral_entropy -=  average_spectral * std::log(average_spectral);
        }
        // calc spectral density divergence
        f_divergence.calc(normalmode_dist_temp, toda_spectrum);
        js_divergece_from_toda = f_divergence.jensen_shannon_divergence();
        h_distance_from_toda = f_divergence.hellinger_distance();
        l2_norm_from_toda = f_divergence.l2_norm();

        f_divergence.calc(normalmode_dist_temp, harmonic_spectrum);
        js_divergece_from_harmonic = f_divergence.jensen_shannon_divergence();
        h_distance_from_harmonic = f_divergence.hellinger_distance();
        l2_norm_from_harmonic = f_divergence.l2_norm();

        f_divergence.calc(normalmode_dist_temp, thermal_spectrum);
        js_divergece_from_thermal = f_divergence.jensen_shannon_divergence();
        h_distance_from_thermal = f_divergence.hellinger_distance();
        l2_norm_from_thermal = f_divergence.l2_norm();

        //input data 
        physical_values["Velocity"][measure_step] << v[observe];
        physical_values["Velocity_total"][measure_step] << v_total;
        physical_values["Displacement"][measure_step] << x[observe];
        physical_values["Displacement_total"][measure_step] << x_total;
        physical_values["Energy_total"][measure_step] <<  hamiltonian.energy(pt,z) / num_particles;
        physical_values["Kinetic_Energy"][measure_step] << hamiltonian.kinetic_energy(pt,z) / num_particles;
        physical_values["Potential_Energy"][measure_step] << hamiltonian.potential_energy(pt,z) / num_particles;
        physical_values["SpectralEntropy"][measure_step] << spectral_entropy;
        physical_values["EffectiveFractionOfModes"][measure_step] << std::exp(spectral_entropy) / num_particles;
        physical_values["JSDivergence_from_Toda"][measure_step] << js_divergece_from_toda;
        physical_values["JSDivergence_from_Harmonic"][measure_step] << js_divergece_from_harmonic;
        physical_values["JSDivergence_from_Thermal"][measure_step] << js_divergece_from_thermal;
        physical_values["L2_from_Toda"][measure_step] << l2_norm_from_toda;
        physical_values["L2_from_Harmonic"][measure_step] << l2_norm_from_harmonic;
        physical_values["L2_from_Thermal"][measure_step] << l2_norm_from_thermal;
        physical_values["HD_from_Toda"][measure_step] << h_distance_from_toda;
        physical_values["HD_from_Harmonic"][measure_step] << h_distance_from_harmonic;
        physical_values["HD_from_Thermal"][measure_step] << h_distance_from_thermal;

        //input normalmode energy
        for(int i = 0; i < num_particles; ++i){
          normalmode_energy[measure_step][i] << normalmode_energy_temp[i];
          normalmode_dist[measure_step][i]   << normalmode_dist_temp[i];
        }
      }

      //time develp
      pt += dt;
      integrator.step(pt, dt, z, hamiltonian);
    }//end time step
  }//end loop

  // correct data each process
  dataput.time_tag() << "[" << process_id <<"th process] start result correct by self" << std::endl;
  for(const auto& key  : physical_value_list){
    for(int i = 0 ; i < N_time_measure; ++i){
      results[key][i] = physical_values[key][i].sum1();
      results["error_" + key][i] = physical_values[key][i].sum2();
    }
  }
  for(int step = 0 ; step < N_time_measure; ++step){
    for(int i = 0 ; i < num_particles; ++i){
      results_2d["Normalmode_Energy"][step][i] = normalmode_energy[step][i].sum1(); 
      results_2d["error_Normalmode_Energy"][step][i] = normalmode_energy[step][i].sum2();
      results_2d["Normalmode_Dist"][step][i] = normalmode_dist[step][i].sum1(); 
      results_2d["error_Normalmode_Dist"][step][i] = normalmode_dist[step][i].sum2();
    }
  }
  dataput.time_tag() << "[" << process_id <<"th process] finish spatial: ready for transfer" << std::endl;


  // correct data to process id 0 
  if(process_id == 0) dataput.time_tag() << "start total result correct" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start normal results" << std::endl;
  std::vector<double> reduced_vector(N_time_measure);
  for(const auto& key  : physical_value_list){
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_time_measure, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_time_measure; ++i) results[key][i] = reduced_vector[i] / N_loop;
    mpi_error = MPI_Reduce(&results["error_" + key][0], &reduced_vector[0], N_time_measure, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_time_measure; ++i) 
      results["error_" + key][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector[i] / N_loop - results[key][i] * results[key][i]) / (N_loop - 1)) : 0;
  }
  if(process_id == 0) dataput.time_tag() << "end normal results" << std::endl;

  std::vector<double> reduced_vector_1(num_particles);
  if(process_id == 0) dataput.time_tag() << "start 2d results" << std::endl;
  for(int step = 0 ; step < N_time_measure; ++step){
    mpi_error = MPI_Reduce(&results_2d["Normalmode_Energy"][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) results_2d["Normalmode_Energy"][step][i] = reduced_vector_1[i] / N_loop; 
    mpi_error = MPI_Reduce(&results_2d["error_Normalmode_Energy"][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) 
      results_2d["error_Normalmode_Energy"][step][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector_1[i] / N_loop
                                           - results_2d["Normalmode_Energy"][step][i] * results_2d["Normalmode_Energy"][step][i])/(N_loop - 1)) : 0; 
    mpi_error = MPI_Reduce(&results_2d["Normalmode_Dist"][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) results_2d["Normalmode_Dist"][step][i] = reduced_vector_1[i] / N_loop; 
    mpi_error = MPI_Reduce(&results_2d["error_Normalmode_Dist"][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) 
      results_2d["error_Normalmode_Dist"][step][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector_1[i] / N_loop
                                           - results_2d["Normalmode_Dist"][step][i] * results_2d["Normalmode_Dist"][step][i])/(N_loop - 1)) : 0; 

  }
  if(process_id == 0) dataput.time_tag() << "end 2d results" << std::endl;

  if(process_id == 0) dataput.time_tag() << "end total result correct" << std::endl;


  // data output by process_id = 0
  if(process_id == 0){
    dataput.time_tag() << "end sampling process" << std::endl;

    //set measure time
    std::unordered_map< std::string, std::vector<double> > domain_time;
    domain_time["time"] = std::vector<double>(N_time_measure);
    double pt = 0;
    double counter = 0;
    for(int step = 0; step < N_time; ++step){
      if((step % N_time_resolve) == 0){
        domain_time["time"][counter] = pt;
        counter += 1;
      }
      pt += dt;
    }

    //set normalmode plot
    std::unordered_map< std::string, std::vector<double> > domain_normalmode;
    domain_normalmode["wave_vector"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_normalmode["wave_vector"][i] = i;
    std::vector<std::string> normalmode_value_list{"Normalmode_Energy","error_Normalmode_Energy",
                                                    "Normalmode_Dist" ,"error_Normalmode_Dist"};

    std::string result_dat("result.dat");
    std::string normalmode_dat("result_normalmode.dat");

    result_dat = result_directory + result_dat;
    normalmode_dat = result_directory + normalmode_dat;

    dataput << "result output to : " 
            << result_dat << " and "
            << normalmode_dat
            << std::endl;

    dataput.output_result(result_dat,"t",t/N_time_measure,N_time_measure,results);
    dataput.output_result_3d(normalmode_dat,domain_time,domain_normalmode,normalmode_value_list,results_2d);

    dataput.time_tag() << " end time " << std::endl; 
  }

  mpi_error = MPI_Finalize();
  return 0;
} //END clXYmodel

