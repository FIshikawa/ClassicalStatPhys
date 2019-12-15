#include <omp.h>
#include <mpi.h>
#include <cstdlib>
#include <memory>
#include <boost/lexical_cast.hpp>
#include <tools/data_recorder.hpp>
#include <tools/auto_correlation_function.hpp>
#include <tools/histogram.hpp>
#include <tools/accumulator.hpp>
#include <tools/time_measure.hpp>
#include <physics/toda_lax_form.hpp>
#include <physics/toda_conservation_fields.hpp>
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

  int Ns = 10; // lenght of chain
  int N_loop = 1; //loop number
  int Ns_observe = Ns; //observed particle 
  int N_time = 1000; //numer of time step
  int N_normalmode = 5; //numer of time step
  int N_time_measure = 10; //measuring of per time or linear
  int k_initial = 0; //start wave vector filled by initialization 
  int n_bin = 1000;
  int order = 0;
  double E_initial = 1.0; // initial energy
  double t = 10; // Whole time
  double t_relax = 10; // Whole time
  double t_sampling = 10; // Whole time
  std::string plot_scale = "linear";
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
  if (argc > input_counter) t_relax =          boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) t_sampling =       boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) E_initial =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) k_initial =        boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_normalmode =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if(Ns < N_normalmode + k_initial){
    std::cerr << "k_initial + N_noramalmode should be lower than Ns" << std::endl;
    std::exit(1);
  } 
  if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time_measure =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) plot_scale =       boost::lexical_cast<std::string>(argv[input_counter]);++input_counter;

  // set private parameters
  Parameters parameters;
  ParameterSet(argc, argv, parameters, input_counter);

  int num_particles = Lattice::set_num_particles(Ns); //nummer of particles

  int N_each = N_loop / N_mpi_parallel; // number of each parallel
  double dt = t / N_time; //time step length
  condition_dat = result_directory + condition_dat;
  tools::TimeMeasure time_measure(t, dt, N_time_measure, plot_scale);
  int N_total_data = time_measure.number_of_total_data();
  int N_time_relax = static_cast<int>(t_relax / dt);
  int N_time_sampling = static_cast<int>(t_sampling / dt);

  tools::DataRecorder dataput(condition_dat); 
  if(process_id == 0){
    dataput.time_tag() << "start time " << std::endl;

    // declare common settings
    dataput << "Declare Configulation " << std::endl 
            << "Explore Non-thermal Fixed Point : Time Develop Sampling " << std::endl
            << "Frist Hamiltonian for relaxation and sampling : " << FirstHamiltonian::name() << std::endl
            << "Secibd Hamiltonian for qunech : " << SecondHamiltonian::name() << std::endl
            << Ensembler::name() << std::endl
            << Integrator::name() << std::endl
            << Lattice::name() << std::endl;

    // declare private settings 
    ParameterDeclare(parameters,dataput);
    dataput << "Number of patricles : num_particles = " << num_particles << std::endl
            << "Length of whole system : Ns = " << Ns << std::endl
            << "Length of observe system : Ns_observe = " << Ns_observe << std::endl
            << "Number of time for measure: N_time_measure = " << N_time_measure << std::endl
            << "Number of result data (time step): N_total_data = " << N_total_data << std::endl
            << "Order of result data : order = " << time_measure.order() << std::endl
            << "Developing time : t = " << t << std::endl
            << "Number of time steps : N_time = " << N_time << std::endl
            << "Relaxation time : t_relax = " << t_relax << std::endl
            << "Number of relaxation time steps : N_time_relax = " << N_time_relax << std::endl
            << "Sampling interval time : t_sampling = " << t_sampling << std::endl
            << "Number of sampling interval time steps : N_time_sampling = " << N_time_sampling << std::endl
            << "Time interbal : dt = " << dt << std::endl
            << "Energy initial : E_initial =" << E_initial << std::endl
            << "Number of non-zero energy normal modes : N_normalmode =" << N_normalmode << std::endl
            << "Start wave vector filled by initialization : k_initial =" << k_initial << std::endl
            << "Number of MPI parallelization : N_mpi_parallel =" << N_mpi_parallel << std::endl
            << "Number of loop : N_loop =" << N_loop << std::endl
            << "Each thread loop : N_each = " << N_each << std::endl 
            << "Plot scale : plot_scale = " << time_measure.plot_scale()<< std::endl
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
                                                  "TodaIntegralTwo_total",
                                                  "TodaIntegralThree_total",
                                                  "TodaIntegralFour_total",
                                                  "TodaIntegralFive_total",
                                                  "VelocityACF",
                                                  "DisplacementACF" ,
                                                  "SpectralEntropy",
                                                  "EffectiveFractionOfModes"
                                                  };
  for(const auto& key : physical_value_list){
    results[key].resize(N_total_data);
    results["error_" + key].resize(N_total_data);
    for(int i = 0 ; i < N_total_data ; ++i){
      results[key][i] = 0.0;
      results["error_"+key][i] = 0.0;
    }
  }

  //set correlatin function 
  std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d{};
  std::vector<std::string> correlation_list{"Velocity_correlation","Displacement_correlation","EnergyField_correlation",
                                            "TodaIntegralThree_correlation","TodaIntegralFour_correlation","TodaIntegralFive_correlation"}; 
  std::vector<std::string> spatial_list{"Velocity","Displacement","EnergyField","TodaIntegralThree","TodaIntegralFour","TodaIntegralFive"}; 
  for(const auto& key : correlation_list){
    results_2d[key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
    results_2d["error_" + key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
  }
   for(const auto& key : spatial_list){
    results_2d[key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
    results_2d["error_" + key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
  }
  results_2d["Normalmode_Energy"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
  results_2d["error_Normalmode_Energy"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
  results_2d["TodaConservations"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
  results_2d["error_TodaConservations"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
  results_2d["TodaEigenvalues"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
  results_2d["error_TodaEigenvalues"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));

  //set histogram maps 
  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_data; 
 	hist_data["Velocity"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(N_loop)),
 	hist_data["Displacement"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(N_loop));
  //end set results and measure vectors

  if(process_id == 0) dataput.time_tag() << "start sampling process: time development " << std::endl;

  //lattice set
  Lattice lattice(Ns);
  int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj)); //interaction table 
  lattice.create_table(pair_table);
  std::vector<double> z(2*num_particles);
  
  //randam generetor set 
  std::size_t seed = 1234;

  //ensember set 
  using RandomGenerator = std::mt19937_64;
  std::minstd_rand seed_gen(seed);
  RandomGenerator mt(seed_gen());
  Ensembler ensembler(num_particles,k_initial,N_normalmode,E_initial);

  //hamiltonian set 
  FirstHamiltonian first_hamiltonian = FirstHamiltonianSet(parameters, num_particles, pair_table, N_adj);
  SecondHamiltonian second_hamiltonian = SecondHamiltonianSet(parameters, num_particles, pair_table, N_adj);
  //integratro set 
  Integrator integrator(2*num_particles);

  //set correlation
  std::unordered_map<std::string, correlation::AutoCorrelationFunction> autocorrelation;
  autocorrelation["VelocityACF"].initialize(N_total_data,N_loop);
  autocorrelation["DisplacementACF"].initialize(N_total_data,N_loop);

  //measure physical values set
  std::unordered_map<std::string, std::vector<tools::Accumulator> > physical_values;
  for(const auto& key : physical_value_list){
    physical_values[key].resize(N_total_data);
    for(int i = 0; i < N_total_data; ++i) physical_values[key][i].reset(); 
  }

  std::unordered_map<std::string, std::vector<std::vector<tools::Accumulator> > > correlation_function, spatial_function;
  for(auto key : spatial_list){
    correlation_function[key].resize(N_total_data);
    spatial_function[key].resize(N_total_data);
    for(int step = 0 ; step < N_total_data ; ++step){
      correlation_function[key][step] = std::vector<tools::Accumulator>(Ns_observe);
      spatial_function[key][step]     = std::vector<tools::Accumulator>(Ns_observe);
      for(int i = 0; i < Ns_observe; ++i){
        correlation_function[key][step][i].reset();
        spatial_function[key][step][i].reset();
      }
    }
  }

  //set normalmode energy
  std::vector<std::vector<tools::Accumulator> > normalmode_energy(N_total_data,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_total_data ; ++step){
    normalmode_energy[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) normalmode_energy[step][i].reset();
  }

  //set toda conservations
  std::vector<std::vector<tools::Accumulator> > toda_conservations(N_total_data,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_total_data ; ++step){
    toda_conservations[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) toda_conservations[step][i].reset();
  }

  //set toda eigenvalues
  std::vector<std::vector<tools::Accumulator> > toda_eigenvalues(N_total_data,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_total_data ; ++step){
    toda_eigenvalues[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) toda_eigenvalues[step][i].reset();
  }

  integrable::TodaLaxForm toda_lax_form(num_particles,0.25,2.0,pair_table,N_adj);
  integrable::TodaConservationFields toda_conservation_fields(num_particles,0.25,2.0,pair_table,N_adj);

  int observe = (num_particles-1) / 2;
  double *x = &z[0];
  double *v = &z[num_particles];
  double pt = 0.0;
  double v_total,x_total,sum_spectral,average_spectral,spectral_entropy;
  std::vector<double> normalmode_energy_temp(num_particles), 
                      toda_conservations_temp(num_particles),
                      toda_eigenvalues_temp(num_particles),
                      toda_integral(5);

  // set initial state
  ensembler.set_initial_state(z,mt);
  N_time_relax += N_time_sampling * (process_id * N_each);
  if(process_id == 0) dataput.time_tag() << "start time development " << std::endl;
  for(int step = 0; step < N_time_relax; ++step) integrator.step(pt, dt, z, first_hamiltonian);
  if(process_id == 0) dataput.time_tag() << "end time development " << std::endl;

  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){

    // relaxation dynamics 
    for(int step = 0; step < N_time_sampling; ++step) integrator.step(pt, dt, z, first_hamiltonian);

    //time develop
    int measure_step = 0;
    time_measure.reset();
    for(int step = 0; step < N_time; ++step){  
      //set total values 
      if(time_measure.check(step)){
        v_total = 0.0;
        x_total = 0.0;
        for(int i = 0; i < num_particles; ++i){
          v_total += v[i] / num_particles;
          x_total += x[i] / num_particles;
        } 
        for(int kind = 1; kind < 6; ++kind){
          toda_integral[kind-1] = 0;
          for(int i = 0; i < num_particles; ++i) toda_integral[kind-1] += toda_conservation_fields(z,i,kind);
        }
        //set normalmode energy
        NormalModeEnergy(z,normalmode_energy_temp);

        sum_spectral = 0;
        spectral_entropy= 0;
        average_spectral= 0;
        for(int i = 0; i < num_particles; ++i) sum_spectral += normalmode_energy_temp[i];
        for(int i = 0; i < num_particles; ++i){
          average_spectral = normalmode_energy_temp[i] / sum_spectral;
          if(average_spectral > 0) spectral_entropy -=  average_spectral * std::log(average_spectral);
        }

        //set toda_conservations
        toda_lax_form.conservations_with_eigenvalues(z, toda_conservations_temp, toda_eigenvalues_temp);

        //input data 
        physical_values["Velocity"][measure_step] << v[observe];
        physical_values["Velocity_total"][measure_step] << v_total;
        autocorrelation["VelocityACF"] << v[observe] - v_total;
        physical_values["Displacement"][measure_step] << x[observe];
        physical_values["Displacement_total"][measure_step] << x_total;
        autocorrelation["DisplacementACF"] << x[observe];
        physical_values["Energy_total"][measure_step] <<  second_hamiltonian.energy(pt,z) / num_particles;
        physical_values["Kinetic_Energy"][measure_step] << second_hamiltonian.kinetic_energy(pt,z) / num_particles;
        physical_values["Potential_Energy"][measure_step] << second_hamiltonian.potential_energy(pt,z) / num_particles;
        physical_values["SpectralEntropy"][measure_step] << spectral_entropy;
        physical_values["EffectiveFractionOfModes"][measure_step] << std::exp(spectral_entropy) / num_particles;
        physical_values["TodaIntegralTwo_total"][measure_step] << toda_integral[1];
        physical_values["TodaIntegralThree_total"][measure_step] << toda_integral[2];
        physical_values["TodaIntegralFour_total"][measure_step] << toda_integral[3];
        physical_values["TodaIntegralFive_total"][measure_step] << toda_integral[4];

      //input hist data 
        hist_data["Velocity"][measure_step][c_loop] = v[observe];
        hist_data["Displacement"][measure_step][c_loop] = x_total;
        //correlation calc
        SpatialMeasurement(correlation_function, spatial_function, lattice, first_hamiltonian, 
                           toda_conservation_fields, Ns_observe, observe, measure_step, z);
        //input normalmode energy
        for(int i = 0; i < num_particles; ++i) {
          normalmode_energy[measure_step][i] << normalmode_energy_temp[i];
          toda_conservations[measure_step][i] << toda_conservations_temp[i];
          toda_eigenvalues[measure_step][i] << toda_eigenvalues_temp[i];
        }
        ++measure_step;
      }

      //time develp
      pt += dt;
      integrator.step(pt, dt, z, second_hamiltonian);
    }//end time step
    if(measure_step != N_total_data){
      std::cerr << "final measure_step : " << measure_step << ", N_total_data : " << N_total_data <<", not same " << std::endl;
      std::exit(112);
    }

  }//end loop

  // correct data each process
  dataput.time_tag() << "[" << process_id <<"th process] start result correct by self" << std::endl;
  for(const auto& key  : physical_value_list){
    for(int i = 0 ; i < N_total_data; ++i){
      results[key][i] = physical_values[key][i].sum1();
      results["error_" + key][i] = physical_values[key][i].sum2();
    }
  }
  for(const auto& element : autocorrelation){
    std::string key = element.first;
    std::vector<double> acf_result = autocorrelation[key].result();
    for(int i = 0 ; i < N_total_data; ++i){
      results[key][i] = acf_result[i] ;
    }      
  }
  for(int step = 0 ; step < N_total_data; ++step){
    for(int i = 0 ; i < num_particles; ++i){
      results_2d["Normalmode_Energy"][step][i] = normalmode_energy[step][i].sum1(); 
      results_2d["error_Normalmode_Energy"][step][i] = normalmode_energy[step][i].sum2();
      results_2d["TodaConservations"][step][i] = toda_conservations[step][i].sum1(); 
      results_2d["error_TodaConservations"][step][i] = toda_conservations[step][i].sum2();
      results_2d["TodaEigenvalues"][step][i] = toda_eigenvalues[step][i].sum1(); 
      results_2d["error_TodaEigenvalues"][step][i] = toda_eigenvalues[step][i].sum2();
    }
    for(int i = 0 ; i < Ns_observe; ++i){
      for(std::string key : {"Displacement","Velocity","EnergyField"}){
        results_2d[key + "_correlation"][step][i] = correlation_function[key][step][i].sum1();
        results_2d["error_" + key + "_correlation"][step][i] = correlation_function[key][step][i].sum2();
        results_2d[key][step][i] = spatial_function[key][step][i].sum1();
        results_2d["error_" + key][step][i] = spatial_function[key][step][i].sum2();
      }
    }
  }
  dataput.time_tag() << "[" << process_id <<"th process] finish correction : ready for transfer" << std::endl;


  // correct data to process id 0 
  if(process_id == 0) dataput.time_tag() << "start total result correct" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start normal results" << std::endl;
  std::vector<double> reduced_vector(N_total_data);
  for(const auto& key  : physical_value_list){
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_total_data; ++i) results[key][i] = reduced_vector[i] / N_loop;
    mpi_error = MPI_Reduce(&results["error_" + key][0], &reduced_vector[0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_total_data; ++i) 
      results["error_" + key][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector[i] / (N_loop -1) - results[key][i] * results[key][i]* (N_loop)/(N_loop-1)) / N_loop) : 0;
  }
  if(process_id == 0) dataput.time_tag() << "end normal results" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start autocorrelation results" << std::endl;
  for(const auto& element : autocorrelation){
    std::string key = element.first;
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_total_data; ++i) results[key][i] = reduced_vector[i] / N_loop ;
  }
  if(process_id == 0) dataput.time_tag() << "end autocorrelation results" << std::endl;

  std::vector<double> reduced_vector_1(num_particles);
  std::vector<double> reduced_vector_2(Ns_observe);
  if(process_id == 0) dataput.time_tag() << "start 2d results" << std::endl;
  for(int step = 0 ; step < N_total_data; ++step){
    // NormalMode
    for(std::string value_name: {"Normalmode_Energy","TodaConservations","TodaEigenvalues"}){
      mpi_error = MPI_Reduce(&results_2d[value_name][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) results_2d[value_name][step][i] = reduced_vector_1[i] / N_loop; 
      mpi_error = MPI_Reduce(&results_2d["error_" + value_name][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) 
        results_2d["error_" + value_name][step][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector_1[i] / (N_loop -1)
                                             - results_2d[value_name][step][i] * results_2d[value_name][step][i] * (N_loop)/(N_loop-1))/N_loop) : 0; 
    }
 

    // correlation values
    for(const auto& key : correlation_list){
      mpi_error = MPI_Reduce(&results_2d[key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) results_2d[key][step][i] = reduced_vector_2[i] / N_loop;
      mpi_error = MPI_Reduce(&results_2d["error_" + key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) 
        results_2d["error_" + key][step][i] = (N_loop > 1) ? std::sqrt((reduced_vector_2[i] / (N_loop-1)
                                    - results_2d[key][step][i] * results_2d[key][step][i] * (N_loop)/(N_loop-1)) / N_loop) : 0;
    }

    // spatial values
    for(const auto& key : spatial_list){
      mpi_error = MPI_Reduce(&results_2d[key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) results_2d[key][step][i] = reduced_vector_2[i] / N_loop;
      mpi_error = MPI_Reduce(&results_2d["error_" + key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) 
        results_2d["error_" + key][step][i] = (N_loop > 1) ? std::sqrt((reduced_vector_2[i] / (N_loop-1)
                                    - results_2d[key][step][i] * results_2d[key][step][i] * (N_loop)/(N_loop-1)) / N_loop) : 0;
    }
  }
  if(process_id == 0) dataput.time_tag() << "end 2d results" << std::endl;

  std::vector<double> hist_vector(N_loop);
  for(auto key : {"Velocity","Displacement"}){
    for(int step = 0 ; step < N_total_data; ++step){
      mpi_error = MPI_Gather(&hist_data[key][step][0], N_each, MPI_DOUBLE, &hist_vector[0], N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for(int i = 0 ; i < N_loop; ++i) hist_data[key][step][i] = hist_vector[i];
    }
  }

  if(process_id == 0) dataput.time_tag() << "end total result correct" << std::endl;


  // data output by process_id = 0
  if(process_id == 0){
    dataput.time_tag() << "end sampling process" << std::endl;

    dataput.time_tag() << "correlation values calc" << std::endl;
    double correlation_error;
    for(const auto& key : spatial_list){
      for(int step = 0 ; step < N_total_data; ++step){
        for(int i = 0; i < Ns_observe; ++i){
          results_2d[key + "_correlation"][step][i] -= results_2d[key][step][i] * results_2d[key][step][0];
          correlation_error = results_2d["error_" + key][step][i] * results_2d[key][step][0]
                              * results_2d["error_" + key][step][i] * results_2d[key][step][0]
                              + results_2d["error_" + key][step][0] * results_2d[key][step][i]
                              * results_2d["error_" + key][step][0] * results_2d[key][step][i]
                              + results_2d["error_" + key + "_correlation"][step][i]
                              * results_2d["error_" + key + "_correlation"][step][i];
          results_2d["error_" + key + "_correlation"][step][i] = std::sqrt(correlation_error);
        }
      }
    }
    dataput.time_tag() << "end correlation values calc" << std::endl;

    dataput.time_tag() << "start stat calc " << std::endl;
    //set Accumulator
    tools::Accumulator accumulator;

    //define histogram class
    std::unordered_map<std::string, tools::Histogram>  histogram;
    histogram["Displacement"].initialize(n_bin);
    histogram["Velocity"].initialize(n_bin);

    dataput.time_tag() << "start stat declaration " << std::endl;
    //set statistical values
    std::unordered_map<std::string, std::vector<double> > stat_values;
    std::vector<std::string> stat_value_list = {
                                                "mean",
                                                "variance",
                                                "skewness",
                                                "kurtosis",
                                                "kurtosis_exess",
                                                "average","error",
                                                "central_moment1",
                                                "central_moment2",
                                                "central_moment3",
                                                "central_moment4"
                                              };

    std::unordered_map<std::string, std::vector< std::vector<double> > > hist_values;
    std::vector<std::string> hist_value_list = {
                                                "hist_Velocity",
                                                "range_Velocity",
                                                "hist_Displacement",
                                                "range_Displacement"
                                                };

    for(auto key  : stat_value_list) stat_values[key].resize(N_total_data);
    for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(n_bin));

    dataput.time_tag() << "end stat declaration " << std::endl;
    //calc histogram and statisitcal values
    for(int i = 0; i < N_total_data; ++i){
      accumulator.reset();   
      for(int j = 0; j < N_loop; ++j){
        accumulator << hist_data["Velocity"][i][j];
      }
      stat_values["mean"][i] = accumulator.mean();
      stat_values["variance"][i] = accumulator.variance();
      stat_values["skewness"][i] = accumulator.skewness();
      stat_values["kurtosis"][i] = accumulator.kurtosis();
      stat_values["kurtosis_exess"][i] = accumulator.kurtosis_excess();
      stat_values["average"][i] = accumulator.average();
      stat_values["error"][i] = accumulator.error();
      stat_values["central_moment1"][i] = accumulator.central_moment1();
      stat_values["central_moment2"][i] = accumulator.central_moment2();
      stat_values["central_moment3"][i] = accumulator.central_moment3();
      stat_values["central_moment4"][i] = accumulator.central_moment4();
      for(auto element : histogram){
        std::string key = element.first; 
        histogram[key](hist_data[key][i]);
        histogram[key].output(hist_values["hist_" + key][i],
                                          hist_values["range_" + key][i]);
      }
    }
    dataput.time_tag() << "end stat calc " << std::endl;

    //set measure time
    std::unordered_map< std::string, std::vector<double> > domain_time;
    domain_time["time"] = std::vector<double>(N_total_data);
    time_measure.reset();
    for(int step = 0; step < N_total_data; ++step) domain_time["time"][step] = time_measure();

    //set normalmode plot
    std::unordered_map< std::string, std::vector<double> > domain_normalmode;
    domain_normalmode["wave_vector"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_normalmode["wave_vector"][i] = i;
    std::vector<std::string> normalmode_value_list{"Normalmode_Energy","error_Normalmode_Energy"};

    //set toda conservations plot
    std::unordered_map< std::string, std::vector<double> > domain_toda_conservations;
    domain_toda_conservations["power"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_toda_conservations["power"][i] = i;
    std::vector<std::string> toda_conservations_list{"TodaConservations","error_TodaConservations"};

    //set toda eigenvalues plot
    std::unordered_map< std::string, std::vector<double> > domain_toda_eigenvalues;
    domain_toda_eigenvalues["number"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_toda_eigenvalues["number"][i] = i;
    std::vector<std::string> toda_eigenvalues_list{"TodaEigenvalues","error_TodaEigenvalues"};

    //add error values of correlation and spatial to string list
    int correlation_length = correlation_list.size();
    for(int i = 0 ; i < correlation_length; ++i){
      correlation_list.push_back("error_" + correlation_list[i]);
    }
    int spatial_length = spatial_list.size();
    for(int i = 0 ; i < spatial_length; ++i){
      spatial_list.push_back("error_" + spatial_list[i]);
    }
    int physical_value_length = physical_value_list.size();
    for(int i =0 ; i < physical_value_length; ++i){
      physical_value_list.push_back("error_" + physical_value_list[i]);
    }

    //set sight number domain
    std::unordered_map< std::string, std::vector<double> > domain_sight;
    domain_sight["site"] = std::vector<double>(Ns_observe);
    for(int i = 0 ; i < Ns_observe; ++i) domain_sight["site"][i] = i;

    std::string result_dat("result.dat");
    std::string hist_dat("result_hist.dat");
    std::string normalmode_dat("result_normalmode.dat");
    std::string toda_conservations_dat("result_toda_conservations.dat");
    std::string toda_eigenvalues_dat("result_toda_eigenvalues.dat");
    std::string stat_dat("result_stat.dat");
    std::string spatial_dat("result_spatial.dat");
    std::string correlation_dat("result_correlation.dat");

    result_dat = result_directory + result_dat;
    hist_dat = result_directory + hist_dat;
    stat_dat = result_directory + stat_dat;
    normalmode_dat = result_directory + normalmode_dat;
    toda_conservations_dat = result_directory + toda_conservations_dat;
    toda_eigenvalues_dat = result_directory + toda_eigenvalues_dat;
    spatial_dat = result_directory + spatial_dat;
    correlation_dat = result_directory + correlation_dat;

    dataput << "result output to : " 
            << result_dat << " , "
            << correlation_dat << " , "
            << spatial_dat << " , "
            << normalmode_dat << " , "
            << toda_conservations_dat << " , "
            << toda_eigenvalues_dat << " , "
            << hist_dat << " and " 
            << stat_dat 
            << std::endl;

    dataput.output_result(result_dat,domain_time,results);
    dataput.output_result_hist(hist_dat,domain_time,hist_value_list,hist_values);
    dataput.output_result_3d(correlation_dat,domain_time,domain_sight,correlation_list,results_2d);
    dataput.output_result_3d(spatial_dat,domain_time,domain_sight,spatial_list,results_2d);
    dataput.output_result_3d(normalmode_dat,domain_time,domain_normalmode,normalmode_value_list,results_2d);
    dataput.output_result_3d(toda_conservations_dat,domain_time,domain_toda_conservations,toda_conservations_list,results_2d);
    dataput.output_result_3d(toda_eigenvalues_dat,domain_time,domain_toda_eigenvalues,toda_eigenvalues_list,results_2d);
    dataput.output_result(stat_dat,domain_time,stat_values);
    dataput.time_tag() << " end time " << std::endl; 
  }
  mpi_error = MPI_Finalize();
  return 0;
} //END clXYmodel

