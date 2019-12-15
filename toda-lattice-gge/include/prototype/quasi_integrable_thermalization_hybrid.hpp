#include <omp.h>
#include <mpi.h>
#include <cstdlib>
#include <memory>
#include <boost/lexical_cast.hpp>
#include <tools/data_recorder.hpp>
#include <tools/auto_correlation_function.hpp>
#include <tools/histogram.hpp>
#include <tools/accumulator.hpp>
#include <physics/normalmode_energy_fftw.hpp>
#include <integrator/yoshida_4th_parallel.hpp>

using Integrator = integrator::Yoshida4thParallel;
// Ensembler, Hamiltonian, Lattice are defined in each src header 

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

  //parameter define
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;
  int Ns = 10; // lenght of chain
  int N_loop = 1; //loop number
  double t = 10; // Whole time
  int N_time = 1000; //numer of time step
  int N_time_resolve = 1; //numer of time step
  int N_omp_parallel = 2; //number of parallel 
  double E_initial = 1.0; // initial energy
  int n_bin = 1000;
  std::string condition_dat = "condi.dat";
  std::string result_directory("./");

 //parameter set 
  int input_counter = 1;
  if (argc > input_counter) result_directory = boost::lexical_cast<std::string>(argv[input_counter]); ++input_counter;
  if (argc > input_counter) Ns =               boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) t =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) J =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) alpha =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  #if  defined(HAMILTONIAN_FPU_GENERALIZED_FIXED_END_PARALLEL_HPP) || defined(HAMILTONIAN_FPU_FIXED_END_PARALLEL_HPP)
  double beta = 1.0; //interaction constant;
  if (argc > input_counter) beta =             boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  #endif
  #if  defined(HAMILTONIAN_FPU_GENERALIZED_FIXED_END_PARALLEL_HPP)
  double gamma = 1.0; //interaction constant;
  if (argc > input_counter) gamma =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;  
  double delta = 1.0; //interaction constant;
  if (argc > input_counter) delta =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;  
  #endif
  if (argc > input_counter) E_initial =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  #if defined(ENSEMBLE_NORMALMODE_ENSEMBLE_FFTW_HPP) || defined(ENSEMBLE_NORMALMODE_ENSEMBLE_HPP) 
  int N_normalmode = 5; //numer of time step
  int k_initial = 0; //start wave vector filled by initialization 
  if (argc > input_counter) k_initial =        boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_normalmode =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if(Ns < N_normalmode + k_initial){
    std::cerr << "k_initial + N_noramalmode should be lower than Ns" << std::endl;
    std::exit(1);
  } 
  #endif
  if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_omp_parallel =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time_resolve =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;

  int num_particles = Lattice::set_num_particles(Ns); //nummer of particles
  int Ns_observe = num_particles / 10; //observed particle 
  int N_time_measure = N_time / N_time_resolve; //number of time span of histogram
  double dt = t / N_time; //time step length
  int N_each = N_loop / N_mpi_parallel; // number of each parallel
  condition_dat = result_directory + condition_dat;

  tools::DataRecorder dataput(condition_dat); 
  if(process_id == 0){
    dataput.time_tag() << "start time " << std::endl;
    dataput << "Declare Configulation " << std::endl 
            << "Explore Non-thermal Fixed Point : FPU model" << std::endl
            <<  Hamiltonian::name() << std::endl
            <<  Ensembler::name() << std::endl
            <<  Integrator::name() << std::endl
            <<  Lattice::name() << std::endl;

   //condition declare
    dataput << "Coupling constant : J = " << J << std::endl
            << "Coupling constant : alpha = " << alpha << std::endl
            #if  defined(HAMILTONIAN_FPU_GENERALIZED_FIXED_END_PARALLEL_HPP) || defined(HAMILTONIAN_FPU_FIXED_END_PARALLEL_HPP)
            << "Coupling constant : beta = " << beta << std::endl
            #endif
            #if  defined(HAMILTONIAN_FPU_GENERALIZED_FIXED_END_PARALLEL_HPP)
            << "Coupling constant : gamma = " << gamma << std::endl
            << "Coupling constant : delta = " << delta << std::endl
            #endif 
            << "Number of patricles : num_particles = " << num_particles << std::endl
            << "Length of whole system : Ns = " << Ns << std::endl
            << "Length of observe system : Ns_observe = " << Ns_observe << std::endl
            << "Number of time resolve for measure: N_time_resolve = " << N_time_resolve << std::endl
            << "Number of time for measure: N_time_resolve = " << N_time_measure << std::endl
            << "Developing ime : t = " << t << std::endl
            << "Number of time steps : N_time = " << N_time << std::endl
            << "Time interbal : dt = " << dt << std::endl
            << "Energy initial : E_initial =" << E_initial << std::endl
            #if defined(ENSEMBLE_NORMALMODE_ENSEMBLE_FFTW_HPP) || defined(ENSEMBLE_NORMALMODE_ENSEMBLE_FFTW_HPP) 
            << "Number of non-zero energy normal modes : N_normalmode =" << N_normalmode << std::endl
            << "Start wave vector filled by initialization : k_initial =" << k_initial << std::endl
            #endif
            << "Number of omp parallelization : N_omp_parallel =" << N_omp_parallel << std::endl
            << "Number of MPI parallelization : N_mpi_parallel =" << N_mpi_parallel << std::endl
            << "Number of loop : N_loop =" << N_loop << std::endl
            << "Each thread loop : N_each = " << N_each << std::endl 
            << "Number of time resolve for hist : N_time_resolve = " << N_time_resolve << std::endl
            << "Result directory : " << result_directory << std::endl;
   //end parameter define 
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
                                                  "VelocityACF",
                                                  "DisplacementACF" ,
                                                  "SpectralEntropy",
                                                  "EffectiveFractionOfModes"
                                                  };
  for(const auto& key : physical_value_list){
    results[key].resize(N_time);
    results["error_" + key].resize(N_time);
    for(int i = 0 ; i < N_time ; ++i){
      results[key][i] = 0.0;
      results["error_"+key][i] = 0.0;
    }
  }

  //set correlatin function 
  std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d{};  
  std::vector<std::string> correlation_list{"Velocity_correlation","Displacement_correlation","NumberOperator_correlation"}; 
  for(const auto& key : correlation_list){
    results_2d[key] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(Ns_observe,0.0));
    results_2d["error_" + key] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(Ns_observe,0.0));
  }
  results_2d["Normalmode_Energy"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));
  results_2d["error_Normalmode_Energy"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));

  //set histogram maps 
  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_data; 
 	hist_data["Velocity"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(N_loop)),
 	hist_data["Displacement"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(N_loop));
  //end set results and measure vectors

  if(process_id == 0) dataput.time_tag() << "start sampling process : time development " << std::endl;

  // set number of thread parallel 
  omp_set_num_threads(N_omp_parallel);

  //lattice set
  Lattice lattice(Ns);
  int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj)); //interaction table 
  lattice.create_table(pair_table);
  std::vector<double> z(2*num_particles);
  #pragma omp parallel for
  for(int i = 0 ; i < num_particles; ++i){
    z[i] = 0.0;
    z[i+num_particles] = 0.0;
  }
  
  //randam generetor set 
  std::size_t seed = 1234 * (1 + process_id);

  //ensember set 
  using RandomGenerator = std::mt19937_64;
  #if defined(ENSEMBLE_NORMALMODE_ENSEMBLE_FFTW_HPP) || defined(ENSEMBLE_NORMALMODE_ENSEMBLE_HPP) 
  std::minstd_rand seed_gen(seed);
  RandomGenerator mt(seed_gen());
  Ensembler ensembler(num_particles,k_initial,N_normalmode,E_initial);
  #elif defined(ENSEMBLE_WATER_BAG_PARALLEL_HPP)
  //set random value for  parallelization
  std::minstd_rand seed_gen(seed);
  std::vector<std::shared_ptr<RandomGenerator> > mt(N_omp_parallel);
  #pragma omp parallel 
  {
  int tid = omp_get_thread_num();
  int num_threads = omp_get_num_threads();
  for (int p = 0; p < num_threads; ++p) {
    if (p == tid)
      mt[p].reset(new RandomGenerator(seed_gen()));
      #pragma omp barrier
  }
  }
  //set random value for  parallelization end 
  Ensembler ensembler(E_initial*2.0, num_particles);
  #endif 

  //hamiltonian set 
  #if defined(HAMILTONIAN_FPU_FIXED_END_PARALLEL_HPP)
  Hamiltonian hamiltonian(num_particles,J,alpha,beta,pair_table,N_adj);
  #elif defined(HAMILTONIAN_FPU_GENERALIZED_FIXED_END_PARALLEL_HPP)
  Hamiltonian hamiltonian(num_particles,J,alpha,beta,gamma,delta,pair_table,N_adj);
  #else
  Hamiltonian hamiltonian(num_particles,J,alpha,pair_table,N_adj);
  #endif

  //integratro set 
  Integrator integrator(2*num_particles);

  //set correlation
  std::unordered_map<std::string, correlation::AutoCorrelationFunction> autocorrelation;
  autocorrelation["VelocityACF"].initialize(N_time,N_loop);
  autocorrelation["DisplacementACF"].initialize(N_time,N_loop);

  //measure physical values set
  std::unordered_map<std::string, std::vector<tools::Accumulator> > physical_values;
  for(const auto& key : physical_value_list){
    physical_values[key].resize(N_time);
  }

  std::unordered_map<std::string, std::vector<std::vector<tools::Accumulator> > > correlation_function;
  correlation_function["Displacement"].resize(N_time_measure);
  correlation_function["Velocity"].resize(N_time_measure);
  correlation_function["NumberOperator"].resize(N_time_measure);
  for(auto key : {"Displacement","Velocity","NumberOperator"}){
    for(int step = 0 ; step < N_time_measure ; ++step) correlation_function[key][step] = std::vector<tools::Accumulator>(Ns_observe);
  }

  //set normalmode energy
  std::vector<std::vector<tools::Accumulator> > normalmode_energy(N_time_measure,std::vector<tools::Accumulator>(num_particles));

  int observe = (num_particles-1) / 2;
  double v_total = 0.0;
  double x_total = 0.0;
  double *x = &z[0];
  double *v = &z[num_particles];
  int target = 0;
  int measure_step = 0;
  double pt = 0.0;
  double number_operator_observe = 0;
  double number_operator_target  = 0;
  double sum_spectral,average_spectral,spectral_entropy;
  std::vector<double> normalmode_energy_temp(num_particles);

  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){

    // set initial state
    ensembler.set_initial_state(z,mt);

    //time develop
    for(int step = 0; step < N_time; ++step){  

      //set total values 
      v_total = 0.0;
      x_total = 0.0;
      for(int i = 0; i < num_particles; ++i){
        v_total += v[i] / num_particles;
        x_total += x[i] / num_particles;
      } 
      //set normalmode energy
      NormalModeEnergyFFTW(z,normalmode_energy_temp);

      sum_spectral = 0;
      spectral_entropy= 0;
      average_spectral= 0;
      for(int i = 0; i < num_particles; ++i) sum_spectral += normalmode_energy_temp[i];
      for(int i = 0; i < num_particles; ++i){
        average_spectral = normalmode_energy_temp[i] / sum_spectral;
        spectral_entropy -=  average_spectral * std::log(average_spectral);
      }

      //input data 
      physical_values["Velocity"][step] << v[observe];
      physical_values["Velocity_total"][step] << v_total;
      autocorrelation["VelocityACF"] << v[observe] - v_total;
      physical_values["Displacement"][step] << x[observe];
      physical_values["Displacement_total"][step] << x_total;
      autocorrelation["DisplacementACF"] << x[observe];
      physical_values["Energy_total"][step] <<  hamiltonian.energy(pt,z) / num_particles;
      physical_values["Kinetic_Energy"][step] << hamiltonian.kinetic_energy(pt,z) / num_particles;
      physical_values["Potential_Energy"][step] << hamiltonian.potential_energy(pt,z) / num_particles;
      physical_values["SpectralEntropy"][step] << spectral_entropy;
      physical_values["EffectiveFractionOfModes"][step] << std::exp(spectral_entropy) / num_particles;

      //input hist data 
      if((step % N_time_resolve) == 0){
        measure_step = step / N_time_resolve;
        hist_data["Velocity"][measure_step][c_loop] = v[observe];
        hist_data["Displacement"][measure_step][c_loop] = x_total;
        //correlation calc
        for(int i = 0; i < Ns_observe; ++i){
          target = lattice.target_number(i,observe);
          correlation_function["Displacement"][measure_step][i] << x[observe] * x[target];
          correlation_function["Velocity"][measure_step][i] << v[observe] * v[target];
          number_operator_observe = (v[observe] * v[observe] + J * x[observe] * x[observe]) / (2*std::sqrt(J));  
          number_operator_target = (v[target] * v[target] + J * x[target] * x[target]) / (2*std::sqrt(J));  
          correlation_function["NumberOperator"][measure_step][i] << number_operator_observe * number_operator_target;
        }
        //input normalmode energy
        for(int i = 0; i < num_particles; ++i) normalmode_energy[measure_step][i] << normalmode_energy_temp[i];
      }

      //time develp
      pt += dt;
      integrator.step(pt, dt, z, hamiltonian);
    }//end time step
  }//end loop

  // correct data each process
  dataput.time_tag() << "[" << process_id <<"th process] start result correct by self" << std::endl;
  for(const auto& key  : physical_value_list){
    for(int i = 0 ; i < N_time; ++i){
      results[key][i] += physical_values[key][i].mean() / N_mpi_parallel;
      results["error_" + key][i] += physical_values[key][i].error() / N_mpi_parallel;
    }
  }
  for(const auto& element : autocorrelation){
    std::string key = element.first;
    std::vector<double> acf_result = autocorrelation[key].result();
    for(int i = 0 ; i < N_time; ++i){
      results[key][i] += acf_result[i] / N_mpi_parallel;
    }      
  }
  for(int step = 0 ; step < N_time_measure; ++step){
    for(int i = 0 ; i < num_particles; ++i){
      results_2d["Normalmode_Energy"][step][i] += normalmode_energy[step][i].mean() / N_mpi_parallel; 
      results_2d["error_Normalmode_Energy"][step][i] += normalmode_energy[step][i].error() / N_mpi_parallel; 
    }
    for(int i = 0 ; i < Ns_observe; ++i){
      for(std::string key : {"Displacement","Velocity","NumberOperator"}){
        results_2d[key + "_correlation"][step][i] += correlation_function[key][step][i].mean() / N_mpi_parallel;
        results_2d["error_" + key + "_correlation"][step][i] += correlation_function[key][step][i].error() / N_mpi_parallel;
      }
    }
  }
  dataput.time_tag() << "[" << process_id <<"th process] finish correction : ready for transfer" << std::endl;


  // correct data to process id 0 
  if(process_id == 0) dataput.time_tag() << "start total result correct" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start normal results" << std::endl;
  std::vector<double> reduced_vector(N_time);
  for(const auto& key  : physical_value_list){
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_time, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_time; ++i) results[key][i] = reduced_vector[i];
    mpi_error = MPI_Reduce(&results["error_" + key][0], &reduced_vector[0], N_time, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_time; ++i) results["error_" + key][i] = reduced_vector[i];
  }
  if(process_id == 0) dataput.time_tag() << "end normal results" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start autocorrelation results" << std::endl;
  for(const auto& element : autocorrelation){
    std::string key = element.first;
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_time, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_time; ++i) results[key][i] = reduced_vector[i] ;
  }
  if(process_id == 0) dataput.time_tag() << "end autocorrelation results" << std::endl;

  std::vector<double> reduced_vector_1(num_particles);
  std::vector<double> reduced_vector_2(Ns_observe);
  if(process_id == 0) dataput.time_tag() << "start 2d results" << std::endl;
  for(int step = 0 ; step < N_time_measure; ++step){
    mpi_error = MPI_Reduce(&results_2d["Normalmode_Energy"][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) results_2d["Normalmode_Energy"][step][i] = reduced_vector_1[i]; 
    mpi_error = MPI_Reduce(&results_2d["error_Normalmode_Energy"][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) results_2d["error_Normalmode_Energy"][step][i] = reduced_vector_1[i]; 
    for(const auto& key : correlation_list){
      mpi_error = MPI_Reduce(&results_2d[key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) results_2d[key][step][i] = reduced_vector_2[i];
      mpi_error = MPI_Reduce(&results_2d["error_" + key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) results_2d["error_" + key][step][i] = reduced_vector_2[i];
    }
  }
  if(process_id == 0) dataput.time_tag() << "end 2d results" << std::endl;

  std::vector<double> hist_vector(N_loop);
  for(auto key : {"Velocity","Displacement"}){
    for(int step = 0 ; step < N_time_measure; ++step){
      mpi_error = MPI_Gather(&hist_data[key][step][0], N_each, MPI_DOUBLE, &hist_vector[0], N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for(int i = 0 ; i < N_loop; ++i) hist_data[key][step][i] = hist_vector[i];
    }
  }

  if(process_id == 0) dataput.time_tag() << "end total result correct" << std::endl;


  // data output by process_id = 0
  if(process_id == 0){
    dataput.time_tag() << "end sampling process" << std::endl;

    dataput.time_tag() << "start stat calc " << std::endl;
    //set Accumulator
    tools::Accumulator accumulator;

    //define histogram class
    std::unordered_map<std::string, tools::Histogram>  histogram;
    histogram["Displacement"].initialize(n_bin);
    histogram["Velocity"].initialize(n_bin);

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

    for(auto key  : stat_value_list) stat_values[key].resize(N_time_measure);
    for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(n_bin));

    //calc histogram and statisitcal values
    for(int i = 0; i < N_time_measure; ++i){
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
    std::vector<std::string> normalmode_value_list{"Normalmode_Energy","error_Normalmode_Energy"};

    //add error values of correlation to string list
    int correlation_length = correlation_list.size();
    for(int i = 0 ; i < correlation_length; ++i){
      correlation_list.push_back("error_" + correlation_list[i]);
    }

    //set sight number domain
    std::unordered_map< std::string, std::vector<double> > domain_sight;
    domain_sight["site"] = std::vector<double>(Ns_observe);
    for(int i = 0 ; i < Ns_observe; ++i) domain_sight["site"][i] = i;

    std::string result_dat("result.dat");
    std::string correlation_dat("result_correlation.dat");
    std::string hist_dat("result_hist.dat");
    std::string normalmode_dat("result_normalmode.dat");
    std::string stat_dat("result_stat.dat");

    result_dat = result_directory + result_dat;
    correlation_dat = result_directory + correlation_dat;
    hist_dat = result_directory + hist_dat;
    stat_dat = result_directory + stat_dat;
    normalmode_dat = result_directory + normalmode_dat;

    dataput << "result output to : " 
            << result_dat << " , "
            << correlation_dat << " , "
            << normalmode_dat << " , "
            << hist_dat << " and " 
            << stat_dat 
            << std::endl;

    dataput.output_result(result_dat,"pt",dt,N_time,results);
    dataput.output_result_hist(hist_dat,domain_time,hist_value_list,hist_values);
    dataput.output_result_3d(correlation_dat,domain_time,domain_sight,correlation_list,results_2d);
    dataput.output_result_3d(normalmode_dat,domain_time,domain_normalmode,normalmode_value_list,results_2d);
    dataput.output_result(stat_dat,domain_time,stat_values);

    dataput.time_tag() << " end time " << std::endl; 
  }

  mpi_error = MPI_Finalize();
  return 0;
} //END clXYmodel

