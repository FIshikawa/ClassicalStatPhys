#include <omp.h>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <tools/data_recorder.hpp>
#include <tools/auto_correlation_function.hpp>
#include <tools/histogram.hpp>
#include <tools/accumulator.hpp>
#include <integrator/yoshida_4th.hpp>

using Integrator = integrator::Yoshida4th;
// Ensembler, Hamiltonian, Lattice are defined in each src header 

int main(int argc, char **argv){
  //parameter define
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;
  double beta = 1.0; //interaction constant;
  double T = 0.0; //temperture 
  int Ns = 10; // lenght of chain
  int N_loop = 2; //loop number
  double t = 10; // Whole time
  int N_time = 1000; //numer of time step
  int N_time_resolve = 1; //numer of time step
  int N_parallel = 2; //number of parallel 
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
  if (argc > input_counter) alpha =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) beta =             boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  #if defined(HAMILTONIAN_VALENCE_FORCE_FIELD_MODEL_TAGGED_EXTERNAL_FIXED_END_HPP) 
  double frequency = 0.0;
  double intensity = 0.0;
  double phase = 0.0;
  if (argc > input_counter) frequency =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) intensity =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) phase =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  #endif
  if (argc > input_counter) T =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_parallel =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time_resolve =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;

  int num_particles = Lattice::set_num_particles(Ns); //nummer of particles
  int Ns_observe = num_particles / 10; //observed particle 
  int N_time_measure = N_time / N_time_resolve; //number of time span of histogram
  double dt = t / N_time; //time step length
  int N_each = N_loop / N_parallel; // number of each parallel
  condition_dat = result_directory + condition_dat;

  tools::DataRecorder dataput(condition_dat); 
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
          << "Coupling constant : beta = " << beta << std::endl
          #if defined(HAMILTONIAN_VALENCE_FORCE_FIELD_MODEL_TAGGED_EXTERNAL_FIXED_END_HPP) 
          << "Intensity of external field : intensity = " << intensity << std::endl
          << "Frequency of external field : frequency = " << frequency << std::endl
          << "phase of external field : phase = " << phase << std::endl
          #endif
          << "Number of patricles : num_particles = " << num_particles << std::endl
          << "Length of whole system : Ns = " << Ns << std::endl
          << "Length of observe system : Ns_observe = " << Ns_observe << std::endl
          << "Number of time resolve for measure: N_time_resolve = " << N_time_resolve << std::endl
          << "Number of time for measure: N_time_resolve = " << N_time_measure << std::endl
          << "Developing ime : t = " << t << std::endl
          << "Number of time steps : N_time = " << N_time << std::endl
          << "Time interbal : dt = " << dt << std::endl
          << "Number of parallelize : N_parallel =" << N_parallel << std::endl
          << "Number of loop : N_loop =" << N_loop << std::endl
          << "Each thread loop : N_each = " << N_each << std::endl 
          << "Number of time resolve for hist : N_time_resolve = " << N_time_resolve << std::endl
          << "Result directory : " << result_directory << std::endl;
 //end parameter define 

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
  for(auto key : physical_value_list){
    results[key].resize(N_time);
    results["error_" + key].resize(N_time);
    for(int i = 0 ; i < N_time ; ++i){
      results[key][i] = 0.0;
      results["error_"+key][i] = 0.0;
    }
  }

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
  for(auto key  : stat_value_list){
    stat_values[key].resize(N_time);
    for(int i = 0; i < N_time; ++i){
      stat_values[key][i] = 0.0;
    }
  }

  //set correlatin function 
  std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d{};  
  std::vector<std::string> correlation_list{"Velocity_correlation","Displace_correlation"}; 
  for(auto key : correlation_list) results_2d[key] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(Ns_observe,0.0));
  results_2d["Normalmode_Energy"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(num_particles,0.0));

  //set histogram maps 
  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_data; 
 	hist_data["Velocity"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(N_loop)),
 	hist_data["Displacement"] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(N_loop));
  //end set results and measure vectors

  dataput.time_tag() << "start sampling process : time development " << std::endl;

  #pragma omp parallel for num_threads(N_parallel) 
  for(int para_number = 0; para_number  < N_parallel ; ++para_number){

  //lattice set
  Lattice lattice(Ns);
  int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj)); //interaction table 
  lattice.create_table(pair_table);
  
  //randam generetor set 
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
 
  //dynamics set 
  std::vector<double> z(2*num_particles);
  Ensembler ensembler(T, num_particles);


  int observe = (num_particles-1) / 2;
  #if defined(HAMILTONIAN_VALENCE_FORCE_FIELD_MODEL_TAGGED_EXTERNAL_FIXED_END_HPP) 
  Hamiltonian  hamiltonian(num_particles,J,alpha,beta,observe,intensity,frequency,phase,pair_table,N_adj);
  #else 
  Hamiltonian hamiltonian(num_particles,J,alpha,beta,pair_table,N_adj); 
  #endif

  Integrator integrator(2*num_particles);

  //set autocorrelation
  std::unordered_map<std::string, correlation::AutoCorrelationFunction> autocorrelation;
  autocorrelation["VelocityACF"].initialize(N_time,N_loop);
  autocorrelation["DisplacementACF"].initialize(N_time,N_loop);

  //measure physical values set
  std::unordered_map<std::string, std::vector<tools::Accumulator> > physical_values;
  for(auto key : physical_value_list){
    physical_values[key].resize(N_time);
  }

  //set correlation
  std::unordered_map<std::string, std::vector<std::vector<tools::Accumulator> > > correlation_function;
  correlation_function["Displacement"].resize(N_time_measure);
  correlation_function["Velocity"].resize(N_time_measure);
  for(auto key : {"Displacement","Velocity"}){
    for(int step = 0 ; step < N_time_measure ; ++step) correlation_function[key][step] = std::vector<tools::Accumulator>(Ns_observe,tools::Accumulator());
  }

  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){
    double v_total = 0.0;
    double x_total = 0.0;
    double *x = &z[0];
    double *v = &z[num_particles];
    int target = 0;
    int measure_step = 0;
    double pt = 0.0;

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


      //input hist data 
      if((step % N_time_resolve) == 0){
        measure_step = step / N_time_resolve;
        hist_data["Velocity"][measure_step][para_number * N_each + c_loop] = v[observe];
        hist_data["Displacement"][measure_step][para_number * N_each + c_loop] = x_total;
        //correlation calc
        for(int i = 0; i < Ns_observe; ++i){
          target = lattice.target_number(i,observe);
          correlation_function["Displacement"][measure_step][i] << x[observe] * x[target];
          correlation_function["Velocity"][measure_step][i] << v[observe] * v[target];
        }
      }

      //time develp
      pt += dt;
      integrator.step(pt, dt, z, hamiltonian);
    }//end time step
  }//end loop

  int tid = omp_get_thread_num();
  int num_threads = omp_get_num_threads();
  {
    #pragma omp critical 
    dataput.time_tag() << "[" << tid <<"th thread] start result correct" << std::endl;
    for(auto key  : physical_value_list){
      for(int i = 0 ; i < N_time; ++i){
        results[key][i] += physical_values[key][i].mean() / N_parallel;
        results["error_" + key][i] += physical_values[key][i].error() / N_parallel;
      }
    }
    for(auto element : autocorrelation){
      std::string key = element.first;
      std::vector<double> acf_result = autocorrelation[key].result();
      for(int i = 0 ; i < N_time; ++i){
        results[key][i] += acf_result[i] / N_parallel;
      }      
    }
    for(int step = 0 ; step < N_time_measure; ++step){
      for(int i = 0 ; i < Ns_observe; ++i){
        results_2d["Displace_correlation"][step][i] += correlation_function["Displacement"][step][i].mean() / N_parallel;
        results_2d["Velocity_correlation"][step][i] += correlation_function["Velocity"][step][i].mean() / N_parallel;
      }
    }
    for(int step = 0 ; step < N_time; ++step){
      stat_values["mean"][step] += physical_values["Velocity"][step].mean()/ N_parallel;
      stat_values["variance"][step] += physical_values["Velocity"][step].variance()/ N_parallel;
      stat_values["skewness"][step] += physical_values["Velocity"][step].skewness()/ N_parallel;
      stat_values["kurtosis"][step] += physical_values["Velocity"][step].kurtosis()/ N_parallel;
      stat_values["kurtosis_exess"][step] += physical_values["Velocity"][step].kurtosis_excess()/ N_parallel;
      stat_values["average"][step] += physical_values["Velocity"][step].average()/ N_parallel;
      stat_values["error"][step] += physical_values["Velocity"][step].error()/ N_parallel;
      stat_values["central_moment1"][step] += physical_values["Velocity"][step].central_moment1()/N_parallel;
      stat_values["central_moment2"][step] += physical_values["Velocity"][step].central_moment2()/N_parallel;
      stat_values["central_moment3"][step] += physical_values["Velocity"][step].central_moment3()/N_parallel;
      stat_values["central_moment4"][step] += physical_values["Velocity"][step].central_moment4()/N_parallel;
    }
    
    dataput.time_tag() << "[" << tid <<"th thread] end result correct" << std::endl;
  }

  }//end parallerize
  dataput.time_tag() << "end sampling process" << std::endl;

  dataput.time_tag() << "start stat calc " << std::endl;
  //set Accumulator

  //define histogram class
  std::unordered_map<std::string, tools::Histogram>  histogram;
  histogram["Displacement"].initialize(n_bin);
  histogram["Velocity"].initialize(n_bin);


  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_values;
  std::vector<std::string> hist_value_list = {
                                              "hist_Velocity",
                                              "range_Velocity",
                                              "hist_Displacement",
                                              "range_Displacement"
                                              };

  for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(N_time_measure, std::vector<double>(n_bin));

  //calc histogram and statisitcal values
  for(int i = 0; i < N_time_measure; ++i){
    for(auto element : histogram){
      std::string key = element.first; 
      histogram[key](hist_data[key][i]);
      histogram[key].output(hist_values["hist_" + key][i],
                                        hist_values["range_" + key][i]);
    }
  }
  dataput.time_tag() << "end stat calc " << std::endl;

  //set stat time 
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

  //set sight number domain
  std::unordered_map< std::string, std::vector<double> > domain_sight;
  domain_sight["site"] = std::vector<double>(Ns_observe);
  for(int i = 0 ; i < Ns_observe; ++i) domain_sight["site"][i] = i;

  std::string result_dat("result.dat");
  std::string correlation_dat("result_correlation.dat");
  std::string hist_dat("result_hist.dat");
  std::string stat_dat("result_stat.dat");

  result_dat = result_directory + result_dat;
  correlation_dat = result_directory + correlation_dat;
  hist_dat = result_directory + hist_dat;
  stat_dat = result_directory + stat_dat;

  dataput << "result output to : " 
          << result_dat << " , "
          << correlation_dat << " , "
          << hist_dat << " and " 
          << stat_dat 
          << std::endl;

  dataput.output_result(result_dat,"pt",dt,N_time,results);
  dataput.output_result(stat_dat,"pt",dt,N_time,stat_values);
  dataput.output_result_hist(hist_dat,domain_time,hist_value_list,hist_values);
  dataput.output_result_3d(correlation_dat,domain_time,domain_sight,correlation_list,results_2d);

  dataput.time_tag() << " end time " << std::endl; 
  return 0;
} //END clXYmodel

