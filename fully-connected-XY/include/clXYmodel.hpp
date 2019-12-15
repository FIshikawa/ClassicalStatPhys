#include <omp.h>
#include <tools/histogram.hpp>
#include <tools/data_recorder.hpp>
#include <tools/auto_correlation_function.hpp>
#include <tools/accumulator.hpp>

//Ensembler, Integrator, Hamiltonian, Lattice are defined in each src header 

int main(int argc, char **argv){

  //parameter define
  double J = 1.0; //std::pow(2.0*M_PI,2.0);
  double T = 0.5; //temperature
  int Ns = 3; // lenght of whole system
  int N_thermalize = 3; // number of thermalization
  int N_loop = 1; //loop number
  int N_parallel = 1; //number of parallel 
  double t = 10; // Whole time
  int N_time = 1000; //numer of time step
  int N_time_resolve = 1; //numer of time step
  int Ns_observe = 3;
  std::string condition_dat = "condi.dat";
  std::string result_dir("./");
  int n_bin = 1000;

 //parameter set 
  int input_counter = 1;
  if (argc > input_counter) result_dir =       boost::lexical_cast<std::string>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_thermalize =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) Ns =               boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) T =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) t =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_parallel =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_time_resolve =   boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) Ns_observe =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;

  // parameter auto set
  int N_each = N_loop/N_parallel;
  int num_particles = Lattice::set_num_particles(Ns);
  double dt = t/N_time;
  int N_hist_length = N_time / N_time_resolve;
  condition_dat = result_dir + condition_dat;

 //parameter set 
  tools::DataRecorder dataput(condition_dat); 
  dataput.time_tag() << "start time " << std::endl;
  dataput << "Declare Configulation " << std::endl 
          << "Equilibration check : classical XY dynamics" << std::endl
          <<  Hamiltonian::name() << std::endl
          <<  Ensembler::name() << std::endl
          <<  Integrator::name() << std::endl
          <<  Lattice::name() << std::endl;

 //condition declare
  dataput << "Coupling constant: J is " << J << std::endl
          << "Expected Temperture: T = " << T << std::endl
          << "Reduced Coupling: J/T = " << J/T << std::endl
          << "Number of patricles :num_particles = " << num_particles << std::endl
          << "Length of whole system :Ns = " << Ns << std::endl
          << "Number of thermalize step : N_thermalize = " << N_thermalize << std::endl
          << "Developing ime : t = " << t << std::endl
          << "Number of time steps : N_time = " << N_time << std::endl
          << "Number of time resolve for hist : N_time_resolve = " << N_time_resolve << std::endl
          << "Time interbal : dt = " << dt << std::endl
          << "Number of bins: n_bin =" << n_bin << std::endl
          << "Number of loop: N_loop =" << N_loop << std::endl
          << "Parallelize : N_parallel = " << N_parallel << std::endl
          << "Each thread loop : N_each = " << N_each << std::endl 
          << "Result into : " << result_dir << std::endl;
 //end parameter define 
 
  //set results and measure vectors
  std::unordered_map<std::string, std::vector<double> > results{};  
  std::vector<std::string> physical_value_list = {
                                                  "Spin",
                                                  "Spin_total",
                                                  "Velocity_total",
                                                  "Velocity",
                                                  "Energy_total",
                                                  "Kinetic_Energy",
                                                  "Potential_Energy",
                                                  "VelocityACF",
                                                  "SpinACF",
                                                  "Spin_totalACF"
                                                };

  for(auto key : physical_value_list){
    results[key].resize(N_time);
    results["error_" + key].resize(N_time);
    for(int i = 0 ; i < N_time ; ++i){
      results[key][i] = 0.0;
      results["error_"+key][i] = 0.0;
    }
  }

  //set correlatin function 
  std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d{};  
  std::vector<std::string> correlation_list{"Velocity_correlation","Spin_correlation"}; 
  for(auto key : correlation_list) results_2d[key] = std::vector< std::vector<double> >(N_time, std::vector<double>(Ns_observe,0.0));

  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_data; 
 	hist_data["Velocity"] = std::vector< std::vector<double> >(N_hist_length, std::vector<double>(N_loop)),
 	hist_data["Spin"] = std::vector< std::vector<double> >(N_hist_length, std::vector<double>(N_loop));
  //end set results and measure vectors

  //start calcluration 
  dataput.time_tag() << "start calculation " << std::endl;
  #pragma omp parallel for num_threads(N_parallel) 
  for(int p = 0; p < N_parallel ; ++p){
  //pararrel tag number def
  int parallel_tag = p;

  //lattice set
  Lattice lattice(Ns);
  int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
  lattice.create_table(pair_table);
  
  //randam genereta set 
  std::size_t seed = p;
  std::mt19937 mt(seed);
 
  //dynamics set 
  std::vector<double> z(2*num_particles);
  Ensembler ensembler(J,T,num_particles);
  Hamiltonian hamiltonian(num_particles,J,pair_table,N_adj); 
  Integrator integrator(2*num_particles);

  //set correlation
  std::unordered_map<std::string, correlation::AutoCorrelationFunction> autocorrelation;
  autocorrelation["VelocityACF"].initialize(N_time,N_loop);
  autocorrelation["SpinACF"].initialize(N_time,N_loop);
  autocorrelation["Spin_totalACF"].initialize(N_time,N_loop);

  //measure physical values set
  std::unordered_map<std::string, std::vector<stat::accumulator> > physical_values;
  for(auto key : physical_value_list){
    physical_values[key].resize(N_time);
  }
  std::unordered_map<std::string, std::vector<std::vector<stat::accumulator> > > correlation_function;
  correlation_function["Spin"].resize(N_time);
  correlation_function["Velocity"].resize(N_time);
  for(auto key : {"Spin","Velocity"}){
    for(int step = 0 ; step < N_time ; ++step) correlation_function[key][step] = std::vector<stat::accumulator>(Ns_observe,stat::accumulator());
  }


  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){
    //time develop
    double pt = 0.0;
    int hist_step = 0;
    int target = 0;
    double *x = &z[0];
    double *v = &z[num_particles];
    double s_total_x = 0.0;
    double s_total_y = 0.0;
    double r_strength = 0.0; 
    double v_total = 0.0;
    int observe = (num_particles-1)/2;//set observe particle

    ensembler.equilibrate_velocity(z,mt);
    for(int counter = 0; counter < N_thermalize; ++counter) ensembler.montecalro(z, counter, hamiltonian, mt);
    for(int step = 0; step < N_time; ++step){  
      v_total = 0.0;
      for(int i = 0; i < num_particles; ++i){ 
        v_total += v[i] / num_particles;
      }
      physical_values["Velocity_total"][step] << v_total;
      autocorrelation["VelocityACF"] << v[observe] - v_total;
      physical_values["Spin"][step] << std::cos(x[observe]);
      autocorrelation["SpinACF"] << std::cos(x[observe]);
      physical_values["Velocity"][step] << v[observe];
      s_total_x = 0.0;
      s_total_y = 0.0;
      for(int i = 0; i < num_particles; ++i){
        s_total_x += std::cos(x[i]) ;
        s_total_y += std::sin(x[i]) ;
      }

      #ifdef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_HPP 
      hamiltonian.set_meanfield(s_total_x, s_total_y); // only Fully-connected
      #endif 

      r_strength = std::sqrt((s_total_x * s_total_x + s_total_y * s_total_y)) / num_particles;
      physical_values["Spin_total"][step] << r_strength;
      autocorrelation["Spin_totalACF"] << r_strength;

      physical_values["Energy_total"][step] <<  hamiltonian.energy(pt,z) / num_particles;
      physical_values["Kinetic_Energy"][step] << hamiltonian.kinetic_energy(pt,z) / num_particles;
      physical_values["Potential_Energy"][step] << hamiltonian.potential_energy(pt,z) / num_particles;

      // correlation calc
      for(int i = 0; i < Ns_observe; ++i){
        target = lattice.target_number(i,observe);
        correlation_function["Spin"][step][i] << std::cos(x[observe] - x[target]);
        correlation_function["Velocity"][step][i] << v[observe] * v[target];
      }
 
      //hist calc
      if((step % N_time_resolve) == 0){
        hist_step = step / N_time_resolve;
        hist_data["Velocity"][hist_step][p * N_each + c_loop] = v[observe];
        hist_data["Spin"][hist_step][p * N_each + c_loop] = r_strength;
      }

      pt += dt;
      integrator.step(pt, dt, z, hamiltonian);

    }//end time step
  }//end loop


  int tid = omp_get_thread_num();
  int num_threads = omp_get_num_threads();
  for (int p = 0; p < num_threads; ++p) {
    #pragma omp critical 
    if (p == tid){
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
      for(int step = 0 ; step < N_time; ++step){
        for(int i = 0 ; i < Ns_observe; ++i){
          results_2d["Spin_correlation"][step][i] += correlation_function["Spin"][step][i].mean() / N_parallel;
          results_2d["Velocity_correlation"][step][i] += correlation_function["Velocity"][step][i].mean() / N_parallel;
        }
      }
    }
  }
  }//end parallelized
  dataput.time_tag() << "end calculation" << std::endl;

  #ifdef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_HPP 
  dataput.time_tag() << "start exact partition " << std::endl;
  //set calc partition function
  exactsolution::PartitionFunctionClassicalXYFullyConnected partition_function(J, num_particles);
  double exact_kinetic_energy, exact_potential_energy, exact_temperture, exact_magnetization;
  exact_temperture = partition_function.T(0.5 * T);
  exact_magnetization = partition_function.M(exact_temperture);
  exact_potential_energy = 0.5 * J * ( 1.0 -  exact_magnetization * exact_magnetization); 
  exact_kinetic_energy = 0.5 * T - exact_potential_energy;
  results["exact_Kinetic_Energy"].resize(N_time);
  results["exact_Potential_Energy"].resize(N_time);
  results["exact_Spin_total"].resize(N_time);
  results["exact_Temperture"].resize(N_time);
  for(int i = 0 ; i < N_time ; ++i){
    results["exact_Kinetic_Energy"][i] = exact_kinetic_energy;
    results["exact_Potential_Energy"][i] = exact_potential_energy;
    results["exact_Spin_total"][i] = exact_magnetization;
    results["exact_Temperture"][i] = exact_temperture;
  }
  dataput.time_tag() << "end exact partition" << std::endl;
  #endif 


  dataput.time_tag() << "start stat calc " << std::endl;
  //set accumulator
  stat::accumulator accum("stat variable");

  //define histogram class
  std::unordered_map<std::string, tools::Histogram>  histogram;
  histogram["Spin"].initialize(n_bin);
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
                                              "hist_Spin",
                                              "range_Spin"
                                              };

  for(auto key  : stat_value_list) stat_values[key].resize(N_hist_length);
  for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(N_hist_length, std::vector<double>(n_bin));

  //calc histogram and statisitcal values
  for(int i = 0; i < N_hist_length; ++i){
    accum.reset();   
    for(int j = 0; j < N_loop; ++j){
      accum << hist_data["Velocity"][i][j];
    }
    stat_values["mean"][i] = accum.mean();
    stat_values["variance"][i] = accum.variance();
    stat_values["skewness"][i] = accum.skewness();
    stat_values["kurtosis"][i] = accum.kurtosis();
    stat_values["kurtosis_exess"][i] = accum.kurtosis_excess();
    stat_values["average"][i] = accum.average();
    stat_values["error"][i] = accum.error();
    stat_values["central_moment1"][i] = accum.central_moment1();
    stat_values["central_moment2"][i] = accum.central_moment2();
    stat_values["central_moment3"][i] = accum.central_moment3();
    stat_values["central_moment4"][i] = accum.central_moment4();
    for(auto element : histogram){
      std::string key = element.first; 
      histogram[key](hist_data[key][i]);
      histogram[key].output(hist_values["hist_" + key][i],
                                        hist_values["range_" + key][i]);
    }
  }
  dataput.time_tag() << "end stat calc " << std::endl;

  //set sight number domain
  std::unordered_map< std::string, std::vector<double> > domain_sight;
  domain_sight["sight"] = std::vector<double>(Ns_observe);
  for(int i = 0 ; i < Ns_observe; ++i) domain_sight["sight"][i] = i;

  //set time domain
  std::unordered_map< std::string, std::vector<double> > domain_time;
  domain_time["time"] = std::vector<double>(N_time);
  for(int i = 0 ; i < N_time; ++i) domain_time["time"][i] = dt * i;

  //set stat time 
  std::unordered_map< std::string, std::vector<double> > domain;
  domain["time"] = std::vector<double>(N_hist_length);
  double pt = 0;
  double counter = 0;
  for(int step = 0; step < N_time; ++step){
    if((step % N_time_resolve) == 0){
      domain["time"][counter] = pt;
      counter += 1;
    }
    pt += dt;
  }

  std::string result_dat("result.dat");
  std::string correlation_dat("result_correlation.dat");
  std::string hist_dat("result_hist.dat");
  std::string stat_dat("result_stat.dat");

  result_dat = result_dir + result_dat;
  correlation_dat = result_dir + correlation_dat;
  hist_dat = result_dir + hist_dat;
  stat_dat = result_dir + stat_dat;

  dataput << "result output to : " 
          << result_dat << " , "
          << correlation_dat << " , "
          << hist_dat << " and " 
          << stat_dat 
          << std::endl;

  dataput.output_result(result_dat,"pt",dt,N_time,results);
  dataput.output_result_3d(correlation_dat,domain_time,domain_sight,correlation_list,results_2d);
  dataput.output_result_hist(hist_dat,domain,hist_value_list,hist_values);
  dataput.output_result(stat_dat,domain,stat_values);

  dataput.time_tag() << " end time " << std::endl; 
  return 0;
} //END clXYmodel

