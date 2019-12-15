#include <omp.h>
#include <tools/histogram.hpp>
#include <tools/data_recorder.hpp>
#include <tools/auto_correlation_function.hpp>
#include <tools/accumulator.hpp>

// Ensembler, Integrator, Hamiltonian, Lattice are defined in each src header 
// If MOLECULAR_DYNAMICAL_EQUILIBRATION is defined, this is MD simulartion for Equilibrium
// In the others, MCMC dynamcis used.

int main(int argc, char **argv){
// std::map<std::string, double> condition_list;
  double J = 1.0;
  double T = 1.0; //start temperature
  double d_T = 0.1; //start temperature
  int Ns = 2; // lenght of whole system
  int N_thermalize = 2; // number of thermalization
  int N_loop = 2; //loop number
  int N_parallel = 2; //number of parallel 
  int N_T = 2; // number of T_list
  int Ns_observe = 2;
  int n_bin = 100;
  std::string condition_dat = "condi.dat";
  std::string result_dir("./");


  //parameter set 
  int input_counter = 1;
  if (argc > input_counter) result_dir =       boost::lexical_cast<std::string>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_thermalize =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) Ns =               boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) T =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) d_T =              boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_T =              boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) N_loop =           boost::lexical_cast<int>(argv[input_counter]);++input_counter;

    #ifdef MOLECULAR_DYNAMICAL_EQUILIBRATION
    double t = 10; // Whole time
    int N_time = 1000; //numer of time step
    if (argc > input_counter) N_time =         boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) t =              boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    double dt = t / N_time;
    #else 
    int N_relax = 2; // relaxation of equilbrium sampling
    if (argc > input_counter) N_relax =        boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    #endif

  if (argc > input_counter) N_parallel =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) Ns_observe =       boost::lexical_cast<int>(argv[input_counter]);++input_counter;
  if (argc > input_counter) n_bin =            boost::lexical_cast<int>(argv[input_counter]);++input_counter;

  // parameter auto set
  int N_each = N_loop/N_parallel;
  int num_particles = Lattice::set_num_particles(Ns);
  condition_dat = result_dir + condition_dat;

  std::vector<double> T_list(N_T,0.0); //temperture span
  for(int i = 0; i < N_T; ++i) T_list[i] = T + (double)i *  d_T;

  //parameter set 
  tools::DataRecorder dataput(condition_dat); 
  dataput.time_tag() << "start time" << std::endl; 
  dataput << "Observe Equilibrium Critical Behavior" << std::endl
          << Hamiltonian::name() << std::endl
          << Ensembler::name() << std::endl
          << Lattice::name() << std::endl;
  
 //condition declare
  dataput << "Coupling constant: J is " << J << std::endl
          << "Temperture: T = " << T << std::endl
          << "Number of Temperture List: N_T = " << N_T << std::endl
          << "Temperture list : T_list = " << T_list[0]; 
  for(int i = 1; i < N_T; ++i) dataput << " , " << T_list[i];
  dataput <<  std::endl;
  dataput << "Number of patricles :num_particles = " << num_particles << std::endl
          << "Length of whole system :Ns = " << Ns << std::endl
          << "Length of observed system :Ns_observe = " << Ns_observe << std::endl
          << "Number of thermalize step : N_thermalize = " << N_thermalize << std::endl
          << "Number of loop: N_loop =" << N_loop << std::endl
          #ifdef MOLECULAR_DYNAMICAL_EQUILIBRATION
          << "Developing ime : t = " << t << std::endl
          << "Number of time steps : N_time = " << N_time << std::endl
          #else 
          << "Interval of relax correlation: N_relax =" << N_relax << std::endl
          #endif
          << "Parallelize : N_parallel = " << N_parallel << std::endl
          << "Each thread loop : N_each = " << N_each << std::endl 
          << "Number of bins: n_bin =" << n_bin << std::endl
          << "Result into : " << result_dir << std::endl;
   
  //set results and measure vectors
  std::unordered_map<std::string, std::vector<double> > results{};  
  std::vector<std::string> physical_value_list = {
                                                  "Spin_total",
                                                  "Square",
                                                  "Quad",
                                                  "Energy_total",
                                                  "Kinetic_Energy",
                                                  "Potential_Energy",
                                                  "Spin_totalACF"
                                                };

  for(auto key : physical_value_list){
    results[key].resize(N_T);
    results["error_" + key].resize(N_T);
    for(int i = 0 ; i < N_T; ++i){
      results[key][i] = 0.0;
      results["error_"+key][i] = 0.0;
    }
  }


  std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d{};  
  std::vector<std::string> correlation_list{"Velocity_correlation","Spin_correlation"}; 
  for(auto key : correlation_list) results_2d[key] = std::vector< std::vector<double> >(N_T, std::vector<double>(Ns_observe,0.0));
  results_2d["Spin_totalACF"] = std::vector< std::vector<double> >(N_T, std::vector<double>(N_each,0.0));

  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_data; 
 	hist_data["Velocity"] = std::vector< std::vector<double> >(N_T, std::vector<double>(N_loop)),
 	hist_data["Spin"] = std::vector< std::vector<double> >(N_T, std::vector<double>(N_loop));
  //end set results and measure vectors

  dataput.time_tag() << "start T change calculation" << std::endl;

  for(int counter_T = 0; counter_T < N_T ; ++counter_T){
    double T_t = T_list[counter_T];
 
    dataput.time_tag() << "start " << counter_T << "th calculation " << std::endl;
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
    std::vector<double> z(2*num_particles,0.0);
    Ensembler ensembler(J,T,num_particles);
    Hamiltonian hamiltonian(num_particles,J,pair_table,N_adj); 

    #ifdef MOLECULAR_DYNAMICAL_EQUILIBRATION
    Integrator integrator(2*num_particles);
    #endif 

    //set correlation
    std::unordered_map<std::string, correlation::AutoCorrelationFunction> autocorrelation;
    autocorrelation["Spin_totalACF"].initialize(N_each,N_loop);
 
    //measure physical values set
    std::unordered_map<std::string, stat::accumulator> physical_values;
    for(auto key : physical_value_list){
      physical_values[key].reset();
    }
    std::unordered_map<std::string, std::vector<stat::accumulator> > correlation_function;
    correlation_function["Spin"] = std::vector<stat::accumulator>(Ns_observe,stat::accumulator());
    correlation_function["Velocity"] = std::vector<stat::accumulator>(Ns_observe,stat::accumulator());

    #ifndef MOLECULAR_DYNAMICAL_EQUILIBRATION
    //thermalize
    ensembler.equilibrate_velocity(z,mt);
    for(int counter = 0; counter < N_thermalize; ++counter) ensembler.montecalro(z,counter, hamiltonian, mt);
    #endif
    
    int target = 0;
    double *x = &z[0];
    double *v = &z[num_particles];
    double s_total_x,s_total_y;
    double r_strength = 0.0; 
    double square = 0.0; 
    int observe = (num_particles-1)/2;//set observe particle

    for(int c_loop=0; c_loop < N_each; ++c_loop){

      //equilibrate step
      #ifdef MOLECULAR_DYNAMICAL_EQUILIBRATION
      for(int counter = 0; counter < num_particles; ++counter) ensembler.montecalro(z,counter, hamiltonian ,mt);
      double pt = 0.0;
      for(int step = 0; step < N_time; ++step){
        #ifdef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_HPP 
        s_total_x = 0.0;
        s_total_y = 0.0;
        for(int i = 0; i < num_particles; ++i){
          s_total_x += std::cos(x[i]) ;
          s_total_y += std::sin(x[i]) ;
        }
        hamiltonian.set_meanfield(s_total_x, s_total_y); // only Fully-connected
        #endif 
        pt += dt;
        integrator.step(pt, dt, z, hamiltonian);
      }//end time step
      #else
      ensembler.equilibrate_velocity(z,mt);
      for(int step = 0; step < N_relax; ++step){
        for(int counter = 0; counter < N_thermalize;){  
          #ifdef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_HPP 
          s_total_x = 0.0;
          s_total_y = 0.0;
          for(int i = 0; i < num_particles; ++i){
            s_total_x += std::cos(x[i]) ;
            s_total_y += std::sin(x[i]) ;
          }
          hamiltonian.set_meanfield(s_total_x, s_total_y); // only Fully-connected
          #endif 
          ensembler.montecalro(z,counter, hamiltonian, mt);
        }
      }
      #endif

      //set total value
      s_total_x = 0.0;
      s_total_y = 0.0;
      for(int i = 0; i < num_particles; ++i){
        s_total_x += std::cos(x[i]) ;
        s_total_y += std::sin(x[i]) ;
      }

      #ifdef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_HPP 
      hamiltonian.set_meanfield(s_total_x, s_total_y); // only Fully-connected
      #endif 

      square = (s_total_x*s_total_x + s_total_y*s_total_y) / (num_particles * num_particles);
      physical_values["Square"] << square;
      physical_values["Quad"] << square * square;

      r_strength = std::sqrt((s_total_x * s_total_x + s_total_y * s_total_y)) / num_particles;
      physical_values["Spin_total"] << r_strength;
      autocorrelation["Spin_totalACF"] << r_strength;

      physical_values["Energy_total"] <<  hamiltonian.energy(0.0,z) / num_particles;
      physical_values["Kinetic_Energy"] << hamiltonian.kinetic_energy(0.0,z) / num_particles;
      physical_values["Potential_Energy"] << hamiltonian.potential_energy(0.0,z) / num_particles;

      hist_data["Velocity"][counter_T][p * N_each + c_loop] = v[observe];
      hist_data["Spin"][counter_T][p * N_each + c_loop] = r_strength;

      // correlation calc
      for(int i = 0; i < Ns_observe; ++i){
        target = lattice.target_number(i,observe);
        correlation_function["Spin"][i] << std::cos(x[observe] - x[target]);
        correlation_function["Velocity"][i] << v[observe] * v[target];
      }
  }//end loop

  int tid = omp_get_thread_num();
  int num_threads = omp_get_num_threads();
  for (int p = 0; p < num_threads; ++p) {
    #pragma omp critical 
    if (p == tid){
      for(auto key  : physical_value_list){
        results[key][counter_T] += physical_values[key].mean() / N_parallel;
        results["error_" + key][counter_T] += physical_values[key].error() / N_parallel;
      }
      for(auto element : autocorrelation){
        std::string key = element.first;
        std::vector<double> acf_result = autocorrelation[key].result();
        for(int i = 0 ; i < N_each; ++i){
          results_2d[key][counter_T][i] += acf_result[i] / N_parallel;
        }
      }
      for(int i = 0 ; i < Ns_observe; ++i){
        results_2d["Spin_correlation"][counter_T][i] += correlation_function["Spin"][i].mean() / N_parallel;
        results_2d["Velocity_correlation"][counter_T][i] += correlation_function["Velocity"][i].mean() / N_parallel;
      }
    }
  }

  }//end parallerize
  dataput.time_tag() << "end " << counter_T << "th calculation " << std::endl;
  }// end T loop
  dataput.time_tag() << "end T change calculation" << std::endl;

  #ifdef HAMILTONIAN_CLASSICAL_XY_FULLYCONNECTED_HPP 
  dataput.time_tag() << "start exact partition " << std::endl;
  //set calc partition function
  exactsolution::PartitionFunctionClassicalXYFullyConnected partition_function(J, num_particles);
  double exact_kinetic_energy, exact_potential_energy, exact_temperture, exact_magnetization;
  results["exact_Kinetic_Energy"].resize(N_T);
  results["exact_Potential_Energy"].resize(N_T);
  results["exact_Spin_total"].resize(N_T);
  results["exact_Temperture"].resize(N_T);
  for(int i = 0 ; i < N_T; ++i){
    exact_magnetization = partition_function.M(T_list[i]);
    exact_potential_energy = 0.5 * J * ( 1.0 -  exact_magnetization * exact_magnetization); 
    exact_kinetic_energy = 0.5 * T_list[i] - exact_potential_energy;
    results["exact_Kinetic_Energy"][i] = exact_kinetic_energy;
    results["exact_Potential_Energy"][i] = exact_potential_energy;
    results["exact_Spin_total"][i] = exact_magnetization;
  }
  dataput.time_tag() << "end exact partition" << std::endl;
  #endif 

  std::unordered_map< std::string, std::vector<double> > domain;
  domain["Temperture"] = std::vector<double>(N_T);
  for(int i = 0; i < N_T; ++i ) domain["Temperture"][i] = T_list[i];

  dataput.time_tag() << "start correlation calc " << std::endl;

  dataput.time_tag() << "end correlation calc " << std::endl;

  //histgoram calc
  dataput.time_tag() << "start stat calc " << std::endl;
  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_values;
  std::vector<std::string> hist_value_list = {
                                              "hist_Velocity",
                                              "range_Velocity",
                                              "hist_Spin",
                                              "range_Spin"
                                              };
  for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(N_T, std::vector<double>(n_bin));
  std::unordered_map<std::string, tools::Histogram>  histogram;
  histogram["Spin"].initialize(n_bin);
  histogram["Velocity"].initialize(n_bin);
  //calc histogram and statisitcal values
  for(int i = 0; i < N_T; ++i){
    for(auto element : histogram){
      std::string key = element.first; 
      histogram[key](hist_data[key][i]);
      histogram[key].output(hist_values["hist_" + key][i],
                                        hist_values["range_" + key][i]);
    }
  }
  dataput.time_tag() << "end stat calc " << std::endl;

  std::string result_dat("result.dat");
  std::string correlation_dat("result_correlation.dat");
  std::string hist_dat("result_hist.dat");

  result_dat = result_dir + result_dat;
  correlation_dat = result_dir + correlation_dat;
  hist_dat = result_dir + hist_dat;

  dataput << "result output to : " 
          << result_dat << " , "
          << hist_dat << " and " 
          << correlation_dat 
          << std::endl;

  dataput.output_result(result_dat,domain,results);
  dataput.output_result_hist(hist_dat,domain,hist_value_list,hist_values);

  dataput.time_tag() << " end time " << std::endl; 
  return 0;
} //END clXYmodel

