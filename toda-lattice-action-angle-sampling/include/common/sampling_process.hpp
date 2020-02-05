#ifndef SAMPLING_PROCESS_HPP
#define SAMPLING_PROCESS_HPP

 using RandomGenerator = std::mt19937_64;

template<typename Settings, typename PhysicalQuanties>
void SamplingProcess(Settings & settings, PhysicalQuanties & physical_quantities){
  //
  std::size_t seed = 1234*(1 + settings.process_id);
  std::minstd_rand seed_gen(seed);
  RandomGenerator mt(seed_gen());

  int num_particles = settings.num_particles;
  int process_id = settings.process_id;
  int N_each = settings.N_each;
  int N_mc = settings.N_mc;
  int N_total_data = settings.N_total_data;

  // defined some classes
  MonteCarloSampler monte_carlo_sampler = settings.monte_carlo_sampler();
  Ensembler ensembler = settings.ensembler();
  Hamiltonian hamiltonian = settings.hamiltonian();
  
  std::vector<double> z(2*num_particles);

  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){
    // set initial state
    ensembler.set_initial_state(z,mt);

    int measure_step = 0;
    for(int step = 0; step < N_mc; ++step){  
      //set total values 
      physical_quantities.measure(z, measure_step);
      //time develp
      monte_carlo_sampler.montecarlo(z, mt);
    }//end time step
    if(measure_step != N_total_data){
      std::cerr << "final measure_step : " << measure_step << ", N_total_data : " << N_total_data <<", not same " << std::endl;
      std::exit(112);
    }
  }//end loop
};

#endif
