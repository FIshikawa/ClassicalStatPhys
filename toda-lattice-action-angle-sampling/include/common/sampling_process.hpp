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
  int N_time = settings.N_time;
  int N_total_data = settings.N_total_data;
  double dt = settings.dt;

  // defined some classes
  MonteCarloSampler monte_carlo_sampler = settings.monte_carlo_sampler();
  Ensembler ensembler = settings.ensembler();
  Hamiltonian hamiltonian = settings.hamiltonian();
  tools::TimeMeasure time_measure = settings.time_measure();
  
  std::vector<double> z(2*num_particles);

  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){
    // set initial state
    ensembler.set_initial_state(z,mt);

    //time develop
    double pt = 0.0;
    int measure_step = 0;
    time_measure.reset();
    for(int step = 0; step < N_time; ++step){  
      //set total values 
      if(time_measure.check(step)) physical_quantities.measure(z, measure_step);
      //time develp
      pt += dt;
      monte_carlo_sampler.montecarlo(z, mt);
    }//end time step
    if(measure_step != N_total_data){
      std::cerr << "final measure_step : " << measure_step << ", N_total_data : " << N_total_data <<", not same " << std::endl;
      std::exit(112);
    }
  }//end loop
};

#endif
