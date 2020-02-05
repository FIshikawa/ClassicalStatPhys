#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <boost/lexical_cast.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <clstatphys/ensemble/hybrid_monte_carlo.hpp>
#include <clstatphys/ensemble/modeoccupancy_ensemble_periodic_boundary_fftw.hpp>
#include <common/settings_common.hpp>

using Ensembler = ensemble::ModeOccupancyEnsemblePeriodicBoundaryFFTW; 
using Hamiltonian = hamiltonian::TodaLattice;
using MonteCarloSamplerOrigin = ensemble::HybridMonteCarlo; 

class MonteCarloSampler: public MonteCarloSamplerOrigin{
public:
  template<class Hamiltonian>
  MonteCarloSampler(
            int num_particles, double temperture, double dt_relax, 
            double relax_time, int total_accept, Hamiltonian hamiltonian)
    : MonteCarloSamplerOrigin(num_particles, temperture, dt_relax, relax_time, total_accept)
  {hamiltonian_ = hamiltonian;}
  
  template<class Rand>
  void montecarlo(std::vector<double> & z, Rand & mt){
    int counter = 0;
    MonteCarloSamplerOrigin::montecarlo(z, counter, hamiltonian_, mt);
  }

private:
  Hamiltonian hamiltonian_;
}; //end MonteCarloSampler definition

struct Settings : public SettingsCommon{
  int N_normalmode = 5; //numer of time step
  int k_initial = 0; //start wave vector filled by initialization 
  double E_initial = 1.0; // initial energy

  int total_accept = 10;
  double relax_time = 10;
  double dt_relax = 0.1;
  double temperture = 1.0;

  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;

  Settings(int argc, char **argv, int & input_counter) 
    : SettingsCommon(argc, argv, input_counter) 
  { 
    set(argc, argv, input_counter);
  }
  Settings() = default;

  inline void set(int argc, char **argv, int & input_counter){
    if (argc > input_counter) E_initial =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) k_initial =        boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) N_normalmode =     boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (Ns < N_normalmode + k_initial){
      std::cerr << "k_initial + N_noramalmode should be lower than Ns" << std::endl;
      std::exit(1);
    } 

    if (argc > input_counter) total_accept       = boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) relax_time         = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) dt_relax           = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) temperture         = boost::lexical_cast<double>(argv[input_counter]);++input_counter;

    if (argc > input_counter) J =     boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) alpha = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  }

  template <class Dataput>
  inline void declare(Dataput & dataput){
    SettingsCommon::declare(dataput);

    dataput <<  "<<System Depend Settings>> " << std::endl
            <<  "  " << Hamiltonian::name() << std::endl
            <<  "  " << Ensembler::name() << std::endl
            <<  "  " << MonteCarloSampler::name() << std::endl

            // declare normalmode ensemble 
            << "  Energy initial : E_initial =" << E_initial << std::endl
            << "  Number of non-zero energy normal modes : N_normalmode =" << N_normalmode << std::endl
            << "  Start wave vector filled by initialization : k_initial =" << k_initial << std::endl

            // declare hybrid monte carlo
            << "  Number of step for fianally accept : total_accept = " << total_accept << std::endl
            << "  Relaxtion time for update : relax_time = " << relax_time << std::endl
            << "  Interbal of time development : dt_relax  = " << dt_relax << std::endl
            << "  Temperture : temperture = " << temperture << std::endl

            // declare hamiltonian 
            << "  Coupling constant : J = " << J << std::endl
            << "  Coupling constant : alpha = " << alpha << std::endl;
  }

  Hamiltonian hamiltonian(){
    Lattice lattice_t(Ns);

    std::vector<std::vector<int> > pair_table(
                                              num_particles,
                                              std::vector<int>(N_adj)
                                              ); 
    lattice_t.create_table(pair_table);

    Hamiltonian hamiltonian_t(
                              num_particles,
                              J,
                              alpha,
                              pair_table,
                              N_adj
                              );
    return hamiltonian_t;
  }

  MonteCarloSampler monte_carlo_sampler(){
    MonteCarloSampler monte_carlo_sampler_t(
                                            num_particles,
                                            temperture,
                                            dt_relax,
                                            relax_time,
                                            total_accept,
                                            hamiltonian()
                                            );
    return monte_carlo_sampler_t;
  }

  Ensembler ensembler(){
    Ensembler ensembler_t(
                          num_particles,
                          k_initial,
                          N_normalmode,
                          E_initial
                          );
    return ensembler_t;
  }

};

#endif
