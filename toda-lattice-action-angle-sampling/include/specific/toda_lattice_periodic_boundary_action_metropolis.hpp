#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <boost/lexical_cast.hpp>
#include <clstatphys/ensemble/modeoccupancy_ensemble_periodic_boundary_fftw.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <clstatphys/ensemble/action_metropolis.hpp>
#include <common/settings_common.hpp>

using MonteCarloSamplerOrigin = ensemble::MetropolisActionToda;
using Ensembler = ensemble::ModeOccupancyEnsemblePeriodicBoundaryFFTW; 
using Hamiltonian = hamiltonian::TodaLattice;

class MonteCarloSampler: public MonteCarloSamplerOrigin{
public:
  MonteCarloSampler(int num, int num_iteration, std::vector<double> betas, 
      double J, double alpha, double dp, double dx) 
    : MonteCarloSamplerOrigin(num,num_iteration,betas,J,alpha,dp,dx), num_(num){} 
  
  template<class Rand>
  void montecarlo(std::vector<double> & z, Rand & mt){
    int counter = 0;
    int counter_t = 0;
    while(counter_t < num_){ 
      MonteCarloSamplerOrigin::montecarlo(z, counter, mt);
      counter_t += 1;
    }
  }
private:
  int num_;
}; //end MonteCarloSampler definition

struct Settings : public SettingsCommon{
  int N_normalmode = 5; //numer of time step
  int k_initial = 0; //start wave vector filled by initialization 
  double E_initial = 1.0; // initial energy
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;
  double dx = 0.01; //interaction constant;
  double dp = 0.01; //interaction constant;

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

    if (argc > input_counter) dp =    boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) dx =    boost::lexical_cast<double>(argv[input_counter]);++input_counter;

    if (argc > input_counter) J =     boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) alpha = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  }

  template <class Dataput>
  inline void declare(Dataput & dataput){
    SettingsCommon::declare(dataput);

    dataput <<  "<<System Depend Settings>> " << std::endl
            <<  "  " << Hamiltonian::name() << std::endl
            <<  "  " << Ensembler::name() << std::endl
            <<  "  " << Lattice::name() << std::endl

            << "  Energy initial : E_initial =" << E_initial << std::endl
            << "  Number of non-zero energy normal modes : N_normalmode =" << N_normalmode << std::endl
            << "  Start wave vector filled by initialization : k_initial =" << k_initial << std::endl
 
            << "  Difference for propose : dx = " << dx << std::endl
            << "  Difference for propose : dp = " << dp << std::endl
 
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

  Ensembler ensembler(){
    Ensembler ensembler_t(
                          num_particles,
                          k_initial,
                          N_normalmode,
                          E_initial
                          );
    return ensembler_t;
  }
  
  MonteCarloSampler monte_carlo_sampler(){
    std::vector<double> betas(num_particles-1,1.0);
    MonteCarloSampler monte_carlo_sampler_t(
                                         num_particles, 
                                         num_iterations, 
                                         betas,
                                         J,
                                         alpha,
                                         dp,
                                         dx
                                         );
    return monte_carlo_sampler_t;
  }

};

#endif
