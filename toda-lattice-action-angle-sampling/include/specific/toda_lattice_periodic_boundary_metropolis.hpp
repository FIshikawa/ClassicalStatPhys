#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <boost/lexical_cast.hpp>
#include <clstatphys/ensemble/hybrid_monte_carlo.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <common/settings_common.hpp>

using EnsemblerOrigin = ensemble::HybridMonteCarlo; 
using Hamiltonian = hamiltonian::TodaLattice;

class Ensembler : public EnsemblerOrigin{
public:
  template<class Hamiltonian>
  Ensembler(
            int num_particles, double temperture, double dt_relax, 
            double relax_time, int total_accept, int initial_relax_time,
            Hamiltonian hamiltonian)
    : EnsemblerOrigin(num_particles, temperture, dt_relax, relax_time, total_accept),
      initial_relax_time_(initial_relax_time) {hamiltonian_ = hamiltonian;}
  
  template<class Rand>
  void set_initial_state(std::vector<double> & z, Rand & mt){
    int counter = 0;
    for(int step = 0; step < initial_relax_time_; ++step) EnsemblerOrigin::montecarlo(z, counter, hamiltonian_, mt);
  }

private:
  int initial_relax_time_;
  Hamiltonian hamiltonian_;
}; //end Ensembler definition

struct Settings : public SettingsCommon{
  int total_accept = 10;
  double initial_relax_time = 10;
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
    if (argc > input_counter) total_accept       = boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) initial_relax_time = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) relax_time         = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) dt_relax           = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) temperture         = boost::lexical_cast<double>(argv[input_counter]);++input_counter;

    if (argc > input_counter) J                  = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) alpha              = boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  }

  template <class Dataput>
  inline void declare(Dataput & dataput){
    SettingsCommon::declare(dataput);

    dataput <<  "<<System Depend Settings>> " << std::endl
            <<  "  " << Hamiltonian::name() << std::endl
            <<  "  " << EnsemblerOrigin::name() << std::endl
            <<  "  " << Lattice::name() << std::endl

            << "  Number of step for fianally accept : total_accept = " << total_accept << std::endl
            << "  Relaxtion time for update : relax_time = " << relax_time << std::endl
            << "  Interbal of time development : dt_relax  = " << dt_relax << std::endl
            << "  Temperture : temperture = " << temperture << std::endl

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
                          temperture,
                          dt_relax,
                          relax_time,
                          total_accept,
                          initial_relax_time,
                          hamiltonian()
                          );
    return ensembler_t;
  }

}; //end Settings definition



#endif
