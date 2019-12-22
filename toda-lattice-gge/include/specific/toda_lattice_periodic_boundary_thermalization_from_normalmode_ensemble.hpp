#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <boost/lexical_cast.hpp>
#include <clstatphys/ensemble/modeoccupancy_ensemble_periodic_boundary_fftw.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <common/settings_common.hpp>

using Ensembler = ensemble::ModeOccupancyEnsemblePeriodicBoundaryFFTW; 
using Hamiltonian = hamiltonian::TodaLattice;

struct Settings : public SettingsCommon{
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;

  Settings(int argc, char **argv, int & input_counter) 
    : SettingsCommon(argc, argv, input_counter) 
  { 
    set(argc, argv, input_counter);
  }
  Settings() = default;

  inline void set(int argc, char **argv, int & input_counter){
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
  
};

#endif
