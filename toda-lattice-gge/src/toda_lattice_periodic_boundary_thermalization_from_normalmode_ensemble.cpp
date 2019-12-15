#include <clstatphys/ensemble/modeoccupancy_ensemble_periodic_boundary_fftw.hpp>
#include <clstatphys/lattice/chain.hpp>
#include <clstatphys/physics/todalattice.hpp>
#include <clstatphys/boost/lexical_cast.hpp>

using Ensembler = ensemble::ModeOccupancyEnsemblePeriodicBoundaryFFTW; 
using Hamiltonian = hamiltonian::TodaLattice;
using Lattice = lattice::Chain;

struct SettingDepend{
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;
  Hamiltonian hamiltonian;
  Ensembler ensembler;
  Lattice lattice;

  template<struct SettingCommon>
  inline void set(int argc, char **argv, SettingCommon & setting_common, int & input_counter){
    if (argc > input_counter) J =     boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) alpha = boost::lexical_cast<double>(argv[input_counter]);++input_counter;

    Lattice lattice_t(parameter.common.Ns);
    lattice = lattice_t;

    std::vector<std::vector<int> > pair_table(
                                              setting_common.num_particles,
                                              std::vector<int>(setting_common.N_adj)
                                              ); 
    lattice.create_table(pair_table);

    Hamiltonian hamiltonian_t(
                              setting_common.common.num_particles,
                              setting_common.depend.J,
                              setting_common.depend.alpha,
                              pair_table,
                              setting_common.common.N_adj
                              );

    hamiltonian = hamiltonian_t;

    Ensembler ensembler_t(
                          setting_common.common.num_particles,
                          setting_common.common.k_initial,
                          setting_common.common.N_normalmode,
                          setting_common.common.E_initial
                          );

    ensembler = ensembler_t;
  }
  
  template <class Dataput>
  inline void declare(Dataput & dataput){
    dataput <<  "System Depend settings " << std::endl
            <<  Hamiltonian::name() << std::endl
            <<  Ensembler::name() << std::endl
            <<  Lattice::name() << std::endl

            << "Coupling constant : J = " << J << std::endl
            << "Coupling constant : alpha = " << alpha << std::endl;
  }
};

struct ResultsDepend{
  inline void set(){};
  inline void output(){};
};

#include<configurations.hpp>
