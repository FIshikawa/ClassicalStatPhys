#include <ensemble/modeoccupancy_ensemble_periodic_boundary_fftw.hpp>
#include <lattice/chain.hpp>
#include <physics/todalattice.hpp>
#include <boost/lexical_cast.hpp>

using Ensembler = ensemble::ModeOccupancyEnsemblePeriodicBoundaryFFTW; 
using Hamiltonian = hamiltonian::TodaLattice;
using Lattice = lattice::Chain;

struct Parameters{
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;
};

template <class Parameters>
inline void ParameterSet(int argc, char **argv, Parameters & parameters, int & input_counter){
  if (argc > input_counter) parameters.J =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.alpha =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;
}

template <class Parameters, class Dataput>
inline void ParameterDeclare(Parameters const parameters, Dataput & dataput){
    dataput << "Coupling constant : J = " << parameters.J << std::endl
            << "Coupling constant : alpha = " << parameters.alpha << std::endl;
}

template <class Parameters,class Table>
inline Hamiltonian HamiltonianSet(Parameters const parameters, int num_particles, Table const pair_table, int const N_adj){
  Hamiltonian hamiltonian(num_particles,parameters.J,parameters.alpha,pair_table,N_adj);
  return hamiltonian;
}

#include <main/quasi_integrable_thermalization_periodic_boundary_conditions.hpp>
#include <main/quasi_integrable_thermalization.hpp>
