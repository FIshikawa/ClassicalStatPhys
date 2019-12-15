#include <ensemble/normalmode_ensemble_fftw.hpp>
#include <lattice/chain_open_boundary.hpp>
#include <physics/fpu_generalized_fixed_end.hpp>
#include <boost/lexical_cast.hpp>

using Ensembler = ensemble::NormalModeEnsembleFFTW; 
using Hamiltonian = hamiltonian::FPUGeneralizedFixedEnd;
using Lattice = lattice::ChainOpenBoundary;

struct Parameters{
  double J = 1.0; //interaction constant;
  double alpha = 1.0; //interaction constant;
  double beta = 1.0; //interaction constant;
  double gamma = 1.0; //interaction constant;
  double delta = 1.0; //interaction constant;
};

template <class Parameters>
inline void ParameterSet(int argc, char **argv, Parameters & parameters, int & input_counter){
  if (argc > input_counter) parameters.J =                boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.alpha =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.beta =             boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.gamma =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;  
  if (argc > input_counter) parameters.delta =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;  
}

template <class Parameters, class Dataput>
inline void ParameterDeclare(Parameters const parameters, Dataput & dataput){
    dataput << "Coupling constant : J = " << parameters.J << std::endl
            << "Coupling constant : alpha = " << parameters.alpha << std::endl
            << "Coupling constant : beta = " << parameters.beta << std::endl
            << "Coupling constant : gamma = " << parameters.gamma << std::endl
            << "Coupling constant : delta = " << parameters.delta << std::endl;
}

template <class Parameters,class Table>
inline Hamiltonian HamiltonianSet(Parameters const parameters, int num_particles, Table const pair_table, int const N_adj){
  Hamiltonian hamiltonian(num_particles,parameters.J,
      parameters.alpha,parameters.beta,parameters.gamma,parameters.delta,pair_table,N_adj);
  return hamiltonian;
}

#include <main/quasi_integrable_thermalization_fixed_end_conditions.hpp>
#include <main/quasi_integrable_thermalization.hpp>
