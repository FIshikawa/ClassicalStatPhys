#include <ensemble/modeoccupancy_ensemble_periodic_boundary_fftw.hpp>
#include <lattice/chain.hpp>
#include <physics/fpu_generalized.hpp>
#include <physics/todalattice.hpp>
#include <boost/lexical_cast.hpp>

using Ensembler = ensemble::ModeOccupancyEnsemblePeriodicBoundaryFFTW; 
using FirstHamiltonian= hamiltonian::TodaLattice;
using SecondHamiltonian = hamiltonian::FPUGeneralized;
using Lattice = lattice::Chain;

struct Parameters{
  double J_fpu = 1.0; //interaction constant;
  double J_toda = 1.0; //interaction constant;
  double alpha_fpu = 1.0; //interaction constant;
  double alpha_toda = 1.0; //interaction constant;
  double beta = 1.0; //interaction constant;
  double gamma = 1.0; //interaction constant;
  double delta = 1.0; //interaction constant;
};

template <class Parameters>
inline void ParameterSet(int argc, char **argv, Parameters & parameters, int & input_counter){
  if (argc > input_counter) parameters.J_fpu =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.alpha_fpu =        boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.beta =             boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.gamma =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;  
  if (argc > input_counter) parameters.delta =            boost::lexical_cast<double>(argv[input_counter]);++input_counter;  
  if (argc > input_counter) parameters.J_toda =           boost::lexical_cast<double>(argv[input_counter]);++input_counter;
  if (argc > input_counter) parameters.alpha_toda =       boost::lexical_cast<double>(argv[input_counter]);++input_counter;

  double diff = 1.0 - parameters.beta /( parameters.alpha_fpu * parameters.alpha_fpu * 2 / 3);
  if(std::abs(diff) < 1e-4)
    parameters.beta = parameters.alpha_fpu * parameters.alpha_fpu * 2 / 3;
  diff = 1.0 - parameters.gamma / (parameters.alpha_fpu * parameters.alpha_fpu * parameters.alpha_fpu / 3);
  if(std::abs(diff) < 1e-4)
    parameters.gamma = parameters.alpha_fpu * parameters.alpha_fpu * parameters.alpha_fpu / 3;
  diff = 1.0 - parameters.delta / (parameters.alpha_fpu * parameters.alpha_fpu * parameters.alpha_fpu * parameters.alpha_fpu * 2 / 15);
  if(std::abs(diff) < 1e-4)
    parameters.delta = parameters.alpha_fpu * parameters.alpha_fpu * parameters.alpha_fpu * parameters.alpha_fpu * 2 / 15;
}

template <class Parameters, class Dataput>
inline void ParameterDeclare(Parameters const parameters, Dataput & dataput){
    dataput << "Coupling constant of fpu : J_fpu = " << parameters.J_fpu << std::endl
            << "Coupling constant of fpu : alpha_fpu = " << parameters.alpha_fpu<< std::endl
            << "Coupling constant : beta = " << parameters.beta << std::endl
            << "Coupling constant : gamma = " << parameters.gamma << std::endl
            << "Coupling constant : delta = " << parameters.delta << std::endl
            << "Coupling constant of toda : J_toda = " << parameters.J_toda << std::endl
            << "Coupling constant of toda : alpha_toda = " << parameters.alpha_toda << std::endl;
}

template <class Parameters,class Table>
inline FirstHamiltonian FirstHamiltonianSet(Parameters const parameters, int num_particles, Table const pair_table, int const N_adj){
  FirstHamiltonian hamiltonian(num_particles,parameters.J_toda,parameters.alpha_toda,pair_table,N_adj);
  return hamiltonian;
}

template <class Parameters,class Table>
inline SecondHamiltonian SecondHamiltonianSet(Parameters const parameters, int num_particles, Table const pair_table, int const N_adj){
  SecondHamiltonian hamiltonian(num_particles,parameters.J_fpu,
      parameters.alpha_fpu,parameters.beta,parameters.gamma,parameters.delta,pair_table,N_adj);
  return hamiltonian;
}

#include <main/quasi_integrable_thermalization_periodic_boundary_conditions.hpp>
#include <main/quasi_integrable_thermalization_sampling.hpp>
