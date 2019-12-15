#include <physics/fpu_fixed_end_parallel.hpp>
#include <lattice/chain_open_boundary.hpp>
#include <ensemble/water_bag_parallel.hpp>

using Ensembler = ensemble::WaterBagParallel; 
using Hamiltonian = hamiltonian::FPUFixedEndParallel;
using Lattice = lattice::ChainOpenBoundary;

#include<main/quasi_integrable_thermalization_hybrid.hpp>
