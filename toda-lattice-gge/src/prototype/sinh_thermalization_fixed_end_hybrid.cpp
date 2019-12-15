#include <physics/sinh_lattice_fixed_end_parallel.hpp>
#include <lattice/chain_open_boundary.hpp>
#include <ensemble/normalmode_ensemble_fftw.hpp>

using Ensembler = ensemble::NormalModeEnsembleFFTW; 
using Hamiltonian = hamiltonian::SinhLatticeFixedEndParallel;
using Lattice = lattice::ChainOpenBoundary;

#include<main/quasi_integrable_thermalization_hybrid.hpp>
