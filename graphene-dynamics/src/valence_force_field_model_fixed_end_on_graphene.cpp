#include <ensemble/water_bag.hpp>
#include <physics/valence_force_field_model_fixed_end.hpp>
#include<lattice/graphene_square_open_boundary.hpp>

using Ensembler = ensemble::WaterBag; 
using Hamiltonian = hamiltonian::ValenceForceFieldModelFixedEnd;
using Lattice = lattice::GrapheneSquareOpenBoundary;

#include<main/quasi_integrable_thermalization_2d.hpp>
