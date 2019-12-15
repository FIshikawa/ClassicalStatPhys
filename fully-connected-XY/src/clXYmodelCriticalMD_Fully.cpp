#define MOLECULAR_DYNAMICAL_EQUILIBRATION
#include <integrator/velocity_velret.hpp>
using Integrator =  integrator::VelocityVelret;
#include <ensemble/water_bag.hpp>
using Ensembler =  ensemble::WaterBag;
#include <physics/classical_xy_fullyconnected.hpp>
#include <physics/partition_function_classical_xy_fullyconnected.hpp>
using Hamiltonian =  hamiltonian::ClassicalXYFullyConnected;
#include <lattice/fullyconnected.hpp>
using Lattice = lattice::FullyConnected;
#include <main/clXYmodelCritical.hpp>

