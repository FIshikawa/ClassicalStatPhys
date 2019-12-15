#define MOLECULAR_DYNAMICAL_EQUILIBRATION
#include <integrator/velocity_velret.hpp>
using Integrator =  integrator::VelocityVelret;
#include <ensemble/water_bag.hpp>
using Ensembler =  ensemble::WaterBag;
#include <physics/classical_xy.hpp>
using Hamiltonian =  hamiltonian::ClassicalXY;
#include <lattice/regularcube.hpp>
using Lattice = lattice::RegularCube;
#include <main/clXYmodelCritical.hpp>
