#define MOLECULAR_DYNAMICAL_EQUILIBRATION
#include <integrator/velocity_velret.hpp>
using Integrator =  integrator::VelocityVelret;
#include <ensemble/water_bag.hpp>
using Ensembler =  ensemble::WaterBag;
#include <physics/classical_xy.hpp>
using Hamiltonian =  hamiltonian::ClassicalXY;
#include <lattice/tesseract.hpp>
using Lattice = lattice::Tesseract;
#include <main/clXYmodelCritical.hpp>
