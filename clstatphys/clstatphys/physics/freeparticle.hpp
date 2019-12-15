#ifndef HAMILTONIAN_FREEPARTICLE_HPP
#define HAMILTONIAN_FREEPARTICLE_HPP

#include <string>
#include <cmath>
#include <vector>
#include <omp.h>

namespace hamiltonian{

class FreeParticle{ // assuming two body interactions 
public:
  // z[0...n-1]: position, z[n...2n-1]: momentum
  static std::string name() { return "Free particle (non-interacted system) Dynamics"; }
  FreeParticle(int num_particles, double J, std::vector<std::vector<int > > table, int N_adj) : n_(num_particles), J_(J), table_(table),Nd_(N_adj){}

  double energy(const double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }
  double potential_energy(const double  t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_ ; ++i) ene += target_potential_energy(i,z,t);
    return ene/2.0;
  }
  double kinetic_energy(const double t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_; ++i) ene += target_kinetic_energy(i,z);
    return ene;
  }

  double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    double ene = 0;
    return ene;
  }

  double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[n_];
    double ene = 0;
    ene += 0.5 * v[l]*v[l];
    return ene;
  }

  // "force" calculation
  void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double *fx = &force[0];
    double *fv = &force[n_];
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 0; i < n_ ; ++i) fv[i] = 0.0;
  }

  int nd()const { return nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  int n_;
  double j_;
  std::vector<std::vector<int> > table_;
  int nd_;
};//end free_particle

} //end namespace

#endif //HAMILTONIAN_FREEPARTICLE_HPP
