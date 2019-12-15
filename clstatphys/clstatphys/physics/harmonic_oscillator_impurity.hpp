#ifndef HAMILTONIAN_HARMONIC_OSCILLATOR_IMPURITY_HPP
#define HAMILTONIAN_HARMONIC_OSCILLATOR_IMPURITY_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class HarmonicOscillatorImpurity {
public:
  static std::string name() { return "Classical Harmonic Oscillator (impurity)"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  HarmonicOscillatorImpurity(int num_particles, double J, double m, std::vector<std::vector<int > > table, int N_adj) : n_(num_particles), J_(J), table_(table),Nd_(N_adj), m_(m), observe_(num_particles - 1 ){}

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
    for(int i = 0; i < Nd_ ; ++i){
      double dx = x[l] - x[table_[l][i]];
      ene += 0.5* J_ * dx * dx;
    }
    return ene;
  }

  double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[n_];
    double ene = 0;
    if(l == observe_) ene += 0.5 * v[l]*v[l] * m_ ;
    else ene += 0.5 * v[l]*v[l];
    return ene;
  }

  // "force" calculation
  void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double *fx = &force[0];
    double *fv = &force[n_];
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 0; i < n_ ; ++i){
      fv[i] = 0.0;
      if(i == observe_){
        for(int d = 0; d < Nd_; ++d){
        fv[i] += -1.0 * J_ / m_  * (x[i] -x[table_[i][d]]);
        }
      }
      else{
        for(int d = 0; d < Nd_; ++d){
        fv[i] += -1.0 * J_ *(x[i] -x[table_[i][d]]);
        }
      }
    }
  }

  int mass()const {return m_;}

  int nd()const { return nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  int n_;
  double j_;
  std::vector<std::vector<int> > table_;
  int nd_;
  double m_;
  int observe_;
};

} //end namespace

#endif //HAMILTONIAN_HARMONIC_OSCILLATOR_IMPURITY_HPP
