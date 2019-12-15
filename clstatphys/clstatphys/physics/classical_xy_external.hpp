#ifndef HAMILTONIAN_CLASSICAL_XY_EXTERNAL_HPP
#define HAMILTONIAN_CLASSICAL_XY_EXTERNAL_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class ClassicalXYExternal {
public:
  // z[0...n-1]: position, z[n...2n-1]: momentum
  static std::string name() { return "Classical XYspin Dynamics with External Field : Periodic and Isotropic"; }
  ClassicalXYExternal(int num_particles, double J, std::vector<std::vector<int > > table, int N_adj, double Inte, double Freq, double D) : n_(num_particles), J_(J), table_(table),Nd_(N_adj),Inte_(Inte),Freq_(Freq),D_(D) {}

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

  double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[n_];
    double ene = 0;
    ene += 0.5 * v[l]*v[l];
    return ene;
  }

  double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    double ene = 0;
    for(int i = 0; i < Nd_ ; ++i){
      double dx = x[l] - x[table_[l][i]];
      ene += J_ *(1.0 -  std::cos(dx));
    }
    ene -= D_ * std::cos(x[l]) ;
    ene += Inte_ * std::cos(Freq_ * t) * std::cos(x[l]);
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
      for(int d = 0; d < Nd_; ++d){
        fv[i] += -1.0 * J_ * std::sin(x[i] -x[table_[i][d]]);
      }
      fv[i] += ( - D_ * std::sin(x[i]) + Inte_ * std::cos( Freq_ * t) * std::sin(x[i]));
    }
  }

  int Nd()const { return Nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  double D_, Freq_, Inte_;
  int n_;
  double J_;
  std::vector<std::vector<int> > table_;
  int Nd_;
};

} //end namespace

#endif // HAMILTONIAN_CLASSICAL_XY_EXTERNAL_HPP 
