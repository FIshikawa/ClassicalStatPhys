#ifndef HAMILTONIAN_CLASSICAL_HEISENBERG_HPP
#define HAMILTONIAN_CLASSICAL_HEISENBERG_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class ClassicalHeisenberg {
public:
  // z[0...n-1]: position, z[n...2n-1]: momentum
  static std::string name() { return "Classical Heisenberg Model"; }
  ClassicalHeisenberg(int num_particles, double J, std::vector<std::vector<int > > table, int N_adj) : n_(num_particles), J_(J), table_(table),Nd_(N_adj){}

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
    const double *y = &z[n_];
    double ene = 0;
    for(int i = 0; i < Nd_ ; ++i){
      double dx = x[l] - x[table_[l][i]];
      double y_n = y[table_[l][i]];
      ene -= J_ *(1 - std::cos(dx) * std::cos(y[l]) * std::cos(y_n) + std::sin(y[l])*std::sin(y_n));
    }
    return ene;
  }
  double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[2*n_];
    const double *w = &z[3*n_];
    double ene = 0;
    ene += 0.5 * v[l]*v[l] + 0.5 * w[l]+w[l];
    return ene;
  }

  // "force" calculation
  void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *y = &z[n_];
    const double *v = &z[2*n_];
    const double *w = &z[3*n_];
    double *fx = &force[0];
    double *fy = &force[n_];
    double *fv = &force[2*n_];
    double *fw = &force[3*n_];
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 0; i < n_ ; ++i) fy[i] = w[i];
    for(int i = 0; i < n_ ; ++i){
      fv[i] = 0.0;
      for(int d = 0; d < Nd_; ++d){
        fv[i] += -1.0 * J_ * std::sin(x[i] -x[table_[i][d]]) * std::cos(y[i]) * std::cos(y[table_[i][d]]);
      }
    }
    for(int i = 0; i < n_ ; ++i){
      fw[i] = 0.0;
      for(int d = 0; d < Nd_; ++d){
        fw[i] += 1.0 * J_ *( - std::cos(x[i] -x[table_[i][d]]) * std::sin(y[i]) * std::cos(y[table_[i][d]]) + std::cos(y[i]) * std::sin(y[table_[i][d]]));
      }
    }
  }

  int Nd()const { return Nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

protected:
  int n_;
  double J_;
  std::vector<std::vector<int> > table_;
  int Nd_;

};

}//end namespace

#endif //HAMILTONIAN_CLASSICAL_HEISENBERG_HPP
