#ifndef HAMILTONIAN_FPU_HPP
#define HAMILTONIAN_FPU_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class FPU {
public:
  static std::string name() { return "Classical Fermi Pasta Ulam Hamiltonian"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  FPU(int num_particles, double J, double alpha, double beta, std::vector<std::vector<int > > table, int N_adj) :  n_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha),beta_(beta){}

  double energy(double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }
  double potential_energy(double  t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_ ; ++i) ene += target_potential_energy(i,z,t);
    return ene/2.0;
  }
  double kinetic_energy(double /* t */, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_; ++i) ene += target_kinetic_energy(i,z);
    return ene;
  }

  double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    double ene = 0;
    for(int i = 0; i < Nd_; ++i){
      double dx = x[l] - x[table_[l][i]];
      ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    }
    for(int i = 0; i < Nd_/2; ++i){
      double dx = x[l] - x[table_[l][i]];
      ene += -1.0 * alpha_ * dx * dx * dx / 3.0;
    }
    for(int i = Nd_/2; i < Nd_ ; ++i){
      double dx = x[l] - x[table_[l][i]];
      ene += alpha_ * dx * dx * dx / 3.0;
    }
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
    for(int i = 0; i < n_ ; ++i){
      fv[i] = 0.0;
      for(int d = 0; d < Nd_; ++d){
        double dx = (x[i] -x[table_[i][d]]);
        fv[i] += -1.0 * (  J_ * dx  +  beta_ * dx * dx * dx);
      }
      for(int d = 0; d < Nd_/2; ++d){
        double dx = (x[i] -x[table_[i][d]]);
        fv[i] += 1.0 * alpha_ * dx * dx;
      }
      for(int d = Nd_/2; d < Nd_ ; ++d){
        double dx = (x[i] -x[table_[i][d]]);
        fv[i] += -1.0 * alpha_ * dx * dx ;
      }
    }
  }

  int Nd()const { return Nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  int n_;
  double J_;
  std::vector<std::vector<int> > table_;
  int Nd_;
  double alpha_;
  double beta_;
};//end FPU

} //end namespace

#endif //HAMILTONIAN_FPU_HPP
