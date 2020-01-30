#ifndef HAMILTONIAN_FPU_FIXED_END_HPP
#define HAMILTONIAN_FPU_FIXED_END_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class FPUFixedEnd { 
public:
  static std::string name() { return "Classical Fermi Pasta Ulam : fixed End : must use Open Boundary Chain Lattice"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  FPUFixedEnd(int num_particles, double J, double alpha, double beta, std::vector<std::vector<int > > table, int N_adj) : n_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha),beta_(beta){}

  double energy(double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }

  // The displacements are defined as r_n = u_n - u_n+1
  // and the Toda potential is defined as r_n^2 + r_n^3 + r_n^4 + r_n^5 + r_n^6 

  double potential_energy(double  t, std::vector<double> const& z) const {
    const double *x = &z[0];
    double ene = 0, dx = 0;
    for(int i = 1; i < n_; ++i){
      dx = x[table_[i][0]] - x[i];
      ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
      ene += alpha_ * dx * dx * dx / 3.0;
    }
    dx = - x[0];
    ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    ene += alpha_ * dx * dx * dx / 3.0;

    dx =   x[n_-1];
    ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    ene += alpha_ * dx * dx * dx / 3.0;

    return ene;
  }

  double kinetic_energy(double /* t */, std::vector<double> const& z) const {
    double ene = 0;
    const double *v = &z[n_];
    for(int i = 0; i < n_; ++i) ene += 0.5 * v[i] * v[i];
    return ene;
  }

  // "force" calculation
  void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double *fx = &force[0];
    double *fv = &force[n_];
    double dx;
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 1; i < n_-1 ; ++i){
      fv[i] = 0.0;

      dx = x[table_[i][0]] - x[i];
      fv[i] += J_ * dx  +  beta_ * dx * dx * dx;
      fv[i] +=  alpha_ * dx * dx;

      dx = x[i] - x[table_[i][1]];
      fv[i] -=  alpha_ * dx * dx;
      fv[i] -= J_ * dx  +  beta_ * dx * dx * dx;
    }
    fv[0] = 0.0;

    dx = - x[0];
    fv[0] += J_ * dx  +  beta_ * dx * dx * dx;
    fv[0] +=  alpha_ * dx * dx;

    dx = x[0] - x[table_[0][1]];
    fv[0] -=  alpha_ * dx * dx;
    fv[0] -= J_ * dx  +  beta_ * dx * dx * dx;

    fv[n_-1] = 0.0;

    dx = x[table_[n_-1][0]] - x[n_-1];
    fv[n_-1] += J_ * dx  +  beta_ * dx * dx * dx;
    fv[n_-1] +=  alpha_ * dx * dx;

    dx = x[n_-1];
    fv[n_-1] -= J_ * dx  +  beta_ * dx * dx * dx;
    fv[n_-1] -=  alpha_ * dx * dx;

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

};//end FPUfixedB

} //end namespace

#endif //HAMILTONIAN_FPU_FIXED_END_HPP
