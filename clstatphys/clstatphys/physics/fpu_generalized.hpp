#ifndef HAMILTONIAN_FPU_GENERALIZED_HPP
#define HAMILTONIAN_FPU_GENERALIZED_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class FPUGeneralized {
public:
  static std::string name() { return "Classical Generalized Fermi Pasta Ulam Hamiltonian"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  FPUGeneralized(int num_particles, double J, double alpha, double beta, double gamma, double delta, 
      std::vector<std::vector<int > > table, int N_adj) :  
    n_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha),beta_(beta),gamma_(gamma),delta_(delta){}
  FPUGeneralized() = default;

  double energy(double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }

  // The displacements are defined as r_n = u_n - u_n+1
  // and the Toda potential is defined as r_n^2 + r_n^3 + r_n^4 + r_n^5 + r_n^6 

  double potential_energy(double  t, std::vector<double> const& z) const {
    const double *x = &z[0];
    double ene = 0, dx = 0;
    for(int i = 0; i < n_; ++i){
      dx = x[table_[i][0]] - x[i];
      ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25 + delta_ * dx * dx * dx * dx * dx * dx / 6.0;
      ene += alpha_ * dx * dx * dx / 3.0 + gamma_ * dx * dx * dx * dx * dx / 5.0;
    }
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
    for(int i = 0; i < n_ ; ++i){
      fv[i] = 0.0;

      dx = x[table_[i][0]] - x[i];
      fv[i] += J_ * dx  +  beta_ * dx * dx * dx + delta_ * dx * dx * dx * dx * dx;
      fv[i] +=  alpha_ * dx * dx + gamma_ * dx * dx * dx * dx;

      dx = x[i] - x[table_[i][1]];
      fv[i] -= J_ * dx  +  beta_ * dx * dx * dx + delta_ * dx * dx * dx * dx * dx;
      fv[i] -=  alpha_ * dx * dx + gamma_ * dx * dx * dx * dx;
    }
  }

  int Nd()const { return Nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  int n_;
  int Nd_;
  double J_;
  double alpha_, beta_, gamma_, delta_;
  std::vector<std::vector<int> > table_;
};//end FPU

} //end namespace

#endif //HAMILTONIAN_FPU_GENERALIZED_HPP
