#ifndef HAMILTONIAN_FPU_FIXED_END_PARALLEL_HPP
#define HAMILTONIAN_FPU_FIXED_END_PARALLEL_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class FPUFixedEndParallel { 
public:
  static std::string name() { return "Classical Fermi Pasta Ulam : fixed End : must use Open Boundary Chain Lattice : Parallel"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  FPUFixedEndParallel(int num_particles, double J, double alpha, double beta, std::vector<std::vector<int > > table, int N_adj) : n_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha),beta_(beta){}

  double energy(double t, std::vector<double> const& z) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double ene = 0;

    #pragma omp parallel for reduction(+:ene)
    for (int l = 0; l < n_; ++l){
      ene += 0.5 * v[l]*v[l] * 2.0;
      for(int i = 0; i < Nd_; ++i){
        double dx = x[l] - x[table_[l][i]];
        ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
      }
      for(int i = 0; i < Nd_/2; ++i){
        double dx = x[table_[l][i]] - x[l];
        ene += alpha_ * dx * dx * dx / 3.0;
      }
      for(int i = Nd_/2; i < Nd_ ; ++i){
        double dx = x[l] - x[table_[l][i]];
        ene += alpha_ * dx * dx * dx / 3.0;
      }
    }
    ene /= 2.0;
    double dx = - x[0];
    ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    ene += alpha_ * dx * dx * dx / 3.0;
    dx = x[n_ - 1];
    ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    ene += alpha_ * dx * dx * dx / 3.0;
    return ene;
  }

  double kinetic_energy(double /* t */, std::vector<double> const& z) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double ene = 0;

    #pragma omp parallel for reduction(+:ene)
    for (int l = 0; l < n_ ; ++l){
      ene += 0.5 * v[l]*v[l];
    }
    return ene;
  }

  double potential_energy(double /* t */, std::vector<double> const& z) const {
    const double *x = &z[0];
    double ene = 0;

    #pragma omp parallel for reduction(+:ene)
    for (int l = 0; l < n_; ++l){
      for(int i = 0; i < Nd_; ++i){
        double dx = x[l] - x[table_[l][i]];
        ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
      }
      for(int i = 0; i < Nd_/2; ++i){
        double dx = x[table_[l][i]] - x[l];
        ene += alpha_ * dx * dx * dx / 3.0;
      }
      for(int i = Nd_/2; i < Nd_ ; ++i){
        double dx = x[l] - x[table_[l][i]];
        ene += alpha_ * dx * dx * dx / 3.0;
      }
    }
    ene /= 2.0;
    double dx = - x[0];
    ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    ene += alpha_ * dx * dx * dx / 3.0;
    dx = x[n_ - 1];
    ene += J_ * dx * dx * 0.5 + beta_ * dx * dx * dx * dx * 0.25;
    ene += alpha_ * dx * dx * dx / 3.0;
    return ene;
  }

  // "force" calculation
  void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double *fx = &force[0];
    double *fv = &force[n_];

    #pragma omp parallel for 
    for(int i = 0; i < n_ ; ++i){ 
      fx[i] = v[i];
      fv[i] = 0.0;
      for(int d = 0; d < Nd_; ++d){
        double dx = (x[i] -x[table_[i][d]]);
        fv[i] -=   J_ * dx  +  beta_ * dx * dx * dx;
      }
      for(int d = 0; d < Nd_/2; ++d){
        double dx = x[table_[i][d]] - x[i];
        fv[i] += alpha_ * dx * dx;
      }
      for(int d = Nd_/2; d < Nd_ ; ++d){
        double dx = x[i] -x[table_[i][d]];
        fv[i] -=  alpha_ * dx * dx ;
      }
    }
    for(auto& i : {0, n_ - 1}) fv[i] -= J_ * x[i] + beta_ * x[i] * x[i] * x[i];
    fv[0] += alpha_ * x[0] * x[0];
    fv[n_ - 1] -= alpha_ * x[n_ - 1] * x[n_ - 1];
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
