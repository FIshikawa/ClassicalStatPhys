#ifndef HAMILTONIAN_VALENCE_FORCE_FIELD_MODEL_TAGGED_EXTERNAL_FIXED_END_HPP
#define HAMILTONIAN_VALENCE_FORCE_FIELD_MODEL_TAGGED_EXTERNAL_FIXED_END_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class ValenceForceFieldModelTaggedExternalFixedEnd{ 
public:
  static std::string name() { return "Classical Fermi Pasta Ulam : fixed End : must use Open Boundary Chain Lattice"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  ValenceForceFieldModelTaggedExternalFixedEnd(int num_particles, double J, double alpha, double beta, int observe, double intensity, double frequency, double phase, std::vector<std::vector<int > > table, int N_adj) 
    : n_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha),
      beta_(beta),intensity_(intensity),frequency_(frequency),phase_(phase), observe_(observe){}

  double energy(double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }
  double potential_energy(double  t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_ ; ++i) ene += target_potential_energy(i,z,t);
    return ene;
  }
  double kinetic_energy(double /* t */, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_; ++i) ene += target_kinetic_energy(i,z);
    return ene;
  }
  double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    double ene = 0;
    double x_nn_sum = 0;
    double dx = 0;
    double dx_2 = 0;
    for(int i = 0; i < Nd_; ++i){
      if(table_[l][i] == l){
        dx = x[l];
        x_nn_sum += 0;
      }
      else{
        dx = x[l] - x[table_[l][i]];
        x_nn_sum += x[table_[l][i]];
      }
      ene +=  0.125 * alpha_ *  dx * dx * dx * dx;
      for(int j = i+1; j < Nd_; ++j){
        if(table_[l][j] == l) dx_2 = x[l];
        else dx_2 = x[l] - x[table_[l][j]];
        ene += beta_ * dx * dx * dx_2 * dx_2;
      }
    }
    ene += J_ * (x_nn_sum - 3.0 * x[l]) * (x_nn_sum - 3.0 * x[l]);
    if(l == observe_) ene -= x[l] * intensity_ * std::sin(frequency_ * t + phase_);
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
    double x_nn_sum;
    double x_nnn_sum;
    double dx;
    double dx_2;
    int site_nn;
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 0; i < n_ ; ++i){
      fv[i] = 0.0;
      x_nn_sum = 0;
      //center calc
      for(int d = 0; d < Nd_; ++d){
        if(table_[i][d] == i) dx = x[i];
        else{
          site_nn = table_[i][d];
          x_nn_sum += x[site_nn];
          dx = x[i] - x[site_nn];
          x_nnn_sum = x[i];
          for(int j = 0 ; j < Nd_; ++j){
            if(table_[site_nn][j] != i){
              if(table_[site_nn][j] == site_nn) dx_2 = x[site_nn];
              else{
                dx_2 = x[site_nn] - x[table_[site_nn][j]];
                x_nnn_sum += x[table_[site_nn][j]];
              }
              fv[i] -= 2.0 * beta_ * dx * dx_2 * dx_2;
            }
          }
          fv[i] -= 2.0 * J_ * (x_nnn_sum - 3.0 * x[site_nn]);
          fv[i] -= 0.5 * alpha_ * dx * dx * dx ;
        }
        fv[i] -= 0.5 * alpha_ * dx * dx * dx ;
        for(int j = d+1; j < Nd_; ++j){
          if(table_[i][j] == i) dx_2 = x[i];
          else dx_2 = x[i] - x[table_[i][j]];
          fv[i] -= 2.0 * beta_ * (dx * dx_2 * dx_2 + dx_2 * dx * dx);
        }
      }
      fv[i] += 6.0 * J_ * (x_nn_sum - 3.0 * x[i]);
      if(i == observe_) fv[i] += intensity_ * std::sin(frequency_ * t + phase_);
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
  double frequency_;
  double intensity_;
  double phase_;
  double observe_;
};//end FPUfixedB

} //end namespace

#endif //HAMILTONIAN_VALENCE_FORCE_FIELD_MODEL_TAGGED_EXTERNAL_FIXED_END_HPP
