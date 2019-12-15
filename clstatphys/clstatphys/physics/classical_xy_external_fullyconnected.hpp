#ifndef CLASSICAL_XY_EXTERNAL_FULLYCONNECTED_HPP
#define CLASSICAL_XY_EXTERNAL_FULLYCONNECTED_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class ClassicalXYExternalFullyConnected {
public:
  // z[0...n-1]: position, z[n...2n-1]: momentum
  static std::string name() { return "Classical XYspin Dynamics Fully Connected with External Field : Periodic and Isotropic"; }
  ClassicalXYExternalFullyConnected(int num_particles, double J, std::vector<std::vector<int > > table, int N_adj, double Inte, double Freq, double D) : n_(num_particles), J_(J / (double)num_particles), table_(table),Nd_(N_adj),Inte_(Inte),Freq_(Freq),D_(D) {s_total_x_ = 0.0; s_total_y_ = 0.0;}

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
    ene =  1.0 * J_ * ( n_ - std::cos(x[l]) * s_total_x_ - std::sin(x[l]) * s_total_y_ );
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
    double s_total_x = 0.0;
    double s_total_y = 0.0;
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 0; i < n_ ; ++i){
      s_total_x += std::cos(x[i]) ;
      s_total_y += std::sin(x[i]) ;
    }

    for(int i = 0; i < n_ ; ++i){
      fv[i] = - 1.0 * J_  *  (  std::sin(x[i]) * s_total_x - std::cos(x[i]) * s_total_y );
      fv[i] += ( - D_ * std::sin(x[i]) +  Inte_ * std::cos( Freq_ * t) * std::sin(x[i]));
    }
  }

  void set_meanfield(double s_total_x, double s_total_y) {
    s_total_x_ = s_total_x ;
    s_total_y_ = s_total_y ;
  }

  double check_meanfield(int i) const {
    if(i == 0) return s_total_x_;
    else return s_total_y_;
  }

  int Nd() const { return Nd_ ;}
  int table(int i, int j) const { return table_[i][j];} 

private:
  int n_;
  double J_;
  std::vector<std::vector<int> > table_;
  int Nd_;
  double s_total_y_;
  double s_total_x_;
  double D_; 
  double Freq_;
  double Inte_;
};//end classicalXYNoeq

} //end namespace

#endif //CLASSICAL_XY_EXTERNAL_FULLYCONNECTED_HPP
