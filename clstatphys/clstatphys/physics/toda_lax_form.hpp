#ifndef TODA_LAX_FORM_HPP
#define TODA_LAX_FORM_HPP

#include <cmath>
#include <vector>
#include <rokko/lapack/syev.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>

namespace integrable{

class TodaLaxForm{
public:
  TodaLaxForm(int num_particles, double J, double alpha, std::vector<std::vector<int > > table, int N_adj) : num_particles_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha){}

  rokko::dlmatrix L_matrix(std::vector<double>& z){
    rokko::dlmatrix L = rokko::dlmatrix::Zero(num_particles_, num_particles_);
    L(0, 0) = toda_momentum(z, 0);
    L(0, 1) = toda_position(z, 0);
    L(0, num_particles_-1) = toda_position(z, num_particles_-1);
    for(int i = 1; i < num_particles_-1; ++i){
      L(i, i)   = toda_momentum(z, i); 
      L(i, i+1) = toda_position(z, i);
      L(i, i-1) = toda_position(z, i-1);
    }
    L(num_particles_-1, num_particles_-1) = toda_momentum(z, num_particles_-1);
    L(num_particles_-1, 0) = toda_position(z, num_particles_-1);
    L(num_particles_-1, num_particles_-2) = toda_position(z, num_particles_-2);
    return L;
  }

  rokko::dlmatrix B_matrix(std::vector<double>& z){
    rokko::dlmatrix B(num_particles_, num_particles_);
    B(0, 1) = -1.0 * toda_position(z, 0);
    B(0, num_particles_-1) = toda_position(z, num_particles_-1);
    for(int i = 1; i < num_particles_-1; ++i){
      B(i, i+1) = toda_position(z, i);
      B(i, i-1) = -1.0 * toda_position(z, i-1);
    }
    B(num_particles_-1, 0) = -1.0 * toda_momentum(z, num_particles_-1);
    B(num_particles_-1, num_particles_-2) = toda_momentum(z, num_particles_-2);
    return B;
  }

  double toda_position(std::vector<double>& z, int const& n){
    const double *x = &z[0];
    return 0.5 * std::sqrt(J_) * alpha_ * std::exp(0.5*alpha_*(x[table_[n][1]]-x[n]));
  }

  double toda_momentum(std::vector<double>& z, int const& n){
    const double *p = &z[num_particles_];
    return 0.5 * alpha_ * p[n];
  }

  void conservations(std::vector<double>& z, std::vector<double>& conservations){
    std::vector<double> w(num_particles_);
    eigenvalues(z,w);
    std::vector<double> w_temp = w;
    int number_of_conservations = conservations.size();
    double total_value = 0;
    for(int i = 0 ; i < number_of_conservations; ++i){
      conservations[i] = 0.0;
      for(int j = 0; j < num_particles_; ++j){
        conservations[i] += w_temp[j];
        w_temp[j] = w_temp[j] * w[j];
      }
    }
  }

  void conservations_with_eigenvalues(
                                      std::vector<double>& z,
                                      std::vector<double>& conservations,
                                      std::vector<double>& eigenvalues_vect
                                      ){
    std::vector<double> w(num_particles_);
    eigenvalues(z,w);
    eigenvalues_vect = w;
    std::vector<double> w_temp = w;
    int number_of_conservations = conservations.size();
    double total_value = 0;
    for(int i = 0 ; i < number_of_conservations; ++i){
      conservations[i] = 0.0;
      for(int j = 0; j < num_particles_; ++j){
        conservations[i] += w_temp[j];
        w_temp[j] = w_temp[j] * w[j];
      }
    }
  }
  void eigenvalues(std::vector<double>& z, std::vector<double>& w){
    rokko::dlmatrix L = L_matrix(z);
    rokko::dlvector w_temp(num_particles_);
    int info = rokko::lapack::syev('V', 'U', L, w_temp);
    for(int i = 0; i < num_particles_; ++i) w[i] = w_temp[i];
  }

private:
  int Nd_;
  int num_particles_;
  double J_;
  double alpha_;
  std::vector<std::vector<int> > table_;
};//end TodaLattice

}// end namespace
#endif
