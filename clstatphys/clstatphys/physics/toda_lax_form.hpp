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
  TodaLaxForm(int num_particles, double J, double alpha, std::string boundary="periodic") 
    : num_particles_(num_particles), J_(J), alpha_(alpha), boundary_(boundary)
  { 
    if(boundary_ == "periodic") matrix_size_ = num_particles_;
    else matrix_size_ = num_particles_ + 1;
  }
  TodaLaxForm() = default;

  rokko::dlmatrix L_matrix(std::vector<double>& z){

    rokko::dlmatrix L = rokko::dlmatrix::Zero(matrix_size_, matrix_size_);
    L(0, 0) = toda_momentum(z, 0);
    L(0, 1) = toda_position(z, 0);
    L(0, matrix_size_-1) = toda_position(z, matrix_size_-1);
    for(int i = 1; i < matrix_size_-1; ++i){
      L(i, i)   = toda_momentum(z, i); 
      L(i, i+1) = toda_position(z, i);
      L(i, i-1) = toda_position(z, i-1);
    }
    L(matrix_size_-1, matrix_size_-1) = toda_momentum(z, matrix_size_-1);
    L(matrix_size_-1, 0) = toda_position(z, matrix_size_-1);
    L(matrix_size_-1, matrix_size_-2) = toda_position(z, matrix_size_-2);
    return L;
  }

  rokko::dlmatrix B_matrix(std::vector<double>& z){
    int matrix_size_;
    rokko::dlmatrix B(matrix_size_, matrix_size_);
    B(0, 1) = -1.0 * toda_position(z, 0);
    B(0, matrix_size_-1) = toda_position(z, matrix_size_-1);
    for(int i = 1; i < matrix_size_-1; ++i){
      B(i, i+1) = toda_position(z, i);
      B(i, i-1) = -1.0 * toda_position(z, i-1);
    }
    B(matrix_size_-1, 0) = -1.0 * toda_momentum(z, matrix_size_-1);
    B(matrix_size_-1, matrix_size_-2) = toda_momentum(z, matrix_size_-2);
    return B;
  }

  // The displacements are defined as r_n = u_n - u_n+1
  // and the Toda potential is defined as J * exp(a * r_n ) -  a * r_n -1

  double toda_position(std::vector<double>& z, int const& n){
    const double *x = &z[0];
    int now_index = n % matrix_size_;
    int target_index = (matrix_size_ + now_index - 1) % matrix_size_;
    double x_now, x_target;

    if(now_index == num_particles_) x_now = 0.0;
    else x_now = x[now_index];

    if(target_index == num_particles_) x_target = 0.0;
    else x_target = x[target_index];

    return 0.5 * std::sqrt(J_) * alpha_ * std::exp(0.5*alpha_*(x_target-x_now));
  }

  double toda_momentum(std::vector<double>& z, int const& n){
    int now_index = n % matrix_size_;
    double p_now;
    const double *p = &z[num_particles_];

    if(now_index == num_particles_) p_now = 0.0;
    else p_now = p[now_index];

    return 0.5 * alpha_ * p_now;
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
  int num_particles_, matrix_size_;
  double J_, alpha_;
  std::string boundary_;
};//end TodaLattice

}// end namespace
#endif
