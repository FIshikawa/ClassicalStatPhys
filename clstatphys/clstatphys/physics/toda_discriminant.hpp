#ifndef TODA_DISCRIMINANT_HPP
#define TODA_DISCRIMINANT_HPP

#include <cmath>
#include <vector>
#include <algorithm>
#include <rokko/lapack/syev.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>

namespace integrable{

class TodaDiscriminant{
public:
  TodaDiscriminant(int num_particles, rokko::dlmatrix L_matrix, std::string boundary="periodic") 
    : num_particles_(num_particles), boundary_(boundary)
  { 
    if(boundary_ == "periodic") matrix_size_ = num_particles_;
    else matrix_size_ = num_particles_ + 1;
    L_matrix_ = L_matrix;
  }
  TodaDiscriminant() = default;

  double operator()(double x) const {
    double phi_now = 0;
    double phi_previous = 1;

    double cos_like = hill_solution(x, phi_now, phi_previous); 

    phi_now = 1;
    phi_previous = 0;
    double sin_like= hill_solution(x, phi_now, phi_previous); 

    return cos_like + sin_like;
  }

  double hill_solution(double x, double phi_now, double phi_previous) const {
    int num_iteration = matrix_size_ - 1;
    if(phi_now > phi_previous) num_iteration += 1; 

    double phi_next;
    int i_next, i_previous;
    for(int i = 0; i < num_iteration; ++i){
      i_next = (matrix_size_ + i + 1) % matrix_size_;
      i_previous = (matrix_size_ + i - 1) % matrix_size_;

      phi_next = -(L_matrix_(i,i_previous) * phi_previous 
                    + L_matrix_(i,i) * phi_now 
                    - x * phi_now) / L_matrix_(i,i_next);
      phi_previous = phi_now;
      phi_now = phi_next;
    }
    return phi_now;
  }

  void total_roots(std::vector<double> & w) const {
    rokko::dlvector w_temp(matrix_size_);

    rokko::dlmatrix L_plus = L_matrix_;

    int info = rokko::lapack::syev('V', 'U', L_plus, w_temp);
    for(int i = 0; i < matrix_size_; ++i) w[i] = w_temp[i];

    rokko::dlmatrix L_minus = L_matrix_;
    L_minus(0,matrix_size_-1) *= -1;
    L_minus(matrix_size_-1,0) *= -1;

    info = rokko::lapack::syev('V', 'U', L_minus, w_temp);
    for(int i = 0; i < matrix_size_; ++i) w[i+matrix_size_] = w_temp[i];
    std::sort(w.begin(), w.end());
  }

  int matrix_size() const { return matrix_size_;}

private:
  int num_particles_, matrix_size_;
  std::string boundary_;
  rokko::dlmatrix L_matrix_;
};//end TodaLattice

}// end namespace
#endif
