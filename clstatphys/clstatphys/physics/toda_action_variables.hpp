#ifndef TODA_ACTION_VARIABLES_HPP
#define TODA_ACTION_VARIABLES_HPP

#include <cmath>
#include <vector>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/tools/simpson.hpp>

void TodaActionVariables(std::vector<double> & action_variables, 
                         integrable::TodaDiscriminant const & discriminant, 
                         int num_iterations=16){

  int matrix_size = discriminant.matrix_size();
  std::vector<double> roots(2*matrix_size);
  discriminant.total_roots(roots);

  auto func = [&discriminant](double x){ 
    if(discriminant(x) > 2) return std::acosh(0.5*discriminant(x));
    else if(discriminant(x) < -2) return std::acosh(-0.5*discriminant(x));
    else return 0.0;
  };

  for(int i = 0; i < matrix_size - 1; ++i){
    action_variables[i] = standards::simpson_1d(func, roots[2*i+1], roots[2*i+2], num_iterations);
  }
}

#endif
