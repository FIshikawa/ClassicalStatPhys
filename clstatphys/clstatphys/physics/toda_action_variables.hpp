#ifndef TODA_ACTION_VARIABLES_HPP
#define TODA_ACTION_VARIABLES_HPP

#include <cmath>
#include <vector>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/tools/simpson.hpp>

void TodaActionVariables(std::vector<double> & action_variables, 
                         integrable::TodaDiscriminant & discriminant, 
                         int num_iterations=16,
                         double epsilon=1e-12){

  int matrix_size = discriminant.matrix_size();
  std::vector<double> roots(2*matrix_size);
  discriminant.total_roots(roots);

  auto func = [&discriminant, &epsilon](double x){ 
    if(discriminant(x) > 0) return std::acosh(0.5*discriminant(x)+epsilon);
    else return std::acosh(-0.5*discriminant(x)+epsilon);
  };
  std::cout << std::scientific;
  std::cout << "0.5 * discriminant : " << 0.5 * discriminant(roots[1]) << std::endl
            << "arcosh at above : " << std::acosh(0.5 * discriminant(roots[1])) << std::endl
            << "arcosh at 1 " << std::acosh(1.0) << std::endl
            << "roots 1 : " << func(roots[1]) << std::endl;

  for(int i = 0; i < matrix_size - 1; ++i){
    action_variables[i] = standards::simpson_1d(func, roots[2*i+1], roots[2*i+2], num_iterations);
  }
}

#endif
