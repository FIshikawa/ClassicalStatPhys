#ifndef TODA_ACTION_VARIABLES_GK_HPP
#define TODA_ACTION_VARIABLES_GK_HPP

#include <cmath>
#include <vector>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace boost::math::quadrature;

void TodaActionVariablesGK(std::vector<double> & action_variables, 
                         integrable::TodaDiscriminant const & discriminant, 
                         int max_depth=5,double relative_error=1e-5){

  int matrix_size = discriminant.matrix_size();
  std::vector<double> roots(2*matrix_size);
  discriminant.total_roots(roots);

  auto func = [&discriminant](double x){ 
    if(discriminant(x) > 2) return std::acosh(0.5*discriminant(x));
    else if(discriminant(x) < -2) return std::acosh(-0.5*discriminant(x));
    else return 0.0;
  };

  for(int i = 0; i < matrix_size - 1; ++i){
    double error;
    if(roots[2*i+1] - roots[2*i+2] == 0)
      action_variables[i] = 0.0;
    else
      action_variables[i] = 2.0 / M_PI 
                          * gauss_kronrod<double, 15>::integrate(
                                                               func, 
                                                               roots[2*i+1], 
                                                               roots[2*i+2], 
                                                               max_depth, 
                                                               relative_error, 
                                                               &error
                                                               );
  }
}

#endif
