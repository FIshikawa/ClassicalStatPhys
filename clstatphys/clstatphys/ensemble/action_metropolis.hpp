#ifndef ENSEMBLE_METROPOLIS_ACTION_TODA_HPP
#define ENSEMBLE_METROPOLIS_ACTION_TODA_HPP

#include <cmath>
#include <string>
#include <random>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/physics/toda_action_variables.hpp>

namespace ensemble{

class MetropolisActionToda{
public:
  static std::string name() { return "Metropolis for sampling action variables of Toda"; }
  MetropolisActionToda(int num, int num_iteration, std::vector<double> betas, 
      double J, double alpha, double dp, double dx) 
    : num_iteration(num_iteration),num_(num),betas_(betas),
      J_(J), alpha_(alpha), dp_(dp), dx_(dx) {} 

  template <class Rand>
  void montecarlo(std::vector<double>& z, int& counter, Rand & mt) const {
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::uniform_int_distribution<> pick_index(0,num_-1);
    std::uniform_real_distribution<> diff(-1.0,1.0);
    std::uniform_real_distribution<> realdist(0.0,1.0);
    double *p = &z[num_]; 
    double *x = &z[0]; 
    int tagged = pick_index(mt);

    double past_p = p[tagged];
    double past_x = x[tagged];

    double new_p = p[tagged] + dp_ * diff(mt);
    double new_x = x[tagged] + dx_ * diff(mt);

    integrable::TodaLaxForm toda_lax_form(num_,J_,alpha_,"periodic");
    std::vector<double> action_past(num_-1);
    rokko::dlmatrix L = toda_lax_form.L_matrix(z);
    integrable::TodaDiscriminant discriminant(num_, L, "periodic");
    TodaActionVariables(action_past, discriminant, num_iteration);

    std::vector<double> action_new(num_-1);
    x[tagged] = new_x;
    p[tagged] = new_p;
    L = toda_lax_form.L_matrix(z);
    discriminant = integrable::TodaDiscriminant(num_, L, "periodic");
    TodaActionVariables(action_new, discriminant, num_iteration);

    double log_likelyhood = 0.0;
    for(int i = 0; i < num_-1; ++i) 
      log_likelyhood -= betas_[i] * (action_new[i] - action_past[i]);
    
    // c.f. dP = std::exp(-1.0/T_ * new_dE)/std::exp(-1.0/T_ * past_dE);
    double dP = std::exp(log_likelyhood);
    double P_accept = realdist(mt);
   
    if(dP > P_accept){
      counter += 1;
    }
    else{
      x[tagged] = past_x;
      p[tagged] = past_p;
    }
  }

private:
  int num_, num_iteration;
  double J_, alpha_, dp_, dx_;
  std::vector<double> betas_;
}; 

} //end namespace

#endif 

