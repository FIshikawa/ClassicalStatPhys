#ifndef TODA_CONSERVATION_FIELDS_HPP
#define TODA_CONSERVATION_FIELDS_HPP

#include <cmath>
#include <vector>

namespace integrable{

  // ref: 
  // Integrals of the Toda lattice
  // M. Hénon
  // Phys. Rev. B 9, 1921 – Published 15 February 1974

class TodaConservationFields{
public:
  TodaConservationFields(int num_particles, double J, double alpha, std::vector<std::vector<int > > table, int N_adj) : num_particles_(num_particles), J_(J), table_(table),Nd_(N_adj),alpha_(alpha){}

  double operator()(std::vector<double>& z, int site, int kind){
    int num_particles = z.size() / 2;
    const double *x = &z[0];
    const double *p = &z[num_particles];

    int next = table_[site][1];
    double dx = x[next] - x[site];
    double V_toda_site = J_ * std::exp( alpha_ * dx);
    dx = x[table_[next][1]] - x[next];
    double V_toda_next = J_ * std::exp( alpha_ * dx);

    double value = 0.0;
    if(kind == 1) value = p[site];
    else if(kind == 2) value = power(p[site],2) * 0.5 + V_toda_site;
    else if(kind == 3) value = power(p[site],3) / 3 + (p[site] + p[next]) * V_toda_site; 
    else if(kind == 4) value = power(p[site],4) / 4 
                                + (power(p[site],2) + p[site] * p[next] + power(p[next],2)) * V_toda_site 
                                  + power(V_toda_site,2) * 0.5 + V_toda_site * V_toda_next;
    else if(kind == 5) value = power(p[site],5) / 5 
                                + (power(p[site],3) + power(p[site],2) * p[next] + p[site] * power(p[next],2) + power(p[next],3) )* V_toda_site
                                  + (p[site] + p[next]) * power(V_toda_site,2) + (p[site] + 2*p[next] + p[table_[next][1]] ) * V_toda_site * V_toda_next;

    else std::cerr << "kind oveer 5, they are not defined" << std::endl; 

    value *= std::pow(0.5*alpha_,kind)*kind;
    return value;
  }

  double power(double a, int b){
    double total = 1.0;
    for(int i = 0; i < b; ++i) total *= a; 
    return total;
  }

private:
  int Nd_;
  int num_particles_;
  double J_;
  double alpha_;
  std::vector<std::vector<int> > table_;
};//end TodaLattice

}// end namespace

#endif //NORMALMODE_ENERGY_HPP
