#ifndef STATDISTANCE_HPP
#define STATDISTANCE_HPP

#include <cmath>
#include <string>
#include <vector>
#include <boost/random.hpp>
#include <boost/math/special_functions/bessel.hpp>

//todo : should change name and implementaion for general case

namespace stat {

class statdistance_equil{
public :
  static std::string name() { return "Statistical Distance from Equilibrium : calc by histogram"; }
  statdistance(unsigned int n_bin, unsigned int N_loop) :  n_bin_(n_bin), N_(N_loop) {KL_ = BC_ = 0.0;}
  
  template <class exactcalc>
  void inputdat(std::vector<double> const& range, std::vector<double> const& hist) const {
    for(int i = 0 ; i < n_bin ; ++i){
    BC_ += std::exp( - 0.5 * (E/T + exact.logZ(T))) / sqrtN_; 
    KL_ += (E/T + exact.logZ(T)) / N_;
  }

  double KLdivergence() const {
    return KL_;
  }
  
  double BhattacharyyaCoefficient() const {
    return BC_;
  }

  double BhattacharyyaDistance() const {
   return -1.0 * std::log(BC_); 
  }

  double HellingerDistance() const {
    return 1 - BC_;
  }

  void reset() const {
    BC_ = KL_ = 0;
  }

private :
  unsigned int n_bin_;
  unsigned int N_;
  double BC_;
  double KL_;
}; //end class

class statdistance_equil{

}; // end class


}
#endif // EXACTCALC_HPP
