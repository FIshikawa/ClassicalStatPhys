#ifndef STATITICAL_DIVERGENCE_HPP
#define STATITICAL_DIVERGENCE_HPP

#include <vector>
#include <cmath>

namespace tools{

class StatisticalDivergence{
public:
  void reset(){
    l2_norm_ = l1_norm_ = hellinger_distance_ = jensen_shannon_divergence_ = 0;
  }
  void calc(std::vector<double> const& p, std::vector<double> const& q){
    int data_size = p.size();
    double product,fraction,middle_point,derivative;
    double l2_norm,l1_norm,b_coefficinet,kl_p,kl_q;
    reset();
    l2_norm = l1_norm = b_coefficinet = kl_p = kl_q = 0.0;
    for(int i = 0; i < data_size; ++i){
      derivative = (p[i] - q[i]) * (p[i] - q[i]);
      l2_norm += derivative;
      derivative = std::sqrt(derivative);
      l1_norm += derivative;
      product = (p[i]*q[i] < 0) ? 0.0 : p[i]*q[i];
      b_coefficinet += std::sqrt(product);
      middle_point = 0.5 * (p[i] + q[i]);
      kl_p += (p[i] > 0) ? p[i] * std::log(p[i] / middle_point) : 0.0;
      kl_q += (q[i] > 0) ? q[i] * std::log(q[i] / middle_point) : 0.0;
    }
    l2_norm_ = std::sqrt(l2_norm);
    l1_norm_ = l1_norm;
    hellinger_distance_ = std::sqrt(1.0 - b_coefficinet); 
    jensen_shannon_divergence_ = 0.5 * (kl_p + kl_q);
  }
  
  double l2_norm(){ return l2_norm_;}
  double l1_norm(){ return l1_norm_;}
  double hellinger_distance() { return hellinger_distance_;}
  double jensen_shannon_divergence() {return jensen_shannon_divergence_;}
private:
  double l2_norm_,l1_norm_,hellinger_distance_,jensen_shannon_divergence_;
}; // end StatisticalDivergence

} // namespace end

#endif //STATITICAL_DIVERGENCE_HPP
