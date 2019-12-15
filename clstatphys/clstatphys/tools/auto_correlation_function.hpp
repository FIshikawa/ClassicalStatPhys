#ifndef AUTO_CORRELATION_FUNCTION_HPP
#define AUTO_CORRELATION_FUNCTION_HPP

#include <string>
#include <vector>
#include <cmath>

namespace correlation{

class AutoCorrelationFunction{
public:
  AutoCorrelationFunction(int dim=1, int Nl=1) : dim_(dim), Nl_(Nl), counter_(0), correlation_(dim,0.0), mean_(dim,0.0){}
  void initialize(int dim, int Nl){
    counter_ = 0;
    dim_ = dim;
    Nl_ = Nl;
    correlation_.resize(dim);
    mean_.resize(dim);
    for(int i = 0 ; i < dim; ++i){
      correlation_[i] = 0.0;
      mean_[i] = 0.0;
    }
  }

  void operator<< (const double value) {
    if(counter_ == 0 || counter_ > dim_ - 1 ){
      value_init_ = value;
      counter_ = 0;
    }
    mean_[counter_] += value / Nl_;
    correlation_[counter_] += value * value_init_ / Nl_;
    counter_++ ;
  }

  std::vector<double> result(){
    std::vector<double> acf_t(dim_, 0.0);
    for(int i = 0; i < dim_ ; ++i) acf_t[i] = correlation_[i] - mean_[counter_] * mean_[0];
    return acf_t;
  }

  void calc( std::vector<double>const& z, std::vector<double>& ACF) const {
    double N = (double)dim_ ; 
    for (int i = 0; i < dim_ ; ++i) ACF[i] += z[i]*z[0]/Nl_ ;
  }

private:
  int dim_, Nl_, counter_;
  double value_init_;
  std::vector<double> correlation_, mean_;
};

} //end namespace

#endif // AUTO_CORRELATION_FUNCTION_HPP

