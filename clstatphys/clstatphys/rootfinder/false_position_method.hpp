#ifndef OPTIMIZATION_FALSE_POSITION_METHOD_HPP
#define OPTIMIZATION_FALSE_POSITION_METHOD_HPP

#include <algorithm>
#include<cmath>
#include<limits>

namespace optimization {

class FalsePositionMethod{
public:
  static std::string name(){return "Optimize Method : False Position Method";}

  template <typename Func>
  int find_zero(Func& f, double x0, double x1,
   double prec = 2.0 * std::numeric_limits<double>::epsilon()) {     
    if (x0 > x1) std::swap(x0, x1);
    double y0 = f(x0); 
    int counter = 1;
    if(std::abs(y0) < prec) {
      zero_ = x0;
      return counter;
    }
    ++counter;
    double y1 = f(x1);
    if(std::abs(y1) < prec){
      zero_ = x1;
      return counter;
    }
    if(y0 * y1 > 0) return -1;

    while(((x1-x0) / std::max(std::abs(x0), std::abs(x1))) > prec &&
          std::min(std::abs(y1), std::abs(y0)) > prec ){
    ++counter;
      double xm = (x0 * y1 - x1 * y0) / (y1 - y0);
      double ym = f(xm);
      if(y0 * ym > 0) {
        x0 = xm; y0 = ym;
      }else if(y1 * ym > 0){
        x1 = xm; y1 = ym;
      }
      else{
        y1 = ym; y0 = ym;
        x1 = xm; x0 = xm;
      }
    }
    zero_ = x0;
    return counter;
  }

  double zero(){ return zero_;}

private:
  double zero_;
}; // end false position

} // end namespace

#endif //OPTIMIZATION_FALSE_POSITION_METHOD_HPP
