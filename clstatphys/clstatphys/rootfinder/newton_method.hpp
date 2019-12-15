#ifndef OPTIMIZATION_NEWTON_METHOD_HPP
#define OPTIMIZATION_NEWTON_METHOD_HPP

#include<cmath>
#include<limits>

namespace optimization {

class NewtonMethod{
public:
  static std::string name(){return "Optimize Method : Newton Method";}
  template <typename Func1, typename Func2>
  int find_zero(Func1& f, Func2& df,  double x,
   double prec = 2.0 * std::numeric_limits<double>::epsilon()) {     
    double y = f(x); 
    int counter = 1;
    if(std::abs(y) < prec) {
      zero_ = x;
      return counter;
    }
    double dy = df(x);
    if(std::abs(dy) < prec) return -1;

    while(std::abs(y) > prec && std::abs(dy) > prec) {
      counter += 1;
      x = x -  y / dy;
      y = f(x);
      dy = df(x);
    }
    if(std::abs(dy) < prec) return -1;
    else{
      zero_ = x;
      return counter;
    }
  } 

  double zero(){ return zero_;}

protected:
  double zero_;
}; // end neweton

}
#endif // OPTIMIZATION_NEWTON_METHOD_HPP
