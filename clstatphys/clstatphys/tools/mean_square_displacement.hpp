#ifndef TOOLS_MEAN_SQUARE_DISPLACEMENT_HPP
#define TOOLS_MEAN_SQUARE_DISPLACEMENT_HPP

#include <string>
#include <vector>
#include <cmath>

namespace correlation{

class MeanSquareDisplacement{
public:
  MeanSquareDisplacement(int dim, int Nl) : dim_(dim), Nl_(Nl),mean_(dim,0.0),error_(dim,0.0) {}
  void calc (std::vector<double>const& x) const {
    double N = (double)dim_;
    for (int i = 0; i< dim_ ; ++i){
      mean_[i] += std::pow(x[i] - x[0] , 2.0)/Nl_ ;
      error_[i] += std::pow(x[i] - x[0] , 4.0)/Nl_ ;
    }
  }

  std::vector<double> result(){ return mean_; }
  std::vector<double> error(){ return error_; }

private:
  int dim_;
  int Nl_;
  std::vector<double> error_, mean_;
};

}//end namespace 

#endif //TOOLS_MEAN_SQUARE_DISPLACEMENT_HPP
