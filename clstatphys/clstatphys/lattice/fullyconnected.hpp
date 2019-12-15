#ifndef LATTICE_FULLYCONNECTED_HPP
#define LATTICE_FULLYCONNECTED_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class FullyConnected {
public:
  static std::string name() { return "Fully_connected";}
  static int set_num_particles(int Ns){return Ns;}
  FullyConnected(unsigned int LatticeSize) : dim_(LatticeSize) {}
  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      for (int j = 0; j < dim_-1; ++j){
        table[i][j] = numberize(i+j+1);
      }
    }
  }

  int numberize(const int i){
    int i_t = (i + dim_);
    i_t %= dim_ ;
    return i_t ;
  }

  int latticize(int d, const int l){ //d is direction(0 ~ ) ex. 0 denotes x, 1 denotes y,~
    int i = l;
    return i;
  }

  int number_adjacent() {return dim_-1;}

private:
  unsigned int dim_;
};

} //end namespace

#endif //LATTICE_FULLYCONNECTED
