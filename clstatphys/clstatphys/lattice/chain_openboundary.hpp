#ifndef LATTICE_CHAIN_OPENBOUNDARY_HPP
#define LATTICE_CHAIN_OPENBOUNDARY_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class ChainOpenBoundary{
public:
  static std::string name() { return "1DChain";}
  static int set_num_particles(int Ns){return Ns;}
  ChainOpenBoundary(unsigned int LatticeSize) : dim_(LatticeSize) {}
  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 1; i < dim_-1; ++i){
      int l = numberize(i);
      int temp = 0;
      for(int d = -1; d < 2; d =d+2){ 
      table[l][temp] = numberize(i+d);
      temp += 1;
      }
    }
    int l = numberize(0);
    table[l][0] = numberize(0 + 1);
    table[l][1] = numberize(0);
    l = numberize(dim_-1);
    table[l][0] = numberize(dim_ -1 - 1);
    table[l][1] = numberize(dim_ -1);

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

  int number_adjacent() {return 2;} 

private:
  unsigned int dim_;
};

}//end namespace

#endif //LATTICE_CHAIN_OPENBOUNDARY_HPP 
