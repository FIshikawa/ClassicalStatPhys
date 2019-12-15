#ifndef LATTICE_CHAIN_HPP
#define LATTICE_CHAIN_HPP

#include <string>
#include <vector>

namespace lattice {

class Chain{
public:
  static std::string name() { return "1DChain";}
  static int set_num_particles(int Ns){return Ns;}
  static int number_adjacent() {return 2;} 

  Chain(int LatticeSize) : dim_(LatticeSize) {}

  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      int l = numberize(i);
      int temp = 0;
      for(int d = -1; d < 2; d =d+2){ 
      table[l][temp] = numberize(i+d);
      temp += 1;
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

  int target_number(const int i, const int observe){
    int sight_x = 0;
    sight_x = latticize(0,observe);
    return  numberize(sight_x + i);
  }

private:
  int dim_;
};

} //end name space

#endif //LATTICE_CHAIN_HPP
