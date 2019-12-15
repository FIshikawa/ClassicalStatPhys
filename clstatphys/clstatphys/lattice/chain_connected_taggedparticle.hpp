#ifndef LATTICE_CHAIN_CONNECTED_TAGGEDPARTICLE_HPP
#define LATTICE_CHAIN_CONNECTED_TAGGEDPARTICLE_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class ChainConnectedTaggedParticle{
public:
  static std::string name() { return "1DChain Lasso (total num_particle must be even)";}
  static int set_num_particles(int Ns){return Ns+1;}
  ChainConnectedTaggedParticle(unsigned int LatticeSize) : dim_(LatticeSize) {}
  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      int l = numberize(i);
      int temp = 0;
      table[l][2] = l;
      for(int d = -1; d < 2; d =d+2){ 
      table[l][temp] = numberize(i+d);
      temp += 1;
      }
    }
   int tagg = dim_;
   int l = numberize((dim_-1)/2);
   table[l][2] = tagg; 
   table[tagg][0] = l;
   table[tagg][1] = tagg;
   table[tagg][2] = tagg;
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

 int number_adjacent() {return 3;} 
private:
  unsigned int dim_;
};

} //end namespace

#endif //LATTICE_CHAIN_CONNECTED_TAGGEDPARTICLE_HPP 
