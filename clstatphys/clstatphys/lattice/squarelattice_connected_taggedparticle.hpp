#ifndef LATTICE_SQUARELATTICE_CONNECTED_TAGGEDPARTICLE_HPP
#define LATTICE_SQUARELATTICE_CONNECTED_TAGGEDPARTICLE_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class SquareLatticeConnectedTaggedParticle {
public:
  static std::string name() { return "SquareLattice lasso (total must be even)";}
  int set_num_particles(int Ns){return Ns*Ns+1;}
  SquareLatticeConnectedTaggedParticle(unsigned int LatticeSize) : dim_(LatticeSize) {}
  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      for (int j = 0; j < dim_; ++j){
          int l = numberize(i,j);
          table[l][4] = l;
          int temp = 0;
          for(int d = -1; d < 2; d =d+2){ 
            table[l][temp] = numberize(i+d,j);
            temp += 1;
            table[l][temp] = numberize(i,j+d);
            temp += 1;
          }
        }
      }
   int l = numberize((dim_-1)/2,(dim_-1)/2);
   int tagg = dim_ * dim_ ;
   table[l][4] = tagg; 
   table[tagg][0] = l;
   for(int temp = 1; temp < 5 ; ++temp) table[tagg][temp] = tagg;
    }

  int numberize(const int i, const int j){
    int i_t = (i + dim_);
    int j_t = (j + dim_);
    i_t %= dim_ ;
    j_t %= dim_ ;
    return i_t + j_t * dim_;
  }

  int latticize(int d, const int l){ //d is direction(0 ~ ) ex. 0 denotes x, 1 denotes y,~
    int j = l/dim_ ;
    int i = (l - dim_ *j);
    if(d==0) return i;
    else return j;
  }

  int number_adjacent() {return 5;}

private:
  unsigned int dim_;
};

}//end namespace

#endif //LATTICE_SQUARELATTICE_CONNECTED_TAGGEDPARTICLE_HPP
