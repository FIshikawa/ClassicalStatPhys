#ifndef LATTICE_REGULARCUBE_CONNECTED_TAGGEDPARTICLE_HPP
#define LATTICE_REGULARCUBE_CONNECTED_TAGGEDPARTICLE_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class RegularCubeConnectedTaggedParticle{
public:
  static std::string name() { return "Regular Cube lasso (total must be even)";}
  static int set_num_particles(int Ns){return Ns*Ns*Ns+1;}
  RegularCubeConnectedTaggedParticle(unsigned int LatticeSize) : dim_(LatticeSize) {}

  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      for (int j = 0; j < dim_; ++j){
        for (int k = 0; k < dim_; ++k){
          int l = numberize(i,j,k);
          int temp = 0;
          table[l][6] = l;
          for(int d = -1; d < 2; d =d+2){ 
            table[l][temp] = numberize(i+d,j,k);
            temp += 1;
            table[l][temp] = numberize(i,j+d,k);
            temp += 1;
            table[l][temp] = numberize(i,j,k+d);
            temp += 1;
          }
        }
      }
    }
   int tagg = dim_ * dim_ * dim_ ;
   int l = numberize((dim_-1)/2,(dim_-1)/2,(dim_-1)/2);
   table[l][6] = tagg; 
   table[tagg][0] = l;
   for(int temp = 1; temp < 7 ; ++temp) table[tagg][temp] = tagg;
   
  }  

  int numberize(const int i, const int j,const int k){
    int i_t = (i + dim_);
    int j_t = (j + dim_);
    int k_t = (k + dim_);
    i_t %= dim_ ;
    j_t %= dim_ ;
    k_t %= dim_ ;
    return i_t + j_t * dim_ + k_t * dim_ * dim_;
  }

  int latticize(int d, const int l){ //d is direction(0 ~ ) ex. 0 denotes x, 1 denotes y,~
    int k = l / (dim_*dim_) ;
    int j = (l - dim_ * dim_* k)/dim_ ;
    int i = (l - dim_ * dim_ * k - dim_ *j);
    if(d==0) return i;
    else if (d==1) return j;
    else return k;
  }

  int number_adjacent() {return 7;}

private:
  unsigned int dim_;
};

} //end namespace

#endif //LATTICE_REGULARCUBE_CONNECTED_TAGGEDPARTICLE_HPP
