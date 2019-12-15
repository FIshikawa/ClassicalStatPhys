#ifndef LATTICE_REGULAR_CUBE_HPP
#define LATTICE_REGULAR_CUBE_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class RegularCube {
public:
  static std::string name() { return "Regular Cube";}
  static int set_num_particles(int Ns){return Ns*Ns*Ns;}
  RegularCube(int LatticeSize) : dim_(LatticeSize) {}

  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      for (int j = 0; j < dim_; ++j){
        for (int k = 0; k < dim_; ++k){
          int l = numberize(i,j,k);
          int temp = 0;
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

  int target_number(const int i, const int observe){
    int sight_x = 0;
    int sight_y = 0;
    int sight_z = 0;
    sight_x = latticize(0,observe);
    sight_y = latticize(1,observe);
    sight_z = latticize(2,observe);
    return numberize(sight_x + i, sight_y, sight_z);
  }

 int number_adjacent() {return 6;}

private:
  int dim_;
};

} //end namespace

#endif //LATTICE_REGULAR_CUBE_HPP 
