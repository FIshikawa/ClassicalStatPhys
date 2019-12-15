#ifndef LATTICE_SQUARELATTICE_HPP
#define LATTICE_SQUARELATTICE_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class SquareLattice {
public:
  static std::string name() { return "SquareLattice";}
  static int set_num_particles(int Ns){return Ns*Ns;}
  SquareLattice(int LatticeSize) : dim_(LatticeSize) {}
  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_; ++i){
      for (int j = 0; j < dim_; ++j){
          int l = numberize(i,j);
          int temp = 0;
          for(int d = -1; d < 2; d =d+2){ 
            table[l][temp] = numberize(i+d,j);
            temp += 1;
            table[l][temp] = numberize(i,j+d);
            temp += 1;
          }
        }
      }
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

  int target_number(const int i, const int observe){
    int sight_x = 0;
    int sight_y = 0;
    sight_x = latticize(0,observe);
    sight_y = latticize(1,observe);
    return  numberize(sight_x + i, sight_y);
  }

  int number_adjacent() {return 4;}

private:
  int dim_;
};//Regular cube end

} //end namespace

#endif //LATTICE_SQUARELATTICE_HPP
