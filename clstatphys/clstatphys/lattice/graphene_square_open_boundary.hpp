#ifndef LATTICE_GRAPHENE_SQUARE_OPEN_BOUNDARY_HPP
#define LATTICE_GRAPHENE_SQUARE_OPEN_BOUNDARY_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class GrapheneSquareOpenBoundary{
public:
  static std::string name() { return "Graphene Square Open Boundary : input Ns is the number of hexagonal :(2*Ns) * (2*Ns+1) square";}
  static int set_num_particles(int Ns){return 2*Ns*(2*Ns+1);}
  GrapheneSquareOpenBoundary(int LatticeSize) : dim_(2*LatticeSize) {}
  void create_table(std::vector<std::vector < int > >& table){
  //for bulk
    for (int i = 0; i < dim_ + 1; ++i){
      for (int j = 0; j < dim_; ++j){
        int l = numberize(i,j);
        if(i == 0) table[l][0] = l;
        else table[l][0] = numberize(i-1,j);
        if(i == dim_) table[l][1] = l;
        else table[l][1] = numberize(i+1,j);
        if((i+j)%2 == 0){
          if(j == dim_ -1 ) table[l][2] = l;
          else table[l][2] = numberize(i,j+1);
        }
        else {
          if(j == 0) table[l][2] = l; 
          else table[l][2] = numberize(i,j-1);
        }
      }
    }
  }

  int numberize(const int i, const int j){
    int i_t = (i + dim_+1);
    int j_t = (j + dim_);
    i_t %= (dim_+1) ;
    j_t %= dim_ ;
    return i_t + j_t * (dim_+1);
  }

  int latticize(int d, const int l){ //d is direction(0 ~ ) ex. 0 denotes x, 1 denotes y,~
    int j = l/(dim_+1);
    int i = (l - (dim_+1) *j);
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

  int number_adjacent() {return 3;}

private:
  int dim_;
};//Regular cube end

} //end namespace

#endif //LLATTICE_GRAPHENE_SQUARE_OPEN_BOUNDARY_HPP
