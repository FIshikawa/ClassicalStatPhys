#ifndef LATTICE_BETHELATTICE_HPP
#define LATTICE_BETHELATTICE_HPP

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace lattice {

class BetheLattice {
public:
  static std::string name() { return "Bethe Lattice";}
  BetheLattice(unsigned int LatticeSize , unsigned int NeighboringNumber) : num_stage(LatticeSize), n_adj(NeighboringNumber) {}
  void create_table(std::vector<std::vector < int > >& table){
    for(int stage = 0; stage < num_stage+1 ; ++stage){
      int num_each_stage = calc_num_each_stage(stage);
      for(int tagg = 0; tagg < num_each_stage ; ++tagg){
        int number  = numberize(stage,tagg);
        if(stage == num_stage){
          for(int i = 0; i < n_adj - 1; ++i){
            table[number][i] = number;
          }
        }
        else if(stage == 0){
          for(int i = 0; i < n_adj; ++i){
            table[number][i] = i+1;
            table[i+1][n_adj-1] = number;
          }
        }
        else{
          for(int i = 0; i < n_adj - 1; ++i){
            int pair_number = numberize(stage+1, tagg*(n_adj-1)+i);
            table[number][i] = pair_number;
            table[pair_number][n_adj-1] = number;
          }
        }
      }
    }
  }

  int numberize(const int stage, const int tagg){// j <= std::pow(n_adj-1,i-1)* n_adj - 1
    if(stage ==0) return 0;
    else return BetheCalc(stage) + 1 + tagg;
  }

  int latticize(int d, const int l){ //d selects stage_number or tagg_number (0, 1) ex. 0 denotes stage, 1 denotes tagg_number
    int stage = num_stage;
    int tagg = 1;
    int result = 0;
    while(tagg != 0){
      tagg = l % BetheCalc(stage);
      --stage;
      if(stage < 2) break;
    }
    if(stage < 2){
      if(l == 0){
        if(d==0) result = 0;
        else result = 0;
      }
      else{
        if(d==0)result = stage;
        else result= l;
      }
    }
   return result;
  }

  int set_num_particles(int Ns){
    if(Ns == 0) return  1;
    else return BetheCalc(Ns+1) + 1;
  }

  int BetheCalc(int stage ){ //N_t => 0 
    int num_t = 0;
    for(int k = 1 ; k <= stage-1; ++k){
     num_t += calc_num_each_stage(k);
    }
    return num_t;
  }

  int calc_num_each_stage(int stage){
    double temp;
    double stage_t = boost::lexical_cast<double>(stage);
    double n_adj_t = boost::lexical_cast<double>(n_adj);
    if(stage == 0) temp = 1;
    else temp = n_adj_t * std::pow(n_adj_t-1, stage_t-1) ;
    return boost::lexical_cast<int>(temp);
 }

  int number_adjacent() {return n_adj;} 

private:
  unsigned int num_stage;
  unsigned int n_adj;

};//Bethe lattice end

} // end namespace

#endif //LATTICE_BETHELATTICE_HPP
