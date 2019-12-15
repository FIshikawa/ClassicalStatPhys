#ifndef SPIN_PLOT_HPP
#define SPIN_PLOT_HPP

#include<vector>
#include<cmath>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<unordered_map>

namespace tools{

 template <class Lattice>
 void SpinPlotXY(const std::string& filename, unsigned int Ns_ObserveSystem, std::vector<double> z, Lattice  lattice) const {
  std::ofstream ofs(filename);
  const double *x = &z[0];
  
  #ifdef SQUARELATTICE_HPP 
    auto lattice_number = [&lattice](int i, int j){return lattice.numberize(i,j)}
  #else
    auto lattice_number = [&lattice](int i, int j){return lattice.numberize(i,j,Ns_/2)}
  #endif 

  }
  for(int i= Ns_/2-Ns_ObserveSystem/2; i < Ns_/2+Ns_ObserveSystem/2; ++i){
    for(int j=Ns_/2-Ns_ObserveSystem/2; j < Ns_/2+Ns_ObserveSystem/2; ++j){
      int l = lattice_number(i,j);
      ofs << i << " " << j << " " << std::cos(x[l]) << " " << std::sin(x[l]) << " " << fmod(x[l]+10000*M_PI,2.0*M_PI)<< std::endl;
     }
   }
 }

 template <class Lattice>
 void SpinPlotHeisenberg(const std::string& filename, unsigned int Ns_ObserveSystem, std::vector<double> z, Lattice  lattice) const {
  std::ofstream ofs(filename.c_str());
  const double *x = &z[0];
  const double *y = &z[Ns_];
  for(int i= Ns_/2-Ns_ObserveSystem/2; i < Ns_/2+Ns_ObserveSystem/2; ++i){
    for(int j=Ns_/2-Ns_ObserveSystem/2; j < Ns_/2+Ns_ObserveSystem/2; ++j){
      for(int k=Ns_/2-Ns_ObserveSystem/2; k < Ns_/2+Ns_ObserveSystem/2; ++k){
        int l = lattice.numberize(i,j,k);
        ofs << i << " " << j << " " << k << " " << std::cos(x[l])*std::sin(y[l]) << " " << std::sin(x[l])*std::sin(y[l]) << " " << std::cos(y[l]) << " " << fmod(y[l]+10000*M_PI,M_PI)/M_PI  << std::endl;
      }
     }
   }
 }
};//end dataputForRegularCube

} // end namespace

#endif //SPIN_PLOT_HPP
