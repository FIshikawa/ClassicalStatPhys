#include<iostream>
#include<cmath>
#include<fstream>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<boost/random.hpp>
#include<boost/random/uniform_int_distribution.hpp>
#include<boost/format.hpp>
#include<map>
#include <rokko/blas.hpp>
#include <rokko/lapack/syev.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>


#include "dataput.hpp"
#include "hamiltonian.hpp"
#include "integrator.hpp"
#include "correlation.hpp"
#include "Lattice.hpp"
#include "calcerror.hpp"
#include "exefftw.hpp"

//typedef Lattice::RegularCubeLasso lattice_t;
//typedef Lattice::BetheLattice lattice_t;
typedef Lattice::BetheLatticelasso lattice_t;
//typedef Lattice::ChainLasso lattice_t;
//typedef Lattice::SquareLatticeLasso lattice_t;
//typedef Lattice::RegularCube lattice_t;
//typedef Lattice::SquareLattice lattice_t;

int main(int argc, char **argv){
 unsigned int Ns = 2 ;
 unsigned int num_particles;
 unsigned int N_adj = 3;// number of adjacent spins

 if (argc > 1) Ns =      boost::lexical_cast<unsigned int>(argv[1]);
 if (argc > 2) N_adj =      boost::lexical_cast<unsigned int>(argv[2]);

 lattice_t lattice(Ns,N_adj);
 //num_particles = lattice_t::set_num_particles(Ns);
 num_particles = lattice.set_num_particles(Ns);
 N_adj = lattice.number_adjacent();
 int counter ; //montecalro counter 

 std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
 lattice.create_table(pair_table);

 std::cout << "Lattice : " << lattice_t::name() << std::endl
           << "Ns is : " << Ns << std::endl
           << "N_adj is : " << boost::lexical_cast<unsigned int>(argv[2]) << std::endl
           << "num_particles is : " << num_particles << std::endl; 

 unsigned int n = num_particles;
 rokko::dlmatrix a(n, n);
 a = rokko::dlmatrix::Identity(n, n);
 for (int i = 0; i < num_particles; ++i) {
   for(int j = 0; j < N_adj; ++j){ 
     int p = pair_table[i][j];
     a(i, p) = -1;
   }
 }

 if(  num_particles < 30 ) std::cout << "Matrix A : " << std::endl << a << std::endl;
 
 return 0;
}

