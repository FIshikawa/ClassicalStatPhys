//square lattice clXYmodel monte carlp
#include<iostream>
#include<cmath>
#include<fstream>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<boost/random.hpp>
#include<boost/random/uniform_int_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

#include "integrator.hpp"
#include "initializer.hpp"
#include "dataput.hpp"
#include "hamiltonian.hpp"
#include "Lattice.hpp"
#include "correlation.hpp"

//typedef integrator::velver integrator_t;
//typedef initializer::equilclHeisenberg initializer_t;
typedef initializer::equilclHeisenbergWolf initializer_t;
//typedef initializer::equilclXYmodel initializer_t;
//typedef initializer::equilclXYmodelWolf initializer_t;
typedef Lattice::RegularCube lattice_t;
//typedef Lattice::SquareLattice lattice_t;
//typedef dataput::dataputForSquareLattice dataput_t;
typedef dataput::dataputForRegularCube_Heisenberg dataput_t;
typedef hamiltonian::classicalHeisenberg hamiltonian_t;
//typedef dataput::dataputForRegularCube dataput_t;
//typedef hamiltonian::classicalXY hamiltonian_t;

int main(int argc, char **argv){
  double J = 1.0;
  int N_relax = 10;
  double T = 0.5; //temperature
  int Ns = 3; // lenght of whole system
  int N_thermalize = 3; // number of thermalization
  int num_particles; //nummer of particles
  double t = 10; // Whole time
  double dt; //size of interval dt = t/N_time
  int N_time = 100; //numer of time step
  int Ns_ObserveSystem = 3;
  std::string condition_dat = "condi.dat";
  std::string result_dir("./test/mapdata/");

 //parameter set 
  if (argc > 1) result_dir =       boost::lexical_cast<std::string>(argv[1]);
  if (argc > 2) N_thermalize =     boost::lexical_cast<int>(argv[2]);
  if (argc > 3) Ns =               boost::lexical_cast<int>(argv[3]);
  if (argc > 4) T =                boost::lexical_cast<double>(argv[4]);
  if (argc > 5) N_time =           boost::lexical_cast<int>(argv[5]);
  if (argc > 6) t =                boost::lexical_cast<double>(argv[6]);
  if (argc > 7) N_relax    =       boost::lexical_cast<int>(argv[7]);
  if (argc > 8) Ns_ObserveSystem = boost::lexical_cast<int>(argv[8]);

 // parameter auto set
  num_particles = lattice_t::set_num_particles(Ns);
  dt = t/(double)N_time;
  condition_dat = result_dir + condition_dat;

 //parameter set 
  std::vector<std::string> experimental_condition(4);
  experimental_condition[0] = "classical Spin models Demonstration with visualization";
  experimental_condition[1] = hamiltonian_t::name();
  experimental_condition[2] = initializer_t::name();
  experimental_condition[3] = lattice_t::name();
  dataput_t::startput(condition_dat.c_str(), experimental_condition);

 //condition declare
  std::ofstream ofs;
  ofs.open(condition_dat.c_str(),std::ios::app);
  if(!ofs)std::cerr<< condition_dat <<"opening failed"<<std::endl;
  ofs  << "Coupling constant: J is " << J << std::endl
       << "Expected Temperture: T = " << T << std::endl
       << "Reduced Coupling: J/T = " << J/T << std::endl
       << "Number of patricles :num_particles = " << num_particles << std::endl
       << "Length of whole system :Ns = " << Ns << std::endl
       << "Number of thermalize step : N_thermalize = " << N_thermalize << std::endl
       << "Number of observed patricles :Ns_ObserveSystem = " << Ns_ObserveSystem << std::endl
       << "Developing ime : t = " << t << std::endl
       << "Number of time steps : N_time = " << N_time << std::endl
       << "Number of relaxation step : N_relax = " << N_relax << std::endl
       << "Time interbal : dt = " << dt << std::endl
       << "Result into : " << result_dir << std::endl;
  ofs.close();

 //dataputer set 
  dataput_t dataputer(Ns,num_particles);
 //lattice set
  lattice_t lattice(Ns);
  unsigned int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();
  int counter = 0.0 ; //montecalro counter 
  double E; //energy
  double K; //kinetic energy
  double U; //potential energy
  double S_total; // total spin


  std::string result_dat("result.dat");
  result_dat = result_dir + result_dat;
  ofs.open(result_dat.c_str());
  //vector set 
  std::vector<double> z(4*num_particles);
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
  lattice.create_table(pair_table);
  //randam genereta set 
  std::size_t seed = 0;
  boost::random::mt19937 mt(seed);
  //initializer set
  initializer_t initializer(J,T,num_particles);
  //initial config set 
  initializer.initset(z,S_total,mt);

  initializer.Momentset(z,T,mt);
  //initial config out put
  //hamiltonian set
  hamiltonian_t ham(num_particles,J,pair_table,N_adj); 

  //set integrator
  //integrator_t integrator(2*num_particles);

  //set dataput
  std::string str(result_dir + "vectorField_");
  std::string str_loop;
  std::string dat(".dat");

 
  ofs << "#output_form pt S_total s_total Energy Kinetic Potential" << std::endl;
  ofs << std::scientific;
  ofs << std::setprecision(10);

  double pt = 0.0;
  for(int step = 0; step < N_time; ++step){
    //integrator.step(pt, dt, z, ham);
    str_loop = boost::lexical_cast<std::string>(step);
    str_loop = str + str_loop + dat;
    dataputer.output_config(str_loop, Ns_ObserveSystem,z,lattice);
    E = ham.energy(pt, z);
    K = ham.kinetic_energy(pt,z);
    U = ham.potential_energy(pt,z);
    ofs << pt << " " <<  S_total << " " <<  S_total/num_particles << " " << E << " " << K << " " << U <<  std::endl;
    for(int counter = 0; counter < N_thermalize; ) initializer.montecalro(z,S_total,counter, ham ,mt);
    pt += dt;
  }

  dataputer.endput(condition_dat.c_str());
  return 0;
}//END clXYmodel
 
 
