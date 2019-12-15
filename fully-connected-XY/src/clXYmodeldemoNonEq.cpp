//square lattice clXYmodel monte carlp
#include<iostream>
#include<cmath>
#include<fstream>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<boost/random.hpp>
#include<boost/random/uniform_int_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

//#include "integrator.hpp"
#include "dissipated_integrator.hpp"
#include "setcondi.hpp"
#include "initializer.hpp"
#include "dataput.hpp"
#include "hamiltonian.hpp"
#include "Lattice.hpp"
#include "correlation.hpp"

//typedef integrator::velver integrator_t;
typedef dissipated_integrator::rk2 integrator_t;
typedef initializer::equilclXYmodelWolf initializer_t;
//typedef initializer::equilLangevin initializer_t;
//typedef initializer::equilclXYmodel initializer_t;
//typedef Lattice::RegularCube lattice_t;
typedef Lattice::SquareLattice lattice_t;
typedef dataput::dataputForSquareLattice dataput_t;
//typedef dataput::dataputForRegularCube dataput_t;
//typedef hamiltonian::classicalXY hamiltonian_t;
//typedef hamiltonian::freeparticle hamiltonian_t;
//typedef hamiltonian::classicalXYNonEq hamiltonian_t;
typedef hamiltonian::classicalXYAnisoNonEq hamiltonian_t;

int main(int argc, char **argv){
 double J = 1.0;
 double E; //energy
 double K; //kinetic energy
 double U; //potential energy
 double T; //temperature
 double S_total; // total spin
 int Ns; // lenght of whole system
 int N_thermalize; // number of thermalization
 int Ns_ObserveSystem; // length of observe system
 int num_particles; //nummer of particles
 int counter ; //montecalro counter 
 double t; // Whole time
 double dt; //size of interval dt = t/N_time
 int N_time; //numer of time step
 unsigned int N_adj;
 double Inte = 0.0; //Intencity of External Force
 double Freq = 0.0; //Frequency of External Force
 double gamma = 1.0; //Friction of Dissipation, satisfying D = gamma * T 
 double D = 1.0; //Anisotropy
 char* condition_dat = "condi.dat";
 char* parameter_dat = "parameter.dat";
 char* result_dat = "result.dat";

//parameter set 

 if (argc > 1) N_thermalize =     boost::lexical_cast<int>(argv[1]);
 if (argc > 2) Ns =               boost::lexical_cast<int>(argv[2]);
 if (argc > 3) Ns_ObserveSystem = boost::lexical_cast<int>(argv[3]);
 if (argc > 4) T =                boost::lexical_cast<double>(argv[4]);
 if (argc > 5) N_time =           boost::lexical_cast<int>(argv[5]);
 if (argc > 6) t =                boost::lexical_cast<double>(argv[6]);
 if (argc > 7)Inte =             boost::lexical_cast<double>(argv[7]);
 if (argc > 8)Freq =             boost::lexical_cast<double>(argv[8]);
 if (argc > 9)D =                boost::lexical_cast<double>(argv[9]);

 // parameter auto set
 num_particles = lattice_t::set_num_particles(Ns);
 dt = t/(double)N_time;

 std::vector<std::string> experimental_condition(5);
 experimental_condition[0] = "Config Movie : Non Equilibrium Dynamics in 3D anisotropic classical XY with periodical external field";
 experimental_condition[1] = hamiltonian_t::name();
 experimental_condition[2] = initializer_t::name();
 experimental_condition[3] = integrator_t::name();
 experimental_condition[4] = lattice_t::name();
 dataput_t::startput(condition_dat, experimental_condition);

 std::ofstream ofs;
 ofs.open(condition_dat,std::ios::app);
 if(!ofs)std::cerr<< condition_dat <<"opening failed"<<std::endl;
 ofs  << "Coupling constant: J is " << J << std::endl
      << "Temperture: T = " << T << std::endl
      << "Reduced Coupling: J/T = " << J/T << std::endl
      << "Number of patricles :num_particles = " << num_particles << std::endl
      << "Length of whole system :Ns = " << Ns << std::endl
      << "Lenght of ObserveSystem :Ns_ObserveSystem = " << Ns_ObserveSystem << std::endl
      << "Number of thermalize step : N_thermalize = " << N_thermalize << std::endl
      << "Developing ime : t = " << t << std::endl
      << "Number of time steps : N_time = " << N_time << std::endl
      << "Time interbal : dt = " << dt << std::endl
      << "Intencity of the external field : Inte = " << Inte << std::endl
      << "Frequency of the external field : Freq = " << Freq << std::endl
      << "Strength of the Anisotropy : D = " << D << std::endl;
 ofs.close();
 //dataputer set 
 dataput_t dataputer(Ns,num_particles,N_time);

 //lattice set
 lattice_t lattice(Ns);
 N_adj = lattice.number_adjacent();
 //vector set 
 std::vector<double> z(2*num_particles);
 std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
 lattice.create_table(pair_table);
 //randam genereta set 
 std::size_t seed = 0;
 boost::random::mt19937 mt(seed);
 //initializer set
 initializer_t initializer(J,T,num_particles);
 //initial config set 
 //initializer.initset(z,S_total,counter,mt);
 initializer.initset(z,counter,mt);
 //initial config out put
 //hamiltonian set
 //hamiltonian_t ham(num_particles,J,pair_table,N_adj,Inte,Freq); 
 hamiltonian_t ham(num_particles,J,D,pair_table,N_adj,Inte,Freq); 
 //set integrator
 //integrator_t integrator(2*num_particles);
 integrator_t integrator(T, gamma, 2*num_particles);
 //set dataput
 std::string str("./mapdata/vectorField_");
 std::string str_loop;
 std::string dat(".dat");

//thermalize
 E = ham.energy(0, z);
 K = ham.kinetic_energy(0,z);
 U = ham.potential_energy(0,z);

 std::cout << "Initial" << " S_total: " <<  S_total << "  s_total: " <<  S_total/num_particles << " E: " << E << " K: " << K << " U: " << U <<  std::endl;
 initializer.Momentset(z,T,mt);
 for(int couter = 0; counter < N_thermalize; ++counter){ 
//  initializer.montecalro(z,S_total,counter, ham ,mt);
  initializer.montecalro(z,counter, ham ,mt);
 }
 double *v = &z[num_particles];
 double v_total = 0;
 for(int i = 0; i < num_particles; ++i) v_total += v[i];
 for(int i = 0; i < num_particles; ++i) v[i] -= v_total;
 //initial energy calc
 E = ham.energy(0, z);
 K = ham.kinetic_energy(0,z);
 U = ham.potential_energy(0,z);
 str_loop = boost::lexical_cast<std::string>("init");
 str_loop = str + str_loop + dat;
 dataputer.output_config(str_loop, Ns_ObserveSystem,z,lattice);

 std::cout << "Thermalized" << " S_total: " <<  S_total << "  s_total: " <<  S_total/num_particles << " E: " << E << " K: " << K << " U: " << U <<  std::endl;

 double pt = 0.0;
 //time step

 ofs.open(result_dat,std::ios::out);
 ofs << "# out put form #1 pt #2 S_total #3 s_total #4 E #5 K #6 U #7 internal E" << std::endl;
 ofs << pt << " " << S_total << " " << S_total/num_particles << " " << E << " " << K << " " << U << " " << K + U << std::endl; 
 //boost::normal_distribution<> normal_dist(0.0,1.0);
 //for(int i = 0; i < num_particles ; ++i) z[i + num_particles] = normal_dist(mt);

 for(int step = 0; step < N_time; ++step){
   //integrator.step(pt, dt, z, ham);
   integrator.step(pt, dt, z, ham, mt);
   E = ham.energy(pt, z);
   K = ham.kinetic_energy(pt,z);
   U = ham.potential_energy(pt,z);
   str_loop = boost::lexical_cast<std::string>(step);
   str_loop = str + str_loop + dat;
   dataputer.output_config(str_loop, Ns_ObserveSystem,z,lattice);
   S_total = 0;
   for(int i = 0; i < num_particles; ++i) S_total += std::cos(z[i]);
   pt += dt;
   ofs << pt << " " << S_total << " " << S_total/num_particles << " " << E << " " << K << " " << U << " " << K + U << std::endl; 
   if(step == N_time/2){ 
     E = ham.energy(pt,z);
     K = ham.kinetic_energy(pt,z);
     U = ham.potential_energy(pt,z);
     std::cout << "Time half steped" << " S_total: " <<  S_total << "  s_total: " <<  S_total/num_particles << " E: " << E << " K: " << K << " U: " << U <<  std::endl;
   }
 
 }//end time develop
 ofs.close();

 E = ham.energy(pt,z);
 K = ham.kinetic_energy(pt,z);
 U = ham.potential_energy(pt,z);

 for(int i = 0; i < num_particles; ++i) S_total += std::cos(z[i]);
 std::cout << "Time steped" << " S_total: " <<  S_total << "  s_total: " <<  S_total/num_particles << " E: " << E << " K: " << K << " U: " << U <<  std::endl;
 //end config out put

 str_loop = boost::lexical_cast<std::string>("fin");
 str_loop = str + str_loop + dat;

 dataputer.output_config(str_loop,Ns_ObserveSystem,z,lattice);

 dataputer.endput(condition_dat);
return 0;
}//END clXYmodel
 
 



   

