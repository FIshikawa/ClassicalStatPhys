#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<memory>
#include<random>
#include<string>
#include<unordered_map>

#include<boost/date_time/posix_time/posix_time.hpp>

#include"WangLandau.hpp"
#include"hamiltonian.hpp"
#include"exactcalc.hpp"
#include"Lattice.hpp"
#include"dataput.hpp"
#include<omp.h>


typedef hamiltonian::classicalXY_FullyConnected_Parallel hamiltonian_t;
typedef metadynamics::WangLandau_Parallelized metadynamics_t;
typedef Lattice::FullyConnected lattice_t; 
typedef dataput::dataput dataput_t;

int main(int argc, char **argv){
  //parameter define
  double J = 1.0; //interaction
  int num_particles = 10; //nummer of particles
  int Ns = 10; // lenght of whole system
  unsigned int dimension = 2 * num_particles;
  double bin = 0.01;
  double E_max = 0.5;
  double E_min = 0.0;
  unsigned int N_bin = int((E_max - E_min) / bin);
  int max_iteration = 1;
  int N_parallel_max = omp_get_num_procs();
  std::string condition_dat = "condi.dat";
  std::string result_dir("./");

 //parameter set 
  if (argc > 1) result_dir =       boost::lexical_cast<std::string>(argv[1]);
  if (argc > 2) Ns =               boost::lexical_cast<unsigned int>(argv[2]);
  if (argc > 3) bin =              boost::lexical_cast<double>(argv[3]);
  if (argc > 4) E_max =            boost::lexical_cast<double>(argv[4]);
  if (argc > 5) E_min =            boost::lexical_cast<double>(argv[5]);
  if (argc > 6) max_iteration =    boost::lexical_cast<unsigned int>(argv[6]);

  // parameter auto set
  num_particles = lattice_t::set_num_particles(Ns);
  condition_dat = result_dir + condition_dat;
  N_bin = int((E_max - E_min) / bin);
  dimension = 2 * num_particles;

  //parameter set 
  dataput_t dataputer(condition_dat);
  boost::posix_time::ptime nowtime = boost::posix_time::second_clock::local_time();
  dataputer << "# start time : " << nowtime << std::endl 
            << "Observe MD relaxed Critical Behavior" << std::endl
            << hamiltonian_t::name() << std::endl
            << metadynamics_t::name() << std::endl
            << lattice_t::name() << std::endl;
 
  //condition declare
  dataputer  << "Coupling constant: J is " << J << std::endl
             << "Number of patricles : num_particles = " << num_particles << std::endl
             << "Length of whole system : Ns = " << Ns << std::endl
             << "Bin width : bin = " << bin << std::endl
             << "Max Energy of search space : E_max = " << E_max << std::endl
             << "Min Energy of search space : E_min = " << E_min << std::endl
             << "Number of bins : N_bin =" << N_bin << std::endl
             << "Max Number of parallel : N_parallel_max = " << N_parallel_max << std::endl
             << "Result into : " << result_dir << std::endl;

  //declare classes 
  metadynamics_t wang_landau(dimension, N_bin, E_max, E_min); 
  lattice_t lattice(Ns);
  int N_adj = lattice.number_adjacent();
  hamiltonian_t  ham(num_particles,J,N_adj);

  //set random value for  parallelization
  unsigned long seed = 1234;
  std::minstd_rand seed_gen(seed);
  std::mt19937 mt(seed_gen());
  std::vector<std::shared_ptr<std::mt19937> > mts(N_parallel_max);
  omp_set_num_threads(N_parallel_max);
  #pragma omp parallel 
  {
  int tid = omp_get_thread_num();
  int num_threads = omp_get_num_threads();
  for (int p = 0; p < num_threads; ++p) {
    if (p == tid)
      mts[p].reset(new std::mt19937(seed_gen()));
      #pragma omp barrier
  }
  }
  //set random value for  parallelization end 

  //set E_calc 
  auto E_calc = [&ham,&num_particles](std::vector<double> const & z){
    double s_total_x = 0.0;
    double s_total_y = 0.0;
    #pragma omp parallel for reduction(+:s_total_x,s_total_y)
    for(int j = 0 ; j < num_particles; ++j){
      s_total_x += std::cos(z[j]);
      s_total_y += std::sin(z[j]);
    }
    ham.set_meanfield(s_total_x, s_total_y); 
    std::cout << " Potential E : " << ham.potential_energy(0.0,z)/(double)num_particles << std::endl;
    return ham.potential_energy(0.0,z)/(double)num_particles;  
  };
  //end E_calc

  auto Proposer = [](std::vector<double> & z, std::vector<std::shared_ptr<std::mt19937> > & mts){
    unsigned int num_particles = z.size() / 2;
    double *x = &z[0];
    double *p = &z[num_particles];
    #pragma omp parallel
    {
    int tid = omp_get_thread_num();
    #pragma omp for
    for(int i = 0; i < num_particles ; ++i){
      std::mt19937& mt = *mts[tid];
      std::uniform_real_distribution<> x_Rand(0.0,2.0*M_PI);
      x[i] = x_Rand(mt);
    }
    }
  };
  //end Proposer 

  //do wang Landau
  wang_landau.montecalro(E_calc,Proposer,max_iteration,mt,mts);

  //output 
  std::unordered_map<std::string, std::vector<double> > physical_values{};  
  std::vector<std::string> physical_value_list = {"bin", "density", "entropy"};
  for(auto key : physical_value_list) physical_values[key] = std::vector<double>(N_bin,0);
  physical_values["bin"] = wang_landau.bins();
  physical_values["density"] = wang_landau.density();
  physical_values["entropy"] = wang_landau.entropy();
 

  std::string result_dat("result.dat");
  result_dat = result_dir + result_dat;
  dataputer.output_result(result_dat, physical_value_list, physical_values);  

  //anounce end time and end total energy to condi.dat
  nowtime = boost::posix_time::second_clock::local_time();
  dataputer << "# end time : " << nowtime << std::endl; 

  return 0;
}
