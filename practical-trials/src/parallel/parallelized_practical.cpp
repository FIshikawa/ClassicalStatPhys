#include<omp.h>
#include<vector>
#include<random>
#include<fstream>
#include<string>
#include<chrono>
#include<iostream>
#include<memory>
#include<unordered_map>
#include<iomanip>

#include "accumulator.hpp"
#include "Lattice.hpp"
#include "initializer.hpp"
#include "integrator.hpp"
#include "hamiltonian.hpp"

int main(int argc, char **argv){

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[test parallelaizeation]\n";

  int N_parallel_max = omp_get_num_procs();
  std::cout << "[limit threads : " << N_parallel_max << "]\n" ;

  //set output 
  std::cout<< "open ouput_file" << std::endl;
  std::ofstream ofs("result.dat",std::ios::out);
  std::ofstream ofsr("result_TimeAtNsMax.dat",std::ios::out);
  std::ofstream ofsrr("result_PhysicalValue.dat",std::ios::out);
  ofs.close();
  ofs.open("result.dat",std::ios::app);
  ofs << std::scientific << std::setprecision(10);
  ofsr << std::scientific << std::setprecision(10);
  ofsrr << std::scientific << std::setprecision(10);
  std::cout << std::scientific << std::setprecision(10);
  ofs << "#output Ns ";
  ofsr << "#output Np time" << std::endl;
  std::vector<std::string> physical_value_list = {"Total_Spin", "Total_Energy", "Total_Velocity"};
  ofsrr << "#output Np ";
  for(auto key : physical_value_list) ofsrr << key << " " ;
  ofsrr << std::endl;

  for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){
    ofs << "time_Np_" << std::to_string(N_parallel) << " ";
  }
  ofs << std::endl;
  //set output end 

  std::cout<< "start Ns iteration" << std::endl;

  int N_iterate_max = 1000;
  int Ns_init = 1e4;
  int Ns_max = 1e6;
  int dNs = (Ns_max - Ns_init)/N_iterate_max;

  for(int N_iterate = 0; N_iterate < N_iterate_max; ++N_iterate){
    int Ns = Ns_init + dNs * N_iterate;

    std::cout<< "create lattice" << std::endl;
    Lattice::FullyConnected lattice(Ns);
    int N_adj = lattice.number_adjacent();
    std::cout << "N_adj : " << N_adj << std::endl;
    int num_particles = lattice.set_num_particles(Ns);
    std::cout << "num_particles : " << num_particles << std::endl;
    std::vector<double> z(2*num_particles);

    //set random value for  parallelization
    unsigned long seed = 1234;
    std::minstd_rand seed_gen(seed);
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

    ofs << Ns;
    for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){

      //target class def
      std::cout<< "generate instance by the classes" << std::endl;
      double J = 1.0;
      double T = 1.0;
      double dt = 0.01;
      double t = 1.0;
      int N_time = (int)(t / dt);
      int N_sample = 10;
      hamiltonian::classicalXY_FullyConnected_Parallel ham(num_particles,J,N_adj);
      integrator::velver_parallel integrator(2*num_particles);
      initializer::WaterBag_Prallelized initializer(J,T,num_particles);
      std::unordered_map<std::string, std::vector<stat::accumulator> > physical_values{};  
      for(auto key : physical_value_list){
        physical_values[key].resize(N_time);
      }
      //target class def end

      omp_set_num_threads(N_parallel);
      
      #pragma omp parallel
      {
      #pragma omp single
      { 
        std::cout << "[in parallel : now threads : " << omp_get_num_threads() << "]\n";
      }
      }

      std::chrono::system_clock::time_point  start, end;
      start = std::chrono::system_clock::now();

      //target behavior 
      double s_total_x;
      double s_total_y;
      double v_total;
      double pt = 0.0;
      double *x = &z[0];
      double *v = &z[num_particles];
      for(int counter = 0; counter < N_sample; ++counter){
        initializer.set_initial_state(z,mts);
        //time develop
        for(int i = 0; i < N_time ; ++i){
          s_total_x = 0.0;
          s_total_y = 0.0;
          #pragma omp parallel for reduction(+:s_total_x,s_total_y)
          for(int j = 0 ; j < num_particles; ++j){
            s_total_x += std::cos(z[j]);
            s_total_y += std::sin(z[j]);
          }
          ham.set_meanfield(s_total_x, s_total_y); 
          v_total = 0.0;
          #pragma omp parallel for reduction(+:v_total)
          for(int j = 0 ; j < num_particles; ++j){
            v_total += v[j];
          }
          physical_values["Total_Spin"][i] << std::sqrt(s_total_x * s_total_x + s_total_y * s_total_y)/num_particles;
          physical_values["Total_Energy"][i] << ham.energy(pt,z) / num_particles;
          physical_values["Total_Velocity"][i] << v_total / num_particles;
          integrator.step(pt,dt,z,ham);
          pt += dt;
        }
        //time develop end
      }
      //target behavior end 

      end = std::chrono::system_clock::now();
      double result_value = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); 

      //result calc 
      for(auto key  : physical_value_list){
        std::cout << "[" << key << "] "
                  << "value : " << physical_values[key][N_time-1].mean() << " +- " << physical_values[key][N_time-1].error() 
                  << std::endl;
      }
      //result calc 

      ofs << " " << result_value; 

      if(N_iterate == N_iterate_max - 1){
        ofsr << N_parallel << " " << result_value << std::endl;
        ofsrr << N_parallel << " " ;
        for(auto key  : physical_value_list) ofsrr << physical_values[key][N_time-1].mean() << " " ; 
        ofsrr << std::endl;
      }
    }
    ofs << std::endl;
  }
  std::cout << "[finish test]\n"; 
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
