#include <mpi.h>
#include <omp.h>

#include <fstream>
#include <string>
#include <chrono>
#include <iostream>
#include <memory>

#include <lattice/fullyconnected.hpp>
#include <ensemble/water_bag_parallel.hpp>
#include <integrator/velocity_velret_parallel.hpp>
#include <physics/classical_xy_fullyconnected_parallel.hpp>

int main(int argc, char **argv){
  int mpierr;
  int process_id;
  int num_process;
  int world_size;
  int world_rank;
  int name_length;
  mpierr = MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  char processor_name_temp[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processor_name_temp, &name_len);
  std::string processor_name(processor_name_temp);

  std::cout << "[test parallelaizeation by hybrid]\n";

  //set parallelized int 
  int N_parallel_max = omp_get_num_procs();
  std::cout << "[limit threads : " << N_parallel_max << "]\n" ;


  //set output 
  std::cout<< "open ouput_file" << std::endl;
  std::ofstream ofs("result.dat",std::ios::out);
  std::ofstream ofsr("result_TimeAtNsMax.dat",std::ios::out);
  ofs.close();
  ofs.open("result.dat",std::ios::app);
  ofs << std::scientific;
  ofsr << std::scientific;
  ofs << "#output Ns ";
  ofsr << "#output Np time" << std::endl;

  for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){
    ofs << "time_Np_" << std::to_string(N_parallel) << " ";
  }
  ofs << std::endl;
  //set output end 

  std::cout<< "start Ns iteration" << std::endl;

  int N_iterate = 10;
  int Ns_init = 10;
  int Ns_max = 100;
  int dNs = (Ns_max - Ns_init)/N_iterate;

  for(int i = 0; i < N_iterate; ++i){
    int Ns = Ns_init + dNs * i;

    std::cout<< "create lattice" << std::endl;
    lattice::FullyConnected lattice(Ns);
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

    //target class def
    std::cout<< "generate instance by the classes" << std::endl;
    double J = 1.0;
    double T = 1.0;
    hamiltonian::ClassicalXYFullyConnectedParallel hamiltonian(num_particles,J,N_adj);
    integrator::VelocityVelretParallel integrator(2*num_particles);
    ensemble::WaterBagPrallel initializer(J,T,num_particles);
    //target class def end

    ofs << Ns;
    for( int N_parallel = 1; N_parallel <= N_parallel_max ; ++N_parallel){
      omp_set_num_threads(N_parallel);
      double return_value = 0.0;
      
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
      double s_total_x = 0.0;
      double s_total_y = 0.0;
      double pt = 0.0;
      double dt = 0.01;
      initializer.set_initial_state(z,mts);
      #pragma omp parallel for reduction(+:s_total_x,s_total_y)
      for(int j = 0 ; j < num_particles; ++j){
        s_total_x += std::cos(z[j]);
        s_total_y += std::sin(z[j]);
      }
      hamiltonian.set_meanfield(s_total_x, s_total_y); 
      integrator.step(pt,dt,z,hamiltonian);
      return_value = hamiltonian.energy(pt,z);
      //target behavior end 
      end = std::chrono::system_clock::now();
      double result_value = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count(); 

      std::cout << " return value : " << return_value << std::endl;
      ofs << " " << result_value; 

      if(i == N_iterate - 1) ofsr << N_parallel << " " << result_value << std::endl;
    }
    ofs << std::endl;
  }
  std::cout << "[finish test]\n"; 
  return 0;
}
