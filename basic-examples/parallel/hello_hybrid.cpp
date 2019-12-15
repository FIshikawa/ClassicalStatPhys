#include <mpi.h>
#include <omp.h>

#include <thread>
#include <iostream>

int main(int argc, char **argv){
  int mpi_error;
  int num_process;
  int process_id;
  int name_length;
  int num_threads;

  mpi_error = MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_process);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processor_name, &name_length);

  #pragma omp parallel
  {
  #pragma omp single
  { 
    num_threads = omp_get_num_threads(); 
  }
  }

  if(process_id == 0){
    #pragma omp parallel
    {
    #pragma omp single
    { 
    std::cout << " num_process : "
              << num_process;
    std::cout << ", in parallel, now threads : " 
              << num_threads 
              << ", omp_get_num_procs : "
              << omp_get_num_procs()
              << std::endl;
    }
    }
  }
  int num_threads_temp;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  mpi_error = MPI_Barrier(MPI_COMM_WORLD);
  for(int prid = 0; prid < num_process ; ++prid){
    if(prid == process_id){
      #pragma omp parallel
      {
      num_threads_temp = omp_get_num_threads();
      int id_thread = omp_get_thread_num();
      for (int thid = 0; thid < num_threads; ++thid) {
        #pragma omp critical 
        if (thid == id_thread){
          std::cout << "Hello world from processor " << processor_name
                <<" : rank " << process_id
                <<" out of " << num_process << " processes" 
                <<" : threads id " << omp_get_thread_num()
                <<" out of " << omp_get_num_threads()
                << std::endl;
        }
      }
      }
    }
    std::this_thread::sleep_for(std::chrono::seconds(1));
    mpi_error = MPI_Barrier(MPI_COMM_WORLD);
  }
  mpi_error = MPI_Finalize();
  return 0;
}
