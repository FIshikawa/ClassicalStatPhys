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
  std::cout << "Hello world from processor " 
            << processor_name
            <<" : rank " 
            << process_id
            <<" out of " 
            << num_process 
            << " processes"
            << std::endl;

  mpi_error = MPI_Finalize();
  return 0;
}
