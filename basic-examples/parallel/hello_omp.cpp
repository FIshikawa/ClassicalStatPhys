#include <omp.h>
#include <iostream>

int main(int argc, char **argv){
  int num_threads;
  int total = 0;
  #pragma omp parallel
  {
  num_threads = omp_get_num_threads();
  int id_thread = omp_get_thread_num();
  for (int p = 0; p < num_threads; ++p) {
    #pragma omp critical 
    if (p == id_thread){
      std::cout << "Hello World from " 
                << id_thread
                << " of "
                << num_threads
                << " threads" << std::endl;
      total += id_thread;
    }
  }
  }
  return 0;
}
