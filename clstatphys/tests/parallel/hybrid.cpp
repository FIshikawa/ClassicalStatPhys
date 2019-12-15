#include <mpi.h>
#include <omp.h>

#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <vector>
#include <thread>

namespace {
class HybridTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    //MPI set up
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
  }
  int mpi_error;
  int num_process;
  int process_id;
  int name_length;
  char** processor_name;
  int num_threads;
};

TEST_F(HybridTest, HelloWorldTest){
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
  if(process_id == 0){
    ASSERT_EQ(num_threads_temp, num_threads);
  } 
  std::this_thread::sleep_for(std::chrono::seconds(1));
  mpi_error = MPI_Barrier(MPI_COMM_WORLD);
}

TEST_F(HybridTest, ReduceTest) {
  std::vector<int> values(num_threads,0);
  #pragma omp parallel for 
  for(int i = 0 ; i < num_threads ; ++i){
    int thread_id = omp_get_thread_num();
    values[i] = thread_id * (process_id + 1);
  }  
  int each_total = 0.0;
  #pragma omp parallel for reduction(+:each_total)
  for(int i = 0 ; i < num_threads ; ++i){
    each_total += values[i];
  }
  int total_value = 0;
  mpi_error = MPI_Reduce(&each_total, &total_value, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(process_id == 0){
    int answer = num_threads * ( num_threads - 1 ) / 2 * num_process * (num_process + 1 ) / 2;
    EXPECT_TRUE(total_value ==  answer);
  }
  mpi_error = MPI_Barrier(MPI_COMM_WORLD);
  std::this_thread::sleep_for(std::chrono::seconds(1));
}

TEST_F(HybridTest, RandomGeneratorTest) {
  //set random value for  parallelization
  unsigned long seed = 1234;
  std::minstd_rand seed_gen(seed);
  std::vector<std::shared_ptr<std::mt19937> > mts(num_threads);
  #pragma omp parallel 
  {
  int thread_id = omp_get_thread_num();
  for (int tag = 0 ; tag < num_process; ++tag){
    if (tag == process_id){
      for (int p = 0; p < num_threads; ++p) {
          if (p == thread_id) mts[p].reset(new std::mt19937(seed_gen()));
          #pragma omp barrier
      }
    }
    else for(int p = 0; p < num_threads; ++p) seed_gen();
  }
  }

  //set random value for  parallelization end 
  std::vector<double> initial_values(num_threads,0.0);
  #pragma omp parallel for
  for(int i = 0 ; i < num_threads ; ++i){
    int thread_id = omp_get_thread_num();
    std::mt19937& mt = *mts[thread_id];
    std::uniform_real_distribution<> dist(0.0, 10.0);
    initial_values[i] = dist(mt);
  }
  std::this_thread::sleep_for(std::chrono::seconds(1));
  mpi_error = MPI_Barrier(MPI_COMM_WORLD);

  // Gather data 
  std::vector<double> initial_values_all(num_threads * num_process,-1.0); 
  mpi_error = MPI_Gather(&initial_values[0], num_threads, MPI_DOUBLE, &initial_values_all[0], num_threads, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //set cout  
  std::cout << std::fixed;
  std::cout << std::setprecision(3);

  //check each data 
  for(int prid = 0; prid < num_process ; ++prid){
    if(prid == process_id){
      for (int thid = 0; thid < num_threads; ++thid) std::cout << initial_values[thid] << " " << std::flush;
    }
    std::this_thread::sleep_for(std::chrono::seconds(1));
    mpi_error = MPI_Barrier(MPI_COMM_WORLD);
  }
  if(process_id == 0){
    std::cout << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    for(int i = 0 ; i < num_threads*num_process ; ++i) std::cout << initial_values_all[i] << " "; std::cout << std::endl;
    for(int i = 0 ; i < num_threads*num_process; ++i){
      EXPECT_TRUE(initial_values_all[i] != -1.0);
    }
    for(int i = 0 ; i < initial_values_all.size() ; ++i){
      for(int j = 0 ; j < num_threads*num_process; ++j){
        if(i != j) {
          EXPECT_TRUE(initial_values_all[i] != initial_values_all[j]);
        }
      }
    }
  }
  mpi_error = MPI_Barrier(MPI_COMM_WORLD);
  std::this_thread::sleep_for(std::chrono::seconds(1));
}

}//end namespace


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int mpi_error;
  mpi_error = MPI_Init(&argc, &argv);
  ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank != 0) delete listeners.Release(listeners.default_result_printer());
  int result = RUN_ALL_TESTS();
  mpi_error = MPI_Finalize();
  return result;
}
