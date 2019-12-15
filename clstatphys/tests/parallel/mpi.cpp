#include <mpi.h>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
class MPITest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_length);
  }
  int mpi_error;
  int num_process;
  int process_id;
  int name_length;
  char** processor_name;
  int N_parallel_max;
};

TEST_F(MPITest, HelloWorldTest) {
  std::cout << "Hello world from processor " << processor_name
            <<" : rank " << process_id
            <<" out of " << num_process << " processes" 
            << std::endl;
}

TEST_F(MPITest, GatherTest) {
  int each_value = process_id;
  std::vector<int> gather_values(num_process);
  mpi_error = MPI_Gather(&each_value, 1, MPI_INT, &gather_values[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(process_id == 0){
    for(int i = 0; i < num_process; ++i){
      EXPECT_EQ(gather_values[i],i);
    }
  }
}

TEST_F(MPITest, ReduceTest) {
  int each_value = process_id;
  int total_value;
  mpi_error = MPI_Reduce(&each_value, &total_value, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(process_id == 0){
    int answer =  num_process * (num_process - 1 ) / 2;
    EXPECT_EQ(total_value,answer);
  }
}

TEST_F(MPITest, ReduceMapTest) {
  std::unordered_map<std::string,int> results;
  results["id"] = process_id;
  int total_value;
  mpi_error = MPI_Reduce(&results["id"], &total_value, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(process_id == 0){
    int answer =  num_process * (num_process - 1 ) / 2;
    EXPECT_EQ(total_value,answer);
  }
}

TEST_F(MPITest, ReduceMapVectorTest) {
  std::unordered_map<std::string,std::vector<int>> results;
  results["id"] = std::vector<int> {process_id,process_id,process_id,process_id};
  std::vector<int> total_vector(4);
  mpi_error = MPI_Reduce(&results["id"][0], &total_vector[0], 4, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(process_id == 0){
    int answer =  num_process * (num_process - 1 ) / 2 * 4;
    int total_value = 0;
    for(int i = 0 ; i < 4; ++i) total_value += total_vector[i];
    EXPECT_EQ(total_value,answer);
  }
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
