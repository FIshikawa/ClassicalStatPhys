#include <iostream>
#include <cstdio>
#include <clstatphys/tools/data_recorder.hpp>

int main() {

  std::string condition_dat = "condition.dat";
  std::cout << "[test for DataRecorder]\n";
  tools::DataRecorder dataput(condition_dat);
  dataput << "Time declare! " 
          << dataput.time() 
          << "Hello! World!" << std::endl
          << "You have two eyes but one mouth!" << std::endl;

  dataput.time_tag() << "check time method" << std::endl;
  std::ifstream check_file(condition_dat,std::ifstream::in);
  if(check_file.good()) std::remove(condition_dat.c_str());
  else std::cout << "Not find  : " << condition_dat << std::endl;
  return 0;
}

