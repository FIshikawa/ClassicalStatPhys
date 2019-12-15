#include <gtest/gtest.h>
#include <cstdio>
#include <clstatphys/tools/data_recorder.hpp>

TEST(DataRecorderTest, BasicTest){
  std::string condition_dat = "condition_basic.dat";
  tools::DataRecorder dataput(condition_dat);
  std::string input_line = "Hello! World!"; 
  dataput << input_line << std::endl;
  dataput.close();
  std::ifstream check_file(condition_dat);
  ASSERT_TRUE(check_file.good());
  std::string read_line;
  std::getline(check_file, read_line);
  EXPECT_TRUE(read_line == input_line);
  check_file.close();
  dataput.initialize_ouputfile();
  dataput.time_tag() << std::endl;
  dataput.close();
  check_file.open(condition_dat);
  std::getline(check_file, read_line);
  EXPECT_FALSE(read_line == input_line);
  if(check_file) std::remove(condition_dat.c_str());
  else std::cerr << "Not find  : " << condition_dat << std::endl;
}

TEST(DataRecorderTest, OutputTest1){
  std::string condition_dat = "condition_test1.dat";
  std::string result_dat = "result_test1.dat";
  tools::DataRecorder dataput(condition_dat);
  dataput << "output test 1" << std::endl;
  dataput.close();
  std::unordered_map<std::string, std::vector<double> > results{};  
  std::vector<std::string> physical_value_list = {
                                                  "Velocity",
                                                  "Energy_total",
                                                  "time"
                                                };
  for(auto key : physical_value_list){
    results[key].resize(10);
    for(int i = 0 ; i < 10; ++i){
      results[key][i] = 1.0E+2;
    }
  }
  dataput.output_result(result_dat,"time",results);
  std::ifstream check_file(result_dat,std::ostream::in);
  std::string read_line;
  std::getline(check_file, read_line);
  EXPECT_TRUE(read_line == "#output_form time Velocity Energy_total " ||
      read_line == "#output_form time Energy_total Velocity " );
  check_file.close();
  for(auto file_name : {condition_dat, result_dat}){
    check_file.open(file_name);
    if(check_file) std::remove(condition_dat.c_str());
    else std::cerr << "Not find  : " << file_name << std::endl;
    check_file.close();
  }
}

TEST(DataRecorderTest, OutputTest2){
  std::string condition_dat = "condition_test2.dat";
  std::string result_dat = "result_test2.dat";
  tools::DataRecorder dataput(condition_dat);
  dataput << "output test 2" << std::endl;
  dataput.close();
  std::unordered_map<std::string, std::vector<double> > results{};  
  std::vector<std::string> physical_value_list = {
                                                  "Velocity",
                                                  "Energy_total"
                                                };
  for(auto key : physical_value_list){
    results[key].resize(10);
    for(int i = 0 ; i < 10; ++i){
      results[key][i] = 1.0E+2;
    }
  }
  dataput.output_result(result_dat,"time",0.01,10,results);
  std::ifstream check_file(result_dat,std::ostream::in);
  std::string read_line;
  std::getline(check_file, read_line);
  EXPECT_TRUE(read_line == "#output_form time Velocity Energy_total " ||
      read_line == "#output_form time Energy_total Velocity " );
  check_file.close();
  for(auto file_name : {condition_dat, result_dat}){
    check_file.open(file_name);
    if(check_file) std::remove(condition_dat.c_str());
    else std::cerr << "Not find  : " << file_name << std::endl;
    check_file.close();
  }
}

TEST(DataRecorderTest, OutputTest3){
  std::string condition_dat = "condition_test3.dat";
  std::string result_dat = "result_test3.dat";
  tools::DataRecorder dataput(condition_dat);
  dataput << "output test 3" << std::endl;
  dataput.close();
  std::unordered_map<std::string, std::vector<double> > results{};  
  std::vector<std::string> physical_value_list = {
                                                  "Velocity",
                                                  "Energy_total"
                                                };
  for(auto key : physical_value_list){
    results[key].resize(10);
    for(int i = 0 ; i < 10; ++i){
      results[key][i] = 1.0E+2;
    }
  }
  std::unordered_map< std::string, std::vector<double> > domain;
  domain["time"] = std::vector<double>(10,1.0);
  dataput.output_result(result_dat,domain,results);
  std::ifstream check_file(result_dat,std::ostream::in);
  std::string read_line;
  std::getline(check_file, read_line);
  EXPECT_TRUE(read_line == "#output_form time Velocity Energy_total " ||
      read_line == "#output_form time Energy_total Velocity " );
  check_file.close();
  for(auto file_name : {condition_dat, result_dat}){
    check_file.open(file_name);
    if(check_file) std::remove(condition_dat.c_str());
    else std::cerr << "Not find  : " << file_name << std::endl;
    check_file.close();
  }
}

TEST(DataRecorderTest, OutputTestHist){
  std::string condition_dat = "condition_hist.dat";
  std::string result_dat = "result_hist.dat";
  tools::DataRecorder dataput(condition_dat);
  dataput << "output hist plot test" << std::endl;
  dataput.close();
  int domain_number = 5;
  int range_number = 10;
  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_values;
  std::vector<std::string> hist_value_list = {
                                              "hist_Velocity",
                                              "range_Velocity",
                                              };
  for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(
      domain_number, std::vector<double>(range_number,1.0));
  std::unordered_map< std::string, std::vector<double> > domain;
  domain["time"] = std::vector<double>(domain_number,1.0);
  dataput.output_result_hist(result_dat,domain,hist_value_list,hist_values);
  std::ifstream check_file(result_dat,std::ostream::in);
  check_file.close();
  for(auto file_name : {condition_dat, result_dat}){
    check_file.open(file_name);
    if(check_file) std::remove(condition_dat.c_str());
    else std::cerr << "Not find  : " << file_name << std::endl;
    check_file.close();
  }
}

TEST(DataRecorderTest, OutputTest3d1){
  std::string condition_dat = "condition_3d_1.dat";
  std::string result_dat = "result_3d_1.dat";
  tools::DataRecorder dataput(condition_dat);
  dataput << "output 3d plot test" << std::endl;
  dataput.close();
  int domain_x_number = 20;
  int domain_y_number = 10;
  std::unordered_map<std::string, std::vector< std::vector<double> > > data_values;
  std::vector<std::string> data_value_list = {
                                              "z",
                                              };
  for(auto key  : data_value_list ) data_values[key] = std::vector< std::vector<double> >(
      domain_y_number, std::vector<double>(domain_x_number,1.0));
  std::unordered_map< std::string, std::vector<double> > domain_x,domain_y;
  domain_x["x"] = std::vector<double>(domain_x_number,1.0);
  domain_y["y"] = std::vector<double>(domain_y_number,1.0);
  for(int i = 0 ; i < domain_x["x"].size(); ++i) domain_x["x"][i] = 0.1 * i; 
  for(int i = 0 ; i < domain_y["y"].size(); ++i) domain_y["y"][i] = 0.1 * i; 
  for(int i = 0 ; i < domain_x["x"].size(); ++i){
    for(int j = 0 ; j < domain_y["y"].size(); ++j){
      data_values["z"][j][i] = 0.1 * i * i * j;
    }
  }
  dataput.output_result_3d(result_dat,domain_y,domain_x,data_value_list,data_values);
  std::ifstream check_file(result_dat,std::ostream::in);
  check_file.close();
  for(auto file_name : {condition_dat, result_dat}){
    check_file.open(file_name);
    if(check_file) std::remove(condition_dat.c_str());
    else std::cerr << "Not find  : " << file_name << std::endl;
    check_file.close();
  }
}

TEST(DataRecorderTest, OutputTest3d2){
  std::string condition_dat = "condition_3d_2.dat";
  std::string result_dat = "result_3d_2.dat";
  tools::DataRecorder dataput(condition_dat);
  dataput << "output 3d plot test" << std::endl;
  dataput.close();
  int domain_x_number = 20;
  int domain_y_number = 10;
  std::unordered_map<std::string, std::vector< std::vector<double> > > data_values;
  std::vector<std::string> data_value_list = {
                                              "z",
                                              "w",
                                              };
  for(auto key  : data_value_list ) data_values[key] = std::vector< std::vector<double> >(
      domain_y_number, std::vector<double>(domain_x_number,1.0));
  std::unordered_map< std::string, std::vector<double> > domain_x,domain_y;
  domain_x["x"] = std::vector<double>(domain_x_number,1.0);
  domain_y["y"] = std::vector<double>(domain_y_number,1.0);
  for(int i = 0 ; i < domain_x["x"].size(); ++i) domain_x["x"][i] = 0.1 * i; 
  for(int i = 0 ; i < domain_y["y"].size(); ++i) domain_y["y"][i] = 0.1 * i; 
  for(int i = 0 ; i < domain_x["x"].size(); ++i){
    for(int j = 0 ; j < domain_y["y"].size(); ++j){
      data_values["z"][j][i] = 0.1 * i * i * j;
      data_values["w"][j][i] = -0.1 * i * j * j;
    }
  }
  dataput.output_result_3d(result_dat,domain_y,domain_x,data_value_list,data_values);
  std::ifstream check_file(result_dat,std::ostream::in);
  check_file.close();
  for(auto file_name : {condition_dat, result_dat}){
    check_file.open(file_name);
    if(check_file) std::remove(condition_dat.c_str());
    else std::cerr << "Not find  : " << file_name << std::endl;
    check_file.close();
  }
}

