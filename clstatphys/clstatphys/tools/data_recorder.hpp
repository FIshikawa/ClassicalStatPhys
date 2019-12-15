#ifndef DATARECORDER_HPP
#define DATARECORDER_HPP

#include<vector>
#include<cmath>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<unordered_map>

namespace tools {

class DataRecorder{
public:
  DataRecorder(const std::string condidat) : condidat_(condidat) {initialize_ouputfile();}

  template <class Sentence>
  DataRecorder& operator<<(const Sentence &sentence) {
    ofs_ << sentence;
    std::cout << sentence;
    return *this;
  }

  DataRecorder& operator<< (std::ostream& (*__pf)(std::ostream&))  { 
    std::ostream &ofs(ofs_);
    __pf(ofs);
    __pf(std::cout); 
    return *this; 
  }

  void initialize_ouputfile(){
    if(ofs_.good()) close();
    ofs_.open(condidat_.c_str(),std::ofstream::out);
    close();
    open();
  }

  void open(){
    ofs_.open(condidat_.c_str(),std::ofstream::app);
  }

  void close(){
    ofs_.close();
  }

  std::string time(){
    std::chrono::system_clock::time_point nowtime = std::chrono::system_clock::now();
    std::time_t nowtime_c = std::chrono::system_clock::to_time_t(nowtime);
    char sentence_temp[100];
    std::strftime(sentence_temp, sizeof(sentence_temp), "%F %T", std::localtime(&nowtime_c));
    std::string sentence = sentence_temp;
    sentence = sentence + " ";
    return sentence;
  } 

  DataRecorder& time_tag(){
    std::chrono::system_clock::time_point nowtime = std::chrono::system_clock::now();
    std::time_t nowtime_c = std::chrono::system_clock::to_time_t(nowtime);
    ofs_ << "[" << std::put_time(std::localtime(&nowtime_c), "%F %T") << "] ";
    std::cout << "[" << std::put_time(std::localtime(&nowtime_c), "%F %T") << "] ";
    return *this;
  }

  void output_result(const std::string filename, std::string pt_name, double dt, int N_t, std::unordered_map<std::string, std::vector<double> >& physical_values){
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific;
    ofs << std::setprecision(10);
    double pt = 0.0;
    ofs << "#output_form " + pt_name + " ";
    std::vector<std::string> physical_value_list;
    for(auto element : physical_values){
      physical_value_list.push_back(element.first);
      ofs << element.first << " " ;
    } 
    ofs << std::endl;
    for(int i = 0; i < N_t; ++i){
      ofs << pt << " " ;
      for(auto key : physical_value_list) ofs << physical_values[key][i] << " " ;
      ofs << std::endl;
      pt += dt;
    }
    ofs.close();
    }

  void output_result(const std::string filename, std::unordered_map<std::string, std::vector<double> > & domain, std::unordered_map<std::string, std::vector<double> > & physical_values){
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific;
    ofs << std::setprecision(10);
    ofs << "#output_form ";
    int dim_domain;
    std::string domain_name;
    for(auto element : domain){
      domain_name = element.first;
      ofs << domain_name << " ";
      dim_domain = element.second.size();
    }
    std::vector<std::string> physical_value_list;
    for(auto element : physical_values){
      physical_value_list.push_back(element.first);
      ofs << element.first << " " ;
    }
    ofs << std::endl;
    for(int i = 0; i < dim_domain; ++i){
      ofs << domain[domain_name][i] << " ";
      for(auto key : physical_value_list) ofs << physical_values[key][i] << " " ;
      ofs << std::endl;
    }
    ofs.close();
    }

  void output_result(const std::string filename, std::string pt_name, std::unordered_map<std::string, std::vector<double> >& physical_values){
    if(physical_values.find(pt_name) != physical_values.end() ){
      std::ofstream ofs(filename.c_str());
      ofs << std::scientific;
      ofs << std::setprecision(10);
      ofs << "#output_form " + pt_name + " ";
      std::vector<std::string> physical_value_list{pt_name};
      for(auto element : physical_values){
        if(element.first != pt_name){
          physical_value_list.push_back(element.first);
          ofs << element.first << " " ;
        }
      } 
      ofs << std::endl;
      int N_t = physical_values[physical_value_list[0]].size();
      for(int i = 0; i < N_t; ++i){
        for(auto key : physical_value_list) ofs << physical_values[key][i] << " " ;
        ofs << std::endl;
      }
      ofs.close();
    }
    else std::cerr << "Not " << pt_name << " is contained in maps" << std::endl;
  }

  void output_result_hist(const std::string filename, std::unordered_map<std::string, std::vector<double> >& domain, std::vector<std::string> physical_value_list, std::unordered_map<std::string, std::vector<std::vector<double> > >& physical_values){
    std::ofstream ofs(filename);
    ofs << std::scientific;
    ofs << std::setprecision(10);
    int dim_domain = 0;
    for(auto element : domain){
      ofs << "#domain_form " << element.first << " ";
      dim_domain = element.second.size();
      for(int i = 0; i < dim_domain ; ++i) ofs << element.second[i] << " ";
    }
    ofs << std::endl;
    ofs << "#output_form ";
    for(int i = 0 ; i < dim_domain; ++i){
      for(auto key : physical_value_list){
        ofs << key << " " ;
      } 
    }
    ofs << std::endl;
    int N_t = physical_values[physical_value_list[0]][0].size();
    for(int i = 0; i < N_t; ++i){
      for(int j = 0 ; j < dim_domain; ++j){
        for(auto key : physical_value_list) ofs << physical_values[key][j][i] << " " ;
      }
      ofs << std::endl;
    }
    ofs.close();
  }

  void output_result_3d(const std::string filename, std::unordered_map<std::string, std::vector<double> >& domain_y, std::unordered_map<std::string, std::vector<double> >& domain_x, std::vector<std::string> physical_value_list, std::unordered_map<std::string, std::vector<std::vector<double> > >& physical_values){
    std::ofstream ofs(filename);
    ofs << std::scientific;
    ofs << std::setprecision(10);
    int dim_domain_y = 0;
    int dim_domain_x = 0;
    std::string domain_x_name;
    for(auto element : domain_y){
      ofs << "#domain_form " << element.first << " ";
      dim_domain_y = element.second.size();
      for(int i = 0; i < dim_domain_y ; ++i) ofs << element.second[i] << " ";
    }
    ofs << std::endl;
    ofs << "#output_form ";
    for(auto element : domain_x){
      ofs << element.first << " ";
      domain_x_name = element.first;
      dim_domain_x = element.second.size();
    }
    for(int i = 0 ; i < dim_domain_y; ++i){
      for(auto key : physical_value_list){
        ofs << key << " " ;
      }
    } 
    ofs << std::endl;
    for(int i = 0 ; i < dim_domain_x; ++i){
      ofs << domain_x[domain_x_name][i] << " ";
      for(int j = 0 ; j < dim_domain_y; ++j){
        for(auto key : physical_value_list) ofs << physical_values[key][j][i] << " " ;
      }
      ofs << std::endl;
    }
    ofs.close();
  }


protected:
 std::string condidat_;
 std::ofstream ofs_;
}; //end dataput class


}//namespace end

#endif // DATARECORDER_HPP

