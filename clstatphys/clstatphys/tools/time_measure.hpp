#ifndef TIME_MEASURE_HPP
#define TIME_MEASURE_HPP

#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

namespace tools {

class TimeMeasure{
public:
  TimeMeasure(double t, double dt, int N_time_measure, std::string scale = "linear", int order = 0) : 
    t_(t), dt_(dt), N_time_measure_(N_time_measure), scale_(scale), order_(order){reset();initialize();}

  void reset(){
    counter_ = 0;
    int order_temp = std::floor(std::log10(t_)+1);
    if(order_ == 0) order_ = order_temp;
  }

  void initialize(){
    check_time_ = std::vector<double>();
    if(scale_ == "log"){
      if(t_ / std::pow(10,order_-1) < dt_){
        std::cerr << "t / 10^(order-1) should be larger than dt " << std::endl;
        std::exit(11);
      }
      double ten_power = std::pow(10,order_-1);
      double t_now_max = t_ / ten_power;
      double now_dt = t_now_max / N_time_measure_; 
      double t_temp = 0.0; 
      for(int order = order_; order > 0; --order){
        ten_power = std::pow(10,order-1);
        t_now_max = t_ / ten_power;
        now_dt = t_now_max / N_time_measure_; 
        if(now_dt/dt_ != std::floor(now_dt/dt_)){
          std::cerr << "order should be divisible by dt" << std::endl;
          std::exit(11);
        }
        check_time_.push_back(t_temp);
        for(int i = static_cast<int>(t_temp/now_dt); i < N_time_measure_-1; ++i){
          t_temp += now_dt;
          check_time_.push_back(t_temp);
        }
        t_temp = t_now_max; 
      }
    }
    else if(scale_ == "linear"){
      double now_dt = t_ / N_time_measure_;
      double t_temp = 0; 
      check_time_.push_back(t_temp);
      for(int i = 0; i < N_time_measure_-1; ++i){
        t_temp += now_dt;
        check_time_.push_back(t_temp);
      }
    }
    else{
      std::cerr << "sclae should be linear or log" << std::endl;
      std::exit(11);
    }
    N_total_data_ = check_time_.size();
  }

  int number_of_total_data(){
    return N_total_data_;
  }

  bool check(int step){
    int now_step = static_cast<int>(std::round(check_time_[counter_] / dt_));
    if(counter_ < N_total_data_ && step == now_step){
      ++counter_;
      return true;
    }
    else return false;
  }

  double now_time(){if(counter_ < N_total_data_) return check_time_[counter_]; else return 0.0;}
  double now_step(){if(counter_ < N_total_data_) 
    return static_cast<int>(std::round(check_time_[counter_]/dt_)); else return 0.0;}

  double operator()(){
    ++counter_;
    return check_time_[counter_-1];
  }
  int order(){return order_;}
  int counter(){return counter_;}
  std::string plot_scale(){return scale_;}

private:
  double t_, dt_;
  std::vector<double> check_time_;
  int counter_, order_, N_time_measure_, N_total_data_;
  std::string scale_;
}; //end dataput class


}//namespace end

#endif // 
