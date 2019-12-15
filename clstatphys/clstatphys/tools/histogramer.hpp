#ifndef TOOLS_HISTOGRAMER_HPP 
#define TOOLS_HISTOGRAMER_HPP 

#include <vector>
#include <numeric>
#include <random>

namespace tools{

class Histogramer{
public:
  Histogramer(int length=1) : length_(length),hist_(length,0),range_(length,0) 
  {v_max = v_min = v_central = v_variance =  dv_ = 0.0;}

  void initialize(int length){
    hist_.resize(length);
    range_.resize(length);
    length_ = length;
    v_max = v_min = v_central = v_variance = dv_ = 0.0;
    for(int i = 0 ; i < length ; ++i){
      hist_[i] = 0.0;
      range_[i] = 0.0;
    }
  }

  void operator()(std::vector<double> const& v,int const num){
    set_range(v,num);
    distribute(v,num);
  }

  void set_range(std::vector<double> const& v, int const  num) {
    v_central = 0.0;
    v_variance = 0.0;
    for( int i = 0 ; i < num; ++i) v_central += v[i] /(double)num;
    for( int i = 0 ; i < num; ++i) v_variance += v[i] * v[i]/(double)num;
    v_variance = ( v_variance - v_central * v_central) / ((double)num - 1) * (double)num;
    v_max = max(v,num);
    v_min = min(v,num);
    if( (v_max - v_central) / (v_min - v_central) < -1) v_min = -1.0 * (v_max - v_central) + v_central;
    else v_max = -1.0 * (v_min - v_central) + v_central;
    dv_ = (v_max - v_min) / ((double)range_.size()-1);
    if( dv_ > 0.1){
      dv_ = 0.1;
      for( int i = 0; i != range_.size(); ++i) range_[i] = dv_ * ((double)i - 0.5 * (double)(range_.size()-1));
    }
    else{
      for( int i = 0; i != range_.size(); ++i) range_[i] = v_min + dv_ * (double)i;
    }
    v_min = range_[0];
    v_max = range_[range_.size()-1];
  }

  void distribute(std::vector<double> const& v, int const num){
    double fractpart, intpart;
    int position = 0;
    for( int i = 0; i < num ; ++i){
      fractpart = std::modf( (v[i] - v_min) / dv_, &intpart); 
      if(fractpart >= 0.5) fractpart = 1.0;
      else fractpart = 0.0;
      position = (int)(intpart + fractpart);
      if(position <= hist_.size()) hist_[position] += 1;
    }
  }

  void output(std::vector<double>& hist, std::vector<double>& range){
    for(int i = 0; i < length_; ++i){
      hist[i] = (double)hist_[i];
      range[i] = range_[i];
    }
  }

  double normalize_const(){
    double normalize_const = 0.0;
    for(int i = 0 ; i < length_; ++i){
      normalize_const += hist_[i] * dv_;
    }
    return normalize_const;
  }

  int size() {return length_;}
  
  double max(std::vector<double> const& v, const int num){
    double temp = 0;
    for(int i = 0; i < num ; ++i){
      if(temp < v[i]) temp = v[i];
    }
    return temp;
  }

  double min(std::vector<double> const& v,const int num){
    double temp = 0.0;
    for(int i = 0; i < num ; ++i){
      if(temp > v[i]) temp = v[i];
    }
    return temp;
  }

  double min(){return v_min;}
  double max(){return v_max;}
  double interval(){return dv_;}
  double central(){return v_central;}
  double variance(){return v_variance;}
  double confirm(){return v_max - v_central + v_min - v_central;}
  double total_count(){return std::accumulate(hist_.begin(),hist_.end(),0);}

  void reset(){
    v_max = v_min = v_central = dv_ = 0.0;
    for(int i = 0; i < length_ ; ++i){
      range_[i] = 0.0;
      hist_[i] = 0;
    }
  }

protected:
  mutable double v_max, v_min, v_central, dv_, v_variance ; 
  mutable std::vector<double> range_ ;
  mutable std::vector<int> hist_ ;
  int length_;

};

} // end namespace stat

#endif //TOOLS_HISTOGRAMER_HPP 
