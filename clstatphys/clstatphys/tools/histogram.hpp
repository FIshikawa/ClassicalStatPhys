#ifndef TOOLS_HISTOGRAM_HPP 
#define TOOLS_HISTOGRAM_HPP 

#include <vector>
#include <numeric>
#include <random>

namespace tools{

class Histogram{
public:
  Histogram(int length=1) : length_(length),hist_(length,0),range_(length,0)
    {max_value_ = 0.0;min_value_ = 0.0;mean_ = 0.0;variance_ = 0.0;bin_width_ = 0.0;}

  void initialize(int length){
    hist_.resize(length);
    range_.resize(length);
    length_ = length;
    max_value_ = min_value_ = mean_ = variance_ = bin_width_ = 0.0;
    for(int i = 0 ; i < length ; ++i){
      hist_[i] = 0.0;
      range_[i] = 0.0;
    }
  }

  void operator()(std::vector<double> const& data){
    int sample_size = data.size();
    set_range(data,sample_size);
    distribute(data,sample_size);
  }

  void set_range(std::vector<double> const& data, const int sample_size) {
    for( int i = 0 ; i < sample_size; ++i) mean_ += data[i] / sample_size;
    for( int i = 0 ; i < sample_size; ++i) variance_ += data[i] * data[i] / sample_size;
    unbiassed_variance_ = ( variance_ - mean_ * mean_) / (sample_size- 1.0) * sample_size;
    max_value_ = max(data,sample_size);
    min_value_ = min(data,sample_size);
    if( (max_value_ - mean_) / (min_value_ - mean_) < -1) min_value_ = -1.0 * (max_value_ - mean_) + mean_;
    else max_value_ = -1.0 * (min_value_ - mean_) + mean_;
    bin_width_ = (max_value_ - min_value_) / ( - 1.0 + range_.size());
    if( bin_width_ > 0.1){
      bin_width_ = 0.1;
      for( int i = 0; i < range_.size(); ++i) range_[i] = bin_width_ * ( i - 0.5 * (range_.size()-1.0));
    }
    else{
      for( int i = 0; i < range_.size(); ++i) range_[i] = min_value_ + bin_width_ * i;
    }
    min_value_ = range_[0];
    max_value_ = range_[range_.size()-1];
  }

  void distribute(std::vector<double> const& data, const int sample_size){
    double fractpart, intpart;
    int position = 0;
    for( int i = 0; i < sample_size ; ++i){
      fractpart = std::modf( (data[i] - min_value_) / bin_width_, &intpart); 
      if(fractpart >= 0.5) fractpart = 1.0;
      else fractpart = 0.0;
      position = static_cast<int>(intpart + fractpart);
      if(position <= hist_.size()) hist_[position] += 1;
    }
  }

  void output(std::vector<double>& hist, std::vector<double>& range){
    for(int i = 0; i < length_; ++i){
      hist[i] = static_cast<double>(hist_[i]);
      range[i] = range_[i];
    }
  }

  double normalize_const(){
    double normalize_const = 0.0;
    for(int i = 0 ; i < length_; ++i){
      normalize_const +=  bin_width_ * hist_[i];
    }
    return normalize_const;
  }

  int size() {return length_;}
  
  double max(std::vector<double> const& data, const int sample_size){
    double temp = 0;
    for(int i = 0; i < sample_size ; ++i){
      if(temp < data[i]) temp = data[i];
    }
    return temp;
  }

  double min(std::vector<double> const& data,const int sample_size){
    double temp = 0.0;
    for(int i = 0; i < sample_size ; ++i){
      if(temp > data[i]) temp = data[i];
    }
    return temp;
  }

  double min(){return min_value_;}
  double max(){return max_value_;}
  double interval(){return bin_width_;}
  double mean(){return mean_;}
  double variance(){return variance_;}
  double unbiassed_variance(){return unbiassed_variance_;}
  double confirm(){return max_value_ - mean_ + min_value_ - mean_;}
  double total_count(){return std::accumulate(hist_.begin(),hist_.end(),0);}

  void reset(){
    max_value_ = min_value_ = mean_ = bin_width_ = 0.0;
    for(int i = 0; i < length_ ; ++i){
      range_[i] = 0.0;
      hist_[i] = 0;
    }
  }

private:
  mutable double max_value_, min_value_, mean_, bin_width_, variance_, unbiassed_variance_; 
  mutable std::vector<double> range_ ;
  mutable std::vector<int> hist_ ;
  int length_;

};

} // end namespace stat

#endif //TOOLS_HISTOGRAM_HPP 
