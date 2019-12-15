#ifndef METADYNAMICS_WANGLANDAU_PARALLELIZED_HPP
#define METADYNAMICS_WANGLANDAU_PARALLELIZED_HPP

#include<random>
#include<vector>
#include<algorithm>
#include<unordered_map>
#include<string>
#include<memory>

namespace metadynamics {

class WangLandauParallelized {
public:
  static std::string name() { return "metadynamics : Wang Landau "; }
  WangLandauParallelized(unsigned int dim, unsigned int N_bin, double E_max, double E_min) : dim_(dim), E_min_(E_min), E_max_(E_max), N_bins_(N_bin),density_(N_bin,0.0),entropy_(N_bin,0.0),bins_(N_bin,0.0){bin_ = (E_max_ - E_min_) / (double)N_bins_;}

  template <typename Func, typename Proposer>
  void montecalro(Func const& f,  Proposer const& proposer, const int max_iteration, std::mt19937 & mt, std::vector<std::shared_ptr<std::mt19937> > & mts, int mc_check = 0){
    if(mc_check == 0) mc_check = N_bins_;
    for(int i = 0 ; i < N_bins_; ++i) bins_[i] = (double)i * bin_ + E_min_;
    std::vector<double> z_current(dim_),z_proposed(dim_);
    std::vector<double> hist(N_bins_,0);
    std::uniform_real_distribution<> accept_rand(0.0,1.0);
    double precision = 1.0;
    int current_bin = propose_bin_z(f,proposer,mts,z_current);
    int mc_step = 0;
    int counter = 0;
    while(mc_step < max_iteration){
     // std::cout << "counter : " << counter << " mc_step : " << mc_step << std::endl;
      int proposed_bin = propose_bin_z(f,proposer,mts,z_proposed);
      if( accept_rand(mt) < std::exp(entropy_[current_bin] - entropy_[proposed_bin]) ){
        z_current = z_proposed;
        current_bin = proposed_bin;
      }
      entropy_[current_bin] += precision;
      hist[current_bin] += 1;
      counter++;
      if(counter % mc_check == 0){  
        //std::cout << " check hist !! " << std::endl;
        double max_hist = *std::max_element(hist.begin(),hist.end());
        double temp = max_hist;
        for(int i = 0 ; i < N_bins_ ; ++i){
          if(hist[i]==0) continue;
          temp = std::min(temp, hist[i]);
        }
        //std::cout << " temp : " << temp << " max_hist : " << max_hist << std::endl; 
        if(1.0  - temp /max_hist  < 1e-2){
          for(int i = 0; i < N_bins_ ; ++i) hist[i] = 0;
          mc_step++;
          std::cout << " mc_step : " << mc_step << std::endl;
          precision = 1.0 / (double)(mc_step+1); 
        }
      }
    }
    set_result();
  }

  template <typename Proposer,typename Func>
  int propose_bin_z(Func const& f, Proposer const& proposer, std::vector<std::shared_ptr<std::mt19937> > & mts, std::vector<double>& z){
    proposer(z,mts);  
    double E = f(z);
    if(E < E_min_ || E > E_max_){
      proposer(z, mts);  
      E = f(z);
    }
    int bin_pos = std::distance(bins_.begin(), std::lower_bound(bins_.begin(), bins_.end(), E));
    return bin_pos;
  }

  void set_result(){
    double max_entropy = *std::max_element(entropy_.begin(),entropy_.end());
    std::vector<double> density(N_bins_);
    double c_normalize = 0.0;
    for(int i = 0 ; i < N_bins_ ; ++i){ 
      if(entropy_[i] == 0) density_[i] = 0.0;
      else density_[i] = std::exp(entropy_[i] - max_entropy);
      c_normalize += density_[i] * bin_;
    }
    for(int i = 0 ; i < N_bins_ ; ++i){ 
      density_[i] /= c_normalize; 
      if(density_[i] > 0.0) entropy_[i] = std::log(density_[i]);
    }
  }

  std::vector<double> density(){return density_;}
  std::vector<double> bins(){ return bins_;}
  std::vector<double> entropy(){return entropy_;}

protected:
  unsigned int N_bins_, dim_;
  double E_min_, E_max_,bin_;
  std::vector<double> entropy_,bins_,density_;
};

} //end namespace

#endif // METADYNAMICS_WANGLANDAU_PARALLELIZED_HPP
