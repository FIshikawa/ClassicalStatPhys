#ifndef METADYNAMICS_WANGLANDAU_2D_HPP
#define METADYNAMICS_WANGLANDAU_2D_HPP

#include<random>
#include<vector>
#include<algorithm>
#include<unordered_map>
#include<string>
#include<memory>

namespace metadynamics {

class WangLandau2D {
public:
  static std::string name() { return "metadynamics : Wang Landau 2D"; }
  WangLandau2D(unsigned int dim, unsigned int N_bin, std::string E_name, double E_max, double E_min, std::string M_name, double M_max, double M_min) : 
  dim_(dim), E_min_(E_min), E_max_(E_max), M_min_(M_min), M_max_(M_max), N_bins_(N_bin), entropy_(N_bin*N_bin,0.0), bins_M_(N_bin,0.0), bins_E_(N_bin,0.0), 
  E_name_(E_name), M_name_(M_name){ bin_E_ = (E_max_ - E_min_) / (double)N_bin ; bin_M_ = (M_max_ - M_min_) / (double)N_bin;}

  template <typename Func, typename Proposer>
  void montecalro(Func const& E_calc, Func const& M_calc, Proposer const& proposer, const int max_iteration, std::mt19937 & mt, int mc_check = 0){
    if(mc_check == 0) mc_check = N_bins_;
    for(int i = 0 ; i < N_bins_; ++i) bins_E_[i] = (double)i * bin_E_ + E_min_;
    for(int i = 0 ; i < N_bins_; ++i) bins_M_[i] = (double)i * bin_M_ + M_min_;
    std::vector<double> z_current(dim_),z_proposed(dim_);
    std::vector<double> hist(N_bins_ * N_bins_ + N_bins_, 0.0);
    std::uniform_real_distribution<> accept_rand(0.0,1.0);
    double precision = 1.0;
    int current_bin = propose_bin_z(E_calc,M_calc,proposer,mt,z_current);
    int mc_step = 0;
    int counter = 0;
    while(mc_step < max_iteration){
     // std::cout << "counter : " << counter << " mc_step : " << mc_step << std::endl;
      int proposed_bin = propose_bin_z(E_calc,M_calc,proposer,mt,z_current);
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
  int propose_bin_z(Func const& E_calc, Func const& M_calc, Proposer const& proposer, std::mt19937 & mt, std::vector<double>& z){
    proposer(z,mt);  
    double E = E_calc(z);
    double M = M_calc(z);
    if(E < E_min_ || E > E_max_){
      proposer(z, mt);  
      E = E_calc(z);
    }
    if(M < M_min_ || M > M_max_){
      proposer(z, mt);  
      M = M_calc(z);
    }
    int bin_pos_E = std::distance(bins_E_.begin(), std::lower_bound(bins_E_.begin(), bins_E_.end(), E));
    int bin_pos_M = std::distance(bins_M_.begin(), std::lower_bound(bins_M_.begin(), bins_M_.end(), M));
    int bin_pos = (N_bins_ -1)* bin_pos_E + bin_pos_M;
    return bin_pos;
  }

  void set_result(){
    density_ = std::vector< std::vector<double> >(N_bins_ , std::vector<double>(N_bins_, 0.0));
    double max_entropy = *std::max_element(entropy_.begin(),entropy_.end());
    std::vector< std::vector<double> > density(N_bins_, std::vector<double>(N_bins_, 0.0));
    double c_normalize = 0.0;
    for(int i = 0 ; i < N_bins_ ; ++i){ 
      for(int j = 0 ; j < N_bins_ ; ++j){ 
        if(entropy_[i*(N_bins_ -1) + j] == 0) density_[i][j] = 0.0;
        else density_[i][j] = std::exp(entropy_[i*(N_bins_ -1) + j] - max_entropy);
        c_normalize += density[i][j] * bin_E_ * bin_M_;
      }
    }
    for(int i = 0 ; i < N_bins_ ; ++i){ 
      for(int j = 0 ; j < N_bins_ ; ++j){ 
        density_[i][j] /= c_normalize; 
        if(density_[i][j] > 0.0) entropy_[i*(N_bins_-1)+j] = std::log(density_[i][j]);
      }
    }
  }

  std::vector< std::vector<double> > density(){
    return density_;
  }

  std::unordered_map<std::string, int > elements_name(){
    std::unordered_map<std::string, int > elements_name;
    elements_name[E_name_] = 0;
    elements_name[M_name_] = 1;
    return elements_name;
  }

  std::unordered_map<std::string, std::vector<double> >bins(){
    std::unordered_map<std::string, std::vector<double> >  bins;
    bins[E_name_] = bins_E_;
    bins[M_name_] = bins_M_;
    return bins;
  }

  std::vector<std::vector<double> > entropy(){
    std::vector< std::vector<double> > entropy(N_bins_, std::vector<double>(N_bins_, 0.0));
    for(int i = 0 ; i < N_bins_ ; ++i){ 
      for(int j = 0 ; j < N_bins_ ; ++j){ 
        entropy[i][j] = entropy_[i*(N_bins_ -1) + j];
      }
    }
    return entropy;
  }

protected:
  unsigned int N_bins_, dim_;
  double E_min_, E_max_, bin_E_, M_min_, M_max_, bin_M_;
  std::string E_name_, M_name_;
  std::vector<double> entropy_, bins_E_, bins_M_;
  std::vector< std::vector<double> > density_;
};//end WangLandau

} //end namespace

#endif //METADYNAMICS_WANGLANDAU_2D_HPP
