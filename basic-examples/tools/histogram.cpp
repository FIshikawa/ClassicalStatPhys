#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<random>
#include<iomanip>

#include<clstatphys/tools/histogram.hpp>

int main() {
  std::cout << "[test Histogram]\n";

  //declare class 
  int n_bin = 100;
  std::cout << "[declare class]" << std::endl;
  tools::Histogram histogram(n_bin);
  //end declare class

  //generate data
  int N_loop = 1E+06;
  std::vector<double> data(N_loop);
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  std::normal_distribution<> dist(0.0, 1.0);
  for(int i = 0; i < N_loop; ++i) data[i] = dist(mt); 
  //end generate data

  //create histogram
  histogram(data);
  //end create histogram

  //output results
  std::cout << "[result put]" << std::endl;
  //set precision
  std::cout << std::scientific;
  std::cout << std::setprecision(15);


  std::vector<double> bins(n_bin), hist(n_bin);
  //output hist and bins
  histogram.output(hist,bins);
  //output normalize const
  double normalize_const = histogram.normalize_const(); 
  auto exact = [](double x) {return 1.0/std::sqrt(2.0 * M_PI) * exp(-0.5*x*x);};
  double c_norm = 0.0;
  double variance = 0.0;
  double mean = 0.0;
  double width = histogram.interval();
  std::cout << "#output bins hist exact" << std::endl;
  for(int i = 0 ; i < n_bin ; ++i){
    std::cout << bins[i] << " " << hist[i]/normalize_const << " " << exact(bins[i]) << std::endl;
    c_norm += hist[i] * width;
    variance += hist[i]/normalize_const * bins[i] * bins[i] * width;
    mean += hist[i]/normalize_const * bins[i] * width;
  }
  std::cout << "# N_loop : 1E+06" << std::endl;
  std::cout << "# n_bin  : 100" << std::endl;
  std::cout << "# bin width : " << width << std::endl;
  std::cout << "# variance " << std::endl ;
  std::cout << "#   from histogram : " << variance << std::endl
            << "#   from count     : " << histogram.variance() << std::endl
            << "#   unbiassed      : " << histogram.unbiassed_variance()<< std::endl 
            << "#   exact          : " << 1.0 << std::endl; 
  std::cout << "# mean " << std::endl
            << "#   from histogram : " <<  mean << std::endl
            << "#   from count     : " <<  histogram.mean() << std::endl
            << "#   exact          : " << 0.0 << std::endl; 
  //end out put results

  std::cout << "[finish test]" << std::endl;

 return 0;
}
