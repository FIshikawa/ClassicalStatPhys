// Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "accumulator.hpp"
#include "tester_dissipated.hpp"
#include "normal_distribution.hpp"
#include "initializer.hpp"
#include "hamiltonian.hpp"
#include "Lattice.hpp"
#include "dissipated_integrator.hpp"
#include "dataput.hpp"
#include "histogramer.hpp"

#include <boost/random.hpp>
#include <iostream>

typedef Lattice::FullyConnected lattice_t;
typedef initializer::equilLangevin initializer_t;
//typedef dissipated_integrator::rk2 integrator_rk2;
//typedef dissipated_integrator::euler integrator_euler;
typedef dissipated_integrator::rk2_cauchy integrator_rk2;
typedef dissipated_integrator::euler_cauchy integrator_euler;
typedef hamiltonian::freeparticle hamiltonian_t;
typedef dataput::dataputForRegularCube dataput_t;
typedef stat::histogramer_particle histogramer_t;

int main(int argc, char **argv){
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  unsigned long mmax = 1024;
  std::size_t seed = 12345;
  int counter = 0;
  unsigned int num_particles = 1024;
  double t = 1.0;
  int dN_t  = 10;
  int N_loop  = 10;
  int N_t  = 10;
  char* condition_dat = "condi.dat";
  std::string result_dirr("./test");

  double J = 0.0;
  double T = 1.0;
  double gamma = 1.0;
  double pt = 0.0;
  unsigned int length = num_particles;

 //parameter set 
 if (argc > 1) result_dirr =      boost::lexical_cast<std::string>(argv[1]);
 if (argc > 2) num_particles =    boost::lexical_cast<int>(argv[2]);
 if (argc > 3) N_loop =           boost::lexical_cast<int>(argv[3]);
 if (argc > 4) dN_t =             boost::lexical_cast<int>(argv[4]);
 if (argc > 5) N_t =              boost::lexical_cast<int>(argv[5]);
 if (argc > 6) length =           boost::lexical_cast<int>(argv[6]);
 if (argc > 7) T =                boost::lexical_cast<double>(argv[7]);
 if (argc > 8) gamma =            boost::lexical_cast<double>(argv[8]);

//parameter set 
 std::vector<std::string> experimental_condition(4);
 experimental_condition[0] = "Dissipation test";
 experimental_condition[1] = hamiltonian_t::name();
 experimental_condition[2] = initializer_t::name();
 experimental_condition[3] = lattice_t::name();
 dataput_t::startput(condition_dat, experimental_condition);

//condition declare
 std::ofstream ofs;
 ofs.open(condition_dat,std::ios::app);
 if(!ofs)std::cerr<< condition_dat <<"opening failed"<<std::endl;
 ofs  << "Temperture: T = " << T << std::endl
      << "Friction Coefficient: gamma = " << gamma << std::endl
      << "Number of patricles :num_particles = " << num_particles << std::endl
      << "Developing ime : t = " << t << std::endl
      << "Start Number of time steps : N_time = " << N_t << std::endl
      << "Number of loop: N_loop =" << N_loop << std::endl
      << "Difference of number of time: dN_t =" << dN_t << std::endl
      << "Result into : " << result_dirr << std::endl;
 ofs.close();

  //dataputer set 
  dataput_t dataputer(1,num_particles);

  boost::mt19937 mt(seed);

  lattice_t lattice(num_particles);
  unsigned int N_adj;// number of adjacent spins
  N_adj = lattice.number_adjacent();

  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
  lattice.create_table(pair_table);

  stat::normal_distribution dist(0.0, T);

  std::vector<double> z(2*num_particles);
  std::vector<double> r(2*num_particles);
  std::vector< std::vector<double> >  
		              result_rk2(20, std::vector<double>(N_loop)),
		              result_euler(20, std::vector<double>(N_loop));

  hamiltonian_t ham(num_particles,J,pair_table,N_adj);

  histogramer_t histogram_euler(length);
  histogramer_t histogram_rk2(length);

  integrator_rk2   rk2(T, gamma, 2*num_particles);
  integrator_euler euler(T, gamma, 2*num_particles);

  initializer_t initializer(J,T,num_particles);

  initializer.coldset(z,counter,mt);
  initializer.coldset(r,counter,mt);
  
  double* v_rk2 = &z[num_particles];
  double* v_euler = &r[num_particles];

  for(int loop = 0; loop < N_loop; ++loop){
    double dt = t/(double)N_t;
    std::cout << "loop : " << loop << " ,dt : " << dt << std::endl;
    for (int i = 0; i < 2*num_particles; ++i) z[i] = 0.0;
    for (int i = 0; i < 2*num_particles; ++i) r[i] = 0.0;
    pt = 0.0;
    for(int step = 0; step < N_t; ++step){
      pt += dt;
      rk2.step(pt, dt, z, ham, mt);
      euler.step(pt, dt, r, ham, mt);
    }
    test4(dist,z,result_rk2,num_particles,loop); 
    test4(dist,r,result_euler,num_particles,loop); 
    N_t = N_t + dN_t;
  }

  histogram_euler(r,num_particles);
  histogram_rk2(z,num_particles);
  std::cout << "[max] rk2 : " << histogram_rk2.max() << " euler : " << histogram_euler.max() << std::endl;
  std::cout << "[min] rk2 : " << histogram_rk2.min() << " euler : " << histogram_euler.min() << std::endl;
  std::cout << "[central] rk2 : " << histogram_rk2.central() << " euler : " << histogram_euler.central() << std::endl;
  std::cout << "[confirm] rk2 : " << histogram_rk2.confirm()  << " euler : " << histogram_euler.confirm() << std::endl;
  std::cout << "[count] rk2 : " << histogram_rk2.total_count()  << " euler : " << histogram_euler.total_count() << std::endl;

  stat::accumulator accum_rk2(integrator_rk2::name());
  stat::accumulator accum_euler(integrator_euler::name());
  for(int i = 0; i < num_particles ; ++i){
    accum_rk2 << z[i + num_particles];
    accum_euler << r[i + num_particles];
  }
  std::cout << "[accumulator central] rk2 : " << accum_rk2.average() << " euler : " << accum_euler.average() << std::endl;
  std::cout << "[accumulator error] rk2 : " << accum_rk2.error() << " euler : " << accum_euler.error() << std::endl;

  double* observable_pointer_rk2[20];
  double* observable_pointer_euler[20];
  std::vector<std::string> observable_list(20);
  
  int tag = 0;
  for(int i = 0; i < 20; ++i){
    observable_pointer_rk2[tag] = &result_rk2[tag][0];
    observable_pointer_euler[tag] = &result_euler[tag][0];
    tag += 1;
  }
  test_result_set(observable_list);

  std::vector<double> hist_euler(length);
  std::vector<double> range_euler(length);
  std::vector<double> hist_rk2(length);
  std::vector<double> range_rk2(length);
  
  histogram_euler.output(hist_euler,range_euler);
  histogram_rk2.output(hist_rk2,range_rk2);
  double* result_hist[4];
  std::vector<std::string> hist_list(4);
  result_hist[0] = &range_rk2[0];
  result_hist[1] = &hist_rk2[0];
  hist_list[0] = "bin of rk2";
  hist_list[1] = "hist of rk2";
  result_hist[2] = &range_euler[0];
  result_hist[3] = &hist_euler[0];
  hist_list[2] = "bin of euler";
  hist_list[3] = "hist of euler";

  std::string result_dat_rk2("/result_rk2.dat");
  std::string result_dat_euler("/result_euler.dat");
  std::string result_dat_hist("/result_hist.dat");
  result_dat_hist = result_dirr + result_dat_hist;
  result_dat_rk2 = result_dirr + result_dat_rk2;
  result_dat_euler = result_dirr + result_dat_euler;
  dataputer.output_result(result_dat_hist.c_str(),condition_dat,"-1bin",1,length,result_hist,hist_list);
  dataputer.output_result(result_dat_rk2.c_str(),condition_dat,"dN_t",dN_t,N_loop,observable_pointer_rk2,observable_list);
  dataputer.output_result(result_dat_euler.c_str(),condition_dat,"dN_t",dN_t,N_loop,observable_pointer_euler,observable_list);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
