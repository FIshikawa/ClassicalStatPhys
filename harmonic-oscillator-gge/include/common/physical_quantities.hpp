#ifndef PHYISCAL_QUANTITIES_HPP
#define PHYISCAL_QUANTITIES_HPP

#include <unordered_map>
#include <clstatphys/tools/accumulator.hpp>
#include <clstatphys/physics/normalmode_energy_periodic_fftw.hpp>

using PhysicalQuantities1d = std::unordered_map<std::string, std::vector<tools::Accumulator> >; 
using PhysicalQuantities2d = std::unordered_map<std::string, std::vector<std::vector<tools::Accumulator> > >;

struct PhysicalQuantities{
  int num_particles;
  int N_total_data;
  int Ns_observe;
  double v_init;
  double x_init;
  std::vector<std::string> quantities_1d_list = {
                                               "Velocity",
                                               "Velocity_total",
                                               "Displacement",
                                               "Displacement_total",
                                               "Energy_total",
                                               "Kinetic_Energy",
                                               "Potential_Energy",
                                               "SpectralEntropy",
                                               "EffectiveFractionOfModes",
                                               "DisplacementACF",
                                               "VelocityACF"
                                               };
  std::vector<std::string> quantities_2d_list = {
                                               "NormalModeEnergy",
                                               "SpaceCorrelation",
                                               "GeneratorCorrelation"
                                               };

  PhysicalQuantities1d quantities_1d;
  PhysicalQuantities2d quantities_2d;

  Hamiltonian hamiltonian;

  template<typename Settings>
  PhysicalQuantities(Settings & settings){
    N_total_data = settings.N_total_data;
    num_particles = settings.num_particles;
    Ns_observe = settings.Ns_observe;

    for(const auto& key : quantities_1d_list){
      quantities_1d[key].resize(N_total_data);
      for(int step = 0; step < N_total_data; ++step) quantities_1d[key][step].reset();
    }
    for(const auto& key : quantities_2d_list){
      quantities_2d[key].resize(N_total_data);
      for(int step = 0 ; step < N_total_data ; ++step){
        quantities_2d[key][step] = std::vector<tools::Accumulator>(Ns_observe);
        for(int i = 0; i < Ns_observe; ++i) quantities_2d[key][step][i].reset();
      }
    }
    hamiltonian = settings.hamiltonian();
  }


  void measure(std::vector<double> z, int & measure_step){
    int observe = (num_particles-1) / 2;
    double *x = &z[0];
    double *v = &z[num_particles];
    double pt = 0.0;
    double v_total,x_total,sum_spectral,average_spectral,spectral_entropy;
    std::vector<double> normalmode_energy_temp(num_particles);

    v_total = 0.0;
    x_total = 0.0;
    for(int i = 0; i < num_particles; ++i){
      v_total += v[i] / num_particles;
      x_total += x[i] / num_particles;
    } 
    //set normalmode energy
    NormalModeEnergyPeriodicFFTW(z,normalmode_energy_temp);

    sum_spectral = 0;
    spectral_entropy= 0;
    average_spectral= 0;
    for(int i = 0; i < num_particles; ++i) sum_spectral += normalmode_energy_temp[i];
    for(int i = 0; i < num_particles; ++i){
      average_spectral = normalmode_energy_temp[i] / sum_spectral;
      if(average_spectral > 0) spectral_entropy -=  average_spectral * std::log(average_spectral);
    }

    //input data 
    quantities_1d["Velocity"][measure_step] << v[observe];
    quantities_1d["Velocity_total"][measure_step] << v_total;
    quantities_1d["Displacement"][measure_step] << x[observe];
    quantities_1d["Displacement_total"][measure_step] << x_total;
    quantities_1d["Energy_total"][measure_step] <<  hamiltonian.energy(pt,z) / num_particles;
    quantities_1d["Kinetic_Energy"][measure_step] << hamiltonian.kinetic_energy(pt,z) / num_particles;
    quantities_1d["Potential_Energy"][measure_step] << hamiltonian.potential_energy(pt,z) / num_particles;
    quantities_1d["SpectralEntropy"][measure_step] << spectral_entropy;
    quantities_1d["EffectiveFractionOfModes"][measure_step] << std::exp(spectral_entropy) / num_particles;


    if(measure_step == 0){
      v_init = v[observe] - v_total;
      x_init = x[observe];
    }
    quantities_1d["VelocityACF"][measure_step] << v_init * (v[observe] - v_total);
    quantities_1d["DisplacementACF"][measure_step] << x_init * x[observe];

    for(int i = 0; i < Ns_observe; ++i) {
      quantities_2d["NormalModeEnergy"][measure_step][i] << normalmode_energy_temp[i];
      quantities_2d["SpaceCorrelation"][measure_step][i] << x[observe] * x[i];
      quantities_2d["GeneratorCorrelation"][measure_step][i] << 0.5 * (v[observe] * v[observe] - x[observe] * x[observe]) *
                                                                0.5 * (v[i] * v[i] - x[i] * x[i]);
    }
    ++measure_step;
  }


  template<typename Settings, typename Dataput>
  void output(Settings & settings, Dataput & dataput){

    std::unordered_map<std::string, std::vector<std::string> > values_name;
    std::unordered_map<std::string, std::string> file_name;

    // normal results define 
    file_name["normal"] = settings.result_directory + "result.dat";
    std::unordered_map<std::string, std::vector<double> > results_normal;
    
    for(const auto& key  : quantities_1d_list){
      values_name["normal"].push_back(key);
      values_name["normal"].push_back("error_" + key);
      results_normal[key].resize(N_total_data);
      results_normal["error_" + key].resize(N_total_data);
      for(int step = 0; step < N_total_data; ++step){
        results_normal[key][step] = quantities_1d[key][step].mean();
        results_normal["error_" + key][step] = quantities_1d[key][step].error();
      }
    }
    std::unordered_map< std::string, std::vector<double> > domain_time;
    domain_time["time"] = std::vector<double>(N_total_data);
    tools::TimeMeasure time_measure = settings.time_measure();
    for(int step = 0; step < N_total_data; ++step) domain_time["time"][step] = time_measure();

    dataput.time_tag() << "start writing : normal" << std::endl;
    dataput.output_result(file_name["normal"], domain_time, results_normal);
    dataput.time_tag() << "finish writing : normal" << std::endl;


    // 2dim resutls set 
    std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d;
    std::unordered_map<std::string, std::unordered_map< std::string, std::vector<double> > > domain_2d;

    for(const auto& key  : quantities_2d_list){
      values_name[key].push_back(key);
      values_name[key].push_back("error_"+key);

      results_2d[key] 
          = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
      results_2d["error_" + key] 
          = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));

      for(int step = 0; step < N_total_data; ++step){
        for(int site = 0; site < Ns_observe; ++site){
          results_2d[key][step][site] = quantities_2d[key][step][site].mean();
          results_2d["error_" + key][step][site] = quantities_2d[key][step][site].error();
        }
      }
    }

    values_name["NormalModeEnergy"] 
      = std::vector<std::string>{"NormalModeEnergy","error_NormalModeEnergy",
                                 "NormalModeEnergyVariance","error_NormalModeEnergyVariance"};
    results_2d["NormalModeEnergyVariance"] 
        = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
    results_2d["error_NormalModeEnergyVariance"] 
        = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));

    for(int step = 0; step < N_total_data; ++step){
      for(int site = 0; site < Ns_observe; ++site){
        results_2d["NormalModeEnergyVariance"][step][site] 
            = quantities_2d["NormalModeEnergy"][step][site].variance();
        results_2d["error_NormalModeEnergyVariance"][step][site] 
          = quantities_2d["NormalModeEnergy"][step][site].variance_error();
      }
    }

    // Eigenvalues
    file_name["SpaceCorrelation"] = settings.result_directory + "result_space_correlation.dat";
    domain_2d["SpaceCorrelation"]["site"] = std::vector<double>(Ns_observe);
    for(int i = 0 ; i < Ns_observe; ++i) domain_2d["SpaceCorrelation"]["site"][i] = i;

    // Normalmodes
    file_name["NormalModeEnergy"] = settings.result_directory + "result_normalmode.dat";
    domain_2d["NormalModeEnergy"]["wave_vector"] = std::vector<double>(Ns_observe);
    for(int i = 0 ; i < Ns_observe; ++i) domain_2d["NormalModeEnergy"]["wave_vector"][i] = i;


    // Toda conservatsion
    file_name["GeneratorCorrelation"] = settings.result_directory + "result_generator_correlation.dat";
    domain_2d["GeneratorCorrelation"]["site"] = std::vector<double>(Ns_observe);
    for(int i = 0 ; i < Ns_observe; ++i) domain_2d["GeneratorCorrelation"]["site"][i] = i;

    // output 2d data 
    for(const auto& key  : quantities_2d_list){
      dataput.time_tag() << "start writing : " << key <<  std::endl;
      dataput.output_result_3d(
                               file_name[key],
                               domain_time,
                               domain_2d[key],
                               values_name[key],
                               results_2d
                               );
      dataput.time_tag() << "finish writing : " << key <<  std::endl;
    }

    dataput << "result output to : " ;
    for(const auto& elements: file_name){
      dataput << elements.second << ", ";
    }
    dataput << std::endl;

  }

}; // end struct define

#endif
