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
                                               "TodaConservations",
                                               "TodaEigenvalues"
                                               };

  PhysicalQuantities1d quantities_1d;
  PhysicalQuantities2d quantities_2d;

  integrable::TodaLaxForm toda_lax_form;
  Hamiltonian hamiltonian;

  template<typename Settings>
  PhysicalQuantities(Settings & settings){
    N_total_data = settings.N_total_data;
    num_particles = settings.num_particles;

    for(const auto& key : quantities_1d_list){
      quantities_1d[key].resize(N_total_data);
      for(int step = 0; step < N_total_data; ++step) quantities_1d[key][step].reset();
    }
    for(const auto& key : quantities_2d_list){
      quantities_2d[key].resize(N_total_data);
      for(int step = 0 ; step < N_total_data ; ++step){
        quantities_2d[key][step] = std::vector<tools::Accumulator>(num_particles);
        for(int i = 0; i < num_particles ; ++i) quantities_2d[key][step][i].reset();
      }
    }
    toda_lax_form = settings.toda_lax_form();
    hamiltonian = settings.hamiltonian();
  }


  void measure(std::vector<double> z, int & measure_step){
    int observe = (num_particles-1) / 2;
    double *x = &z[0];
    double *v = &z[num_particles];
    double pt = 0.0;
    double v_total,x_total,sum_spectral,average_spectral,spectral_entropy;
    std::vector<double> normalmode_energy_temp(num_particles), 
                        toda_conservations_temp(num_particles),
                        toda_eigenvalues_temp(num_particles);
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
  
    //set toda_conservations
    toda_lax_form.conservations_with_eigenvalues(z, toda_conservations_temp, toda_eigenvalues_temp);
  
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

    for(int i = 0; i < num_particles; ++i) {
      quantities_2d["NormalModeEnergy"][measure_step][i] << normalmode_energy_temp[i];
      quantities_2d["TodaConservations"][measure_step][i] << toda_conservations_temp[i];
      quantities_2d["TodaEigenvalues"][measure_step][i] << toda_eigenvalues_temp[i];
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
      results_normal[key].resize(N_total_data);
      results_normal["error_" + key].resize(N_total_data);
      for(int step = 0; step < N_total_data; ++step){
        results_normal[key][step] = quantities_1d[key][step].mean();
        results_normal["error_" + key][step] = quantities_1d[key][step].error();
      }
      values_name["normal"].push_back(key);
      values_name["normal"].push_back("error_" + key);
    }
    std::unordered_map< std::string, std::vector<double> > domain_time;
    domain_time["time"] = std::vector<double>(N_total_data);
    tools::TimeMeasure time_measure = settings.time_measure();
    for(int step = 0; step < N_total_data; ++step) domain_time["time"][step] = time_measure();

    dataput.output_result(file_name["normal"], domain_time, results_normal);


    // 2dim resutls set 
    std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d;
    std::unordered_map<std::string, std::unordered_map< std::string, std::vector<double> > > domain_2d;

    for(const auto& key  : quantities_2d_list){
      results_2d[key] 
          = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
      results_2d["error_" + key] 
          = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
      for(int step = 0; step < N_total_data; ++step){
        for(int site = 0; site < num_particles; ++site){
          results_2d[key][step][site] = quantities_2d[key][step][site].mean();
          results_2d["error_" + key][step][site] = quantities_2d[key][step][site].error();
          values_name[key].push_back(key);
          values_name[key].push_back("error_"+key);
        }
      }
    }

    values_name["TodaConservations"] 
      = std::vector<std::string>{"TodaConservations","error_TodaConservations",
                                 "TodaConservationsVariance","error_TodaConservationsVariance"};
    results_2d["TodaConservationsVariance"] 
        = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
    results_2d["error_TodaConservationsVariance"] 
        = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));

    for(int step = 0; step < N_total_data; ++step){
      for(int site = 0; site < num_particles; ++site){
        results_2d["TodaConservationsVariance"][step][site] 
            = quantities_2d["TodaConservations"][step][site].variance();
        results_2d["error_TodaConservationsVariance"][step][site] 
          = quantities_2d["TodaConservations"][step][site].variance_error();
      }
    }

    // Eigenvalues
    file_name["TodaEigenvalues"] = settings.result_directory + "result_toda_eigenvalues.dat";
    domain_2d["TodaEigenvalues"]["number"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_2d["TodaEigenvalues"]["number"][i] = i;

    // Normalmodes
    file_name["NormalModeEnergy"] = settings.result_directory + "result_normalmode.dat";
    domain_2d["NormalModeEnergy"]["wave_vector"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_2d["NormalModeEnergy"]["wave_vector"][i] = i;


    // Toda conservatsion
    file_name["TodaConservations"] = settings.result_directory + "result_toda_conservations.dat";
    domain_2d["TodaConservations"]["power"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_2d["TodaConservations"]["power"][i] = i;

    // output 2d data 
    for(const auto& key  : quantities_2d_list){
      dataput.output_result_3d(
                               file_name[key],
                               domain_time,
                               domain_2d[key],
                               values_name[key],
                               results_2d
                               );
    }

    dataput << "result output to : " ;
    for(const auto& elements: file_name){
      dataput << elements.second << ", ";
    }
    dataput << std::endl;

  }

}; // end struct define

#endif
