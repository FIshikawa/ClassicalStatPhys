#ifndef PHYISCAL_QUANTITIES_HPP
#define PHYISCAL_QUANTITIES_HPP

#include <unordered_map>
#include <algorithm>
#include <clstatphys/tools/accumulator.hpp>
#include <clstatphys/physics/normalmode_energy_periodic_fftw.hpp>
#include <clstatphys/physics/toda_discriminant.hpp>
#include <clstatphys/physics/toda_action_variables.hpp>

using PhysicalQuantities1d = std::unordered_map<std::string, std::vector<tools::Accumulator> >; 
using PhysicalQuantities2d = std::unordered_map<std::string, std::vector<std::vector<tools::Accumulator> > >;

struct PhysicalQuantities{
  int num_particles;
  int N_total_data;
  int num_iterations;
  double total_energy_init, spectral_entropy_init;
  std::vector<double> conservations_init, eigenvalues_init, action_variables_init;
  std::vector<std::string> quantities_1d_list = {
                                               "Energy_total",
                                               "Kinetic_Energy",
                                               "Potential_Energy",
                                               "EffectiveFractionOfModes",
                                               "SpectralEntropy",
                                               "Energy_totalACF",
                                               "SpectralEntropyACF",
                                               };
  std::vector<std::string> quantities_2d_list = {
                                               "NormalModeEnergy",
                                               "TodaConservations",
                                               "TodaEigenvalues",
                                               "TodaActionVariables",
                                               "TodaConservationsACF",
                                               "TodaEigenvaluesACF",
                                               "TodaActionVariablesACF",
                                               };

  PhysicalQuantities1d quantities_1d;
  PhysicalQuantities2d quantities_2d;

  integrable::TodaLaxForm toda_lax_form;
  Hamiltonian hamiltonian;

  template<typename Settings>
  PhysicalQuantities(Settings & settings){
    N_total_data = settings.N_total_data;
    num_particles = settings.num_particles;
    num_iterations = settings.num_iterations;

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
    
    // init initial vector for ACF
    conservations_init = std::vector<double>(num_particles,0.0);
    eigenvalues_init= std::vector<double>(num_particles,0.0);
    action_variables_init= std::vector<double>(num_particles,0.0);

    toda_lax_form = settings.toda_lax_form();
    hamiltonian = settings.hamiltonian();
  }


  void measure(std::vector<double> const & z, int & measure_step){
    double pt = 0.0;
    double v_total,x_total,sum_spectral,average_spectral,spectral_entropy,total_energy;
    std::vector<double> normalmode_energy_temp(num_particles,0.0), 
                        toda_action_variables_temp(num_particles,0.0),
                        toda_conservations_temp(num_particles,0.0),
                        toda_eigenvalues_temp(num_particles,0.0);
    //set normalmode energy
    NormalModeEnergyPeriodicFFTW(z,normalmode_energy_temp);

    //set action variables 
    rokko::dlmatrix L = toda_lax_form.L_matrix(z);
    integrable::TodaDiscriminant discriminant(num_particles, L, "periodic");
    TodaActionVariables(toda_action_variables_temp, discriminant, num_iterations);
  
    //calc spectral enetropy
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
    total_energy = hamiltonian.energy(pt,z) / num_particles;

    quantities_1d["Energy_total"][measure_step] << total_energy;
    quantities_1d["Kinetic_Energy"][measure_step] << hamiltonian.kinetic_energy(pt,z) / num_particles;
    quantities_1d["Potential_Energy"][measure_step] << hamiltonian.potential_energy(pt,z) / num_particles;
    quantities_1d["SpectralEntropy"][measure_step] << spectral_entropy;
    quantities_1d["EffectiveFractionOfModes"][measure_step] << std::exp(spectral_entropy) / num_particles;
  
  
    if(measure_step == 0){
      total_energy_init = total_energy;
      spectral_entropy_init = spectral_entropy;
      for(int i = 0; i < num_particles; ++i){
        conservations_init[i] = toda_conservations_temp[i];
        eigenvalues_init[i] = toda_eigenvalues_temp[i];
        action_variables_init[i] = toda_action_variables_temp[i];
      }
    }
    quantities_1d["Energy_totalACF"][measure_step] << total_energy_init * total_energy;
    quantities_1d["SpectralEntropyACF"][measure_step] << spectral_entropy_init * spectral_entropy;

    for(int i = 0; i < num_particles; ++i) {
      quantities_2d["NormalModeEnergy"][measure_step][i] << normalmode_energy_temp[i];
      quantities_2d["TodaConservations"][measure_step][i] << toda_conservations_temp[i];
      quantities_2d["TodaEigenvalues"][measure_step][i] << toda_eigenvalues_temp[i];
      quantities_2d["TodaActionVariables"][measure_step][i] << toda_action_variables_temp[i];
      quantities_2d["TodaConservationsACF"][measure_step][i] 
                                      << conservations_init[i] * toda_conservations_temp[i];
      quantities_2d["TodaEigenvaluesACF"][measure_step][i] 
                                      << eigenvalues_init[i] * toda_eigenvalues_temp[i];
      quantities_2d["TodaActionVariablesACF"][measure_step][i] 
                                      << action_variables_init[i] * toda_action_variables_temp[i];
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
    for(int step = 0; step < N_total_data; ++step) domain_time["time"][step] = step;

    // set ACF 
    for(const auto& key  : quantities_1d_list){
      if(std::find(quantities_1d_list.begin(), quantities_1d_list.end(), key + "ACF") != quantities_1d_list.end()){
        for(int step = 0; step < N_total_data; ++step){
          results_normal[key+"ACF"][step] -= results_normal[key][0] * results_normal[key][step]; 
          results_normal[key+"ACF"][step] /= std::sqrt(quantities_1d[key][0].variance())
                                                * std::sqrt(quantities_1d[key][step].variance());
          results_normal["error_"+key+"ACF"][step] = 0.0; 
        }
      }
    }

  dataput.output_result(file_name["normal"], domain_time, results_normal);


    // 2dim resutls set 
    std::unordered_map<std::string, std::vector< std::vector<double> > > results_2d;
    std::unordered_map<std::string, std::unordered_map< std::string, std::vector<double> > > domain_2d;

    for(const auto& key  : quantities_2d_list){
      values_name[key].push_back(key);
      values_name[key].push_back("error_"+key);

      results_2d[key] 
          = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
      results_2d["error_" + key] 
          = std::vector< std::vector<double> >(N_total_data, std::vector<double>(num_particles,0.0));
      for(int step = 0; step < N_total_data; ++step){
        for(int site = 0; site < num_particles; ++site){
          results_2d[key][step][site] = quantities_2d[key][step][site].mean();
          results_2d["error_" + key][step][site] = quantities_2d[key][step][site].error();
        }
      }
    }

    for(const std::string & key : {"TodaConservations","TodaEigenvalues","TodaActionVariables"}){
      values_name[key] 
        = std::vector<std::string>{key,"error_" + key,
                                   key + "ACF","error_" + key + "ACF"};
      // set ACF 
      for(int step = 0; step < N_total_data; ++step){
        for(int site = 0; site < num_particles; ++site){
          results_2d[key+"ACF"][step][site] -= results_2d[key][0][site] * results_2d[key][step][site]; 
          results_2d[key+"ACF"][step][site] /= std::sqrt(quantities_2d[key][0][site].variance())
                                                * std::sqrt(quantities_2d[key][step][site].variance());
          results_2d["error_"+key+"ACF"][step][site] = 0.0;
        }
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

    // Toda action variables
    file_name["TodaActionVariables"] = settings.result_directory + "result_toda_action_variables.dat";
    domain_2d["TodaActionVariables"]["number"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_2d["TodaActionVariables"]["number"][i] = i;


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
