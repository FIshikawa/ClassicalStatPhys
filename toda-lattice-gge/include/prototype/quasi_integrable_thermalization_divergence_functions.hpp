#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

void SetHarmonicDist(std::vector<double> & harmonic_spectrum,int N_normalmode, int k_initial){
  int num_particles = harmonic_spectrum.size();

  for(int i = k_initial ; i < k_initial + N_normalmode/2; ++i){
    harmonic_spectrum[i] = 1.0 / N_normalmode;
  } 
  for(int i = num_particles - N_normalmode/2 - (k_initial-1) ; i < num_particles - (k_initial-1); ++i){
    harmonic_spectrum[i] = 1.0 / N_normalmode;
  } 
}

int SetTodaDist(std::vector<double> & toda_spectrum, const double E_initial){
  int num_particles = toda_spectrum.size();
  if(num_particles < 1024){
    std::cerr << "num_particles should be 1024" << std::endl;
    std::exit(1);
  }
  std::string referece_name;
  if(E_initial == 0.1 * 1024) referece_name = "toda_spectrum_E01.dat";
  else if(E_initial == 0.01 * 1024) referece_name = "toda_spectrum_E001.dat";
  else if(E_initial == 0.001 * 1024) referece_name = "toda_spectrum_E0001.dat";
  else if(E_initial == 0.0025 * 1024) referece_name = "toda_spectrum_E00025.dat";
  else if(E_initial == 0.005 * 1024) referece_name = "toda_spectrum_E0005.dat";
  else if(E_initial == 0.05 * 1024) referece_name = "toda_spectrum_E005.dat";
  else if(E_initial == 0.075 * 1024) referece_name = "toda_spectrum_E0075.dat";
  else{
    std::cerr << E_initial << " is unexpected " << std::endl;
    std::exit(1124);
  }
  referece_name = "./reference_data/" + referece_name;
  std::ifstream toda_reference(referece_name);
  if(!toda_reference){
    std::cerr << "failed open" << std::endl;
    std::exit(1124);
  }
  std::string data_temp;
  int counter = 0;
  while(std::getline(toda_reference,data_temp)){
    toda_spectrum[counter] = std::stod(data_temp);
    counter += 1;
  }
  return counter;
}
