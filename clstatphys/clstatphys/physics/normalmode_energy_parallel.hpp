#ifndef NORMALMODE_ENERGY_PARALLEL_HPP
#define NORMALMODE_ENERGY_PARALLEL_HPP

void NormalModeEnergyParallel(const std::vector<double>& z, std::vector<double>& Ek){
  //energy calc
  int num_particles = z.size() / 2;
  double C0 = std::sqrt(2.0 / (double)num_particles);

  #pragma omp parallel for 
  for(int i = 0 ; i < num_particles ; ++i){
    double Pk = 0.0;
    double Qk = 0.0;
    for(int j = 0 ; j < num_particles ; ++j){
      Pk +=  C0 * std::sin(M_PI * (double)(i * j) / (double)num_particles) * z[num_particles + j]; 
      Qk +=  C0 * std::sin(M_PI * (double)(i * j) / (double)num_particles) * z[j]; 
    }
    double omega = 2.0 * std::sin(0.5 * M_PI * (double)i  / (double)num_particles);
    Ek[i] = 0.5 * ( Pk * Pk + omega * omega * Qk * Qk);
  }
}

#endif //NORMALMODE_ENERGY_PARALLEL_HPP
