#include <physics/normalmode_energy_fftw.hpp>

inline void NormalModeEnergy(std::vector<double> & z, std::vector<double> & normalmode_energy_temp){
  NormalModeEnergyFFTW(z, normalmode_energy_temp);
}

template <class CorrelationFunction, class SpatialFunction, class Parameters>
inline void SpatialMeasurement(CorrelationFunction & correlation_function, SpatialFunction & spatial_function, 
    Lattice & lattice, Parameters const & parameters,int Ns_observe,int observe, int measure_step, std::vector<double> & z){
  double number_operator_observe = 0;
  double number_operator_target  = 0;
  int target = 0;
  int J = parameters.J;
  int num_particles = z.size()/2;
  double *x = &z[0];
  double *v = &z[num_particles];
  for(int i = 0; i < Ns_observe; ++i){
    target = lattice.target_number(i,observe);
    spatial_function["Displacement"][measure_step][i] << x[observe];
    spatial_function["Velocity"][measure_step][i] << v[observe];
    correlation_function["Displacement"][measure_step][i] << x[observe] * x[target];
    correlation_function["Velocity"][measure_step][i] << v[observe] * v[target];
    number_operator_observe = (v[observe] * v[observe] + J * x[observe] * x[observe]) / (2*std::sqrt(J));  
    number_operator_target = (v[target] * v[target] + J * x[target] * x[target]) / (2*std::sqrt(J));  
    spatial_function["NumberOperator"][measure_step][i] << number_operator_target;
    correlation_function["NumberOperator"][measure_step][i] << number_operator_observe * number_operator_target;
  }
}


