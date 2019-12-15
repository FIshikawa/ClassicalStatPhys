#include <physics/normalmode_energy_periodic_fftw.hpp>
#include <iostream>

inline void NormalModeEnergy(std::vector<double> & z, std::vector<double> & normalmode_energy_temp){
  NormalModeEnergyPeriodicFFTW(z, normalmode_energy_temp);
}

template <class TodaConservationField, class Hamiltonian, class CorrelationFunction, class SpatialFunction>
inline void SpatialMeasurement(CorrelationFunction & correlation_function, SpatialFunction & spatial_function, 
    Lattice & lattice, Hamiltonian & hamiltonian, TodaConservationField & toda_conservation_fields, 
    int Ns_observe,int observe, int measure_step, std::vector<double> & z){

  int target = 0;
  int num_particles = z.size()/2;
  double energy_field;
  double energy_field_observe, energy_field_target;
  double x_value, v_value, x_correlation, v_correlation, energy_field_correlation;
  double *x = &z[0];
  double *v = &z[num_particles];
  std::vector<double> toda_integral_values(5), toda_integral_correlation(5);
  for(int i = 0; i < Ns_observe; ++i){
    energy_field = 0;
    energy_field_observe = 0;
    energy_field_target = 0;
    energy_field_correlation = 0;
    x_value = 0;
    v_value = 0;
    x_correlation = 0;
    v_correlation = 0;
    for(int kind = 1; kind < 6; ++kind){
      toda_integral_values[kind-1] = 0.0;
      toda_integral_correlation[kind-1] = 0.0;
    }
    for(int observe_temp = 0; observe_temp < num_particles; ++observe_temp){
      target = lattice.target_number(i,observe_temp);
      energy_field_observe = hamiltonian.energy_field(observe_temp, z, 0.0);  
      energy_field_target = hamiltonian.energy_field(target,z, 0.0); 
      x_value += x[observe_temp];
      v_value += v[observe_temp];
      x_correlation += x[observe_temp] * x[target];
      v_correlation += v[observe_temp] * v[target];
      energy_field += energy_field_observe;
      energy_field_correlation += energy_field_observe * energy_field_target;
      for(int kind = 1; kind < 6; ++kind){
        double tagged_value = toda_conservation_fields(z,observe_temp,kind);
        double target_value = toda_conservation_fields(z,target,kind);
        toda_integral_values[kind-1] += tagged_value;
        toda_integral_correlation[kind-1] += tagged_value * target_value;
      }
    }
    x_value /= num_particles;
    v_value /= num_particles;
    x_correlation /= num_particles;
    v_correlation /= num_particles;
    energy_field /= num_particles;
    energy_field_correlation /= num_particles;
    for(int kind = 1; kind < 6; ++kind){
      toda_integral_values[kind-1] /= num_particles;
      toda_integral_correlation[kind-1] /= num_particles;
    }

    spatial_function["Displacement"][measure_step][i] << x_value;
    spatial_function["Velocity"][measure_step][i] << v_value;
    correlation_function["Displacement"][measure_step][i] << x_correlation;
    correlation_function["Velocity"][measure_step][i] << v_correlation;
    spatial_function["EnergyField"][measure_step][i] << energy_field;
    correlation_function["EnergyField"][measure_step][i] << energy_field_correlation;
    spatial_function["TodaIntegralThree"][measure_step][i] << toda_integral_values[2];
    correlation_function["TodaIntegralThree"][measure_step][i] << toda_integral_correlation[2];
    spatial_function["TodaIntegralFour"][measure_step][i] << toda_integral_values[3];
    correlation_function["TodaIntegralFour"][measure_step][i] << toda_integral_correlation[3];
    spatial_function["TodaIntegralFive"][measure_step][i] << toda_integral_values[4];
    correlation_function["TodaIntegralFive"][measure_step][i] << toda_integral_correlation[4];
  }
}

