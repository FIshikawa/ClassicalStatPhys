void SamplingProcess(){
  int observe = (num_particles-1) / 2;
  double *x = &z[0];
  double *v = &z[num_particles];
  double pt = 0.0;
  double v_total,x_total,sum_spectral,average_spectral,spectral_entropy;
  std::vector<double> normalmode_energy_temp(num_particles), 
                      toda_conservations_temp(num_particles),
                      toda_eigenvalues_temp(num_particles);

  if(process_id == 0) dataput.time_tag() << "start loop: sampling" << std::endl;
  //loop calc
  for(int c_loop=0; c_loop < N_each; ++c_loop){

    // set initial state
    ensembler.set_initial_state(z,mt);

    //time develop
    int measure_step = 0;
    time_measure.reset();
    for(int step = 0; step < N_time; ++step){  
      //set total values 
      if(time_measure.check(step)){
        v_total = 0.0;
        x_total = 0.0;
        for(int i = 0; i < num_particles; ++i){
          v_total += v[i] / num_particles;
          x_total += x[i] / num_particles;
        } 
        //set normalmode energy
        NormalModeEnergy(z,normalmode_energy_temp);

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
        physical_values["Velocity"][measure_step] << v[observe];
        physical_values["Velocity_total"][measure_step] << v_total;
        autocorrelation["VelocityACF"] << v[observe] - v_total;
        physical_values["Displacement"][measure_step] << x[observe];
        physical_values["Displacement_total"][measure_step] << x_total;
        autocorrelation["DisplacementACF"] << x[observe];
        physical_values["Energy_total"][measure_step] <<  hamiltonian.energy(pt,z) / num_particles;
        physical_values["Kinetic_Energy"][measure_step] << hamiltonian.kinetic_energy(pt,z) / num_particles;
        physical_values["Potential_Energy"][measure_step] << hamiltonian.potential_energy(pt,z) / num_particles;
        physical_values["SpectralEntropy"][measure_step] << spectral_entropy;
        physical_values["EffectiveFractionOfModes"][measure_step] << std::exp(spectral_entropy) / num_particles;

      //input hist data 
        hist_data["Velocity"][measure_step][c_loop] = v[observe];
        hist_data["Displacement"][measure_step][c_loop] = x_total;
        //correlation calc
        //SpatialMeasurement(correlation_function, spatial_function, lattice, parameters, Ns_observe, observe, measure_step, z);
        //input normalmode energy
        for(int i = 0; i < num_particles; ++i) {
          normalmode_energy[measure_step][i] << normalmode_energy_temp[i];
          toda_conservations[measure_step][i] << toda_conservations_temp[i];
          toda_eigenvalues[measure_step][i] << toda_eigenvalues_temp[i];
        }
        ++measure_step;
      }

      //time develp
      pt += dt;
      integrator.step(pt, dt, z, hamiltonian);
    }//end time step
    if(measure_step != N_total_data){
      std::cerr << "final measure_step : " << measure_step << ", N_total_data : " << N_total_data <<", not same " << std::endl;
      std::exit(112);
    }

  }//end loop
}
