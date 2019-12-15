void Finalize(){

  // correct data each process
  dataput.time_tag() << "[" << process_id <<"th process] start result correct by self" << std::endl;
  for(const auto& key  : physical_value_list){
    for(int i = 0 ; i < N_total_data; ++i){
      results[key][i] = physical_values[key][i].sum1();
      results["error_" + key][i] = physical_values[key][i].sum2();
    }
  }
  for(const auto& element : autocorrelation){
    std::string key = element.first;
    std::vector<double> acf_result = autocorrelation[key].result();
    for(int i = 0 ; i < N_total_data; ++i){
      results[key][i] = acf_result[i] ;
    }      
  }
  for(int step = 0 ; step < N_total_data; ++step){
    for(int i = 0 ; i < num_particles; ++i){
      results_2d["Normalmode_Energy"][step][i] = normalmode_energy[step][i].sum1(); 
      results_2d["error_Normalmode_Energy"][step][i] = normalmode_energy[step][i].sum2();
      results_2d["TodaConservations"][step][i] = toda_conservations[step][i].sum1(); 
      results_2d["error_TodaConservations"][step][i] = toda_conservations[step][i].sum2();
      results_2d["TodaConservationsVariance"][step][i] = toda_conservations[step][i].sum2(); 
      results_2d["error_TodaConservationsVariance"][step][i] = toda_conservations[step][i].sum4(); 
      results_2d["TodaEigenvalues"][step][i] = toda_eigenvalues[step][i].sum1(); 
      results_2d["error_TodaEigenvalues"][step][i] = toda_eigenvalues[step][i].sum2();
      results_2d["TodaEigenvaluesVariance"][step][i] = toda_eigenvalues[step][i].sum2(); 
      results_2d["error_TodaEigenvaluesVariance"][step][i] = toda_eigenvalues[step][i].sum4();
    }
//    for(int i = 0 ; i < Ns_observe; ++i){
//      for(std::string key : {"Displacement","Velocity","NumberOperator"}){
//        results_2d[key + "_correlation"][step][i] = correlation_function[key][step][i].sum1();
//        results_2d["error_" + key + "_correlation"][step][i] = correlation_function[key][step][i].sum2();
//        results_2d[key][step][i] = spatial_function[key][step][i].sum1();
//        results_2d["error_" + key][step][i] = spatial_function[key][step][i].sum2();
//      }
//    }
  }
  dataput.time_tag() << "[" << process_id <<"th process] finish correction : ready for transfer" << std::endl;


  // correct data to process id 0 
  if(process_id == 0) dataput.time_tag() << "start total result correct" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start normal results" << std::endl;
  std::vector<double> reduced_vector(N_total_data);
  for(const auto& key  : physical_value_list){
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_total_data; ++i) results[key][i] = reduced_vector[i] / N_loop;
    mpi_error = MPI_Reduce(&results["error_" + key][0], &reduced_vector[0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_total_data; ++i) 
      results["error_" + key][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector[i] / (N_loop - 1) - results[key][i] * results[key][i]* (N_loop)/(N_loop-1)) / N_loop) : 0;
  }
  if(process_id == 0) dataput.time_tag() << "end normal results" << std::endl;

  if(process_id == 0) dataput.time_tag() << "start autocorrelation results" << std::endl;
  for(const auto& element : autocorrelation){
    std::string key = element.first;
    mpi_error = MPI_Reduce(&results[key][0], &reduced_vector[0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(process_id == 0) for(int i = 0 ; i < N_total_data; ++i) results[key][i] = reduced_vector[i] / N_loop ;
  }
  if(process_id == 0) dataput.time_tag() << "end autocorrelation results" << std::endl;

  std::vector<double> reduced_vector_1(num_particles);
  std::vector<double> reduced_vector_2(Ns_observe);
  if(process_id == 0) dataput.time_tag() << "start 2d results" << std::endl;
  for(int step = 0 ; step < N_total_data; ++step){
    // NormalMode
    for(std::string value_name: {"Normalmode_Energy","TodaConservations","TodaEigenvalues","TodaConservationsVariance","TodaEigenvaluesVariance"}){
      mpi_error = MPI_Reduce(&results_2d[value_name][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) results_2d[value_name][step][i] = reduced_vector_1[i] / N_loop;
      mpi_error = MPI_Reduce(&results_2d["error_" + value_name][step][0], &reduced_vector_1[0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(process_id == 0) for(int i = 0 ; i < num_particles; ++i) 
        results_2d["error_" + value_name][step][i] = (N_loop > 1 ) ? std::sqrt((reduced_vector_1[i] / (N_loop-1)
                                             - results_2d[value_name][step][i] * results_2d[value_name][step][i]* (N_loop)/(N_loop-1))/N_loop) : 0; 
    }
 

//    // correlation values
//    for(const auto& key : correlation_list){
//      mpi_error = MPI_Reduce(&results_2d[key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) results_2d[key][step][i] = reduced_vector_2[i] / N_loop;
//      mpi_error = MPI_Reduce(&results_2d["error_" + key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) 
//        results_2d["error_" + key][step][i] = (N_loop > 1) ? std::sqrt((reduced_vector_2[i] / (N_loop -1)
//                                    - results_2d[key][step][i] * results_2d[key][step][i]) / N_loop ) : 0;
//    }
//
//    // spatial values
//    for(const auto& key : spatial_list){
//      mpi_error = MPI_Reduce(&results_2d[key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) results_2d[key][step][i] = reduced_vector_2[i] / N_loop;
//      mpi_error = MPI_Reduce(&results_2d["error_" + key][step][0], &reduced_vector_2[0], Ns_observe, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//      if(process_id == 0) for(int i = 0 ; i < Ns_observe; ++i) 
//        results_2d["error_" + key][step][i] = (N_loop > 1) ? std::sqrt((reduced_vector_2[i] / (N_loop -1)
//                                    - results_2d[key][step][i] * results_2d[key][step][i]) / N_loop ) : 0;
//    }
  }
  if(process_id == 0) dataput.time_tag() << "end 2d results" << std::endl;

  std::vector<double> hist_vector(N_loop);
  for(auto key : {"Velocity","Displacement"}){
    for(int step = 0 ; step < N_total_data; ++step){
      mpi_error = MPI_Gather(&hist_data[key][step][0], N_each, MPI_DOUBLE, &hist_vector[0], N_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for(int i = 0 ; i < N_loop; ++i) hist_data[key][step][i] = hist_vector[i];
    }
  }

  if(process_id == 0) dataput.time_tag() << "end total result correct" << std::endl;

}
