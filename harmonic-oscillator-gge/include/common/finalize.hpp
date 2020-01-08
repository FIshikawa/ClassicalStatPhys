
template<typename Settings, typename PhysicalQuanties>
void Finalize(Settings & settings, PhysicalQuanties & physical_quantities){
  int num_particles = settings.num_particles;
  int N_total_data  = settings.N_total_data;
  int process_id    = settings.process_id;
  int mpi_error = 0;

  std::unordered_map<std::string, std::vector<double> > corrections;
  std::unordered_map<std::string, std::vector<double> > reductions;
  std::vector<std::string> values_name{"sum1","sum2","sum3","sum4","count"};
  for(const auto& key : values_name){
    corrections[key] = std::vector<double>(N_total_data, 0.0);
    reductions[key] = std::vector<double>(N_total_data, 0.0);
  }

  for(const auto& key  : physical_quantities.quantities_1d_list){
    for(int step = 0; step < N_total_data; ++step){
      corrections["sum1"][step] = physical_quantities.quantities_1d[key][step].sum1();
      corrections["sum2"][step] = physical_quantities.quantities_1d[key][step].sum2();
      corrections["sum3"][step] = physical_quantities.quantities_1d[key][step].sum3();
      corrections["sum4"][step] = physical_quantities.quantities_1d[key][step].sum4();
      corrections["count"][step] = static_cast<double>(physical_quantities.quantities_1d[key][step].count());
    }
    for(const auto& name: values_name)
      mpi_error = MPI_Reduce(&corrections[name][0], &reductions[name][0], N_total_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(process_id == 0){
      for(int step = 0; step < N_total_data; ++step){
        physical_quantities.quantities_1d[key][step].reset();
        physical_quantities.quantities_1d[key][step].add(reductions["sum1"][step], reductions["sum2"][step], 
                                                         reductions["sum3"][step], reductions["sum4"][step],
                                                         static_cast<long>(reductions["count"][step]));
      }
    }
  }

  for(const auto& key : values_name){
    corrections[key] = std::vector<double>(num_particles, 0.0);
    reductions[key] = std::vector<double>(num_particles, 0.0);
  }

  for(const auto& key  : physical_quantities.quantities_2d_list){
    for(int step = 0; step < N_total_data; ++step){
      for(int site = 0; site < num_particles; ++site){
        corrections["sum1"][site] = physical_quantities.quantities_2d[key][step][site].sum1();
        corrections["sum2"][site] = physical_quantities.quantities_2d[key][step][site].sum2();
        corrections["sum3"][site] = physical_quantities.quantities_2d[key][step][site].sum3();
        corrections["sum4"][site] = physical_quantities.quantities_2d[key][step][site].sum4();
        corrections["count"][site] = static_cast<double>(physical_quantities.quantities_2d[key][step][site].count());
      }
      for(const auto& name: values_name)
        mpi_error = MPI_Reduce(&corrections[name][0], &reductions[name][0], num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if(process_id == 0){
        for(int site = 0; site < num_particles; ++site){
          physical_quantities.quantities_2d[key][step][site].reset();
          physical_quantities.quantities_2d[key][step][site].add(reductions["sum1"][site], reductions["sum2"][site], 
                                                           reductions["sum3"][site], reductions["sum4"][site],
                                                           static_cast<long>(reductions["count"][site]));
        }
      }
    }
  }
};
