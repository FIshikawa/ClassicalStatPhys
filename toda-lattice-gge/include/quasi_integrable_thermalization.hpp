#include <tools/auto_correlation_function.hpp>
#include <tools/histogram.hpp>
#include <tools/accumulator.hpp>
#include <physics/toda_lax_form.hpp>


int main(int argc, char **argv){

  //set correlatin function 
//  std::vector<std::string> correlation_list{"Velocity_correlation","Displacement_correlation","NumberOperator_correlation"}; 
//  std::vector<std::string> spatial_list{"Velocity","Displacement","NumberOperator"}; 
//  for(const auto& key : correlation_list){
//    results_2d[key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
//    results_2d["error_" + key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
//  }
//   for(const auto& key : spatial_list){
//    results_2d[key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
//    results_2d["error_" + key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(Ns_observe,0.0));
//  }

  //set histogram maps 
  std::unordered_map<std::string, std::vector< std::vector<double> > > hist_data; 
 	hist_data["Velocity"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(N_loop)),
 	hist_data["Displacement"] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(N_loop));
  //end set results and measure vectors

  if(process_id == 0) dataput.time_tag() << "start sampling process : time development " << std::endl;



  //set correlation
  std::unordered_map<std::string, correlation::AutoCorrelationFunction> autocorrelation;
  autocorrelation["VelocityACF"].initialize(N_total_data,N_loop);
  autocorrelation["DisplacementACF"].initialize(N_total_data,N_loop);

  //measure physical values set
  std::unordered_map<std::string, std::vector<tools::Accumulator> > physical_values;
  for(const auto& key : physical_value_list){
    physical_values[key].resize(N_total_data);
    for(int i = 0; i < N_total_data; ++i) physical_values[key][i].reset(); 
  }


  //set normalmode energy
  std::vector<std::vector<tools::Accumulator> > normalmode_energy(N_total_data,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_total_data ; ++step){
    normalmode_energy[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) normalmode_energy[step][i].reset();
  }

  //set toda conservations
  std::vector<std::vector<tools::Accumulator> > toda_conservations(N_total_data,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_total_data ; ++step){
    toda_conservations[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) toda_conservations[step][i].reset();
  }

  //set toda eigenvalues
  std::vector<std::vector<tools::Accumulator> > toda_eigenvalues(N_total_data,std::vector<tools::Accumulator>(num_particles));
  for(int step = 0 ; step < N_total_data ; ++step){
    toda_eigenvalues[step] = std::vector<tools::Accumulator>(num_particles);
    for(int i = 0; i < num_particles ; ++i) toda_eigenvalues[step][i].reset();
  }




  // data output by process_id = 0
  if(process_id == 0){
    dataput.time_tag() << "end sampling process" << std::endl;

//    dataput.time_tag() << "correlation values calc" << std::endl;
//    double correlation_error;
//    for(const auto& key : spatial_list){
//      for(int step = 0 ; step < N_total_data; ++step){
//        for(int i = 0; i < Ns_observe; ++i){
//          results_2d[key + "_correlation"][step][i] -= results_2d[key][step][i] * results_2d[key][step][0];
//          correlation_error = results_2d["error_" + key][step][i] * results_2d[key][step][0]
//                              * results_2d["error_" + key][step][i] * results_2d[key][step][0]
//                              + results_2d["error_" + key][step][0] * results_2d[key][step][i]
//                              * results_2d["error_" + key][step][0] * results_2d[key][step][i]
//                              + results_2d["error_" + key + "_correlation"][step][i]
//                              * results_2d["error_" + key + "_correlation"][step][i];
//          results_2d["error_" + key + "_correlation"][step][i] = std::sqrt(correlation_error);
//        }
//      }
//    }
//    dataput.time_tag() << "end correlation values calc" << std::endl;

    dataput.time_tag() << "start stat calc " << std::endl;
    //set Accumulator
    tools::Accumulator accumulator;

    //define histogram class
    std::unordered_map<std::string, tools::Histogram>  histogram;
    histogram["Displacement"].initialize(n_bin);
    histogram["Velocity"].initialize(n_bin);

    dataput.time_tag() << "start stat declaration " << std::endl;
    //set statistical values
    std::unordered_map<std::string, std::vector<double> > stat_values;
    std::vector<std::string> stat_value_list = {
                                                "mean",
                                                "variance",
                                                "skewness",
                                                "kurtosis",
                                                "kurtosis_exess",
                                                "average","error",
                                                "central_moment1",
                                                "central_moment2",
                                                "central_moment3",
                                                "central_moment4"
                                              };

    std::unordered_map<std::string, std::vector< std::vector<double> > > hist_values;
    std::vector<std::string> hist_value_list = {
                                                "hist_Velocity",
                                                "range_Velocity",
                                                "hist_Displacement",
                                                "range_Displacement"
                                                };

    for(auto key  : stat_value_list) stat_values[key].resize(N_total_data);
    for(auto key  : hist_value_list) hist_values[key] = std::vector< std::vector<double> >(N_total_data, std::vector<double>(n_bin));

    dataput.time_tag() << "end stat declaration " << std::endl;
    //calc histogram and statisitcal values
    for(int i = 0; i < N_total_data; ++i){
      accumulator.reset();   
      for(int j = 0; j < N_loop; ++j){
        accumulator << hist_data["Velocity"][i][j];
      }
      stat_values["mean"][i] = accumulator.mean();
      stat_values["variance"][i] = accumulator.variance();
      stat_values["skewness"][i] = accumulator.skewness();
      stat_values["kurtosis"][i] = accumulator.kurtosis();
      stat_values["kurtosis_exess"][i] = accumulator.kurtosis_excess();
      stat_values["average"][i] = accumulator.average();
      stat_values["error"][i] = accumulator.error();
      stat_values["central_moment1"][i] = accumulator.central_moment1();
      stat_values["central_moment2"][i] = accumulator.central_moment2();
      stat_values["central_moment3"][i] = accumulator.central_moment3();
      stat_values["central_moment4"][i] = accumulator.central_moment4();
      for(auto element : histogram){
        std::string key = element.first; 
        histogram[key](hist_data[key][i]);
        histogram[key].output(hist_values["hist_" + key][i],
                                          hist_values["range_" + key][i]);
      }
    }
    dataput.time_tag() << "end stat calc " << std::endl;

    //set measure time
    std::unordered_map< std::string, std::vector<double> > domain_time;
    domain_time["time"] = std::vector<double>(N_total_data);
    time_measure.reset();
    for(int step = 0; step < N_total_data; ++step) domain_time["time"][step] = time_measure();

    //set normalmode plot
    std::unordered_map< std::string, std::vector<double> > domain_normalmode;
    domain_normalmode["wave_vector"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_normalmode["wave_vector"][i] = i;
    std::vector<std::string> normalmode_value_list{"Normalmode_Energy","error_Normalmode_Energy"};

    //set toda conservations plot
    std::unordered_map< std::string, std::vector<double> > domain_toda_conservations;
    domain_toda_conservations["power"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_toda_conservations["power"][i] = i;
    std::vector<std::string> toda_conservations_list{"TodaConservations","error_TodaConservations",
                                                      "TodaConservationsVariance","error_TodaConservationsVariance"};

    //set toda eigenvalues plot
    std::unordered_map< std::string, std::vector<double> > domain_toda_eigenvalues;
    domain_toda_eigenvalues["number"] = std::vector<double>(num_particles);
    for(int i = 0 ; i < num_particles ; ++i) domain_toda_eigenvalues["number"][i] = i;
    std::vector<std::string> toda_eigenvalues_list{"TodaEigenvalues","error_TodaEigenvalues",
                                                    "TodaEigenvaluesVariance","error_TodaEigenvaluesVariance"};

//    //add error values of correlation and spatial to string list
//    int correlation_length = correlation_list.size();
//    for(int i = 0 ; i < correlation_length; ++i){
//      correlation_list.push_back("error_" + correlation_list[i]);
//    }
//    int spatial_length = spatial_list.size();
//    for(int i = 0 ; i < spatial_length; ++i){
//      spatial_list.push_back("error_" + spatial_list[i]);
//    }
    int physical_value_length = physical_value_list.size();
    for(int i =0 ; i < physical_value_length; ++i){
      physical_value_list.push_back("error_" + physical_value_list[i]);
    }

//    //set sight number domain
//    std::unordered_map< std::string, std::vector<double> > domain_sight;
//    domain_sight["site"] = std::vector<double>(Ns_observe);
//    for(int i = 0 ; i < Ns_observe; ++i) domain_sight["site"][i] = i;

    std::string result_dat("result.dat");
    std::string hist_dat("result_hist.dat");
    std::string normalmode_dat("result_normalmode.dat");
    std::string toda_conservations_dat("result_toda_conservations.dat");
    std::string toda_eigenvalues_dat("result_toda_eigenvalues.dat");
    std::string stat_dat("result_stat.dat");
//    std::string spatial_dat("result_spatial.dat");
//    std::string correlation_dat("result_correlation.dat");

    result_dat = result_directory + result_dat;
    hist_dat = result_directory + hist_dat;
    stat_dat = result_directory + stat_dat;
    normalmode_dat = result_directory + normalmode_dat;
    toda_conservations_dat = result_directory + toda_conservations_dat;
    toda_eigenvalues_dat = result_directory + toda_eigenvalues_dat;
//    spatial_dat = result_directory + spatial_dat;
//    correlation_dat = result_directory + correlation_dat;

    dataput << "result output to : " 
            << result_dat << " , "
//            << correlation_dat << " , "
//            << spatial_dat << " , "
            << normalmode_dat << " , "
            << toda_conservations_dat << " , "
            << toda_eigenvalues_dat << " , "
            << hist_dat << " and " 
            << stat_dat 
            << std::endl;

    dataput.output_result(result_dat,domain_time,results);
    dataput.output_result_hist(hist_dat,domain_time,hist_value_list,hist_values);
//    dataput.output_result_3d(correlation_dat,domain_time,domain_sight,correlation_list,results_2d);
//    dataput.output_result_3d(spatial_dat,domain_time,domain_sight,spatial_list,results_2d);
    dataput.output_result_3d(normalmode_dat,domain_time,domain_normalmode,normalmode_value_list,results_2d);
    dataput.output_result_3d(toda_conservations_dat,domain_time,domain_toda_conservations,toda_conservations_list,results_2d);
    dataput.output_result_3d(toda_eigenvalues_dat,domain_time,domain_toda_eigenvalues,toda_eigenvalues_list,results_2d);
    dataput.output_result(stat_dat,domain_time,stat_values);
    dataput.time_tag() << " end time " << std::endl; 
  }
  mpi_error = MPI_Finalize();
  return 0;
} //END clXYmodel

