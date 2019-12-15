#include<iostream>
#include<cmath>
#include<fstream>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<boost/random.hpp>
#include<boost/random/uniform_int_distribution.hpp>
#include<string>

#include "initializer.hpp"
#include "dataput.hpp"
#include "hamiltonian.hpp"
#include "integrator.hpp"
#include "correlation.hpp"
#include "Lattice.hpp"
#include "calcerror.hpp"
#include "exefftw.hpp"
#include "histogramer.hpp"
#include "dissipated_integrator.hpp"
#include "accumulator.hpp"
#include "exactcalc.hpp"

//typedef Lattice::RegularCubeLasso lattice_t;
//typedef Lattice::ChainLasso lattice_t;
//typedef Lattice::SquareLatticeLasso lattice_t;
//typedef Lattice::RegularCube lattice_t;
typedef Lattice::ChainOpenBoundary lattice_t;
//typedef Lattice::SquareLattice lattice_t;
//typedef initializer::equilclXYmodel_impurity initializer_t;
//typedef initializer::equilclXYmodelWolf initializer_t;
//typedef initializer::equilclXYmodel initializer_t;
typedef initializer::steadyCauchy initializer_t;
//typedef initializer::HarmonicChainOB initializer_t;
//typedef initializer::equilLangevin initializer_t;
//typedef integrator::velver integrator_t;
typedef integrator::rk4 integrator_t;
//typedef integrator::euler integrator_t;
//typedef dissipated_integrator::rk2 integrator_t;
//typedef dissipated_integrator::rk2_cauchy integrator_t;
//typedef dissipated_integrator::euler_cauchy integrator_t;
typedef dataput::dataputForRegularCube dataput_t;
//typedef dataput::dataputForSquareLattice dataput_t;
typedef hamiltonian::classicalXY hamiltonian_t;
//typedef hamiltonian::classicalXY_impurity hamiltonian_t;
//typedef hamiltonian::HarmonicOscillator_impurity hamiltonian_t;
//typedef hamiltonian::classicalXYNonEq hamiltonian_t;
//typedef hamiltonian::classicalXYAnisoNonEq hamiltonian_t;
//typedef hamiltonian::HarmonicOscillator hamiltonian_t;
typedef stat::histogramer histogramer_t;
//typedef exactcalc::partition_function_harmonic_chain calcpartition_t;
typedef exactcalc::partition_function_clXY_chain calcpartition_t;



int main(int argc, char **argv){
 double J = 1.0; //std::pow(2.0*M_PI,2.0);
 double T = 0.5; //temperature
 double gamma = 1.0; //gamma
 int Ns = 1; // lenght of whole system
 int N_thermalize = 100; // number of thermalization
 int num_particles; //nummer of particles
 int N_loop = 1; //loop number
 int N_each; //each loop number
 int N_parallel = 1; //number of parallel 
 double t = 10; // Whole time
 double dt; //size of interval dt = t/N_time
 int N_time = 1000; //numer of time step
 int N_time_resolve = 1; //numer of time step
 char* condition_dat = "condi.dat";
 std::string result_dirr("./test");
 unsigned int n_bin = 1000;
 double strength = (double)num_particles; // range of the uniform random 
 int N_T_range = 100; // number of T range
 double T_max = 2.0; //max of T range
 double T_min = 0.05; //min of T range
 double dT; // inteval of T range
 //unsigned int n_bin = N_loop;


//parameter set 
 if (argc > 1) result_dirr =      boost::lexical_cast<std::string>(argv[1]);
 if (argc > 2) N_thermalize =     boost::lexical_cast<int>(argv[2]);
 if (argc > 3) Ns =               boost::lexical_cast<int>(argv[3]);
 if (argc > 4) T =                boost::lexical_cast<double>(argv[4]);
 if (argc > 5) gamma =            boost::lexical_cast<double>(argv[5]);
 if (argc > 6) N_time =           boost::lexical_cast<int>(argv[6]);
 if (argc > 7) t =                boost::lexical_cast<double>(argv[7]);
 if (argc > 8) n_bin =            boost::lexical_cast<unsigned int>(argv[8]);
 if (argc > 9) N_loop =           boost::lexical_cast<int>(argv[9]);
 if (argc > 10) N_parallel =      boost::lexical_cast<int>(argv[10]);
 if (argc > 11) N_time_resolve =  boost::lexical_cast<int>(argv[11]);
 if (argc > 12) strength  =       boost::lexical_cast<double>(argv[12]);
 if (argc > 13) N_T_range   =     boost::lexical_cast<int>(argv[13]);
 if (argc > 14) T_max   =         boost::lexical_cast<double>(argv[14]);
 if (argc > 15) T_min   =         boost::lexical_cast<double>(argv[15]);

// parameter auto set
 N_each = N_loop/N_parallel;
 num_particles = lattice_t::set_num_particles(Ns);
 dt = t/(double)N_time;
 int N_hist_length = N_time / N_time_resolve;
 strength = (double)num_particles;
 dT = (T_max - T_min) / (double)(N_T_range-1);

//parameter set 
 std::vector<std::string> experimental_condition(6);
 experimental_condition[0] = "Equilibration check : classical XY dynamics";
 experimental_condition[1] = hamiltonian_t::name();
 experimental_condition[2] = initializer_t::name();
 experimental_condition[3] = integrator_t::name();
 experimental_condition[4] = lattice_t::name();
 experimental_condition[5] = calcpartition_t::name();
 dataput_t::startput(condition_dat, experimental_condition);

//condition declare
 std::ofstream ofs;
 ofs.open(condition_dat,std::ios::app);
 if(!ofs)std::cerr<< condition_dat <<"opening failed"<<std::endl;
 ofs  << "Coupling constant: J is " << J << std::endl
      << "gamma : gamma = " << gamma << std::endl
      << "Temperture: T = " << T << std::endl
      << "Reduced Coupling: J/T = " << J/T << std::endl
      << "Number of patricles :num_particles = " << num_particles << std::endl
      << "Length of whole system :Ns = " << Ns << std::endl
      << "Number of thermalize step : N_thermalize = " << N_thermalize << std::endl
      << "Developing ime : t = " << t << std::endl
      << "Number of time steps : N_time = " << N_time << std::endl
      << "Number of time resolve for hist : N_time_resolve = " << N_time_resolve << std::endl
      << "Time interbal : dt = " << dt << std::endl
      << "Number of bins: n_bin =" << n_bin << std::endl
      << "Number of loop: N_loop =" << N_loop << std::endl
      << "Parallelize : N_parallel = " << N_parallel << std::endl
      << "Each thread loop : N_each = " << N_each << std::endl 
      << "Strength of the uniform random : strength = " << strength << std::endl 
      << "Number of T range : N_T_range = " << N_T_range << std::endl 
      << "Max T range : T_max = " << T_max << std::endl 
      << "min T range : T_min = " << T_min << std::endl 
      << "Interval of T range : dT = " << dT << std::endl
      << "Result into : " << result_dirr << std::endl;
 ofs.close();

//dataputer set 
 dataput_t dataputer(Ns,num_particles);
 
 std::vector< std::vector<double> >  
		              S(N_parallel, std::vector<double>(N_time)),
		              V_total(N_parallel, std::vector<double>(N_time)),
		              eS(N_parallel, std::vector<double>(N_time)),
		              V(N_parallel, std::vector<double>(N_time)),
		              eV(N_parallel, std::vector<double>(N_time)),
		              E(N_parallel, std::vector<double>(N_time)),
		              eE(N_parallel, std::vector<double>(N_time)),
		              K(N_parallel, std::vector<double>(N_time)),
		              eK(N_parallel, std::vector<double>(N_time)),
		              U(N_parallel, std::vector<double>(N_time)),
		              eU(N_parallel, std::vector<double>(N_time)),
		              S_total(N_parallel, std::vector<double>(N_time)),
		              eS_total(N_parallel, std::vector<double>(N_time)),
		              S_totalACF(N_parallel, std::vector<double>(N_time)),
		              eS_totalACF(N_parallel, std::vector<double>(N_time)),
		              VACF(N_parallel, std::vector<double>(N_time)),
		              SACF(N_parallel, std::vector<double>(N_time)),
		              eVACF(N_parallel, std::vector<double>(N_time)),
		              eSACF(N_parallel, std::vector<double>(N_time));
 std::vector< std::vector< std::vector<double> > >
                  SSCF(N_parallel, std::vector<std::vector<double> >(Ns, std::vector<double>(N_time,0)));
 std::vector< std::vector<double> >  
		              hist_data_v(N_hist_length, std::vector<double>(N_loop)),
		              hist_data_S(N_hist_length, std::vector<double>(N_loop));

 std::vector< std::vector<double> >
                  distdistance(N_parallel,std::vector<double>(N_T_range,0.0));

 
 #pragma omp parallel for num_threads(N_parallel) 
 for(int p = 0; p < N_parallel ; ++p){

 //pararrel tag number def
 int parallel_tag = p;

 //lattice set
 lattice_t lattice(Ns);
 unsigned int N_adj;// number of adjacent spins
 N_adj = lattice.number_adjacent();
 int counter = 0.0 ; //montecalro counter 

 
 //set observe particle
 unsigned int observe = num_particles-1;

 //vector set 
 std::vector<double> z(2*num_particles);

 //set oberve vect
 std::vector<double> moment(N_time),spin(N_time),momentError(N_time),spinError(N_time),s_total(N_time),v_total(N_time),s_totalError(N_time);

 //table set 
 std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));
 lattice.create_table(pair_table);

 
 //randam genereta set 
 std::size_t seed = p;

 //seed = p;
 boost::random::mt19937 mt(seed);

 //initializer set
 initializer_t initializer(J,T,num_particles);

 //initial config set 
 initializer.initset(z,mt);

 //hamiltonian set
 hamiltonian_t ham(num_particles,J,pair_table,N_adj); 
 //initial energy calc

 //set calc partition function
 calcpartition_t calc_partition(J, num_particles);
 //set correlator
 correlation::acf acf_t(N_time,N_loop);
 //integrator set
 //integrator_t  integrator(T, gamma, 2*num_particles);
 integrator_t  integrator(2*num_particles);
 
 
 // calc each partition function
 //thermalize
 double pt = 0.0;
 int hist_step = 0;
 double v_total_init = 0.0;

 //loop calc
 for(int c_loop=0; c_loop < N_each; ++c_loop){
   //time develop
   int sight_x = 0;
   int sight_y = 0;
   int sight_z = 0;
   int target = 0;
   double *x = &z[0];
   double *v = &z[num_particles];
   v_total_init = 0.0;
   initializer.hotset(z,counter,mt,strength);
   initializer.Momentset(z,T,mt);
   for(int counter = 0; counter < N_thermalize; ++counter) initializer.montecalro(z, counter, ham, mt);
   for(int i = 0; i < num_particles; ++i){
     v_total_init += v[i];
   }
   v_total_init = v_total_init / num_particles;
   for(int i = 0; i < num_particles; ++i) v[i] -= v_total_init;
   for(int step = 0; step < N_time; ++step){  
     v_total[step] = v[0]/num_particles;
     for(int i = 1; i < num_particles; ++i) v_total[step] += v[i]/num_particles;
     V_total[p][step] = v_total[step];
     moment[step] = v[observe] - v_total[step] ;
     spin[step] = std::cos(x[observe]);
     momentError[step] = v[observe]*v[observe];
     spinError[step] = std::cos(x[observe])*std::cos(x[observe]);
     V[p][step] += v[observe]/N_loop;
     eV[p][step] += momentError[step]/N_loop;
     S[p][step] += spin[step]/N_loop;
     eS[p][step] += spinError[step]/N_loop;
     s_total[step] = std::cos(x[0])/num_particles;
     for(int i = 1; i < num_particles; ++i) s_total[step] += std::cos(x[i])/num_particles;
     S_total[p][step] += s_total[step]/N_loop;
     eS_total[p][step] += std::pow(s_total[step],2.0)/N_loop;
     s_totalError[step] = std::pow(s_total[step],2.0);
     E[p][step] += ham.energy(pt,z) / N_loop;
     eE[p][step] += ham.energy(pt,z) * ham.energy(pt,z) / N_loop;
     K[p][step] += ham.kinetic_energy(pt,z) / N_loop;
     eK[p][step] += ham.kinetic_energy(pt,z) * ham.kinetic_energy(pt,z) / N_loop;
     U[p][step] += ham.potential_energy(pt,z) / N_loop;
     eU[p][step] += ham.potential_energy(pt,z) * ham.potential_energy(pt,z) / N_loop;
     for(int i = 0; i < Ns; ++i){
      sight_x = lattice.latticize(0,observe);
      //sight_y = lattice.latticize(1,observe);
      //sight_z = lattice.latticize(2,observe);
      //target = lattice.numberize(sight_x + i, sight_y);
      target = lattice.numberize(sight_x + i);
      SSCF[p][i][step] += std::cos(x[observe] - x[target])/N_loop;
     }
     pt += dt;
     if((step % N_time_resolve) == 0){
      hist_step = step / N_time_resolve;
      hist_data_v[hist_step][p * N_each + c_loop] = v[observe];
      hist_data_S[hist_step][p * N_each + c_loop] = s_total[step]; 
     }
     if(step == 0){
      for(int i = 0; i < N_T_range; ++i){
        double T_t = T_min + (double)i * dT;       
        double P_value = 0.0;
        P_value = ham.energy(pt,z)/T_t + calc_partition.logZ(T_t);
        P_value = P_value / N_loop;
        distdistance[p][i] += P_value;
      }
     }
     integrator.step(pt, dt, z, ham);
   }//end time step

   //calc ACF
   acf_t.calc(moment,VACF[p]);
   acf_t.calc(spin,SACF[p]);
   acf_t.calc(momentError,eVACF[p]);
   acf_t.calc(spinError,eSACF[p]);
   acf_t.calc(s_total,S_totalACF[p]);
   acf_t.calc(s_totalError,eS_totalACF[p]);
   
   //thermalize
 }//end loop

}//end parallerize

 //sum for parallelization
  for(int step = 0; step < N_time;++step){
    for(int p = 1; p < N_parallel; ++p){
      V[0][step] += V[p][step];
      eV[0][step] += eV[p][step];
      S[0][step] += S[p][step];
      eS[0][step] += eS[p][step];
      E[0][step] += E[p][step];
      eE[0][step] += eE[p][step];
      K[0][step] += K[p][step];
      eK[0][step] += eK[p][step];
      U[0][step] += U[p][step];
      eU[0][step] += eU[p][step];
      VACF[0][step]+= VACF[p][step];
      SACF[0][step] += SACF[p][step];
      eVACF[0][step] += eVACF[p][step];
      eSACF[0][step]+= eSACF[p][step];
      V_total[0][step] += V_total[p][step];
      S_total[0][step]+= S_total[p][step];
      eS_total[0][step] += eS_total[p][step];
      S_totalACF[0][step] += S_totalACF[p][step];
      eS_totalACF[0][step] += eS_totalACF[p][step];
      for(int j = 0; j < Ns; ++j) SSCF[0][j][step] += SSCF[p][j][step];
    }
  }

 calcerror::serrorForVec calcerror_t(N_time,N_loop);
 calcerror_t.calc(eV[0],V[0]);
 calcerror_t.calc(eS[0],S[0]);
 calcerror_t.calc(eE[0],E[0]);
 calcerror_t.calc(eK[0],K[0]);
 calcerror_t.calc(eU[0],U[0]);
 calcerror_t.calc(eVACF[0],VACF[0]);
 calcerror_t.calc(eSACF[0],SACF[0]);
 calcerror_t.calc(eS_total[0],S_total[0]);
 calcerror_t.calc(eS_totalACF[0],S_totalACF[0]);

 std::vector<double> fVACF(N_time),ifVACF(N_time),fSACF(N_time),ifSACF(N_time);//,fS_totalACF(N_time),ifS_totalACF(N_time);
 exefftw::fftranser ft_t(N_time);
 ft_t.trans(VACF[0],fVACF,ifVACF);
 ft_t.trans(SACF[0],fSACF,ifSACF);
// ft_t.trans(S_totalACF[0],fS_totalACF,ifS_totalACF);

 //dataput
 //list for output
 double* observable_pointer[23];
 std::vector<std::string> observable_list(23);
 double* observable_pointer_correlation[Ns];
 std::vector<std::string> observable_list_correlation(Ns);

 int tag = 0;
 observable_pointer[tag] = &V[0][0];
 observable_list[tag] = "V";
 tag += 1;
 observable_pointer[tag] = &eV[0][0];
 observable_list[tag] = "eV";
 tag += 1;
 observable_pointer[tag] = &S[0][0];
 observable_list[tag] = "S";
 tag += 1;
 observable_pointer[tag] = &eS[0][0];
 observable_list[tag] = "eS";
 tag += 1;
 observable_pointer[tag] = &VACF[0][0];
 observable_list[tag] = "VACF";
 tag += 1;
 observable_pointer[tag] = &eVACF[0][0];
 observable_list[tag] = "eVACF";
 tag += 1;
 observable_pointer[tag] = &SACF[0][0];
 observable_list[tag] = "SACF";
 tag += 1;
 observable_pointer[tag] = &eSACF[0][0];
 observable_list[tag] = "eSACF";
 tag += 1;
 observable_pointer[tag] = &fVACF[0];
 observable_list[tag] = "fVACF";
 tag += 1;
 observable_pointer[tag] = &ifVACF[0];
 observable_list[tag] = "ifVACF";
 tag += 1;
 observable_pointer[tag] = &fSACF[0];
 observable_list[tag] = "fSACF";
 tag += 1;
 observable_pointer[tag] = &ifSACF[0];
 observable_list[tag] = "ifSACF";
 tag += 1;
 observable_pointer[tag] = &V_total[0][0];
 observable_list[tag] = "V_total";
 tag += 1;
 observable_pointer[tag] = &S_total[0][0];
 observable_list[tag] = "S_total";
 tag += 1;
 observable_pointer[tag] = &eS_total[0][0];
 observable_list[tag] = "eS_total";
 tag += 1;
 observable_pointer[tag] = &S_totalACF[0][0];
 observable_list[tag] = "S_totalACF";
 tag += 1;
 observable_pointer[tag] = &eS_totalACF[0][0];
 observable_list[tag] = "eS_totalACF";
 tag += 1;
 observable_pointer[tag] = &E[0][0];
 observable_list[tag] = "E";
 tag += 1;
 observable_pointer[tag] = &eE[0][0];
 observable_list[tag] = "eE";
 tag += 1;
 observable_pointer[tag] = &K[0][0];
 observable_list[tag] = "K";
 tag += 1;
 observable_pointer[tag] = &eK[0][0];
 observable_list[tag] = "eK";
 tag += 1;
 observable_pointer[tag] = &U[0][0];
 observable_list[tag] = "U";
 tag += 1;
 observable_pointer[tag] = &eU[0][0];
 observable_list[tag] = "eU";
 tag += 1;
 for(int i = 0; i < Ns; ++i){
  observable_list_correlation[i] = " ";
  observable_pointer_correlation[i] = &SSCF[0][i][0];
 }

 std::vector< std::vector<double> >  
		              hist_v(N_hist_length, std::vector<double>(n_bin)),
		              hist_S(N_hist_length, std::vector<double>(n_bin)),
		              range_v(N_hist_length, std::vector<double>(n_bin)),
		              range_S(N_hist_length, std::vector<double>(n_bin));

 std::vector< std::vector<double> > stat_value(11, std::vector<double>(N_hist_length));
 stat::accumulator accum("stat variable");
 
 double* observable_pointer_stat[11];
 std::vector<std::string> observable_list_stat(11);

 observable_list_stat[0] = "mean";
 observable_list_stat[1] = "variance";
 observable_list_stat[2] = "skewness";
 observable_list_stat[3] = "kurtosis";
 observable_list_stat[4] = "kurtosis_excess";
 observable_list_stat[5] = "average";
 observable_list_stat[6] = "error";
 observable_list_stat[7] = "central_moment1";
 observable_list_stat[8] = "central_moment2";
 observable_list_stat[9] = "central_moment3";
 observable_list_stat[10] = "central_moment4";

 double* observable_pointer_hist[4*N_hist_length];
 std::vector<std::string> observable_list_hist(4*N_hist_length);
 tag = 0;
 int label = 0;
 int label_t = 0;
 int tag_t = 2*N_hist_length;
 //histogram 
 histogramer_t histogram_v(n_bin);
 histogramer_t histogram_S(n_bin);
 for(int i = 0; i < N_hist_length; ++i){
  for(int j = 0; j < N_loop; ++j){
    accum << hist_data_v[i][j];
  }

  stat_value[0][i] = accum.mean();
  stat_value[1][i] = accum.variance();
  stat_value[2][i] = accum.skewness();
  stat_value[3][i] = accum.kurtosis();
  stat_value[4][i] = accum.kurtosis_excess();
  stat_value[5][i] = accum.average();
  stat_value[6][i] = accum.error();
  stat_value[7][i] = accum.central_moment1();
  stat_value[8][i] = accum.central_moment2();
  stat_value[9][i] = accum.central_moment3();
  stat_value[10][i] = accum.central_moment4();
  accum.reset();   
  histogram_v(hist_data_v[i],N_loop);
  histogram_v.output(hist_v[i],range_v[i]);
  histogram_S(hist_data_S[i],N_loop);
  histogram_S.output(hist_S[i],range_S[i]);
  histogram_v.reset();
  histogram_S.reset();
  observable_pointer_hist[tag] = &range_v[label][0];
  observable_list_hist[tag] = "bin_v of each time";
  tag += 1;
  observable_pointer_hist[tag] = &hist_v[label][0];
  observable_list_hist[tag] = "hist_v of each time";
  tag += 1;
  label += 1;
  observable_pointer_hist[tag_t] = &range_S[label_t][0];
  observable_list_hist[tag_t] = "bin_S of each time";
  tag_t += 1;
  observable_pointer_hist[tag_t] = &hist_S[label_t][0];
  observable_list_hist[tag_t] = "hist_S of each time";
  tag_t +=1;
  label_t += 1;
 }


 for(int i = 0; i < 11; ++i){
  observable_pointer_stat[i] = &stat_value[i][0];
 } 

 double* observable_pointer_dist[2];
 std::vector<std::string> observable_list_dist(2);
 
 observable_list_dist[0] = "T range";
 observable_list_dist[1] = "dist distance";

 std::vector<double> T_range(N_T_range,0.0);
 for(int i = 0; i < N_T_range; ++i){
  T_range[i] = T_min + dT * (double)i;
  for(int j = 1; j < N_parallel; ++j){
    distdistance[0][i] += distdistance[j][i];
  }
 }

 observable_pointer_dist[0] = &T_range[0];
 observable_pointer_dist[1] = &distdistance[0][0];
 


 std::string result_dat("/result.dat");
 std::string correlation_dat("/result_correlation.dat");
 std::string hist_dat("/result_hist.dat");
 std::string stat_dat("/result_stat.dat");
 std::string dist_dat("/result_dist.dat");
 result_dat = result_dirr + result_dat;
 correlation_dat = result_dirr + correlation_dat;
 hist_dat = result_dirr + hist_dat;
 stat_dat = result_dirr + stat_dat;
 dist_dat = result_dirr + dist_dat;

 std::cout << "result output to : " << result_dat << " , " << correlation_dat << " , " << hist_dat << " , " << stat_dat << " and " << dist_dat <<  std::endl;

 dataputer.output_result(result_dat.c_str(),condition_dat,"pt",dt,N_time,observable_pointer,observable_list);
 dataputer.output_result(correlation_dat.c_str(),condition_dat,"pt",dt,N_time,observable_pointer_correlation,observable_list_correlation);
 dataputer.output_result(hist_dat.c_str(),condition_dat,"-1bin",1,n_bin,observable_pointer_hist,observable_list_hist);
 dataputer.output_result(stat_dat.c_str(),condition_dat,"output_t",1,N_hist_length,observable_pointer_stat,observable_list_stat);
 dataputer.output_result(dist_dat.c_str(),condition_dat,"-1temperture",1,N_T_range,observable_pointer_dist,observable_list_dist);

 dataputer.endput(condition_dat);

return 0;
} //END clXYmodel

