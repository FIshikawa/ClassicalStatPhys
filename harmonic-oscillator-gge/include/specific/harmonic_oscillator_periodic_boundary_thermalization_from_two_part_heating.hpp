#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <boost/lexical_cast.hpp>
#include <clstatphys/physics/harmonic_oscillator.hpp>
#include <common/settings_common.hpp>
#include <common/two_part_heating.hpp>

using Ensembler = ensemble::TwoPartHeating; 
using Hamiltonian = hamiltonian::HarmonicOscillator;

struct Settings : public SettingsCommon{
  double J = 1.0; //interaction constant;
  double heating_temperture = 2.0;
  int left1=0;
  int right1=1;
  int left2=2;
  int right2=3;

  Settings(int argc, char **argv, int & input_counter) 
    : SettingsCommon(argc, argv, input_counter) 
  { 
    set(argc, argv, input_counter);
  }
  Settings() = default;

  inline void set(int argc, char **argv, int & input_counter){
    if (argc > input_counter) J                   =  boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) heating_temperture  =  boost::lexical_cast<double>(argv[input_counter]);++input_counter;
    if (argc > input_counter) left1               =  boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) right1              =  boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) left2               =  boost::lexical_cast<int>(argv[input_counter]);++input_counter;
    if (argc > input_counter) right2              =  boost::lexical_cast<int>(argv[input_counter]);++input_counter;

    if( (right1 - left1 + 1) + (right2 - left2 + 1) > num_particles){
      std::cerr << "num_heating_particles must be lower than num_particles" << std::endl;
      std::exit(1);
    }
  }

  template <class Dataput>
  inline void declare(Dataput & dataput){
    SettingsCommon::declare(dataput);

    dataput <<  "<<System Depend Settings>> " << std::endl
            <<  "  " << Hamiltonian::name() << std::endl
            <<  "  " << Ensembler::name() << std::endl
            <<  "  " << Lattice::name() << std::endl
            << "  Coupling constant : J = " << J << std::endl
            << "  Heating temperture : heating_temperture = " << heating_temperture << std::endl
            << "  Heating fisrt part left : left1 = " << left1 << std::endl
            << "  Heating fisrt part right : right1= " << right1 << std::endl
            << "  Heating second part left : left2 = " << left2 << std::endl
            << "  Heating second part right : right2 = " << right2<< std::endl;
  }

  Hamiltonian hamiltonian(){
    Lattice lattice_t(Ns);

    std::vector<std::vector<int> > pair_table(
                                              num_particles,
                                              std::vector<int>(N_adj)
                                              ); 
    lattice_t.create_table(pair_table);

    Hamiltonian hamiltonian_t(
                              num_particles,
                              J,
                              pair_table,
                              N_adj
                              );
    return hamiltonian_t;
  }

  Ensembler ensembler(){
    Ensembler ensembler_t(
                          left1,
                          right1,
                          left2,
                          right2,
                          heating_temperture,
                          temperture,
                          num_particles
                          );
    return ensembler_t;
  }
  
};

#endif
