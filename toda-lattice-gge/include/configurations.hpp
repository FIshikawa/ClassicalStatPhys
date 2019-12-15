#include<settings.hpp>
#include<results.hpp>
#include<finalize.hpp>
#include<sampling_process.hpp>

struct Parameters{
  ParametersCommon common;
  ParametersDepend depend;

  void declare(tools::DataRecorder & dataput){
    common.declare(dataput);
    depend.declare(dataput);
  }

  void set(int args, char **argv){
    int input_counter = 0;
    common.set(args, argv, input_counter);
    depend.set(args, argv, input_counter);
  }
}

struct Modules{
  Modules(Parameters parameters) : 
    common(parameters), depend(parameters) {}

  ModulesCommon common;
  ModulesDepend depend;
}

#include<frame.hpp>
