// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: InputConfig struct definition
//
// -----------------------------------------------------------------------------

#ifndef INPUT_CONFIG_H
#define INPUT_CONFIG_H

#include <string>

struct MeshConfig
{
  double x_left;
  double x_right;
  int n_elements;
};

struct FEMConfig
{
  int order;
};

struct PDEConfig
{
  std::string type;
  double advection_speed;
};

struct FluxConfig
{
  std::string type;
  double alpha;
};

struct TimeIntegratorConfig
{
  std::string type;
};

struct TimeStepControllerConfig
{
  std::string type;
  double dt;
  double cfl;
};

struct OutputConfig
{
  std::string directory;
  std::string prefix;
  double output_dt;
  bool write_initial = false;
  bool write_final = false;
};

struct InitialConditionConfig
{
  std::string expression;
};

struct RunConfig
{
  double final_time;
};

struct InputConfig
{
  MeshConfig                mesh;
  FEMConfig                 fem;
  PDEConfig                 pde;
  FluxConfig                flux;
  TimeIntegratorConfig      time_integrator;
  TimeStepControllerConfig  time_step_controller;
  OutputConfig              output;
  InitialConditionConfig    initial_condition;
  RunConfig                 run;
};

#endif
