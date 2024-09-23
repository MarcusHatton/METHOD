// Navier-Stokes-like simulations for
// initial testing of sub-grid closure
// from coviariant filtering scheme

#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "SubgridNS.h"
//#include "boundaryConds.h"
#include "parallelBoundaryConds.h"
// #include "rkSplit.h"
// #include "backwardsRK.h"
// #include "RKPlus.h"
#include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "parallelEnv.h"
//#include "serialEnv.h"
//#include "serialSaveDataHDF5.h"
#include "parallelSaveDataHDF5.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  // int nx(65536);
  // int nx(32768);
  int nx(800);
  int ny(800);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(20.0);
  double cfl(0.4);
  // double gamma(0.001);
  // double sigma(0.001);
  // These parameters work with IMEX SSP2; given that tau_q << dt,
  // we should not expect them to work with the explicit solver, and indeed
  // it fails very quickly.
  // Note it's the ratio that matters to it being stable with IMEX.
  // Need gamma/sigma <~ 1 or so, o/w wavespeed will be too
  // big, and things fail. At moderate resolution you can get away with
  // a bigger ratio, but once resolution gets above around 1024 this needs
  // to be compatible with the CFL limit.
  //
  // With really steep initial data there can be minor Gibbs oscillation
  // effects, but even at crazy resolutions (65k) these are small provided
  // the CFL limit is met.
  bool output(false);
  int nreports(10);
  std::string output_dir = "SubGrid/2d/KHRandom/";

  ParallelEnv env(&argc, &argv, 8, 5, 1);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  data_args.gamma = 4.0/3.0;
  data_args.reportItersPeriod = 2000;
  // These of course should no longer be used, but can leave them for now as 'scaling parameters'
  // Perhaps to control the filtering length effect!
  const std::vector<double> toy_params           { {1.0, 5.0e-1,  1.0, 5.0e-1,  1.0, 5.0e-1, 4.0} }; 
  // scale ratio currently does nothing - it is set from the ratio of resolutions to the 800x800 calibration one. e.g. s_r=4.0 for 200x200 subgrid sim
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi", "scale_ratio"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  NS model(&data, false);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  //ISKHInstabilitySingleFluid init(&data, 1);
  KHRandomInstabilitySingleFluid init(&data); // optional:magnetic fields, seed

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  // RK2B timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, output_dir+"dp_"+std::to_string(nx)+"x"+std::to_string(ny)+"x"+std::to_string(nz)+"_0", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    ParallelSaveDataHDF5 save_in_loop(&data, &env, output_dir+"dp_"+std::to_string(nx)+"x"+std::to_string(ny)+"x"+std::to_string(nz)+"_"+std::to_string(n+1), ParallelSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  return 0;

}
