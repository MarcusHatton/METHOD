// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "IS.h"
#include "boundaryConds.h"
// #include "rkSplit.h"
// #include "backwardsRK.h"
#include "RKPlus.h"
// #include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include "serialSaveDataHDF5.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  // int nx(65536);
  // int nx(32768);
  int nx(400);
  int ny(0);
  int nz(0);
  double xmin(-2.0);
  double xmax(3.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(2.0);
  double cfl(0.1);
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
  bool output(true);
  int nreports(5);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  data_args.gamma = 5.0/3.0;
  const std::vector<double> toy_params           { {1.0e-15, 1.0e-1,  1.0e-2, 5.0e-2,  1.0e-15, 1.0e-1} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);
  data_args.reportItersPeriod = 1;

  Data data(data_args, &env);

  // Choose particulars of simulation
  IS model(&data); // bool for alt_C2P

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  Outflow bcs(&data);
  // Periodic bcs(&data);

  Simulation sim(&data, &env);

  IS_Shocktube_1D init(&data, 0); //direction given by second arg (int)
  // Blob2dToyQ init(&data);
  //ISKHInstabilitySingleFluid init(&data, 1);
  //Shocktube_Chab21 init(&data);  

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
  // SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  RK2B timeInt(&data, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&data, &env, "1d/bulk/ds_"+std::to_string(nx)+"_0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/bulk/ds_"+std::to_string(nx)+"_"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  return 0;

}
