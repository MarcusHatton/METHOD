// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_q.h"
//#include "boundaryConds.h"
#include "parallelBoundaryConds.h"
// #include "rkSplit.h"
// #include "backwardsRK.h"
#include "SSP2.h"
#include "fluxVectorSplitting.h"
<<<<<<< HEAD
#include "serialEnv.h"
// #include "serialSaveDataHDF5.h"
#include "serialSaveData.h"
=======
//#include "serialEnv.h"
//#include "parallelEnv.h"
#include "platformEnv.h"
//#include "serialSaveDataHDF5.h"
#include "parallelSaveDataHDF5.h"
>>>>>>> 33bc5c1490bc7a4e27ac9c21bfd06e9b9d8e09a1
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  // int nx(65536);
  // int nx(32768);
  int nx(2048);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
<<<<<<< HEAD
  double endTime(5.0);
  double cfl(0.4);
=======
  double endTime(0.1);
  double cfl(0.1);
>>>>>>> 33bc5c1490bc7a4e27ac9c21bfd06e9b9d8e09a1
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
<<<<<<< HEAD
  int nreports(5);
=======
  int nreports(100);
>>>>>>> 33bc5c1490bc7a4e27ac9c21bfd06e9b9d8e09a1

  //SerialEnv env(&argc, &argv, 1, 1, 1);
  ParallelEnv env(&argc, &argv, 2, 1, 1);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  const std::vector<double> toy_params { {1.0e-3, 1.0e-2} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q"};
  const int n_toy_params(2);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  ToyQ model(&data);
  //ToyQFunctional model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  //Outflow bcs(&data);
  ParallelOutflow bcs(&data, &env);
  //Periodic bcs(&data);

  Simulation sim(&data, &env);

  BlobToyQ init(&data, 1.0); // turn on or off initial flux
  // Blob2dToyQ init(&data);

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

<<<<<<< HEAD
  //SerialSaveDataHDF5 save(&data, &env, "1d/data_1em4_serial0", SerialSaveDataHDF5::OUTPUT_ALL);  
  //SerialSaveData save(&data, &env, "1d/data_1em4_serial0", SerialSaveData::saveAll);  
=======
  //SerialSaveDataHDF5 save(&data, &env, "1d/Initial_Flux_Test/data_1em4_serial_0", SerialSaveDataHDF5::OUTPUT_ALL);
  ParallelSaveDataHDF5 save(&data, &env, "1d/Initial_Flux_Test/data_1em4_serial_0", ParallelSaveDataHDF5::OUTPUT_ALL);
  //SerialSaveDataHDF5 save(&data, &env, "1d/data_1em4_serial0", SerialSaveDataHDF5::OUTPUT_ALL);
>>>>>>> 33bc5c1490bc7a4e27ac9c21bfd06e9b9d8e09a1

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, nullptr);

  //save.saveAll();

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
<<<<<<< HEAD
    //SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/data_1em4_serial"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    //SerialSaveData save_in_loop(&data, &env, "1d/data_1em4_serial"+std::to_string(n+1), SerialSaveData::saveAll);
=======
    //SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/Initial_Flux_Test/data_1em4_serial_"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    ParallelSaveDataHDF5 save_in_loop(&data, &env, "1d/Initial_Flux_Test/data_1em4_serial_"+std::to_string(n+1), ParallelSaveDataHDF5::OUTPUT_ALL);
    //SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/data_1em4_serial_"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
>>>>>>> 33bc5c1490bc7a4e27ac9c21bfd06e9b9d8e09a1
    sim.evolve(output);
    //save_in_loop.saveAll();
  }

  return 0;

}
