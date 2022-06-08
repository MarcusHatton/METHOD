/**
 * @Author: Alex James Wright <alex>
 * @Date:   2019-09-30T15:33:00+01:00
 * @Email:  alex.j.wright2@gmail.com
 * @Last modified by:   alex
 * @Last modified time: 2020-02-05T15:48:38+00:00
 * @License: MIT
 */

// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "boundaryConds.h"
#include "rkSplit2ndOrder.h"
#include "fluxVectorSplitting.h"
#include "hybrid.h"
#include "serialSaveDataHDF5.h"
#include "serialEnv.h"
#include "weno.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {


  // Set up domain
  int Ng(4);
  int nx(400);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.4);
  double cfl(0.1);
  // double cp(1.0);
  // int frameSkip(40);
  bool output(false);
  int safety(-1);

  // HYBRID PARAMETERS
  double tau(1.0e-2);
  double tauCrossOver(1.0e-3);
  double tauSpan(5.0e-4);
  bool functionalTau(true);
  bool useDEIFY(true);
  int nreports(10);

  SerialEnv env(&argc, &argv, 1, 1, 1);


  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.gamma = 5.0/3.0;
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  data_args.sfunctionalTau(functionalTau);
  data_args.sTau(tau);
  const std::vector<double> toy_params           { {1.0e-15, 1.0e-1,  5.0e-3, 1.0e-1,  1.0e-15, 1.0e-1} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);
  Data data(data_args, &env);

  // Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
  //           cfl, Ng, xi, tau, cp, frameSkip);

  // Choose particulars of simulation
  Hybrid model(&data, tauCrossOver, tauSpan, useDEIFY);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  model.setupDEIFY(&fluxMethod);

  Outflow bcs(&data);

  Simulation sim(&data, &env);

  IS_Shocktube_1D_Para init(&data, 0); //direction given by second arg (int)

  RKSplit2 timeInt(&data, &model, &bcs, &fluxMethod, NULL);
  // RK2B timeInt(&data, &model, &bcs, &fluxMethod);

  // SerialSaveData save(&data, &env, 1);
  SerialSaveDataHDF5 save(&data, &env, "1d/bulk/data_serial_0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  // sim.evolve(output, safety);

  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/bulk/data_serial_TIx_"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  return 0;

}
