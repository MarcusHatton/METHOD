#include "gtest/gtest.h"
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "SSP2.h"
#include "SSP3.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "weno.h"
#include "wenoUpwinds.h"
#include "serialEnv.h"
#include "serialSaveData.h"
#include <cstdlib>
#include <cstdio>

TEST(SSP2, IMEX2BenchmarkForParallelCode)
{
  /*
    Run a resistive MHD problem to test the IMEX2 scheme. This test will run
    the simulation and save the output data in the TestData directory, ready
    to be compared with the parallel output.
  */
  double sigma(0);

  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(64, 16, 0, 0, 1, 0, 1, 0, 1, 0.05, &env,
            0.5, 4, 4.0/3.0, sigma);

  // Choose particulars of simulation
  SRRMHD model(&data);
  Weno3 weno(&data);
  FVS fluxMethod(&data, &weno, &model);
  Outflow bcs(&data);
  Simulation sim(&data, &env);
  BrioWuSingleFluid init(&data);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  SerialSaveData save(&data, &env);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "SSP2");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}


TEST(SSP3, IMEX3BenchmarkForParallelCode)
{
  /*
    Run a resistive MHD problem to test the IMEX3 scheme. This test will run
    the simulation and save the output data in the TestData directory, ready
    to be compared with the parallel output.
  */
  double sigma(0);

  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(64, 16, 0, 0, 1, 0, 1, 0, 1, 0.05, &env,
            0.5, 4, 4.0/3.0, sigma);

  // Choose particulars of simulation
  SRRMHD model(&data);
  Weno3 weno(&data);
  FVS fluxMethod(&data, &weno, &model);
  Outflow bcs(&data);
  Simulation sim(&data, &env);
  BrioWuSingleFluid init(&data);
  SSP3 timeInt(&data, &model, &bcs, &fluxMethod);
  SerialSaveData save(&data, &env);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "SSP3");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
