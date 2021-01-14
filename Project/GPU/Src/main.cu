// Serial main
#include "simData.h"
#include "parallelCheckpointArgs.h"
#include "simulation.h"
#include "initFunc.h"
#include "initFuncFromCheckpoint.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "parallelBoundaryConds.h"
#include "rkSplit.h"
#include "SSP2.h"
#include "parallelSaveDataHDF5.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <omp.h>

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))

using namespace std;

int main(int argc, char *argv[]) {


  // Set up domain
  int Ng(4);
  int nx(64);
  int ny(16);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(-1.5);
  double zmax(1.5);
  //double endTime(0.0005);
  double endTime(0.01);
  double cfl(0.1);
  double gamma(4.0/3.0);
  double sigma(0);
  bool output(true);
  int safety(180);
  int nxRanks(2);
  int nyRanks(2);
  int nzRanks(1);

  char * ptr(0);
  //! Overwrite any variables that have been passed in as main() arguments
  for (int i(0); i < argc; i++) {
    if (strcmp(argv[i], "sigma") == 0) {
      sigma = (double)strtol(argv[i+1], &ptr, 10);
    }
  }

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  const char* filename = "data_t0.checkpoint.hdf5";

  ParallelCheckpointArgs checkpointArgs(filename, &env);
  checkpointArgs.endTime=endTime;

  Data data(checkpointArgs, &env);

  //Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            //cfl, Ng, gamma, sigma);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  ParallelFlow bcs(&data, &env);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluid init(&data, 1);

  RK2 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_parallel", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime()- startTime);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
