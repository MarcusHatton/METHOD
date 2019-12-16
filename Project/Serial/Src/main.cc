/**
 * @Author: Alex James Wright <alex>
 * @Date:   2019-09-30T15:33:00+01:00
 * @Email:  alex.j.wright2@gmail.com
 * @Last modified by:   alex
 * @Last modified time: 2019-10-11T13:26:33+01:00
 * @License: MIT
 */



// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "SSP2.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "REGIME.h"
#include "hybrid.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))

using namespace std;

int main(int argc, char *argv[]) {


  const double MU(1000);
  // Set up domain
  int Ng(4);
  int nx(256);
  int ny(256);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1);
  double ymax(1);
  double zmin(-1.0);
  double zmax(1.0);
  double endTime(3.0);
  double cfl(0.4);
  double gamma(4.0/3.0);
  double sigma(500);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(40);
  bool output(true);
  int safety(frameSkip);
  bool functionalSigma(true);
  double gam(0.2);
  double sigmaCrossOver(600);
  double sigmaSpan(500);
  bool useREGIME(false);


  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip,
            functionalSigma, gam);

  // Choose particulars of simulation
  // SRRMHD model(&data);
  Hybrid model(&data, sigmaCrossOver, sigmaSpan, useREGIME);

  FVS fluxMethod(&data, &model);

  model.setSubgridModel(&fluxMethod);

  Simulation sim(&data);

  KHInstabilitySingleFluid init(&data, 1);

  Flow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod, NULL);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  sim.evolve(output, safety);
  // sim.updateTime();
  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
