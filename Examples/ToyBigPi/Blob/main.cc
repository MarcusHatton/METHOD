// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_Pi.h"
#include "boundaryConds.h"
//#include "rkSplit.h"
#include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include "serialSaveDataHDF5.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(2048);
  int ny(0);
  int nz(0);
  double xmin(-1.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(0.0);
  double endTime(1.0);
  double cfl(0.1);
  /*
  double gamma(0.01); // eta (shear vis)
  double sigma(0.1); // tau_pi any smaller than 0.01 is unstable...
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  */
  bool output(false);
  int nreports(10);

  SerialEnv env(&argc, &argv, 1, 1, 1);

//  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
//            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip, reportItersPeriod);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  const std::vector<double> toy_params { {0.15, 2.0} };
  const std::vector<std::string> toy_param_names = {"eta", "tau_pi"};
  const int n_toy_params(2);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  ToyPi model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  Outflow bcs(&data);

  Simulation sim(&data, &env);

  BlobToyPi init(&data);

//  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  
  SerialSaveDataHDF5 save(&data, &env, "1d/data_serial0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
//  double startTime(omp_get_wtime());

  // Run until end time and save results
  // sim.evolve(output);
  save.saveAll();

//  double timeTaken(omp_get_wtime() - startTime);

//  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);
//  printf("\nCompleted %d iterations.\n", data.iters);

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/data_serial"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }
  
  return 0;

}
