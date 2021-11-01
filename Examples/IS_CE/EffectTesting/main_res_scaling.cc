// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "ISCE.h"
#include "boundaryConds.h"
// #include "rkSplit.h"
// #include "backwardsRK.h"
#include "RKPlus.h"
// #include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "DEIFY.h"
#include "serialEnv.h"
#include "serialSaveDataHDF5.h"
#include "weno.h"
#include <cstring>
#include "sys/stat.h"

using namespace std;

int main(int argc, char *argv[]) {

  int nxs[] = {1000,2000,4000,8000,16000};
  int nx = 0;

  for(int i=0; i<5; i++) {
    nx = nxs[i];
    cout << nx << std::endl;
    std::string dirpath = "./1d/shear/res/"+std::to_string(nx);
    mkdir(dirpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
  // Set up domain
  int Ng(4);
  //int nx(3200);
  int ny(0);
  int nz(0);
  double xmin(-10.0);
  double xmax(10.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.8);
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
  bool output(false);
  int nreports(50);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  data_args.gamma = 5.0/3.0;
  const std::vector<double> toy_params           { {1.0e-15, 1.0e-1,  1.0e-15, 1.0e-1,  1.0e-3, 1.0e-1} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  ISCE model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  DEIFY ModelExtension(&data, &fluxMethod);

  Outflow bcs(&data);
  // Periodic bcs(&data);

  Simulation sim(&data, &env);

  //ISCE_Shocktube_1D_Para init(&data, 0); //direction given by second arg (int)
  //ISCE_Shocktube_1D_Perp init(&data, 0); // para = v aligned with dir, 
					 // perp = v2 always non-trivial one
  // Blob2dToyQ init(&data);
  //ISKHInstabilitySingleFluid init(&data, 1);
  //Shocktube_Chab21 init(&data);  
  IS_ShearTest init(&data);
  //IS_BulkHeatTest init(&data);

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
  // SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  RK2B timeInt(&data, &model, &bcs, &fluxMethod, &ModelExtension);

  SerialSaveDataHDF5 save(&data, &env, "1d/shear/res/"+std::to_string(nx)+"/data_serial_TIx_0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/shear/res/"+std::to_string(nx)+"/data_serial_TIx_"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  } // end of for-loop

  return 0;

}
