#include "simulation.h"
#include "IS.h"
#include "boundaryConds.h"
#include "serialEnv.h"
#include "initFunc.h"
#include "simData.h"

#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(99);
  int ny(73);
  int nz(5);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(2.0);
 // bool output(false);
 // int nreports(5);
  double endTime=0.0;
  SerialEnv env(&argc, &argv, 1, 1, 1);
  double tol=1e-14;

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sNg(Ng);
  const std::vector<double> toy_params { {5.0e-2, 2.0e-1, 2.0e-1, 2.0e-1, 0, 2.0e-1} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);
  Data * d(&data);

  // Choose particulars of simulation
  IS model(&data);

  // Outflow bcs(&data);
  Periodic bcs(&data);
  Simulation sim(&data, &env);

//  IS_Shocktube_1D init(&data, 0);
//  d->gamma = 4.0/3.0;
//  ISKHInstabilitySingleFluid init(&data, 0);
  IS_C2PStressTest init(&data);

  model.primsToAll(data.cons, data.prims, data.aux);

  // Copy initial primitives to compare to after C2P
  double * init_prims = (double *) malloc(sizeof(double)*6*d->Nx*d->Ny*d->Nz);
  memcpy(init_prims, data.prims, sizeof(double)*6*d->Nx*d->Ny*d->Nz);

  // Add noise to the prims so the C2P must work to make them correct
  srand (time(NULL));

  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

        for(int prim_n(0); prim_n < 6; prim_n++){
          int noise_factor = rand() % 10; // Between 0 and 9
          // create some randomness in whether the noise is positive or negative
          if (noise_factor % 2 == 0) {
          data.prims[ID(prim_n, i, j, k)] += data.prims[ID(prim_n, i, j, k)]*(noise_factor/100);
          } else {
          data.prims[ID(prim_n, i, j, k)] -= data.prims[ID(prim_n, i, j, k)]*(noise_factor/100);
          }
        }
      
      }
    }    
  }  

  model.getPrimitiveVars(data.cons, data.prims, data.aux);

  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
        
        for(int prim_n(0); prim_n < 6; prim_n++){
          // horrible hack workaround for abs not working
          if ( (data.prims[ID(prim_n, i, j, k)] - init_prims[ID(prim_n, i, j, k)]) > tol
              || (data.prims[ID(prim_n, i, j, k)] - init_prims[ID(prim_n, i, j, k)]) < -tol ) {
            printf("diff (C2P - init): %.17g \n", data.prims[ID(prim_n, i, j, k)] - init_prims[ID(prim_n, i, j, k)]);
            printf("pos (prim_n,i,j,k): %i, %i, %i, %i \n", prim_n, i, j, k);
          } 
        }
      
      }
    }    
  }

  free(init_prims);

  return 0;

}
