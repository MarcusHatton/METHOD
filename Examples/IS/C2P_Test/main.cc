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
  int nx(50);
  int ny(20);
  int nz(20);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
 // bool output(false);
 // int nreports(5);
  double endTime=0.0;
  SerialEnv env(&argc, &argv, 1, 1, 1);
  double tol=1e-10;

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

  IS_Shocktube_1D init(&data, 0);
  model.primsToAll(data.cons, data.prims, data.aux);

  double * init_prims = (double *) malloc(sizeof(double)*6*d->Nx*d->Ny*d->Nz);
  memcpy(init_prims, data.prims, sizeof(double)*6*d->Nx*d->Ny*d->Nz);

  model.getPrimitiveVars(data.cons, data.prims, data.aux);

  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
        
        for(int prim_n(0); prim_n < 6; prim_n++){
          if ( (data.prims[ID(prim_n, i, j, k)] - init_prims[ID(prim_n, i, j, k)]) > tol
              || (data.prims[ID(prim_n, i, j, k)] - init_prims[ID(prim_n, i, j, k)]) < -tol ) {
            printf("diff (C2P - init): %.17g \n", data.prims[ID(prim_n, i, j, k)] - init_prims[ID(prim_n, i, j, k)]);
            printf("pos (i,j,k): %i, %i, %i \n", i, j, k);
          }
        }
      
      }
    }    
  }

  free(init_prims);

  return 0;

}
