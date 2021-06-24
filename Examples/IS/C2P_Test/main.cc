#include "IS.h"
#include "boundaryConds.h"
#include "serialSaveDataHDF5.h"
#include "serialEnv.h"

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
  bool output(false);
  int nreports(5);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  const std::vector<double> toy_params { {5.0e-2, 2.0e-1, 2.0e-1, 2.0e-1, 0, 2.0e-1} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  IS model(&data);

  // Outflow bcs(&data);
  Periodic bcs(&data);

  IS_Shocktube_1D init(&data);

  auto init_prims = data->prims;

  model.getPrimitiveVars(&data.cons, &data.prims, &data.aux);

  return 0;

}
