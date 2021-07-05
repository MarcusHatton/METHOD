// Serial main
#include "simData.h"
#include "simulation.h"
#include "serialEnv.h"
#include "boundaryConds.h"
#include "IS.h"
#include <cstring>
#include <math.h>
#include <cmath>

using namespace std;


int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  // int nx(65536);
  // int nx(32768);
  int nx(20);
  int ny(20);
  int nz(20);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.0001);
  double cfl(0.4);
//  bool output(false);
//  int nreports(1);

  double nxRanks(1);
  double nyRanks(1);
  double nzRanks(1);

  SerialEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  // const std::vector<double> toy_params { {0.0e-2, 2.0e-1, 2.0e-2, 2.0e-1, 0, 2.0e-1} };
  const std::vector<double> toy_params { {0.0, 2.0e-1, 0.0, 2.0e-1, 0, 2.0e-1} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q", "zeta", "tau_Pi", "eta", "tau_pi"};
  const int n_toy_params(6);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  IS model(&data, false); // alt_C2P?
  Outflow bcs(&data);
  Simulation sim(&data, &env);
  
  // Set up data

  Data * d(model.data);
  double * orig_prims;
  orig_prims = new double[d->Ntot * d->Nprims]();

  double v_base = 1e-6;
  double v_max = 2.5e-1;
  double n_base = 1e-8;
  double n_max = 1e-1;
  double p_base = 1e-6;
  double p_max = 1e-1;

  double dn = pow(n_max/n_base,1/nx);
  double dv = pow(v_max/v_base,1/nz);
  double dp = pow(p_max/p_base,1/ny);
  
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        orig_prims[ID(5, i, j, k)] = n_base * pow(dn, i-Ng); // n
        orig_prims[ID(3, i, j, k)] = p_base * pow(dp, j-Ng); // p
        orig_prims[ID(4, i, j, k)] = orig_prims[ID(5, i, j, k)] + orig_prims[ID(3, i, j, k)]/(d->gamma-1); // rho
        orig_prims[ID(0, i, j, k)] = v_base * pow(dv, k-Ng); // v1
        orig_prims[ID(1, i, j, k)] = -2*v_base * pow(dv, k-Ng); // v2
        orig_prims[ID(2, i, j, k)] = -3*v_base * pow(dv, k-Ng); // v3
        orig_prims[ID(9, i, j, k)] = 1e-2 * n_base * pow(dn/1.5, i-Ng); // Pi
        orig_prims[ID(6, i, j, k)] = 1e-2 * v_base * pow(dv/1.2, k-Ng); // q1
        orig_prims[ID(7, i, j, k)] = 2e-2 * v_base * pow(dv/1.2, k-Ng); // q2
        orig_prims[ID(8, i, j, k)] = -3e-2 * v_base * pow(dv/1.2, k-Ng); // q3
        orig_prims[ID(10, i, j, k)] = 1e-2 * n_base * pow(dn/1.4, i-Ng); // pi11
        orig_prims[ID(11, i, j, k)] = 2e-2 * n_base * pow(dn/1.4, i-Ng); // pi12
        orig_prims[ID(12, i, j, k)] = 3e-2 * n_base * pow(dn/1.4, i-Ng); // pi13
        orig_prims[ID(13, i, j, k)] =-1e-2 * n_base * pow(dn/1.4, i-Ng); // pi22
        orig_prims[ID(14, i, j, k)] =-2e-2 * n_base * pow(dn/1.4, i-Ng); // pi23
        orig_prims[ID(15, i, j, k)] =-3e-2 * n_base * pow(dn/1.4, i-Ng); // pi33
        
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int v(0); v < d->Nprims; v++) {
          d->prims[ID(v, i, j, k)] = orig_prims[ID(v, i, j, k)];
        }
      }
    }
  }

  //P2C
  model.primsToAll(d->cons, d->prims, d->aux);

  //Perturb
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int v(0); v < d->Nprims; v++) {
          d->prims[ID(v, i, j, k)] *= 0.99;
        }
      }
    }
  }

  //C2P
  model.getPrimitiveVars(d->cons, d->prims, d->aux);

  //Cross-check
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int v(0); v < d->Nprims; v++) {
          if ( abs(d->prims[ID(v, i, j, k)] - orig_prims[ID(v, i, j, k)]) > 1e-5 * abs(orig_prims[ID(v, i, j, k)]) ) {
//	            || ((d->prims[ID(v, i, j, k)] - orig_prims[ID(v, i, j, k)]) > 1e-5 * orig_prims[ID(v, i, j, k)]) ) {
            cout << "Fails " << v << " " << i << " " << j << " " << k << endl;
            cout << "Origs ";
            for (int vz(0); vz < d->Nprims; vz++) {
              cout << orig_prims[ID(vz, i, j, k)] << " ";
            }
            cout << endl;
            cout << "Prims ";
            for (int vz(0); vz < d->Nprims; vz++) {
              cout << d->prims[ID(vz, i, j, k)] << " ";
            }
            cout << endl;
            cout << "Diffs ";
            for (int vz(0); vz < d->Nprims; vz++) {
              cout << orig_prims[ID(vz, i, j, k)] - d->prims[ID(vz, i, j, k)] << " ";
            }
            cout << endl;
            exit(0);
          }
        }
      }
    }
  }

  return 0;

}
