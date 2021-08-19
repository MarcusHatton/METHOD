#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "ISCE.h"
#include "cminpack.h"
#include "wenoUpwinds.h"

ISCE::ISCE() : Model()
{
  this->Ncons = 5;
  this->Nprims = 16;
  this->Naux = 59;
}

ISCE::ISCE(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 5;
  this->Nprims = (this->data)->Nprims = 16;
  this->Naux = (this->data)->Naux = 59;

  // Solutions for C2P all cells
  solution = (double *) malloc(sizeof(double)*data->Nx*data->Ny*data->Nz);

  smartGuesses = 0;
  
  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("S1");
  this->data->consLabels.push_back("S2");  this->data->consLabels.push_back("S3");
  this->data->consLabels.push_back("Tau");

  // 0
  this->data->primsLabels.push_back("v1");   this->data->primsLabels.push_back("v2");
  this->data->primsLabels.push_back("v3");
  // 3
  this->data->primsLabels.push_back("p");   this->data->primsLabels.push_back("rho");
  this->data->primsLabels.push_back("n");   
  // 6
  this->data->primsLabels.push_back("q1");  this->data->primsLabels.push_back("q2");
  this->data->primsLabels.push_back("q3");
  // 9
  this->data->primsLabels.push_back("Pi");  
  // 10
  this->data->primsLabels.push_back("pi11");   this->data->primsLabels.push_back("pi12");
  this->data->primsLabels.push_back("pi13");  this->data->primsLabels.push_back("pi22");
  this->data->primsLabels.push_back("pi23");  this->data->primsLabels.push_back("pi33");

  // 0
  this->data->auxLabels.push_back("h");     this->data->auxLabels.push_back("T");
  this->data->auxLabels.push_back("e");     this->data->auxLabels.push_back("W");
  // 4
  this->data->auxLabels.push_back("q0");    this->data->auxLabels.push_back("qv");
  this->data->auxLabels.push_back("pi00");  this->data->auxLabels.push_back("pi01");
  this->data->auxLabels.push_back("pi02");  this->data->auxLabels.push_back("pi03");
  this->data->auxLabels.push_back("Theta"); this->data->auxLabels.push_back("vsqrd");
  // 12
  this->data->auxLabels.push_back("q1NS");  this->data->auxLabels.push_back("q2NS");
  this->data->auxLabels.push_back("q3NS");
  // 15
  this->data->auxLabels.push_back("PiNS");    
  // 16
  this->data->auxLabels.push_back("pi11NS"); this->data->auxLabels.push_back("pi12NS");
  this->data->auxLabels.push_back("pi13NS"); this->data->auxLabels.push_back("pi22NS");
  this->data->auxLabels.push_back("pi23NS"); this->data->auxLabels.push_back("pi33NS");
  // 22
  this->data->auxLabels.push_back("q1LO");  this->data->auxLabels.push_back("q2LO");
  this->data->auxLabels.push_back("qLO");
  // 25
  this->data->auxLabels.push_back("PiLO");    
  // 26
  this->data->auxLabels.push_back("pi11LO"); this->data->auxLabels.push_back("pi12LO");
  this->data->auxLabels.push_back("pi13LO"); this->data->auxLabels.push_back("pi22LO");
  this->data->auxLabels.push_back("pi23LO"); this->data->auxLabels.push_back("pi33LO");
  // 32
  this->data->auxLabels.push_back("a1");     this->data->auxLabels.push_back("a2");   
  this->data->auxLabels.push_back("a3");

  // 35
  this->data->auxLabels.push_back("dtp");  this->data->auxLabels.push_back("dtrho");
  this->data->auxLabels.push_back("dtn");
  // 38
  this->data->auxLabels.push_back("dtv1");
  this->data->auxLabels.push_back("dtv2");  this->data->auxLabels.push_back("dtv3");
  // 41
  this->data->auxLabels.push_back("dtW");   this->data->auxLabels.push_back("dtT"); 
  // 43
  this->data->auxLabels.push_back("dtq1NS");  this->data->auxLabels.push_back("dtq2NS");
  this->data->auxLabels.push_back("dtq3NS");
  // 46
  this->data->auxLabels.push_back("dtPiNS");    
  // 47
  this->data->auxLabels.push_back("dtpi11NS"); this->data->auxLabels.push_back("dtpi12NS");
  this->data->auxLabels.push_back("dtpi13NS"); this->data->auxLabels.push_back("dtpi22NS");
  this->data->auxLabels.push_back("dtpi23NS"); this->data->auxLabels.push_back("dtpi33NS");
  // 53
  this->data->auxLabels.push_back("dtD"); this->data->auxLabels.push_back("dtS1");
  this->data->auxLabels.push_back("dtS2"); this->data->auxLabels.push_back("dtS3");
  this->data->auxLabels.push_back("dtTau");  this->data->auxLabels.push_back("dtE"); 
 

}

ISCE::~ISCE()
{
  free(solution);
}

void ISCE::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  // printf("ToyQ model does not implement sourceTermSingleCell\n");
  // exit(1);

  Data * d(this->data);

  // D
  source[0] = 0.0;
  // S1,2,3
  source[1] = 0.0; 
  source[2] = 0.0;
  source[3] = 0.0; 
  // Tau
  source[4] = 0.0;
  
}

void ISCE::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int k(0); k < this->data->Nz; k++) {

        // D
        source[ID(D, i, j, k)] = 0.0;
        // S1,2,3
        source[ID(S1, i, j, k)] = 0.0;
        source[ID(S2, i, j, k)] = 0.0;
        source[ID(S3, i, j, k)] = 0.0;
        // Tau
        source[ID(Tau, i, j, k)] = 0.0;        
      }
    }
  }

}

//! Residual function to minimize in the format required by cminpack
/*!
    I know, its a horrible layout, alas we're stuck with it.

    Parameters
    ----------
    p : pointer to struct
      Struct contains additional arguments that are required (if any)
    n : int
      Size of system
    x : pointer to double
      The array containing the initial guess
    fvec : pointer to double
      The array containing the solution
    iflag : int
      Error flag
*/
int ISCEresidual(void *ptr, int n, const double *x, double *fvec, int iflag)
{

//  Data * d(this->data);
  
  // Retrieve additional arguments
  ISCE::Args * args = (ISCE::Args*) ptr;

  // Values must make sense
  // Think this should be taken out for now - need new sensible values
  /*
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  */
  
  double E_rf = args->Tau_rf + args->D_rf;
  double vsqrd_rf = (args->S1_rf*args->S1_rf + args->S2_rf*args->S2_rf + args->S3_rf*args->S3_rf)/((E_rf + x[0])*(E_rf + x[0]));
  double W_rf(1 / sqrt(1 - vsqrd_rf));
  double n_rf(args->D_rf / W_rf);
  double rho_plus_p_rf = ((E_rf + x[0])/(W_rf*W_rf));
  double v1_rf = args->S1_rf/(rho_plus_p_rf*W_rf*W_rf);
  double v2_rf = args->S2_rf/(rho_plus_p_rf*W_rf*W_rf);
  double v3_rf = args->S3_rf/(rho_plus_p_rf*W_rf*W_rf);
  double p_rf = (rho_plus_p_rf - n_rf)*((args->gamma-1)/args->gamma);
  double rho_rf = rho_plus_p_rf - p_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 1 || n_rf < 0 || abs(v1_rf) >= 1 || abs(v2_rf) >= 1 || abs(v3_rf) >= 1 || vsqrd_rf >= 1) {
    printf("EEK");
    fvec[0] = 1e6;
    return 0;
  }

  fvec[0] = p_rf - x[0];

  return 0;
}

void ISCE::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                      // Additional arguments structure
  const int sys_size(1);                     // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1e-5;          // Tolerance of rootfinder
  const int lwa = 8;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Set additional args for rootfind
  args.D_rf = cons[Cons::D];
  args.S1_rf = cons[Cons::S1];
  args.S2_rf = cons[Cons::S2];
  args.S3_rf = cons[Cons::S3];
  args.Tau_rf = cons[Cons::Tau];
  args.gamma = d->gamma;
  
  sol[0] = prims[Prims::p]; // Guess the pressure

  // Solve residual = 0
  info = __cminpack_func__(hybrd1) (&ISCEresidual, &args, sys_size, sol, res,
                                    tol, wa, lwa);
  // If root find fails, add failed cell to the list
  if (info!=1) {
    //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
    throw std::runtime_error("C2P could not converge.\n");
  }
  aux[Aux::vsqrd] = (cons[Cons::S1]*cons[Cons::S1] + cons[Cons::S2]*cons[Cons::S2] 
                      + cons[Cons::S3]*cons[Cons::S3] - sol[3])
                      /((cons[Cons::Tau] + cons[Cons::D] + sol[0])*(cons[Cons::Tau]  + cons[Cons::D] + sol[0]));
  aux[Aux::W] = 1 / sqrt((1-aux[Aux::vsqrd]));
  prims[Prims::n] = cons[Cons::D] / aux[Aux::W];
  double rho_plus_p = (cons[Cons::Tau] + cons[Cons::D] + sol[0])/(aux[Aux::W]*aux[Aux::W]);
  prims[Prims::v1] = cons[Cons::S1]/(rho_plus_p*aux[Aux::W]*aux[Aux::W]);
  prims[Prims::v2] = cons[Cons::S2]/(rho_plus_p*aux[Aux::W]*aux[Aux::W]);
  prims[Prims::v3] = cons[Cons::S3]/(rho_plus_p*aux[Aux::W]*aux[Aux::W]);
  prims[Prims::p] = (rho_plus_p - prims[Prims::n])*((d->gamma-1)/d->gamma);
  prims[Prims::rho] = rho_plus_p - prims[Prims::p];
  
  // Repeating the ones here that depend on v1,v2,v3...
  aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
  aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
        
  aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
  aux[Aux::T] = prims[Prims::p] / prims[Prims::n];     
   
}

void ISCE::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int sys_size(1);              // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1e-5;          // Tolerance of rootfinder
  const int lwa = 8;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array
  std::vector<Failed> fails;          // Vector of failed structs. Stores location of failed cons2prims cells.
  
  /* - fix n->sys_size if ever used
  int maxfev(50);
  int ml(n);
  int mu(n);
  double epsfcn(0.1);
  double diag[n];
  int mode(1);
  double factor(100);
  int nprint(-1);
  int nfev(0);
  double fjac[n][n];
  int ldfjac(n);
  int lr(10);
  double r[lr];
  double qtf[n];
  double wa1[n];                     // Work array
  double wa2[n];                     // Work array
  double wa3[n];                     // Work array
  double wa4[n];                     // Work array
  */
  
  // Y1-3,U,Z11-33
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
            
        // Set additional args for rootfind
        args.D_rf = cons[ID(Cons::D, i, j, k)];
        args.S1_rf = cons[ID(Cons::S1, i, j, k)];
        args.S2_rf = cons[ID(Cons::S2, i, j, k)];
        args.S3_rf = cons[ID(Cons::S3, i, j, k)];
        args.Tau_rf = cons[ID(Cons::Tau, i, j, k)];
        args.gamma = d->gamma;
        
        // Guess the pressure
        sol[0] = prims[ID(Prims::p, i, j, k)];

        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&ISCEresidual, &args, sys_size, sol, res,
                                          tol, wa, lwa);

        // Another solver not in use currently  
        // info = __cminpack_func__(hybrd) (&ISCEresidual, &args, sys_size, sol, res,
        //                                  tol, maxfev, ml, mu, epsfcn, &diag[0], mode, factor, nprint, &nfev, &fjac[0][0], ldfjac, &r[0], lr, &qtf[0], &wa1[0], &wa2[0], &wa3[0], &wa4[0]);                                  

        // If root find fails, add failed cell to the list
        if (info!=1) {
          printf("%i info\n",info);
          printf("(%i, %i, %i) failed\n", i, j, k);
          printf("(%g, %g, %g, %g) res\n", res[0], res[1], res[2], res[3]);
          printf("(%g, %g, %g, %g) sol\n", sol[0], sol[1], sol[2], sol[3]);
          std::cout << "Prims ";
          for (int prim_count(0); prim_count < d->Nprims; prim_count++) {
            std::cout << d->prims[ID(prim_count, i, j, k)] << " ";
          }
          std::cout << std::endl;
          std::cout << "Cons  ";
          for (int con_count(0); con_count < d->Ncons; con_count++) {
            std::cout << d->cons[ID(con_count, i, j, k)] << " ";
          }
          std::cout << std::endl;
          std::cout << "Aux   ";
          for (int aux_count(0); aux_count < d->Naux; aux_count++) {
            std::cout << d->aux[ID(aux_count, i, j, k)] << " ";
          }
          std::cout << std::endl;
          //printf("(%f, %f, %f, %f, %f) prims\n",  prims[ID(Prims::p, i, j, k)], prims[ID(Prims::Pi, i, j, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::v1, i, j, k)], prims[ID(Prims::q1, i, j, k)]);
          //printf("(%f, %f, %f, %f) aux  \n",  aux[ID(Aux::W, i, j, k)], aux[ID(Aux::qv, i, j, k)], aux[ID(Aux::pi00, i, j, k)], aux[ID(Aux::pi01, i, j, k)]);
          //printf("(%f, %f, %f, %f) cons \n",  cons[ID(Cons::Y1, i, j, k)], cons[ID(Cons::D, i, j, k)], cons[ID(Cons::S1, i, j, k)], cons[ID(Cons::Tau, i, j, k)]);
          exit(0);
          Failed fail = {i, j, k};
          fails.push_back(fail);
        }
        else {
          // Now have the correct values for Chi, Sigma123
          solution[ID(0, i, j, k)] = sol[0];
        }      
      
      } // End k-loop
    } // End j-loop
  } // End i-loop
  
  /*
  
  // ################################## Smart guessing ########################### //
  // Are there any failures?
  if (fails.size() > 0) {
    int x, y, z;
    // Loop through any failed cells and try again, using the mean of successfull
    // surrounding cells solutions as an initial estimate
    for (Failed fail : fails) {
      x = fail.x;
      y = fail.y;
      z = fail.z;
      // Vector to contain successful neighbours
      std::vector<Failed> neighbours;
      if (x > 0) neighbours.push_back(Failed {x-1, y, z});
      if (y > 0) neighbours.push_back(Failed {x, y-1, z});
      if (z > 0) neighbours.push_back(Failed {x, y, z-1});
      if (x < d->Nx - 1) neighbours.push_back(Failed {x+1, y, z});
      if (y < d->Ny - 1) neighbours.push_back(Failed {x, y+1, z});
      if (z < d->Nz - 1) neighbours.push_back(Failed {x, y, z+1});

      sol[0] = 0;
      sol[1] = 0;
      sol[2] = 0;
      sol[3] = 0;
      for (Failed neighbour : neighbours) {
        sol[0] += solution[ID(0, neighbour.x, neighbour.y, neighbour.z)];
        sol[1] += solution[ID(1, neighbour.x, neighbour.y, neighbour.z)];
        sol[2] += solution[ID(2, neighbour.x, neighbour.y, neighbour.z)];
        sol[3] += solution[ID(3, neighbour.x, neighbour.y, neighbour.z)];
      }
      sol[0] /= neighbours.size();
      sol[1] /= neighbours.size();
      sol[2] /= neighbours.size();
      sol[3] /= neighbours.size();
      // Solve residual = 0
      info = __cminpack_func__(hybrd1) (&ISresidual, &args, n, sol, res,
                                        tol, wa, lwa);
      if (info != 1) {
        printf("Smart guessing did not work, exiting\n");
        printf("(%i) \n",info);
        printf("(%d, %d, %d) failed\n", fail.x, fail.y, fail.z);
        std::exit(1);
      } else {
        smartGuesses++;
        printf("Smart guessing worked!\n");
        solution[ID(0, x, y, z)] = sol[0];
        solution[ID(1, x, y, z)] = sol[1];
        solution[ID(2, x, y, z)] = sol[2];
        solution[ID(3, x, y, z)] = sol[3];
      } // else
    } // for

  } // if

  */

    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {

        // C2P Scheme as outlined in HP/FYR
        aux[ID(Aux::vsqrd, i, j, k)] = (cons[ID(Cons::S1, i, j, k)]*cons[ID(Cons::S1, i, j, k)] 
                                  + cons[ID(Cons::S2, i, j, k)]*cons[ID(Cons::S2, i, j, k)]
                                  + cons[ID(Cons::S3, i, j, k)]*cons[ID(Cons::S3, i, j, k)])
                                  /((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])*(cons[ID(Cons::Tau, i, j, k)] 
                                  + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)]));
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt((1-aux[ID(Aux::vsqrd, i, j, k)]));
        prims[ID(Prims::n, i, j, k)] = cons[ID(Cons::D, i, j, k)] / aux[ID(Aux::W, i, j, k)];
        double rho_plus_p = ((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])
                                            /(aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]));
        prims[ID(Prims::v1, i, j, k)] = cons[ID(Cons::S1, i, j, k)]/(rho_plus_p*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);
        prims[ID(Prims::v2, i, j, k)] = cons[ID(Cons::S2, i, j, k)]/(rho_plus_p*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);
        prims[ID(Prims::v3, i, j, k)] = cons[ID(Cons::S3, i, j, k)]/(rho_plus_p*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);
        prims[ID(Prims::p, i, j, k)] = (rho_plus_p - prims[ID(Prims::n, i, j, k)])*((d->gamma-1)/d->gamma);
        prims[ID(Prims::rho, i, j, k)] = rho_plus_p - prims[ID(Prims::p, i, j, k)];
       
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];

      } // End k-loop
    } // End j-loop
  } // End i-loop  
        
}

void ISCE::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // printf("Calling primsToAll\n");

  // W, q_kv^k, pi^0_0
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(Aux::vsqrd, i, j, k)] = prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - aux[ID(Aux::vsqrd, i, j, k)] );
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        prims[ID(Prims::rho, i, j, k)] = prims[ID(Prims::n, i, j, k)]*(1+aux[ID(Aux::e, i, j, k)]);
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
        aux[ID(Aux::h, i, j, k)] = 1 + aux[ID(Aux::e, i, j, k)] + prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
      }
    }
  }

  // Conserveds are now EulerSR form (no dissipation)
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(Cons::D, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)];
        // S1,2,3
        cons[ID(Cons::S1, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)]; 
        cons[ID(Cons::S2, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)]; 
        cons[ID(Cons::S3, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)];
        // Tau
        cons[ID(Cons::Tau, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] 
        - (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]);
      }  
    }
  }

/*
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(Cons::D, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)];
        
        // S1,2,3
        cons[ID(Cons::S1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2 + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)]
        + tau_q*(aux[ID(Aux::q1LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2 + aux[ID(Aux::q1LO, i, j, k)] + aux[ID(Aux::q2LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)]
        + tau_Pi*aux[ID(Aux::PiLO, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]
        + tau_pi*(aux[ID(Aux::pi11LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi12LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13LO, i, j, k)]*prims[ID(Prims::v3, i, j, k)]);

        cons[ID(Cons::S2, i, j, k)] = ux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2 + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi21NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)]
        + tau_q*(aux[ID(Aux::q1LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q2LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2 + aux[ID(Aux::q2LO, i, j, k)] + aux[ID(Aux::q3LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)]
        + tau_Pi*aux[ID(Aux::PiLO, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]
        + tau_pi*(aux[ID(Aux::pi21LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi22LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23LO, i, j, k)]*prims[ID(Prims::v3, i, j, k)]);

        cons[ID(Cons::S3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2 + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)] + aux[ID(Aux::pi31NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi32NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)]
        + tau_q*(aux[ID(Aux::q1LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::q2LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::q3LO, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2 + aux[ID(Aux::q3LO, i, j, k)])*aux[ID(Aux::W, i, j, k)]
        + tau_Pi*aux[ID(Aux::PiLO, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]
        + tau_pi*(aux[ID(Aux::pi31LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi32LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi33LO, i, j, k)]*prims[ID(Prims::v3, i, j, k)]);

        // E
        cons[ID(Cons::Tau, i, j, k)] = aaux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)]) - aux[ID(Aux::PiNS, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] - prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]
        + 2*tau_q*(aux[ID(Aux::q1LO, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2LO, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3LO, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)]
        + tau_Pi*(sqr(aux[ID(Aux::W, i, j, k)]) - 1)*aux[ID(Aux::PiLO, i, j, k)]
        + tau_pi*(aux[ID(Aux::pi11LO, i, j, k)] + aux[ID(Aux::pi22LO, i, j, k)] + aux[ID(Aux::pi33LO, i, j, k)]);

        // Y1-3
        cons[ID(Cons::Y1, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::q1, i, j, k)];
        cons[ID(Cons::Y2, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::q2, i, j, k)];
        cons[ID(Cons::Y3, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::q3, i, j, k)];
        //printf("%f, %f, %f, q123\n", prims[ID(Prims::q1, i, j, k)], cons[ID(Cons::Y1, i, j, k)]);
        // U
        cons[ID(Cons::U, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::Pi, i, j, k)];
        // Z11-33
        cons[ID(Cons::Z11, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::pi11, i, j, k)];
        cons[ID(Cons::Z12, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::pi12, i, j, k)]; 
        cons[ID(Cons::Z13, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::pi13, i, j, k)];
        cons[ID(Cons::Z22, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::pi22, i, j, k)];        
        cons[ID(Cons::Z23, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::pi23, i, j, k)];
        cons[ID(Cons::Z33, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::pi33, i, j, k)];
      }  
    }
  }
*/

}

void ISCE::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Dv
        f[ID(0, i, j, k)] = cons[ID(Cons::D, i, j, k)]*prims[ID(dir, i, j, k)];
        // Sv + ..
        for (int nvar(0); nvar < 3; nvar++) {
          f[ID(1+nvar, i, j, k)] = cons[ID(Cons::S1+nvar, i, j, k)]*prims[ID(dir, i, j, k)]; // + ( prims[ID(Prims::q1+dir, i, j, k)] * prims[ID(Prims::v1+nvar, i, j, k)]  
            // - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+nvar, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)] ) * aux[ID(Aux::W, i, j, k)];
          // (p+Pi)delta_ij
          if (dir == nvar) {
            f[ID(1+nvar, i, j, k)] += (prims[ID(Prims::p, i, j, k)]); // + prims[ID(Prims::Pi, i, j, k)]);
          }
        }
        /*
        //  pi^i_j  
        if (dir == 0) {
          for (int nvar(0); nvar < 3; nvar++) {
            f[ID(1+nvar, i, j, k)] += prims[ID(Prims::pi11+nvar, i, j, k)];
          }
        } else if (dir == 1) {
          f[ID(1, i, j, k)] += prims[ID(Prims::pi12, i, j, k)];
          f[ID(2, i, j, k)] += prims[ID(Prims::pi22, i, j, k)];
          f[ID(3, i, j, k)] += prims[ID(Prims::pi23, i, j, k)];
        } else if (dir == 2) {
          f[ID(1, i, j, k)] += prims[ID(Prims::pi13, i, j, k)];
          f[ID(2, i, j, k)] += prims[ID(Prims::pi23, i, j, k)];
          f[ID(3, i, j, k)] += prims[ID(Prims::pi33, i, j, k)];
        } else {
          throw std::runtime_error("Flux direction is not 0, 1 or 2");
        }
        */

        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Cons::Tau, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * prims[ID(dir, i, j, k)];
          // + (prims[ID(Prims::q1+dir, i, j, k)] - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)])*aux[ID(Aux::W, i, j, k)]
          // + aux[ID(Aux::pi01+dir, i, j, k)];
      } // End k loop
    } // End j loop
  } // End i loop
}

