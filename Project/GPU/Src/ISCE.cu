//! Special relativistic dissipative hydro model - ISCE
/*!
    This script contains the function definitions for the ISCE model. The form
  of the quations has been taken from Dionysopoulou and we use a divergence cleaning method
  taken from Muddle.
    For detailed documentation about the methods contained herein, see ISCE.h
  and model.h.
*/

#include <stdexcept>
#include <cstdio>
#include "ISCE.h"
#include "cminpack.h"
#include "cudaErrorCheck.h"

#define TOL 1.0e-12
#define EPS 1.0e-4
#define MAXITER 50

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDCons(var, idx, jdx, kdx) ( (var) + (idx)*(d->Ncons)*(d->Nz)*(d->Ny) + (jdx)*(d->Ncons)*(d->Nz) + (kdx)*(d->Ncons)  )
#define IDPrims(var, idx, jdx, kdx) ( (var) + (idx)*(d->Nprims)*(d->Nz)*(d->Ny) + (jdx)*(d->Nprims)*(d->Nz) + (kdx)*(d->Nprims)  )
#define IDAux(var, idx, jdx, kdx) ( (var) + (idx)*(d->Naux)*(d->Nz)*(d->Ny) + (jdx)*(d->Naux)*(d->Nz) + (kdx)*(d->Naux)  )

// C2P residual and rootfinder (Serial)
static double residual(const double, const double, const double, const double, double);
static int newton(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma, int i, int j, int k);

// C2P residual and rootfinder (Parallel)
__device__
static double residualParallel(const double Z, const double StildeSq, const double D, const double tauTilde, double gamma);
__device__
static int newtonParallel(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma);
__global__
static void getPrimitiveVarsParallel(double *cons, double *prims, double *aux, double *guess, int stream, double gamma, double sigma, int Ncons, int Nprims, int Naux, int origWidth, int streamWidth);

ISCE::ISCE() : Model()
{
  modType_t = ModelType::ISCE;
  this->Ncons = 5;
  this->Nprims = 16;
  this->Naux = 59;

  cudaHostAlloc((void **)&singleCons, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singlePrims, sizeof(double) * this->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleAux, sizeof(double) * this->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleSource, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);
}

ISCE::~ISCE()
{
  // Free up
  cudaFreeHost(singleCons);
  cudaFreeHost(singlePrims);
  cudaFreeHost(singleAux);
  cudaFreeHost(singleSource);

  delete c2pArgs;
}

ISCE::ISCE(Data * data) : Model(data)
{
  modType_t = ModelType::ISCE;
  this->Ncons = (this->data)->Ncons = 5;
  this->Nprims = (this->data)->Nprims = 16;
  this->Naux = (this->data)->Naux = 59;

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

  // Single cell work arrays
  cudaHostAlloc((void **)&singleCons, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singlePrims, sizeof(double) * this->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleAux, sizeof(double) * this->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleSource, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);

  c2pArgs = new C2PArgs(this->data);
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

void ISCE::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  // Syntax
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

void ISCE::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // Syntax
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
  
  aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
  aux[Aux::T] = prims[Prims::p] / prims[Prims::n];     
  aux[Aux::h] = 1 + aux[Aux::e] + prims[Prims::p] / prims[Prims::n];

}

void ISCE::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // First need to copy data to the device
  // A single cell requires all cons variables and aux10 to start the guessing
  // Rearrange data into host arrays ready for copying
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++) {
          c2pArgs->cons_h[IDCons(var, i, j, k)] = cons[ID(var, i, j, k)];
        }
        c2pArgs->guess_h[ID(0, i, j, k)] = aux[ID(10, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < c2pArgs->Nstreams; i++) {
    // Which cell is at the left bound?
    int lcell(i * c2pArgs->streamWidth);
    // Which cell is at the right bound?
    int rcell(lcell + c2pArgs->streamWidth);
    if (rcell > d->Ncells) rcell = d->Ncells;
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));

    // Send stream's data
    gpuErrchk( cudaMemcpyAsync(c2pArgs->cons_d[i], c2pArgs->cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, c2pArgs->stream[i]) );
    gpuErrchk( cudaMemcpyAsync(c2pArgs->guess_d[i], c2pArgs->guess_h + lcell, inMemsize, cudaMemcpyHostToDevice, c2pArgs->stream[i]) );

    // Call kernel and operate on data
    getPrimitiveVarsParallel <<< c2pArgs->bpg, c2pArgs->tpb,
        c2pArgs->tpb * c2pArgs->cellMem, c2pArgs->stream[i] >>> (c2pArgs->cons_d[i],
        c2pArgs->prims_d[i], c2pArgs->aux_d[i], c2pArgs->guess_d[i], i, d->gamma, d->sigma, d->Ncons,
        d->Nprims, d->Naux, c2pArgs->streamWidth, width);


    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(c2pArgs->prims_h + lcell*d->Nprims, c2pArgs->prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, c2pArgs->stream[i]) );
    gpuErrchk( cudaMemcpyAsync(c2pArgs->aux_h + lcell*d->Naux, c2pArgs->aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, c2pArgs->stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );

  // Rearrange data back into arrays
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        for (int var(0); var < d->Nprims; var++) {
          prims[ID(var, i, j, k)] = c2pArgs->prims_h[IDPrims(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          aux[ID(var, i, j, k)] = c2pArgs->aux_h[IDAux(var, i, j, k)];
        }
      }
    }
  }
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
}

static double residual(const double Z, const double StildeSq, const double D, const double tauTilde, double gamma)
{
  // Decalre variables
  double vsq, W, rho, h, p, resid;

  vsq = StildeSq / (Z * Z);

  // Sanity check
  if (vsq >= 1.0 || Z < 0) return 1.0e6;

  // Continue
  W = 1 / sqrt(1 - vsq);
  rho = D / W;
  h = Z / (rho * W * W);
  p = (gamma - 1) * (h - rho) / gamma;

  // Second sanity check
  if (rho < 0 || p < 0 || W < 1 || h < 0) return 1.0e6;

  // Values are physical, compute residual
  resid = (1 - (gamma - 1) / (W * W * gamma)) * Z + ((gamma - 1) /
          (W * gamma) - 1) * D - tauTilde;

  return resid;

}


static int newton(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma, int i, int j, int k)
{
  // Rootfind data
  double bestX;
  double x0(*Z);
  double eps(EPS);
  double x1(x0 + eps);
  double tol(TOL);
  double x2;
  double bestF;
  double f0(residual(x0, StildeSq, D, tauTilde, gamma));
  double f1(residual(x1, StildeSq, D, tauTilde, gamma));
  int iter;
  int maxiter(MAXITER);
  int found(0);

  // If root can not be found return the best so far
  bestX = x0; bestF = f0;
  for (iter=0; iter<maxiter; iter++) {
    if (fabs(f0) < tol) {
      *Z = x0;
      found = 1;
      break;
    }

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x1 = x0;
    x0 = x2;
    f1 = f0;
    f0 = residual(x0, StildeSq, D, tauTilde, gamma);
    if (fabs(f0) < fabs(bestF)) {
      bestX = x0;
      bestF = f0;
    }
  }

  if (!found) {
    // Store result of Z=rho*h*W**2
    *Z = bestX;
    char s[200];
    sprintf(s, "C2P could not converge in cell (%d, %d, %d)\n", i, j, k);

    throw std::runtime_error(s);
  }
  return 1;
}

// /*!
//     This is the device version of the getPrimitiveVars that takes a streams data
//     and computes the rest of the prims and aux vars. This is called when
//     ISCE::getPrimitiveVars is required, i.e. all cells need to be found.
// */
__global__
static void getPrimitiveVarsParallel(double *streamCons, double *streamPrims, double *streamAux, double *guess, int stream, double gamma, double sigma, int Ncons, int Nprims, int Naux, int origWidth, int streamWidth)
{
  // First need thread indicies
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  // const int gID(lID + stream * origWidth);        //!< global index (in domain)
  // Allocate shared memory
  extern __shared__ double sharedArray [];
  double * cons = &sharedArray[tID * (Ncons + Nprims + Naux)];
  double * prims = &cons[Ncons];
  double * aux = &prims[Nprims];

  if (lID < streamWidth) {

    // Load conserved vector into shared memory, and the initial guess
    for (int i(0); i < Ncons; i++) cons[i] = streamCons[lID * Ncons + i];
    aux[10] = guess[lID];

    // Set Bx/y/z and Ex/y/z field in prims
    prims[5] = cons[5]; prims[6] = cons[6]; prims[7] = cons[7];
    prims[8] = cons[8]; prims[9] = cons[9]; prims[10] = cons[10];

    // Bsq, Esq
    aux[7] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
    aux[8] = cons[8] * cons[8] + cons[9] * cons[9] + cons[10] * cons[10];

    // Sbarx, Sbary, Sbarz
    aux[12] = cons[1] - (cons[9] * cons[7] - cons[10] * cons[6]);
    aux[13] = cons[2] - (cons[10] * cons[5] - cons[8] * cons[7]);
    aux[14] = cons[3] - (cons[8] * cons[6] - cons[9] * cons[5]);
    // Sbarsq, tauBar
    aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
    aux[16] = cons[4] - 0.5 * (aux[7] + aux[8]);

    // Solve
    newtonParallel(&aux[10], aux[15], cons[0], aux[16], gamma);

    // vsq
    aux[9] = aux[15] / (aux[10] * aux[10]);

    // W
    aux[1] = 1.0 / sqrt(1 - aux[9]);
    // rho
    prims[0] = cons[0] / aux[1];
    // h
    aux[0] = aux[10] / (prims[0] * aux[1] * aux[1]);
    // e
    aux[2] = (aux[0] - 1) / gamma;
    // c
    aux[3] = sqrt((aux[2] * gamma * (gamma - 1)) / aux[0]);
    // p
    prims[4] = prims[0] * aux[2] * (gamma - 1);
    // vx, vy, vz
    prims[1] = aux[12] / aux[10];
    prims[2] = aux[13] / aux[10];
    prims[3] = aux[14] / aux[10];
    // vE
    aux[11] = prims[1] * cons[8] + prims[2] * cons[9] + prims[3] * cons[10];
    // Jx, Jy, Jz
    aux[4] = cons[13] * prims[1] + aux[1] * sigma * (cons[8] + (prims[2] * cons[7] -
             prims[3] * cons[6]) - aux[11] * prims[1]);
    aux[5] = cons[13] * prims[2] + aux[1] * sigma * (cons[9] + (prims[3] * cons[5] -
             prims[1] * cons[7]) - aux[11] * prims[2]);
    aux[6] = cons[13] * prims[3] + aux[1] * sigma * (cons[10] + (prims[1] * cons[6] -
             prims[2] * cons[5]) - aux[11] * prims[3]);

  }

  // Copy data back from shared memory into device arrays
  for (int i(0); i < Nprims; i++) streamPrims[lID * Nprims + i] = prims[i];
  for (int i(0); i < Naux; i++) streamAux[lID * Naux + i] = aux[i];

}

__device__
static int newtonParallel(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma)
{
  // Rootfind data
  double x0(*Z);
  double x1(x0 + EPS);
  double x2;
  double f0(residualParallel(x0, StildeSq, D, tauTilde, gamma));
  double f1(residualParallel(x1, StildeSq, D, tauTilde, gamma));
  int iter;

  for (iter=0; iter<MAXITER; iter++) {
    if (fabs(f0) < TOL) {
      *Z = x0;
      return 1;
    }

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x1 = x0;
    x0 = x2;
    f1 = f0;
    f0 = residualParallel(x0, StildeSq, D, tauTilde, gamma);
  }

  return 0;
}

__device__
static double residualParallel(const double Z, const double StildeSq, const double D, const double tauTilde, double gamma)
{
  // Declare variables
  double  W;

  // Sanity check
  if (Z < 0) return 1.0e6;

  // Continue
  W = 1 / sqrt(1 - (StildeSq / (Z * Z)));

  // Values are physical, compute residual
  return (1 - (gamma - 1) / (W * W * gamma)) * Z + ((gamma - 1) /
          (W * gamma) - 1) * D - tauTilde;

}


__device__
void ISCE_D::getPrimitiveVarsSingleCell(double * cons, double * prims, double * aux)
{

    // Set Bx/y/z and Ex/y/z field in prims
    prims[5] = cons[5]; prims[6] = cons[6]; prims[7] = cons[7];
    prims[8] = cons[8]; prims[9] = cons[9]; prims[10] = cons[10];

    // Bsq, Esq
    aux[7] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
    aux[8] = cons[8] * cons[8] + cons[9] * cons[9] + cons[10] * cons[10];

    // Sbarx, Sbary, Sbarz
    aux[12] = cons[1] - (cons[9] * cons[7] - cons[10] * cons[6]);
    aux[13] = cons[2] - (cons[10] * cons[5] - cons[8] * cons[7]);
    aux[14] = cons[3] - (cons[8] * cons[6] - cons[9] * cons[5]);
    // Sbarsq, tauBar
    aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
    aux[16] = cons[4] - 0.5 * (aux[7] + aux[8]);


    // Solve
    newtonParallel(&aux[10], aux[15], cons[0], aux[16], args->gamma);

    // vsq
    aux[9] = aux[15] / (aux[10] * aux[10]);
    // W
    aux[1] = 1.0 / sqrt(1 - aux[9]);
    // rho
    prims[0] = cons[0] / aux[1];
    // h
    aux[0] = aux[10] / (prims[0] * aux[1] * aux[1]);
    // e
    aux[2] = (aux[0] - 1) / args->gamma;
    // c
    aux[3] = sqrt((aux[2] * args->gamma * (args->gamma - 1)) / aux[0]);
    // p
    prims[4] = prims[0] * aux[2] * (args->gamma - 1);
    // vx, vy, vz
    prims[1] = aux[12] / aux[10];
    prims[2] = aux[13] / aux[10];
    prims[3] = aux[14] / aux[10];
    // vE
    aux[11] = prims[1] * cons[8] + prims[2] * cons[9] + prims[3] * cons[10];
    // Jx, Jy, Jz
    aux[4] = cons[13] * prims[1] + aux[1] * args->sigma * (cons[8] + (prims[2] * cons[7] -
             prims[3] * cons[6]) - aux[11] * prims[1]);
    aux[5] = cons[13] * prims[2] + aux[1] * args->sigma * (cons[9] + (prims[3] * cons[5] -
             prims[1] * cons[7]) - aux[11] * prims[2]);
    aux[6] = cons[13] * prims[3] + aux[1] * args->sigma * (cons[10] + (prims[1] * cons[6] -
             prims[2] * cons[5]) - aux[11] * prims[3]);
}

__device__
void ISCE_D::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source)
{
  // D
  source[0] = 0.0;
  // S1,2,3
  source[1] = 0.0; 
  source[2] = 0.0;
  source[3] = 0.0; 
  // Tau
  source[4] = 0.0;
}
