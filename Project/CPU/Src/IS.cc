#include "IS.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

// enums to save looking up numbering of C/P/As when using ID accessor.

enum Cons { D, S1, S2, S3, Tau, Y1, Y2, Y3, U, Z11, Z12, Z13, Z22, Z23, Z33 };
enum Prims { v1, v2, v3, p, rho, n, q1, q2, q3, Pi, pi11, pi12, pi13, pi22, pi23, pi33 };
enum Aux { h, T, e, W, q0, qv, pi00, pi01, pi02, pi03, q1NS, q2NS, q3NS, PiNS, 
           pi11NS, pi12NS, pi13NS, pi22NS, pi23NS, pi33NS, Theta, dv1dt, 
           dv2dt, dv3dt, a1, a2, a3, vsqrd, dWdt, rho_plus_p };


IS::IS() : Model()
{
  this->Ncons = 15;
  this->Nprims = 16;
  this->Naux = 30;
}

IS::IS(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 15;
  this->Nprims = (this->data)->Nprims = 16;
  this->Naux = (this->data)->Naux = 30;

  // Solutions for C2P all cells
  solution = (double *) malloc(sizeof(double)*4*data->Nx*data->Ny*data->Nz);

  smartGuesses = 0;
  
  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("S1");
  this->data->consLabels.push_back("S2");  this->data->consLabels.push_back("S3");
  this->data->consLabels.push_back("Tau");  this->data->consLabels.push_back("Y1");
  this->data->consLabels.push_back("Y2");  this->data->consLabels.push_back("Y3");
  this->data->consLabels.push_back("U");  this->data->consLabels.push_back("Z11");
  this->data->consLabels.push_back("Z12");  this->data->consLabels.push_back("Z13");
  this->data->consLabels.push_back("Z22");  this->data->consLabels.push_back("Z23");
  this->data->consLabels.push_back("Z33");

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
  // 10
  this->data->auxLabels.push_back("q1NS");  this->data->auxLabels.push_back("q2NS");
  this->data->auxLabels.push_back("q3NS");
  // 13
  this->data->auxLabels.push_back("PiNS");    
  // 14
  this->data->auxLabels.push_back("pi11NS"); this->data->auxLabels.push_back("pi12NS");
  this->data->auxLabels.push_back("pi13NS"); this->data->auxLabels.push_back("pi22NS");
  this->data->auxLabels.push_back("pi23NS"); this->data->auxLabels.push_back("pi33NS");
  // 20
  this->data->auxLabels.push_back("Theta");  this->data->auxLabels.push_back("dv1dt");
  this->data->auxLabels.push_back("dv2dt");  this->data->auxLabels.push_back("dv3dt");
  this->data->auxLabels.push_back("a1");     this->data->auxLabels.push_back("a2");   
  this->data->auxLabels.push_back("a3");     this->data->auxLabels.push_back("vsqrd");
  this->data->auxLabels.push_back("dWdt");   this->data->auxLabels.push_back("rho_plus_p");
}

IS::~IS()
{
  free(solution);
}


void IS::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  // printf("ToyQ model does not implement sourceTermSingleCell\n");
  // exit(1);

//  double kappa = this->data->optionalSimArgs[0];
  double tau_q = this->data->optionalSimArgs[1];
//  double zeta = this->data->optionalSimArgs[2];
  double tau_Pi = this->data->optionalSimArgs[3];
//  double eta = this->data->optionalSimArgs[4];
  double tau_pi = this->data->optionalSimArgs[5];

  // D
  source[0] = 0.0;
  // S1,2,3
  source[1] = 0.0; 
  source[2] = 0.0;
  source[3] = 0.0; 
  // Tau
  source[4] = 0.0;
  // Y1,2,3
  source[5] = (n / tau_q) * (aux[q1NS] - prims[q1]);
  source[6] = (n / tau_q) * (aux[q2NS] - prims[q2]);
  source[7] = (n / tau_q) * (aux[q3NS] - prims[q3]);
  // U
  source[8] = (n / tau_Pi) * (aux[PiNS] - prims[Pi]);
  // Z11,12,13,22,23,33
  source[9] = (n / tau_pi) * (aux[pi11NS] - prims[pi11]);
  source[10] = (n / tau_pi) * (aux[pi12NS] - prims[pi12]);
  source[11] = (n / tau_pi) * (aux[pi13NS] - prims[pi13]);
  source[12] = (n / tau_pi) * (aux[pi22NS] - prims[pi22]);
  source[13] = (n / tau_pi) * (aux[pi23NS] - prims[pi23]);
  source[14] = (n / tau_pi) * (aux[pi33NS] - prims[pi33]);
  
}

void IS::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

//  double kappa = this->data->optionalSimArgs[0];
  double tau_q = this->data->optionalSimArgs[1];
//  double zeta = this->data->optionalSimArgs[2];
  double tau_Pi = this->data->optionalSimArgs[3];
//  double eta = this->data->optionalSimArgs[4];
  double tau_pi = this->data->optionalSimArgs[5];

/* 

  // q_j,NS 10
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
          aux[ID(q1NS, i, j, k)] = -kappa* ( (aux[ID(T, i+1, j, k)] - aux[ID(T, i-1, j, k)])/(2*d->dx) );
          aux[ID(q2NS, i, j, k)] = -kappa* ( (aux[ID(T, i, j+1, k)] - aux[ID(T, i, j-1, k)])/(2*d->dy) );
          aux[ID(q3NS, i, j, k)] = -kappa* ( (aux[ID(T, i, j, k+1)] - aux[ID(T, i, j, k-1)])/(2*d->dz) );

          if(aux[ID(q1NS, i, j, k)] < -10 || aux[ID(q1NS, i, j, k)] > 10) {
            printf("(%i, %i, %i) ijk\n", i, j, k);
            printf("(%f, %f, %f) \n", aux[ID(T, i+1, j, k)], aux[ID(T, i-1, j, k)], aux[ID(q1NS, i, j, k)]);
          }
//          aux[ID(q1NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (log(aux[ID(T, i+1, j, k)]) - log(aux[ID(T, i-1, j, k)]))/(2*d->dx) );
//          aux[ID(q2NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (log(aux[ID(T, i, j+1, k)]) - log(aux[ID(T, i, j-1, k)]))/(2*d->dy) );
//          aux[ID(q3NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (log(aux[ID(T, i, j, k+1)]) - log(aux[ID(T, i, j, k-1)]))/(2*d->dz) );

      }
    }
  }  
  // Theta 20 then Pi,NS 13 
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(Theta, i, j, k)] = aux[ID(dWdt, i, j, k)] + (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(W, i, j+1, k)]*prims[ID(v2, i, j+1, k)] - aux[ID(W, i, j, k)]*prims[ID(v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v3, i, j, k+1)] - aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k-1)])/(2*d->dz);
        // Pi,NS = -zeta*Theta
        aux[ID(PiNS, i, j, k)] = -zeta * aux[ID(Theta, i, j, k)];
      }
    }
  }  
  // pi^l_j,NS 14
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // 11
        aux[ID(pi11NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 12
        aux[ID(pi12NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v2, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(W, i, j+1, k)]*prims[ID(v1, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v1, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 13
        aux[ID(pi13NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v3, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v3, i-1, j, k)])/(2*d->dx)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 22
        aux[ID(pi22NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v1, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v1, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 23
        aux[ID(pi23NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v3, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v2, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 33
        aux[ID(pi33NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
      }
    }
  }  

*/

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
        // Y1,2,3
        source[ID(Y1, i, j, k)] = (n / tau_q) * (aux[ID(q1NS, i, j, k)] - prims[ID(q1, i, j, k)]);
        source[ID(Y2, i, j, k)] = (n / tau_q) * (aux[ID(q1NS, i, j, k)] - prims[ID(q2, i, j, k)]);
        source[ID(Y3, i, j, k)] = (n / tau_q) * (aux[ID(q1NS, i, j, k)] - prims[ID(q3, i, j, k)]);        
        // U
        source[ID(U, i, j, k)] = (n / tau_Pi) * (aux[ID(PiNS, i, j, k)] - prims[ID(Pi, i, j, k)]);
        // Z11,12,13,22,23,33
        source[ID(Z11, i, j, k)] = (n / tau_pi) * (aux[ID(pi11NS, i, j, k)] - prims[ID(pi11, i, j, k)]);
        source[ID(Z12, i, j, k)] = (n / tau_pi) * (aux[ID(pi12NS, i, j, k)] - prims[ID(pi12, i, j, k)]);
        source[ID(Z13, i, j, k)] = (n / tau_pi) * (aux[ID(pi13NS, i, j, k)] - prims[ID(pi13, i, j, k)]);
        source[ID(Z22, i, j, k)] = (n / tau_pi) * (aux[ID(pi22NS, i, j, k)] - prims[ID(pi22, i, j, k)]); 
        source[ID(Z23, i, j, k)] = (n / tau_pi) * (aux[ID(pi23NS, i, j, k)] - prims[ID(pi23, i, j, k)]);
        source[ID(Z33, i, j, k)] = (n / tau_pi) * (aux[ID(pi33NS, i, j, k)] - prims[ID(pi33, i, j, k)]);                       
        
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
int ISresidual(void *ptr, int n, const double *x, double *fvec, int iflag)
{

//  Data * d(this->data);
  
  // Retrieve additional arguments
  IS::Args * args = (IS::Args*) ptr;

  // Values must make sense
  // Think this should be taken out for now - need new sensible values
  /*
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  */

  double vsqrd_rf = ((args->S1_rf - x[1])*(args->S1_rf - x[1]) + (args->S2_rf - x[2])*(args->S2_rf - x[2]) + (args->S3_rf - x[3])*(args->S3_rf - x[3]))/((args->Tau_rf + x[0])*(args->Tau_rf + x[0]));
  double W_rf(1 / sqrt(1 - vsqrd_rf));
  double n_rf(args->D_rf / W_rf);
  double rho_plus_p_rf = (args->Tau_rf + x[0])/(W_rf*W_rf) - args->Pi_rf;
  double v1_rf = (args->S1_rf - x[1])/((rho_plus_p_rf + args->Pi_rf)*W_rf*W_rf);
  double v2_rf = (args->S2_rf - x[2])/((rho_plus_p_rf + args->Pi_rf)*W_rf*W_rf);
  double v3_rf = (args->S3_rf - x[3])/((rho_plus_p_rf + args->Pi_rf)*W_rf*W_rf);
  double p_rf = (rho_plus_p_rf - n_rf)*((args->gamma-1)/args->gamma);
  double rho_rf = rho_plus_p_rf - p_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 0 || v1_rf >= 1 || v2_rf >= 1 || v3_rf >= 1) {
    printf("EEK");
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  
  double pi00_rf = args->pi11_rf + args->pi22_rf + args->pi33_rf;
  double qv_rf = args->q1_rf*v1_rf + args->q2_rf*v2_rf + args->q3_rf*v3_rf;
  double pi01_rf = args->pi11_rf*v1_rf + args->pi12_rf*v2_rf + args->pi13_rf*v3_rf; // dbl check sign on orthogonality relation
  double pi02_rf = args->pi12_rf*v1_rf + args->pi22_rf*v2_rf + args->pi23_rf*v3_rf;
  double pi03_rf = args->pi13_rf*v1_rf + args->pi23_rf*v2_rf + args->pi33_rf*v3_rf;

  fvec[0] = p_rf + args->Pi_rf + n_rf*W_rf - 2*qv_rf*W_rf - pi00_rf - x[0];
  fvec[1] = (args->q1_rf + qv_rf*v1_rf)*W_rf + pi01_rf - x[1];
  fvec[2] = (args->q2_rf + qv_rf*v2_rf)*W_rf + pi02_rf - x[2];
  fvec[3] = (args->q3_rf + qv_rf*v3_rf)*W_rf + pi03_rf - x[3];

//  printf("fvec[0]");

  return 0;
}

void IS::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  Data * d(this->data);

  // Do what we can first before root-find

  // Y1-3,U,Z11-33
  for (int ndissvar(0); ndissvar < 10; ndissvar++) {
    prims[q1+ndissvar] = cons[Y1+ndissvar] / cons[D];
  }

  aux[pi00] = prims[pi11] + prims[pi22] + prims[pi33]; // this one can be done here fine
  // what about these? need them in the guesses...
  aux[qv] = prims[q1]*prims[v1] + prims[q2]*prims[v2] + prims[q3]*prims[v3];
  aux[pi01] = prims[pi11]*prims[v1] + prims[pi12]*prims[v2] + prims[pi13]*prims[v3]; // dbl check sign on orthogonality relation
  aux[pi02] = prims[pi12]*prims[v1] + prims[pi22]*prims[v2] + prims[pi23]*prims[v3]; // dbl check sign on orthogonality relation
  aux[pi03] = prims[pi13]*prims[v1] + prims[pi23]*prims[v2] + prims[pi33]*prims[v3]; // dbl check sign on orthogonality relation

  // Hybrd1 set-up
  Args args;                      // Additional arguments structure
  const int n(4);                     // Size of system
  double sol[4];                      // Guess and solution vector
  double res[4];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.4e-8;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Set additional args for rootfind
  args.D_rf = cons[D];
  args.S1_rf = cons[S1];
  args.S2_rf = cons[S2];
  args.S3_rf = cons[S3];
  args.Tau_rf = cons[Tau];
  args.q1_rf = prims[q1];
  args.q2_rf = prims[q2];
  args.q3_rf = prims[q3];
  args.Pi_rf = prims[Pi];
  args.pi11_rf = aux[pi11];
  args.pi12_rf = aux[pi12];
  args.pi13_rf = aux[pi13];
  args.pi22_rf = aux[pi22];
  args.pi23_rf = aux[pi23];
  args.pi33_rf = aux[pi33];
  args.gamma = d->gamma;
  
  sol[0] = prims[p] + prims[Pi] + prims[n]*aux[W] - 2*aux[qv]*aux[W] - aux[pi00];
  sol[1] = (prims[q1] + aux[qv]*prims[v1])*aux[W] + aux[pi01];
  sol[2] = (prims[q2] + aux[qv]*prims[v2])*aux[W] + aux[pi02];
  sol[3] = (prims[q3] + aux[qv]*prims[v3])*aux[W] + aux[pi03];

  // Solve residual = 0
  info = __cminpack_func__(hybrd1) (&ISresidual, &args, n, sol, res,
                                    tol, wa, lwa);
  // If root find fails, add failed cell to the list
  if (info!=1) {
    //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
    throw std::runtime_error("C2P could not converge.\n");
  }
  aux[vsqrd] = ((cons[S1] - sol[1])*(cons[S1] - sol[1]) + (cons[S2] - sol[2])*(cons[S2] - sol[2]) + (cons[S3] - sol[3])*(cons[S3] - sol[3]))/((cons[Tau] + sol[0])*(cons[Tau] + sol[0]));
  aux[W] = 1 / (1-aux[vsqrd]);
  prims[n] = cons[D] / aux[W];
  aux[rho_plus_p] = (cons[Tau] + sol[0])/(aux[W]*aux[W]) - prims[Pi];
  prims[v1] = (cons[S1] - sol[1])/((aux[rho_plus_p] + prims[Pi])*aux[W]*aux[W]);
  prims[v2] = (cons[S2] - sol[2])/((aux[rho_plus_p] + prims[Pi])*aux[W]*aux[W]);  
  prims[v3] = (cons[S3] - sol[3])/((aux[rho_plus_p] + prims[Pi])*aux[W]*aux[W]);  
  prims[p] = (aux[rho_plus_p] - prims[n])*((d->gamma-1)/d->gamma);
  prims[rho] = aux[rho_plus_p] - prims[p];
  
  // Repeating the ones here that depend on v1,v2,v3...
  aux[qv] = prims[q1]*prims[v1] + prims[q2]*prims[v2] + prims[q3]*prims[v3];
  aux[pi01] = prims[pi11]*prims[v1] + prims[pi12]*prims[v2] + prims[pi13]*prims[v3]; // dbl check sign on orthogonality relation
  aux[pi02] = prims[pi12]*prims[v1] + prims[pi22]*prims[v2] + prims[pi23]*prims[v3]; // dbl check sign on orthogonality relation
  aux[pi03] = prims[pi13]*prims[v1] + prims[pi23]*prims[v2] + prims[pi33]*prims[v3]; // dbl check sign on orthogonality relation
        
  aux[e] = prims[p] / (prims[n]*(d->gamma-1));
  aux[T] = prims[p] / prims[n];     
   
}

void IS::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(4);                     // Size of system
  double sol[4];                      // Guess and solution vector
  double res[4];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.4e-8;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array
  std::vector<Failed> fails;          // Vector of failed structs. Stores location of failed cons2prims cells.
  
  /*
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
  /*
  for (int i(d->is -1); i < d->ie +1; i++) {
    for (int j(d->js -1); j < d->je +1; j++) {
      for (int k(d->ks -1); k < d->ke +1; k++) {
  */      
        for (int ndissvar(0); ndissvar < 10; ndissvar++) {
          prims[ID(q1+ndissvar, i, j, k)] = cons[ID(Y1+ndissvar, i, j, k)] / cons[ID(D, i, j, k)];
        }
        //aux[ID(W, i, j, k)] = 1 / ( 1 - (prims[ID(v1, i, j, k)]**2 + prims[ID(v2, i, j, k)]**2 + prims[ID(v3, i, j, k)]**2) )**0.5; // Point in re-calcing this here?
        aux[ID(qv, i, j, k)] = (prims[ID(q1, i, j, k)] * prims[ID(v1, i, j, k)]) 
                               + (prims[ID(q2, i, j, k)] * prims[ID(v2, i, j, k)]) 
                               + (prims[ID(q3, i, j, k)] * prims[ID(v3, i, j, k)]);
        aux[ID(pi00, i, j, k)] = prims[ID(pi11, i, j, k)] + prims[ID(pi22, i, j, k)] + prims[ID(pi33, i, j, k)];
        aux[ID(pi01, i, j, k)] = prims[ID(pi11, i, j, k)]*prims[ID(v1, i, j, k)] 
                                 + prims[ID(pi12, i, j, k)]*prims[ID(v2, i, j, k)] 
                                 + prims[ID(pi13, i, j, k)]*prims[ID(v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(pi02, i, j, k)] = prims[ID(pi12, i, j, k)]*prims[ID(v1, i, j, k)] 
                                 + prims[ID(pi22, i, j, k)]*prims[ID(v2, i, j, k)] 
                                 + prims[ID(pi23, i, j, k)]*prims[ID(v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(pi03, i, j, k)] = prims[ID(pi13, i, j, k)]*prims[ID(v1, i, j, k)] 
                                 + prims[ID(pi23, i, j, k)]*prims[ID(v2, i, j, k)] 
                                 + prims[ID(pi33, i, j, k)]*prims[ID(v3, i, j, k)]; // dbl check sign on orthogonality relation
            
        // Set additional args for rootfind
        args.D_rf = cons[ID(D, i, j, k)];
        args.S1_rf = cons[ID(S1, i, j, k)];
        args.S2_rf = cons[ID(S2, i, j, k)];
        args.S3_rf = cons[ID(S3, i, j, k)];
        args.Tau_rf = cons[ID(Tau, i, j, k)];
        args.q1_rf = prims[ID(q1, i, j, k)];
        args.q2_rf = prims[ID(q2, i, j, k)];
        args.q3_rf = prims[ID(q3, i, j, k)];
        args.Pi_rf = prims[ID(Pi, i, j, k)];
        args.pi11_rf = prims[ID(pi11, i, j, k)];
        args.pi12_rf = prims[ID(pi12, i, j, k)];
        args.pi13_rf = prims[ID(pi13, i, j, k)];
        args.pi22_rf = prims[ID(pi22, i, j, k)];
        args.pi23_rf = prims[ID(pi23, i, j, k)];
        args.pi33_rf = prims[ID(pi33, i, j, k)];
        args.gamma = d->gamma;
        
        sol[0] = prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)] + prims[ID(n, i, j, k)]*aux[ID(W, i, j, k)] 
                 - 2*aux[ID(qv, i, j, k)]*aux[ID(W, i, j, k)] - aux[ID(pi00, i, j, k)];
        sol[1] = (prims[ID(q1, i, j, k)] + aux[ID(qv, i, j, k)]*prims[ID(v1, i, j, k)])*aux[ID(W, i, j, k)] + aux[ID(pi01, i, j, k)];
        sol[2] = (prims[ID(q2, i, j, k)] + aux[ID(qv, i, j, k)]*prims[ID(v2, i, j, k)])*aux[ID(W, i, j, k)] + aux[ID(pi02, i, j, k)];
        sol[3] = (prims[ID(q3, i, j, k)] + aux[ID(qv, i, j, k)]*prims[ID(v3, i, j, k)])*aux[ID(W, i, j, k)] + aux[ID(pi03, i, j, k)];

        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&ISresidual, &args, n, sol, res,
                                          tol, wa, lwa);
//        info = __cminpack_func__(hybrd) (&ISresidual, &args, n, sol, res,
//                                          tol, maxfev, ml, mu, epsfcn, &diag[0], mode, factor, nprint, &nfev, &fjac[0][0], ldfjac, &r[0], lr, &qtf[0], &wa1[0], &wa2[0], &wa3[0], &wa4[0]);        
                                                                   
        // If root find fails, add failed cell to the list
        if (info!=1) {
          printf("(%f, %f) Y1,q1\n",  cons[ID(Y1, i, j, k)], prims[ID(q1, i, j, k)]);
          printf("(%f, %f, %f, %f) sol\n", sol[0],sol[1],sol[2],sol[3]);
          printf("%i info\n",info);
          printf("(%i, %i, %i) failed\n", i, j, k);
          printf("(%f, %f, %f, %f) res\n", res[0], res[1], res[2], res[3]);
          printf("(%f, %f, %f, %f) sol\n", sol[0], sol[1], sol[2], sol[3]);
          printf("(%f, %f, %f, %f, %f) prims\n",  prims[ID(p, i, j, k)], prims[ID(Pi, i, j, k)], prims[ID(n, i, j, k)], prims[ID(v1, i, j, k)], prims[ID(q1, i, j, k)]);
          printf("(%f, %f, %f, %f) aux\n",  aux[ID(W, i, j, k)], aux[ID(qv, i, j, k)], aux[ID(pi00, i, j, k)], aux[ID(pi01, i, j, k)]);
          printf("(%f, %f, %f, %f) cons\n",  cons[ID(Y1, i, j, k)], cons[ID(D, i, j, k)], cons[ID(S1, i, j, k)], cons[ID(Tau, i, j, k)]);
          exit(0);
          Failed fail = {i, j, k};
          fails.push_back(fail);
        }
        else {
          // Now have the correct values for Chi, Sigma123
          solution[ID(0, i, j, k)] = sol[0];
          solution[ID(1, i, j, k)] = sol[1];
          solution[ID(2, i, j, k)] = sol[2];
          solution[ID(3, i, j, k)] = sol[3];
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
/*
  for (int i(d->is -1); i < d->ie +1; i++) {
    for (int j(d->js -1); j < d->je +1; j++) {
      for (int k(d->ks -1); k < d->ke +1; k++) {
*/
        // C2P Scheme as outlined in HP/FYR
        aux[ID(vsqrd, i, j, k)] = ((cons[ID(S1, i, j, k)] - solution[ID(1, i, j, k)])*(cons[ID(S1, i, j, k)] - solution[ID(1, i, j, k)]) 
                                  + (cons[ID(S2, i, j, k)] - solution[ID(2, i, j, k)])*(cons[ID(S2, i, j, k)] - solution[ID(2, i, j, k)])
                                  + (cons[ID(S3, i, j, k)] - solution[ID(3, i, j, k)])*(cons[ID(S3, i, j, k)] - solution[ID(3, i, j, k)]))
                                  /((cons[ID(Tau, i, j, k)] + solution[ID(0, i, j, k)])*(cons[ID(Tau, i, j, k)] + solution[ID(0, i, j, k)]));
        aux[ID(W, i, j, k)] = 1 / (1-aux[ID(vsqrd, i, j, k)]);
        prims[ID(n, i, j, k)] = cons[ID(D, i, j, k)] / aux[ID(W, i, j, k)];
        aux[ID(rho_plus_p, i, j, k)] = (cons[ID(Tau, i, j, k)] + solution[ID(0, i, j, k)])/(aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)]) - prims[ID(Pi, i, j, k)];
        prims[ID(v1, i, j, k)] = (cons[ID(S1, i, j, k)] - solution[ID(1, i, j, k)])/((aux[ID(rho_plus_p, i, j, k)] 
                                 + prims[ID(Pi, i, j, k)])*aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)]);
        prims[ID(v2, i, j, k)] = (cons[ID(S2, i, j, k)] - solution[ID(2, i, j, k)])/((aux[ID(rho_plus_p, i, j, k)] 
                                 + prims[ID(Pi, i, j, k)])*aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)]);  
        prims[ID(v3, i, j, k)] = (cons[ID(S3, i, j, k)] - solution[ID(3, i, j, k)])/((aux[ID(rho_plus_p, i, j, k)] 
                                 + prims[ID(Pi, i, j, k)])*aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)]);  
        prims[ID(p, i, j, k)] = (aux[ID(rho_plus_p, i, j, k)] - prims[ID(n, i, j, k)])*((d->gamma-1)/d->gamma);
        prims[ID(rho, i, j, k)] = aux[ID(rho_plus_p, i, j, k)] - prims[ID(p, i, j, k)];

        // Again, repeating this here once the correct values for v1,v2,v3 have been set...
        aux[ID(qv, i, j, k)] = (prims[ID(q1, i, j, k)] * prims[ID(v1, i, j, k)]) + (prims[ID(q2, i, j, k)] * prims[ID(v2, i, j, k)]) 
                               + (prims[ID(q3, i, j, k)] * prims[ID(v3, i, j, k)]);
        // aux[ID(pi00, i, j, k)] = prims[ID(pi11, i, j, k)] + prims[ID(pi22, i, j, k)] + prims[ID(pi33, i, j, k)]; // Should be unncessary here
        aux[ID(pi01, i, j, k)] = prims[ID(pi11, i, j, k)]*prims[ID(v1, i, j, k)] 
                                 + prims[ID(pi12, i, j, k)]*prims[ID(v2, i, j, k)] 
                                 + prims[ID(pi13, i, j, k)]*prims[ID(v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(pi02, i, j, k)] = prims[ID(pi12, i, j, k)]*prims[ID(v1, i, j, k)] 
                                 + prims[ID(pi22, i, j, k)]*prims[ID(v2, i, j, k)] 
                                 + prims[ID(pi23, i, j, k)]*prims[ID(v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(pi03, i, j, k)] = prims[ID(pi13, i, j, k)]*prims[ID(v1, i, j, k)] 
                                 + prims[ID(pi23, i, j, k)]*prims[ID(v2, i, j, k)] 
                                 + prims[ID(pi33, i, j, k)]*prims[ID(v3, i, j, k)]; // dbl check sign on orthogonality relation
        
        aux[ID(e, i, j, k)] = prims[ID(p, i, j, k)] / (prims[ID(n, i, j, k)]*(d->gamma-1));
        aux[ID(T, i, j, k)] = prims[ID(p, i, j, k)] / prims[ID(n, i, j, k)];

      } // End k-loop
    } // End j-loop
  } // End i-loop  

        
  // Recompute NS variables here? Let's try it



  double kappa = this->data->optionalSimArgs[0];
//  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
//  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
//  double tau_pi = this->data->optionalSimArgs[5];
        
  // q_j,NS 10
  for (int i(d->is +1); i < d->ie -1; i++) {
    for (int j(d->js +1); j < d->je -1; j++) {
      for (int k(d->ks +1); k < d->ke -1; k++) {
        aux[ID(q1NS, i, j, k)] = -kappa* ( (aux[ID(T, i+1, j, k)] - aux[ID(T, i-1, j, k)])/(2*d->dx) );
        aux[ID(q2NS, i, j, k)] = -kappa* ( (aux[ID(T, i, j+1, k)] - aux[ID(T, i, j-1, k)])/(2*d->dy) );
        aux[ID(q3NS, i, j, k)] = -kappa* ( (aux[ID(T, i, j, k+1)] - aux[ID(T, i, j, k-1)])/(2*d->dz) );
        
        /*
        printf("(%i, %i, %i) ijk\n", i, j, k);
        printf("(%f, %f, %f) \n", aux[ID(T, i+1, j, k)], aux[ID(T, i-1, j, k)], aux[ID(q1NS, i, j, k)]);
        printf("(%f, %f, %f, %f) \n", prims[ID(p, i+1, j, k)], prims[ID(n, i+1, j, k)], prims[ID(p, i-1, j, k)], prims[ID(n, i-1, j, k)]);
        */
//          aux[ID(q1NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (log(aux[ID(T, i+1, j, k)]) - log(aux[ID(T, i-1, j, k)]))/(2*d->dx) );
//          aux[ID(q2NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (log(aux[ID(T, i, j+1, k)]) - log(aux[ID(T, i, j-1, k)]))/(2*d->dy) );
//          aux[ID(q3NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (log(aux[ID(T, i, j, k+1)]) - log(aux[ID(T, i, j, k-1)]))/(2*d->dz) );
 
        // Theta 20 then Pi,NS 13 
        aux[ID(Theta, i, j, k)] = aux[ID(dWdt, i, j, k)] + (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(W, i, j+1, k)]*prims[ID(v2, i, j+1, k)] - aux[ID(W, i, j, k)]*prims[ID(v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v3, i, j, k+1)] - aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k-1)])/(2*d->dz);
        // Pi,NS = -zeta*Theta
        aux[ID(PiNS, i, j, k)] = -zeta * aux[ID(Theta, i, j, k)];
  
        // pi^l_j,NS 14
        // 11
        aux[ID(pi11NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 12
        aux[ID(pi12NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v2, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(W, i, j+1, k)]*prims[ID(v1, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v1, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 13
        aux[ID(pi13NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v3, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v3, i-1, j, k)])/(2*d->dx)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 22
        aux[ID(pi22NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v1, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v1, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 23
        aux[ID(pi23NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v3, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v2, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 33
        aux[ID(pi33NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
      }
    }
  }  
  

  
}

void IS::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // printf("Calling primsToAll\n");

  // W, q_kv^k, pi^0_0
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(vsqrd, i, j, k)] = prims[ID(v1, i, j, k)]*prims[ID(v1, i, j, k)] 
                                  + prims[ID(v2, i, j, k)]*prims[ID(v2, i, j, k)] 
                                  + prims[ID(v3, i, j, k)]*prims[ID(v3, i, j, k)];
        aux[ID(W, i, j, k)] = 1 / sqrt( 1 - aux[ID(vsqrd, i, j, k)] );
        aux[ID(qv, i, j, k)] = (prims[ID(q1, i, j, k)] * prims[ID(v1, i, j, k)]) + (prims[ID(q2, i, j, k)] * prims[ID(v2, i, j, k)]) 
                               + (prims[ID(q3, i, j, k)] * prims[ID(v3, i, j, k)]);
        aux[ID(pi00, i, j, k)] = prims[ID(pi11, i, j, k)] + prims[ID(pi22, i, j, k)] + prims[ID(pi33, i, j, k)];

        aux[ID(e, i, j, k)] = prims[ID(p, i, j, k)] / (prims[ID(n, i, j, k)]*(d->gamma-1));
        prims[ID(rho, i, j, k)] = prims[ID(n, i, j, k)]*(1+aux[ID(e, i, j, k)]);
        aux[ID(T, i, j, k)] = prims[ID(p, i, j, k)] / prims[ID(n, i, j, k)];
        aux[ID(h, i, j, k)] = 1 + aux[ID(e, i, j, k)] + prims[ID(p, i, j, k)] / prims[ID(n, i, j, k)];
      }
    }
  }

  double kappa = this->data->optionalSimArgs[0];
//  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
//  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
//  double tau_pi = this->data->optionalSimArgs[5];
  
  // q_j,NS 10
  for (int i(1); i < d->Nx-1; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(1); k < d->Nz-1; k++) {
          aux[ID(q1NS, i, j, k)] = -kappa* ( (aux[ID(T, i+1, j, k)] - aux[ID(T, i-1, j, k)])/(2*d->dx) );
          aux[ID(q2NS, i, j, k)] = -kappa* ( (aux[ID(T, i, j+1, k)] - aux[ID(T, i, j-1, k)])/(2*d->dy) );
          aux[ID(q3NS, i, j, k)] = -kappa* ( (aux[ID(T, i, j, k+1)] - aux[ID(T, i, j, k-1)])/(2*d->dz) );
      }
    }
  }  
  // Theta 20 then Pi,NS 13 
  for (int i(1); i < d->Nx-1; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        aux[ID(Theta, i, j, k)] = aux[ID(dWdt, i, j, k)] + (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(W, i, j+1, k)]*prims[ID(v2, i, j+1, k)] - aux[ID(W, i, j, k)]*prims[ID(v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v3, i, j, k+1)] - aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k-1)])/(2*d->dz);
        // Pi,NS = -zeta*Theta
        aux[ID(PiNS, i, j, k)] = -zeta * aux[ID(Theta, i, j, k)];
      }
    }
  }  
  // pi^l_j,NS 14
  for (int i(1); i < d->Nx-1; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        // 11
        aux[ID(pi11NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 12
        aux[ID(pi12NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v2, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(W, i, j+1, k)]*prims[ID(v1, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v1, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 13
        aux[ID(pi13NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v3, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v3, i-1, j, k)])/(2*d->dx)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 22
        aux[ID(pi22NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v1, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v1, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 23
        aux[ID(pi23NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v3, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v2, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 33
        aux[ID(pi33NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
      }
    }
  }  

  // pi^0_j
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) { // Minus signs here (?)
         aux[ID(pi01, i, j, k)] = prims[ID(pi11, i, j, k)]*prims[ID(v1, i, j, k)] 
                                  + prims[ID(pi12, i, j, k)]*prims[ID(v2, i, j, k)] 
                                  + prims[ID(pi13, i, j, k)]*prims[ID(v3, i, j, k)];
         aux[ID(pi02, i, j, k)] = prims[ID(pi12, i, j, k)]*prims[ID(v1, i, j, k)] 
                                  + prims[ID(pi22, i, j, k)]*prims[ID(v2, i, j, k)] 
                                  + prims[ID(pi23, i, j, k)]*prims[ID(v3, i, j, k)];
         aux[ID(pi03, i, j, k)] = prims[ID(pi13, i, j, k)]*prims[ID(v1, i, j, k)] 
                                  + prims[ID(pi23, i, j, k)]*prims[ID(v2, i, j, k)] 
                                  + prims[ID(pi33, i, j, k)]*prims[ID(v3, i, j, k)];
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(D, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)];
        // S1,2,3
        cons[ID(S1, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)] * prims[ID(v1, i, j, k)] 
          + (prims[ID(q1, i, j, k)] + aux[ID(qv, i, j, k)] * prims[ID(v1, i, j, k)]) * aux[ID(W, i, j, k)] + aux[ID(pi01, i, j, k)];
        cons[ID(S2, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)] * prims[ID(v2, i, j, k)] 
          + (prims[ID(q2, i, j, k)] + aux[ID(qv, i, j, k)] * prims[ID(v2, i, j, k)]) * aux[ID(W, i, j, k)] + aux[ID(pi02, i, j, k)];
        cons[ID(S3, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)] * prims[ID(v3, i, j, k)] 
          + (prims[ID(q3, i, j, k)] + aux[ID(qv, i, j, k)] * prims[ID(v3, i, j, k)]) * aux[ID(W, i, j, k)] + aux[ID(pi03, i, j, k)];
        // Tau
        cons[ID(Tau, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]*aux[ID(W, i, j, k)] 
        - (prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)] + prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)]) 
        + 2*aux[ID(qv, i, j, k)]*aux[ID(W, i, j, k)] + aux[ID(pi00, i, j, k)];
        // Y1-3
        cons[ID(Y1, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(q1, i, j, k)];
        cons[ID(Y2, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(q2, i, j, k)];
        cons[ID(Y3, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(q3, i, j, k)];
        //printf("%f, %f, %f, q123\n", prims[ID(q1, i, j, k)], cons[ID(Y1, i, j, k)]);
        // U
        cons[ID(U, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(Pi, i, j, k)];
        // Z11-33
        cons[ID(Z11, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(pi11, i, j, k)];
        cons[ID(Z12, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(pi12, i, j, k)]; 
        cons[ID(Z13, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(pi13, i, j, k)];
        cons[ID(Z22, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(pi22, i, j, k)];        
        cons[ID(Z23, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(pi23, i, j, k)];
        cons[ID(Z33, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(pi33, i, j, k)];
      }  
    }
  }

/*

  for (int i(1); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(v1, i+1, j, k)]-prims[ID(v1, i-1, j, k)])/(2*d->dx);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(1, i, j, k)] = (prims[ID(v1, i, j+1, k)]-prims[ID(v1, i, j-1, k)])/(2*d->dy);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        aux[ID(2, i, j, k)] = (prims[ID(v1, i, j, k+1)]-prims[ID(v1, i, j, k-1)])/(2*d->dz);
      }
    }
  }
*/

}

void IS::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Dv
        f[ID(0, i, j, k)] = cons[ID(D, i, j, k)]*prims[ID(dir, i, j, k)];
        // Sv + ..
        for (int nvar(0); nvar < 3; nvar++) {
          f[ID(1+nvar, i, j, k)] = cons[ID(S1+nvar, i, j, k)]*prims[ID(dir, i, j, k)] + ( prims[ID(q1+dir, i, j, k)] * prims[ID(v1+nvar, i, j, k)]  
            - aux[ID(qv, i, j, k)]*prims[ID(v1+nvar, i, j, k)]*prims[ID(v1+dir, i, j, k)] ) * aux[ID(W, i, j, k)];
          // (p+Pi)delta_ij
          if (dir == nvar) {
            f[ID(1+nvar, i, j, k)] += (prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]);
          }
        }
        //  pi^i_j  
        if (dir == 0) {
          for (int nvar(0); nvar < 3; nvar++) {
            f[ID(1+nvar, i, j, k)] += prims[ID(pi11+nvar, i, j, k)];
          }
        } else if (dir == 1) {
          f[ID(1, i, j, k)] += prims[ID(pi12, i, j, k)];
          f[ID(2, i, j, k)] += prims[ID(pi22, i, j, k)];
          f[ID(3, i, j, k)] += prims[ID(pi23, i, j, k)];
        } else if (dir == 2) {
          f[ID(1, i, j, k)] += prims[ID(pi13, i, j, k)];
          f[ID(2, i, j, k)] += prims[ID(pi23, i, j, k)];
          f[ID(3, i, j, k)] += prims[ID(pi33, i, j, k)];
        } else {
          throw std::runtime_error("Flux direction is not 0, 1 or 2");
        }

        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Tau, i, j, k)] + prims[ID(p, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (prims[ID(q1+dir, i, j, k)] - aux[ID(qv, i, j, k)]*prims[ID(v1+dir, i, j, k)])*aux[ID(W, i, j, k)]
          + aux[ID(pi01+dir, i, j, k)];
        // Y1-3,U,Z11-33 *v
        for (int ndissvar(0); ndissvar < 10; ndissvar++) {                
          f[ID(Y1+ndissvar, i, j, k)] = cons[ID(Y1+ndissvar, i, j, k)]*prims[ID(dir, i, j, k)];
        }
      } // End k loop
    } // End j loop
  } // End i loop
}



/* Model with functional kappa dependence */

/*

ToyQFunctional::ToyQFunctional() : ToyQ()
{
}

ToyQFunctional::ToyQFunctional(Data * data) : ToyQ(data)
{
}

ToyQFunctional::~ToyQFunctional()
{
}

double kappa_of_T(double T, double kappa_0) {
  // return kappa_0 / (0.1 + T + T*T);
  // return kappa_0 / (1.0 + 1e-2*T);
  double kT = kappa_0 * T;
  return kT * T / (kT * kT + 0.25); // Andreas' Slides (bulk viscosity!!!)
}

double tau_q_of_T(double T, double tau_q_0) {
  // return tau_q_0 / (0.1 + 0.5 * T + T*T);
  // return tau_q_0 / (1.0 + 1e-3*T);
  double tT = tau_q_0 * T;
  return tT * T / (tT * tT + 0.25);
}

void ToyQFunctional::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  
  double kappa_0 = this->data->optionalSimArgs[0];
  double tau_q_0 = this->data->optionalSimArgs[1];

  source[0] = 0.0;
  for (int dir(0); dir < 3; dir++) {
    source[1+dir] = -(kappa_of_T(cons[0], kappa_0) * aux[dir] + prims[1+dir]) / tau_q_of_T(cons[0], tau_q_0);
  }
}

void ToyQFunctional::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  double kappa_0 = d->optionalSimArgs[0]; 
  double tau_q_0 = d->optionalSimArgs[1];

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        source[ID(0, i, j, k)] = 0.0;
        for (int dir(0); dir < 3; dir++) {
          source[ID(1+dir, i, j, k)] = -(kappa_of_T(cons[ID(0, i, j, k)], kappa_0) * aux[ID(dir, i, j, k)] +
                                         prims[ID(v2+dir, i, j, k)]) / tau_q_of_T(cons[ID(0, i, j, k)], tau_q_0);
        }
      }
    }
  }
}

*/

