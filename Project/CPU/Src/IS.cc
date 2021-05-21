#include "IS.h"
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
           pi11NS, pi12NS, pi13NS, pi22NS, pi23NS, pi33NS, Theta, a0, a1, a2, a3 };  


IS::IS() : Model()
{
  this->Ncons = ;
  this->Nprims = ;
  this->Naux = ;
}

IS::IS(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = ;
  this->Nprims = (this->data)->Nprims = ;
  this->Naux = (this->data)->Naux = ;

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
  this->data->primsLabels.push_back("q0");    this->data->primsLabels.push_back("qv");
  this->data->primsLabels.push_back("pi00");  this->data->primsLabels.push_back("pi01");
  this->data->primsLabels.push_back("pi02");  this->data->primsLabels.push_back("pi03");
  // 10
  this->data->primsLabels.push_back("q1NS");  this->data->primsLabels.push_back("q2NS");
  this->data->primsLabels.push_back("q3NS");
  // 13
  this->data->primsLabels.push_back("PiNS");    
  // 14
  this->data->primsLabels.push_back("pi11NS"); this->data->primsLabels.push_back("pi12NS");
  this->data->primsLabels.push_back("pi13NS"); this->data->primsLabels.push_back("pi22NS");
  this->data->primsLabels.push_back("pi23NS"); this->data->primsLabels.push_back("pi33NS");
  // 20
  this->data->auxLabels.push_back("Theta");  this->data->auxLabels.push_back("a0");
  this->data->auxLabels.push_back("a1");     this->data->auxLabels.push_back("a2");   
  this->data->auxLabels.push_back("a3");    
}

IS::~IS()
{
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
  source[12] = (n / tau_pi) * (aux[pi22NS] - prims[pi22);
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

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
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
        source[ID(8, i, j, k)] = (n / tau_Pi) * (aux[ID(PiNS, i, j, k)] - prims[ID(Pi, i, j, k)]);
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
int ISresidual(void *p, int n, const double *x, double *fvec, int iflag)
{

  // Retrieve additional arguments
  Args * args = (Args*) p;

  // Values must make sense
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  double Bsq(args->Bsq);
  double Ssq(args->Ssq);
  double BS(args->BS);
  double W(1 / sqrt(1 - x[0]));
  double rho(args->D / W);
  double h(x[1] / (rho * W * W));
  double pr((h - 1) * rho * (args->g - 1) / args->g);
  if (pr < 0 || rho < 0 || h < 0 || W < 1) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  // Values should be OK
  fvec[0] = (x[1] + Bsq) * (x[1] + Bsq) * x[0] - (2 * x[1] + Bsq) * BS * BS / (x[1] * x[1]) - Ssq;
  fvec[1] = x[1] + Bsq - pr - Bsq / (2 * W * W) - BS * BS / (2 * x[1] * x[1]) - args->D - args->tau;

  return 0;
}

void IS::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  // Y1-3,U,Z11-33
  for (int ndissvar(0); ndissvar < 10; nvar++) {
    prims[q1+ndissvar] = cons[Y1+ndissvar] / cons[D];
  }

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(4);                     // Size of system
  double sol[4];                      // Guess and solution vector
  double res[4];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.4e-8;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Set additional args for rootfind
  args.Tau_rf = cons[Tau];
  arg.Pi_rf = prims[Pi];
  

  sol[0] = prims[p] + prims[Pi] + prims[n]*aux[W] - 2*aux[qv]*aux[W] - aux[pi00];
  sol[1] = (prims[q1] + aux[qv]*prims[v1])*aux[W] + aux[pi01];
  sol[2] = (prims[q2] + aux[qv]*prims[v2])*aux[W] + aux[pi02];
  sol[3] = (prims[q3] + aux[qv]*prims[v3])*aux[W] + aux[pi03];

  // Solve residual = 0
  info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                    tol, wa, lwa);
  // If root find fails, add failed cell to the list
  if (info!=1) {
    //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
    throw std::runtime_error("C2P could not converge.\n");
  }
  









  // Note that this freezes the auxilliary variables - they are not computed in
  // this function as they are non-local.
}

void IS::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);
  
  // Y1-3,U,Z11-33
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int ndissvar(0); ndissvar < 10; nvar++) {
          prims[ID(ndissvar+8, i+1, j, k)] = cons[ID(ndissvar+5, i+1, j, k)] / cons[0];
        }
      }
    }
  }

/*
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(v1, i+1, j, k)]-prims[ID(v1, i-1, j, k)])/(2*d->dx);
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(1, i, j, k)] = (prims[ID(v1, i, j+1, k)]-prims[ID(v1, i, j-1, k)])/(2*d->dy);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(2, i, j, k)] = (prims[ID(v1, i, j, k+1)]-prims[ID(v1, i, j, k-1)])/(2*d->dz);
          }
        }
      }
    }
  }
*/
}

void IS::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // printf("Calling primsToAll\n");

  // W
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(W, i, j, k)] = 1 / ( 1 - (prims[ID(v1, i, j, k)]**2 + prims[ID(v2, i, j, k)]**2 + prims[ID(v3, i, j, k)]**2) )**0.5;
      }
    }
  }

  // q_k v^k
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(qv, i, j, k)] = (prims[ID(q1, i, j, k)] * prims[ID(v1, i, j, k)]) + (prims[ID(q2, i, j, k)] * prims[ID(v2, i, j, k)]) + (prims[ID(q3, i, j, k)] * prims[ID(v3, i, j, k)]);
      }
    }
  }
  
  // pi^0_0
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(pi00, i, j, k)] = prims[ID(pi11, i, j, k)] + prims[ID(pi22, i, j, k)] + prims[ID(pi33, i, j, k)];
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
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
          aux[ID(q1NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (ln(aux[ID(T, i+1, j, k)]) - ln(aux[ID(T, i-1, j, k)]))/(2*d->dx) );
          aux[ID(q2NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (ln(aux[ID(T, i, j+1, k)]) - ln(aux[ID(T, i, j-1, k)]))/(2*d->dy) );
          aux[ID(q3NS, i, j, k)] = -kappa*aux[ID(T, i, j, k)] * ( (ln(aux[ID(T, i, j, k+1)]) - ln(aux[ID(T, i, j, k-1)]))/(2*d->dz) );
      }
    }
  }  
  // Theta 20 then Pi,NS 13 
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(Theta, i, j, k)] = ( dtu0 + (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(W, i, j+1, k)]*prims[ID(v2, i, j+1, k)] - aux[ID(W, i, j, k)]*prims[ID(v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v3, i, j, k+1)] - aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k-1)])/(2*d->dz) );
        // Pi,NS = -zeta*Theta
        aux[ID(PiNS, i, j, k)] = -zeta * aux[ID(Theta, i, j, k)];
      }
    }
  }  
  // pi^l_j,NS 14
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // 11
        aux[ID(pi11NS, i, j, k)] = -2*eta*( (aux[ID(W, i+1, j, k)]*prims[ID(v1, i+1, j, k)] - aux[ID(W, i-1, j, k)]*prims[ID(v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v1, i, j, k)])**2)*aux[ID(Theta, i, j, k)] );
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
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])**2)*aux[ID(Theta, i, j, k)] );
        // 23
        aux[ID(pi23NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j+1, k)]*prims[ID(v3, i, j+1, k)] - aux[ID(W, i, j-1, k)]*prims[ID(v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(W, i, j, k+1)]*prims[ID(v2, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(W, i, j, k)]*prims[ID(v2, i, j, k)])*(aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)]))*aux[ID(Theta, i, j, k)] );
        // 33
        aux[ID(pi33NS, i, j, k)] = -2*eta*( (aux[ID(W, i, j, k+1)]*prims[ID(v1, i, j, k+1)] - aux[ID(W, i, j, k-1)]*prims[ID(v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(W, i, j, k)]*prims[ID(v3, i, j, k)])**2)*aux[ID(Theta, i, j, k)] );
      }
    }
  }  

  // pi^0_j
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) { // Minus signs here (?)
         aux[ID(pi01, i, j, k)] = prims[ID(pi11, i, j, k)]*prims[ID(v1, i, j, k)] + prims[ID(pi12, i, j, k)]*prims[ID(v2, i, j, k)] + prims[ID(pi13, i, j, k)]*prims[ID(v3, i, j, k)];
         aux[ID(pi02, i, j, k)] = prims[ID(pi12, i, j, k)]*prims[ID(v1, i, j, k)] + prims[ID(pi22, i, j, k)]*prims[ID(v2, i, j, k)] + prims[ID(pi23, i, j, k)]*prims[ID(v3, i, j, k)];
         aux[ID(pi03, i, j, k)] = prims[ID(pi13, i, j, k)]*prims[ID(v1, i, j, k)] + prims[ID(pi23, i, j, k)]*prims[ID(v2, i, j, k)] + prims[ID(pi33, i, j, k)]*prims[ID(v3, i, j, k)];
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(D, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)];
        // S1,2,3
        cons[ID(S1, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]**2 * prims[ID(v1, i, j, k)] 
          + (prims[ID(q1, i, j, k)] + aux[ID(qv, i, j, k)] * prims[ID(v1, i, j, k)]) * aux[ID(W, i, j, k)] + aux[ID(pi01, i, j, k)];
        cons[ID(S2, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]**2 * prims[ID(v2, i, j, k)] 
          + (prims[ID(q2, i, j, k)] + aux[ID(qv, i, j, k)] * prims[ID(v2, i, j, k)]) * aux[ID(W, i, j, k)] + aux[ID(pi02, i, j, k)];
        cons[ID(S3, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]**2 * prims[ID(v3, i, j, k)] 
          + (prims[ID(q3, i, j, k)] + aux[ID(qv, i, j, k)] * prims[ID(v3, i, j, k)]) * aux[ID(W, i, j, k)] + aux[ID(pi03, i, j, k)];
        // Tau
        cons[ID(Tau, i, j, k)] = (prims[ID(rho, i, j, k)] + prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]) * aux[ID(W, i, j, k)]**2 
        - (prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)] + prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)]) 
        + 2*aux[ID(qv, i, j, k)]*aux[ID(W, i, j, k)] + aux[ID(pi00, i, j, k)];
        // Y1-3
        cons[ID(Y1, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(q1, i, j, k)];
        cons[ID(Y2, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(q2, i, j, k)];
        cons[ID(Y3, i, j, k)] = prims[ID(n, i, j, k)] * aux[ID(W, i, j, k)] * prims[ID(q3, i, j, k)];
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
          f[ID(nvar+1, i, j, k)] = cons[ID(S1+nvar, i, j, k)]*prims[ID(dir, i, j, k)] + ( prims[ID(q1+dir, i, j, k)] * prims[ID(v1+nvar, i, j, k)]  
            - aux[ID(qv, i, j, k)]*prims[ID(v1+nvar, i, j, k)]*prims[ID(v1+dir, i, j, k)] ) * aux[ID(W, i, j, k)] + prims[ID(pi11+nvar+dir, i, j, k)]pi_ij ; // need to add pi21,pi31,pi32 for this to work 
          // (p+Pi)delta_ij
          if (dir == nvar-1) {
            f[ID(nvar, i, j, k)] += (prims[ID(p, i, j, k)] + prims[ID(Pi, i, j, k)]);
          }
        }
        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Tau, i, j, k)] + prims[ID(p, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (prims[ID(q1+dir, i, j, k)] - aux[ID(qv, i, j, k)]*prims[ID(v1+dir, i, j, k)]) + aux[ID(pi01+dir, i, j, k)];
        // Y1-3,U,Z11-33 *v
        for (int nvar(5); nvar < 15; nvar++) {                
          f[ID(nvar, i, j, k)] = cons[ID(nvar, i, j, k)]*prims[ID(dir, i, j, k)];
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

