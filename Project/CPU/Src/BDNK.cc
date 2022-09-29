#include "BDNK.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

template<typename T>
T sqr(T x) { return ((x) * (x)); }

IS::IS() : Model()
{
  this->Ncons = 5;
  this->Nprims = 6;
  this->Naux = 39;
}

IS::IS(Data * data, bool alt_C2P=false) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 5;
  this->Nprims = (this->data)->Nprims = 6;
  this->Naux = (this->data)->Naux = 39;

  // Solutions for C2P all cells
  solution = (double *) malloc(sizeof(double)*4*data->Nx*data->Ny*data->Nz);

  // Vector for storing variable at previous time-step...
  // the 7 here is for the 7 time-deriv variables currently needed... should be automated really not hard-set
  prev_vars = (double *) malloc(sizeof(double)*7*data->Nx*data->Ny*data->Nz); 

  smartGuesses = 0;
  
  alternative_C2P = alt_C2P;
  
  // 0  
  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("S1");
  this->data->consLabels.push_back("S2");  this->data->consLabels.push_back("S3");
  this->data->consLabels.push_back("Tau");
  // 5

  // 0
  this->data->primsLabels.push_back("v1");   this->data->primsLabels.push_back("v2");
  this->data->primsLabels.push_back("v3");
  // 3
  this->data->primsLabels.push_back("p");   this->data->primsLabels.push_back("rho");
  this->data->primsLabels.push_back("n");
  // 6


  // 0
  this->data->auxLabels.push_back("A");      this->data->auxLabels.push_back("Theta");
  this->data->auxLabels.push_back("Pi");
  // 3
  this->data->auxLabels.push_back("q0");     this->data->auxLabels.push_back("q1");
  this->data->auxLabels.push_back("q2");     this->data->auxLabels.push_back("q3");   
  this->data->auxLabels.push_back("qv");
  // 8
  this->data->auxLabels.push_back("pi11");   this->data->auxLabels.push_back("pi12");
  this->data->auxLabels.push_back("pi13");   this->data->auxLabels.push_back("pi22");
  this->data->auxLabels.push_back("pi23");   this->data->auxLabels.push_back("pi33");
  this->data->auxLabels.push_back("pi00");   this->data->auxLabels.push_back("pi01");
  this->data->auxLabels.push_back("pi02");   this->data->auxLabels.push_back("pi03");
  // 18
  this->data->auxLabels.push_back("h");      this->data->auxLabels.push_back("T");
  this->data->auxLabels.push_back("e");      this->data->auxLabels.push_back("W");
  // 22
  this->data->auxLabels.push_back("dpdt");   this->data->auxLabels.push_back("drhodt");
  this->data->auxLabels.push_back("dndt");   this->data->auxLabels.push_back("dv1dt");
  this->data->auxLabels.push_back("dv2dt");  this->data->auxLabels.push_back("dv3dt");
  this->data->auxLabels.push_back("dWdt");
  // 29
  this->data->auxLabels.push_back("a1");     this->data->auxLabels.push_back("a2");   
  this->data->auxLabels.push_back("a3");     this->data->auxLabels.push_back("vsqrd");
  this->data->auxLabels.push_back("rho_plus_p");
  // 34 - spatially calculated time derivatives for comparison purposes!
  this->data->auxLabels.push_back("dtv1");  this->data->auxLabels.push_back("dtv2");
  this->data->auxLabels.push_back("dtv3");
  this->data->auxLabels.push_back("dtp");  this->data->auxLabels.push_back("dtrho");
  // 39
}

IS::~IS()
{
  free(solution);
}

// First and second order "minmod" functions for slope-limiting
double minmodGradFO(double im1, double i, double ip1, double dX) {

  double FDGrad = (-1.0*i + 1*ip1)/dX;
  double BDGrad = (1.0*i - 1*im1)/dX;
  if ( (FDGrad < 0 && BDGrad > 0) || (FDGrad > 0 && BDGrad < 0) ) {
    return 0;
  } else {
    return abs(FDGrad) < abs(BDGrad) ? FDGrad : BDGrad;
  }
}

double minmodGradSO(double im2, double im1, double i, double ip1, double ip2, double dX) {
  
  double FDGrad = (-1.5*i + 2*ip1 - 0.5*ip2)/dX;
  double BDGrad = (1.5*i - 2*im1 + 0.5*im2)/dX;
  if ( (FDGrad < 0 && BDGrad > 0) || (FDGrad > 0 && BDGrad < 0) ) {
    return 0;
  } else {
    return abs(FDGrad) < abs(BDGrad) ? FDGrad : BDGrad;
  }
}



void IS::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
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

void IS::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  // double kappa = this->data->optionalSimArgs[0];
  // double tau_q = this->data->optionalSimArgs[1];
  // double zeta = this->data->optionalSimArgs[2];
  // double tau_Pi = this->data->optionalSimArgs[3];
  // double eta = this->data->optionalSimArgs[4];
  // double tau_epsilon = this->data->optionalSimArgs[5];

/*

  double dxT;
  double dyT;
  double dzT;
  
  double dxux;
  double dyuy;
  double dzuz;

  double dxuy;
  double dxuz;
  double dyux;
  double dyuz;
  double dzux;
  double dzuy;

  // q_j,NS 10
  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

        dxT = minmodGradFO(aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i+1, j, k)], d->dx);
        dyT = minmodGradFO(aux[ID(Aux::T, i, j-1, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j+1, k)], d->dy);
        dzT = minmodGradFO(aux[ID(Aux::T, i, j, k-1)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j, k+1)], d->dz);
        
        dxux = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)], d->dx);
        dyuy = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v2, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)], d->dy);
        dzuz = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);                          

        dxuy = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)], d->dx);
        dxuz = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)], d->dx);
       
        dyux = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)], d->dy);
        dyuz = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)], d->dy);

        dzux = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)], d->dz);  
        dzuy = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);  

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
  
  double E_rf = args->Tau_rf + args->D_rf; 
  double vsqrd_rf = ((args->S1_rf - x[1])*(args->S1_rf - x[1]) + (args->S2_rf - x[2])*(args->S2_rf - x[2]) + (args->S3_rf - x[3])*(args->S3_rf - x[3]))/((E_rf + x[0])*(E_rf + x[0]));
  double W_rf(1 / sqrt(1 - vsqrd_rf));
  double n_rf(args->D_rf / W_rf);
  double rho_plus_p_rf = ((E_rf + x[0])/(W_rf*W_rf)) - (args->Pi_rf + args->A_rf);
  double v1_rf = (args->S1_rf - x[1])/((rho_plus_p_rf + args->Pi_rf + args->A_rf)*W_rf*W_rf);
  double v2_rf = (args->S2_rf - x[2])/((rho_plus_p_rf + args->Pi_rf + args->A_rf)*W_rf*W_rf);
  double v3_rf = (args->S3_rf - x[3])/((rho_plus_p_rf + args->Pi_rf + args->A_rf)*W_rf*W_rf);
  double p_rf = (rho_plus_p_rf - n_rf)*((args->gamma-1)/args->gamma);
  double rho_rf = rho_plus_p_rf - p_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 1 || n_rf < 0 || abs(v1_rf) >= 1 || abs(v2_rf) >= 1 || abs(v3_rf) >= 1 || vsqrd_rf >= 1) {
    printf("EEK");
    fvec[0] = fvec[1] = fvec[2] = fvec[3] = 1e6;
    return 0;
  }
  
  double pi00_rf = args->pi11_rf + args->pi22_rf + args->pi33_rf;
  double qv_rf = args->q1_rf*v1_rf + args->q2_rf*v2_rf + args->q3_rf*v3_rf;
  double pi01_rf = args->pi11_rf*v1_rf + args->pi12_rf*v2_rf + args->pi13_rf*v3_rf; // dbl check sign on orthogonality relation
  double pi02_rf = args->pi12_rf*v1_rf + args->pi22_rf*v2_rf + args->pi23_rf*v3_rf;
  double pi03_rf = args->pi13_rf*v1_rf + args->pi23_rf*v2_rf + args->pi33_rf*v3_rf;

  fvec[0] = p_rf + args->Pi_rf + args->A_rf - 2*qv_rf*W_rf - pi00_rf - x[0];
  fvec[1] = (args->q1_rf + qv_rf*v1_rf)*W_rf + pi01_rf - x[1];
  fvec[2] = (args->q2_rf + qv_rf*v2_rf)*W_rf + pi02_rf - x[2];
  fvec[3] = (args->q3_rf + qv_rf*v3_rf)*W_rf + pi03_rf - x[3];

  return 0;
}

int ISAlternativeResidual(void *ptr, int n, const double *x, double *fvec, int iflag)
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
  
  double vsqrd_rf = x[0]*x[0]*((args->S1_rf - x[1])*(args->S1_rf - x[1]) + (args->S2_rf - x[2])*(args->S2_rf - x[2]) + (args->S3_rf - x[3])*(args->S3_rf - x[3]))/((args->D_rf)*(args->D_rf));
  double W_rf(1 / sqrt(1 - vsqrd_rf));
  double n_rf(args->D_rf / W_rf);
  double v1_rf = x[0]*(args->S1_rf - x[1])/args->D_rf;
  double v2_rf = x[0]*(args->S2_rf - x[2])/args->D_rf;
  double v3_rf = x[0]*(args->S3_rf - x[3])/args->D_rf;
  double pi00_rf = args->pi11_rf + args->pi22_rf + args->pi33_rf;
  double qv_rf = args->q1_rf*v1_rf + args->q2_rf*v2_rf + args->q3_rf*v3_rf;
  double p_rf = args->D_rf*((1/x[0]) -1) - args->Pi_rf + 2*qv_rf*W_rf + pi00_rf - args->Tau_rf;
  double rho_rf = n_rf + p_rf/(args->gamma-1);
  double H_rf = 1 + (p_rf*(args->gamma/(args->gamma-1)) + args->Pi_rf)/n_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 1 || n_rf < 0 || abs(v1_rf) >= 1 || abs(v2_rf) >= 1 || abs(v3_rf) >= 1 || vsqrd_rf >= 1) {
    printf("EEK");
    fvec[0] = fvec[1] = fvec[2] = fvec[3] = 1e6;
    return 0;
  }
  
  double pi01_rf = args->pi11_rf*v1_rf + args->pi12_rf*v2_rf + args->pi13_rf*v3_rf; // dbl check sign on orthogonality relation
  double pi02_rf = args->pi12_rf*v1_rf + args->pi22_rf*v2_rf + args->pi23_rf*v3_rf;
  double pi03_rf = args->pi13_rf*v1_rf + args->pi23_rf*v2_rf + args->pi33_rf*v3_rf;

  fvec[0] = 1/(W_rf*H_rf) - x[0];
  fvec[1] = (args->q1_rf + qv_rf*v1_rf)*W_rf + pi01_rf - x[1];
  fvec[2] = (args->q2_rf + qv_rf*v2_rf)*W_rf + pi02_rf - x[2];
  fvec[3] = (args->q3_rf + qv_rf*v3_rf)*W_rf + pi03_rf - x[3];

  return 0;
}




void IS::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  Data * d(this->data);

  // Do what we can first before root-find

  aux[Aux::pi00] = aux[Aux::pi11] + aux[Aux::pi22] + aux[Aux::pi33]; // this one can be done here fine
  // what about these? need them in the guesses...
  aux[Aux::qv] = aux[Aux::q1]*prims[Prims::v1] + aux[Aux::q2]*prims[Prims::v2] + aux[Aux::q3]*prims[Prims::v3];
  aux[Aux::pi01] = aux[Aux::pi11]*prims[Prims::v1] + aux[Aux::pi12]*prims[Prims::v2] + aux[Aux::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi02] = aux[Aux::pi12]*prims[Prims::v1] + aux[Aux::pi22]*prims[Prims::v2] + aux[Aux::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi03] = aux[Aux::pi13]*prims[Prims::v1] + aux[Aux::pi23]*prims[Prims::v2] + aux[Aux::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation

  // Hybrd1 set-up
  Args args;                      // Additional arguments structure
  const int sys_size(4);                     // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1e-6;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Set additional args for rootfind
  args.D_rf = cons[Cons::D];
  args.S1_rf = cons[Cons::S1];
  args.S2_rf = cons[Cons::S2];
  args.S3_rf = cons[Cons::S3];
  args.Tau_rf = cons[Cons::Tau];
  args.A_rf = aux[Aux::A];
  args.q1_rf = aux[Aux::q1];
  args.q2_rf = aux[Aux::q2];
  args.q3_rf = aux[Aux::q3];
  args.Pi_rf = aux[Aux::Pi];
  args.pi11_rf = aux[Aux::pi11];
  args.pi12_rf = aux[Aux::pi12];
  args.pi13_rf = aux[Aux::pi13];
  args.pi22_rf = aux[Aux::pi22];
  args.pi23_rf = aux[Aux::pi23];
  args.pi33_rf = aux[Aux::pi33];
  args.gamma = d->gamma;
  
  if (alternative_C2P) {
  
    sol[0] = 1/(aux[Aux::W]*(1 + (prims[Prims::p]*(d->gamma/(d->gamma-1)) + aux[Aux::Pi])/prims[Prims::n]));
    sol[1] = (aux[Aux::q1] + aux[Aux::qv]*prims[Prims::v1])*aux[Aux::W] + aux[Aux::pi01];
    sol[2] = (aux[Aux::q2] + aux[Aux::qv]*prims[Prims::v2])*aux[Aux::W] + aux[Aux::pi02];
    sol[3] = (aux[Aux::q3] + aux[Aux::qv]*prims[Prims::v3])*aux[Aux::W] + aux[Aux::pi03];
  
    // Solve residual = 0
    info = __cminpack_func__(hybrd1) (&ISAlternativeResidual, &args, sys_size, sol, res,
                                      tol, wa, lwa);
    // If root find fails, add failed cell to the list
    if (info!=1) {
      //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
      throw std::runtime_error("C2P could not converge.\n");
    }
    aux[Aux::vsqrd] = sol[0]*sol[0]*((cons[Cons::S1] - sol[1])*(cons[Cons::S1] - sol[1]) + (cons[Cons::S2] - sol[2])*(cons[Cons::S2] - sol[2]) 
                      + (cons[Cons::S3] - sol[3])*(cons[Cons::S3] - sol[3]))/(cons[Cons::D]*cons[Cons::D]);
    aux[Aux::W] = (1 / sqrt(1 - aux[Aux::vsqrd]));
    prims[Prims::n] = cons[Cons::D] / aux[Aux::W];
    prims[Prims::v1] = sol[0]*(cons[Cons::S1] - sol[1])/cons[Cons::D];
    prims[Prims::v2] = sol[0]*(cons[Cons::S2] - sol[2])/cons[Cons::D];
    prims[Prims::v3] = sol[0]*(cons[Cons::S3] - sol[3])/cons[Cons::D];
    aux[Aux::pi00] = aux[Aux::pi11] + aux[Aux::pi22] + aux[Aux::pi33]; // not sure we need this again here tbh
    aux[Aux::qv] = aux[Aux::q1]*prims[Prims::v1] + aux[Aux::q2]*prims[Prims::v2] + aux[Aux::q3]*prims[Prims::v3];
    prims[Prims::p] = cons[Cons::D]*(1/sol[0] -1) - aux[Aux::Pi] + 2*aux[Aux::qv]*aux[Aux::W] + aux[Aux::pi00] - cons[Cons::Tau];
    prims[Prims::rho] = prims[Prims::n] + prims[Prims::p]/(d->gamma-1);
    
    // Repeating the ones here that depend on v1,v2,v3...
    aux[Aux::pi01] = aux[Aux::pi11]*prims[Prims::v1] + aux[Aux::pi12]*prims[Prims::v2] + aux[Aux::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi02] = aux[Aux::pi12]*prims[Prims::v1] + aux[Aux::pi22]*prims[Prims::v2] + aux[Aux::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi03] = aux[Aux::pi13]*prims[Prims::v1] + aux[Aux::pi23]*prims[Prims::v2] + aux[Aux::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
          
    aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
    aux[Aux::T] = prims[Prims::p] / prims[Prims::n]; 
    
  
  } else {
  
    sol[0] = prims[Prims::p] + aux[Aux::Pi] + aux[Aux::A] - 2*aux[Aux::qv]*aux[Aux::W] - aux[Aux::pi00];
    sol[1] = (aux[Aux::q1] + aux[Aux::qv]*prims[Prims::v1])*aux[Aux::W] + aux[Aux::pi01];
    sol[2] = (aux[Aux::q2] + aux[Aux::qv]*prims[Prims::v2])*aux[Aux::W] + aux[Aux::pi02];
    sol[3] = (aux[Aux::q3] + aux[Aux::qv]*prims[Prims::v3])*aux[Aux::W] + aux[Aux::pi03];
  
    // Solve residual = 0
    info = __cminpack_func__(hybrd1) (&ISresidual, &args, sys_size, sol, res,
                                      tol, wa, lwa);
    // If root find fails, add failed cell to the list
    if (info!=1) {
      //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
      throw std::runtime_error("C2P could not converge.\n");
    }
    aux[Aux::vsqrd] = ((cons[Cons::S1] - sol[1])*(cons[Cons::S1] - sol[1]) + (cons[Cons::S2] - sol[2])*(cons[Cons::S2] - sol[2]) 
                       + (cons[Cons::S3] - sol[3])*(cons[Cons::S3] - sol[3]))
                       /((cons[Cons::Tau] + cons[Cons::D] + sol[0])*(cons[Cons::Tau]  + cons[Cons::D] + sol[0]));
    aux[Aux::W] = 1 / sqrt((1-aux[Aux::vsqrd]));
    prims[Prims::n] = cons[Cons::D] / aux[Aux::W];
    aux[Aux::rho_plus_p] = (cons[Cons::Tau] + cons[Cons::D] + sol[0])/(aux[Aux::W]*aux[Aux::W]) - (aux[Aux::Pi] + aux[Aux::A]);
    prims[Prims::v1] = (cons[Cons::S1] - sol[1])/((aux[Aux::rho_plus_p] + aux[Aux::Pi] + aux[Aux::A])*aux[Aux::W]*aux[Aux::W]);
    prims[Prims::v2] = (cons[Cons::S2] - sol[2])/((aux[Aux::rho_plus_p] + aux[Aux::Pi] + aux[Aux::A])*aux[Aux::W]*aux[Aux::W]);  
    prims[Prims::v3] = (cons[Cons::S3] - sol[3])/((aux[Aux::rho_plus_p] + aux[Aux::Pi] + aux[Aux::A])*aux[Aux::W]*aux[Aux::W]);  
    prims[Prims::p] = (aux[Aux::rho_plus_p] - prims[Prims::n])*((d->gamma-1)/d->gamma);
    prims[Prims::rho] = aux[Aux::rho_plus_p] - prims[Prims::p];
    
    // Repeating the ones here that depend on v1,v2,v3...
    aux[Aux::qv] = aux[Aux::q1]*prims[Prims::v1] + aux[Aux::q2]*prims[Prims::v2] + aux[Aux::q3]*prims[Prims::v3];
    aux[Aux::pi01] = aux[Aux::pi11]*prims[Prims::v1] + aux[Aux::pi12]*prims[Prims::v2] + aux[Aux::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi02] = aux[Aux::pi12]*prims[Prims::v1] + aux[Aux::pi22]*prims[Prims::v2] + aux[Aux::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi03] = aux[Aux::pi13]*prims[Prims::v1] + aux[Aux::pi23]*prims[Prims::v2] + aux[Aux::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
          
    aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
    aux[Aux::T] = prims[Prims::p] / prims[Prims::n];     
 
  }
   
}

void IS::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int sys_size(4);              // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1e-6;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
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

  // std::cout << d->is << "\t" << d->ie << "\t" << d->js << "\t" << d->je << std::endl;

  // for (int i(d->is); i < d->ie; i++) {
  //   for (int j(d->js); j < d->je; j++) {
  //     for (int k(d->ks); k < d->ke; k++) {

  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

  // for (int i(d->is-1); i < d->ie+1; i++) {
  //   for (int j(d->js-1); j < d->je+1; j++) {
  //     for (int k(d->ks-1); k < d->ke+1; k++) {

        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - (prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
          + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)]
          + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)]) ); // Point in re-calcing this here?
        aux[ID(Aux::qv, i, j, k)] = (aux[ID(Aux::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) 
                               + (aux[ID(Aux::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (aux[ID(Aux::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
        aux[ID(Aux::pi00, i, j, k)] = aux[ID(Aux::pi11, i, j, k)] + aux[ID(Aux::pi22, i, j, k)] + aux[ID(Aux::pi33, i, j, k)];
        aux[ID(Aux::pi01, i, j, k)] = aux[ID(Aux::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + aux[ID(Aux::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + aux[ID(Aux::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi02, i, j, k)] = aux[ID(Aux::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + aux[ID(Aux::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + aux[ID(Aux::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi03, i, j, k)] = aux[ID(Aux::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + aux[ID(Aux::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + aux[ID(Aux::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation

        sol[0] = prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)] - 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] - aux[ID(Aux::pi00, i, j, k)];
        sol[1] = (aux[ID(Aux::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
        sol[2] = (aux[ID(Aux::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
        sol[3] = (aux[ID(Aux::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];

        // if (i==3 && j ==0 && k==0)
        //   std::cout << cons[ID(Cons::S1, i, j, k)] << "\t" << prims[ID(Prims::rho, i, j, k)] << "\t" << prims[ID(Prims::p, i, j, k)] << "\t" 
        //   << aux[ID(Aux::Pi, i, j, k)] << "\t" << aux[ID(Aux::A, i, j, k)] << "\t" << aux[ID(Aux::q1, i, j, k)] << "\t" << aux[ID(Aux::qv, i, j, k)] << 
        //  "\t" <<  prims[ID(Prims::v1, i, j, k)] << "\t" << aux[ID(Aux::pi01, i, j, k)] << "\t" << aux[ID(Aux::W, i, j, k)] << std::endl;

        // Set additional args for rootfind
        args.D_rf = cons[ID(Cons::D, i, j, k)];
        args.S1_rf = cons[ID(Cons::S1, i, j, k)];
        args.S2_rf = cons[ID(Cons::S2, i, j, k)];
        args.S3_rf = cons[ID(Cons::S3, i, j, k)];
        args.Tau_rf = cons[ID(Cons::Tau, i, j, k)];
        args.A_rf = aux[ID(Aux::A, i, j, k)];
        args.q1_rf = aux[ID(Aux::q1, i, j, k)];
        args.q2_rf = aux[ID(Aux::q2, i, j, k)];
        args.q3_rf = aux[ID(Aux::q3, i, j, k)];
        args.Pi_rf = aux[ID(Aux::Pi, i, j, k)];
        args.pi11_rf = aux[ID(Aux::pi11, i, j, k)];
        args.pi12_rf = aux[ID(Aux::pi12, i, j, k)];
        args.pi13_rf = aux[ID(Aux::pi13, i, j, k)];
        args.pi22_rf = aux[ID(Aux::pi22, i, j, k)];
        args.pi23_rf = aux[ID(Aux::pi23, i, j, k)];
        args.pi33_rf = aux[ID(Aux::pi33, i, j, k)];
        args.gamma = d->gamma;
 
        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&ISresidual, &args, sys_size, sol, res,
                                          tol, wa, lwa);
//        info = __cminpack_func__(hybrd) (&ISresidual, &args, sys_size, sol, res,
//                                          tol, maxfev, ml, mu, epsfcn, &diag[0], mode, factor, nprint, &nfev, &fjac[0][0], ldfjac, &r[0], lr, &qtf[0], &wa1[0], &wa2[0], &wa3[0], &wa4[0]);        

                                                                   
        // If root find fails, add failed cell to the list
        if (info!=1) {
          printf("%i info\n",info);
          printf("(%i, %i, %i) failed\n", i, j, k);
          printf("(%g, %g, %g, %g) res\n", res[0], res[1], res[2], res[3]);
          printf("(%g, %g, %g, %g) sol\n", sol[0], sol[1], sol[2], sol[3]);
          std::cout << "Prims ";
          for (int vz(0); vz < d->Nprims; vz++) {
            std::cout << d->prims[ID(vz, i, j, k)] << " ";
          }
          std::cout << std::endl;
          std::cout << "Cons  ";
          for (int vz(0); vz < d->Ncons; vz++) {
            std::cout << d->cons[ID(vz, i, j, k)] << " ";
          }
          std::cout << std::endl;
          std::cout << "Aux   ";
          for (int vz(0); vz < d->Naux; vz++) {
            std::cout << d->aux[ID(vz, i, j, k)] << " ";
          }
          std::cout << std::endl;
          //printf("(%f, %f, %f, %f, %f) prims\n",  prims[ID(Prims::p, i, j, k)], aux[ID(Aux::Pi, i, j, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::v1, i, j, k)], aux[ID(Aux::q1, i, j, k)]);
          //printf("(%f, %f, %f, %f) aux  \n",  aux[ID(Aux::W, i, j, k)], aux[ID(Aux::qv, i, j, k)], aux[ID(Aux::pi00, i, j, k)], aux[ID(Aux::pi01, i, j, k)]);
          //printf("(%f, %f, %f, %f) cons \n",  cons[ID(Cons::Y1, i, j, k)], cons[ID(Cons::D, i, j, k)], cons[ID(Cons::S1, i, j, k)], cons[ID(Cons::Tau, i, j, k)]);
          exit(0);
          Failed fail = {i, j, k};
          fails.push_back(fail);
        } else {
          // Now have the correct values for Chi, Sigma123
          solution[ID(0, i, j, k)] = sol[0];
          solution[ID(1, i, j, k)] = sol[1];
          solution[ID(2, i, j, k)] = sol[2];
          solution[ID(3, i, j, k)] = sol[3];
        }      
      
      } // End k-loop
    } // End j-loop
  } // End i-loop
  
  // for (int i(d->is-1); i < d->ie+1; i++) {
  //   for (int j(d->js-1); j < d->je+1; j++) {
  //     for (int k(d->ks-1); k < d->ke+1; k++) {

  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

  // for (int i(d->is); i < d->ie; i++) {
  //   for (int j(d->js); j < d->je; j++) {
  //     for (int k(d->ks); k < d->ke; k++) {

        // C2P Scheme as outlined in HP/FYR
        aux[ID(Aux::vsqrd, i, j, k)] = ((cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])*(cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)]) 
                                  + (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])*(cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])
                                  + (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])*(cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)]))
                                  /((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])*(cons[ID(Cons::Tau, i, j, k)] 
                                  + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)]));
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt((1-aux[ID(Aux::vsqrd, i, j, k)]));
        prims[ID(Prims::n, i, j, k)] = cons[ID(Cons::D, i, j, k)] / aux[ID(Aux::W, i, j, k)];
        aux[ID(Aux::rho_plus_p, i, j, k)] = ((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])
                                            /(aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)])) - (aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]);
        prims[ID(Prims::v1, i, j, k)] = (cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                  + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);
        prims[ID(Prims::v2, i, j, k)] = (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                  + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);  
        prims[ID(Prims::v3, i, j, k)] = (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                  + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);  
        prims[ID(Prims::p, i, j, k)] = (aux[ID(Aux::rho_plus_p, i, j, k)] - prims[ID(Prims::n, i, j, k)])*((d->gamma-1)/d->gamma);
        prims[ID(Prims::rho, i, j, k)] = aux[ID(Aux::rho_plus_p, i, j, k)] - prims[ID(Prims::p, i, j, k)];
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];

      } // End k-loop
    } // End j-loop
  } // End i-loop  

  // double kappa = d->optionalSimArgs[0];
  // double tau_q = d->optionalSimArgs[1];
  // double zeta = d->optionalSimArgs[2];
  // double tau_Pi = d->optionalSimArgs[3];
  // double eta = d->optionalSimArgs[4];
  // double tau_epsilon = d->optionalSimArgs[5];

  double beta_epsilon;
  double beta_n;  

  double kappa { d->optionalSimArgs[0] };
  double tau_q;
  double zeta { d->optionalSimArgs[2] };
  double eta;
  double tau_epsilon;
  double tau_Pi;
  
  double dxrho;
  double dyrho;
  double dzrho;

  double dxn;
  double dyn;
  double dzn;

  double dxux;
  double dyuy;
  double dzuz;

  double dxuy;
  double dxuz;
  double dyux;
  double dyuz;
  double dzux;
  double dzuy;

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

          dxrho = minmodGradFO(prims[ID(Prims::rho, i-1, j, k)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i+1, j, k)], d->dx);
          dyrho = minmodGradFO(prims[ID(Prims::rho, i, j-1, k)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i, j+1, k)], d->dy);
          dzrho = minmodGradFO(prims[ID(Prims::rho, i, j, k-1)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i, j, k+1)], d->dz);
          
          dxn = minmodGradFO(prims[ID(Prims::n, i-1, j, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::n, i+1, j, k)], d->dx);
          dyn = minmodGradFO(prims[ID(Prims::n, i, j-1, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::n, i, j+1, k)], d->dy);
          dzn = minmodGradFO(prims[ID(Prims::n, i, j, k-1)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::n, i, j, k+1)], d->dz);

          dxux = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                              aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)], d->dx);
          dyuy = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v2, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                              aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)], d->dy);
          dzuz = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                              aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);                          

          dxuy = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                              aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)], d->dx);
          dxuz = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                              aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)], d->dx);
        
          dyux = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                              aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)], d->dy);
          dyuz = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                              aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)], d->dy);

          dzux = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                              aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)], d->dz);  
          dzuy = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                              aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);  

          tau_q = (3.0/4.0)*d->optionalSimArgs[1]*pow(prims[ID(Prims::rho, i, j, k)],0.25);
          eta = d->optionalSimArgs[3]*pow(prims[ID(Prims::rho, i, j, k)],0.25);
          tau_epsilon = (3.0/4.0)*d->optionalSimArgs[4]*pow(prims[ID(Prims::rho, i, j, k)],0.25);
          tau_Pi = tau_epsilon/3.0;
          beta_epsilon = tau_q*(d->gamma -1) - (d->gamma -1)*(kappa*prims[ID(Prims::rho, i, j, k)]*
            (prims[ID(Prims::rho, i, j, k)]+prims[ID(Prims::p, i, j, k)]))/(prims[ID(Prims::n, i, j, k)]**2 * prims[ID(Prims::p, i, j, k)]);

          aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + dxux + dyuy + dzuz;

          aux[ID(Aux::A, i, j, k)] = -tau_epsilon * ( aux[ID(Aux::W, i, j, k)]*(aux[ID(Aux::drhodt, i, j, k)] + 
          prims[ID(Prims::v1, i, j, k)]*dxrho + 
          prims[ID(Prims::v2, i, j, k)]*dyrho + 
          prims[ID(Prims::v3, i, j, k)]*dzrho ) +
          (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*aux[ID(Aux::Theta, i, j, k)]  );

          aux[ID(Aux::Pi, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)] + (tau_Pi/tau_epsilon)*aux[ID(Aux::A, i, j, k)];

          // if (i==4 && j ==0 && k==0)
          //   std::cout << zeta << "\t" << aux[ID(Aux::Pi, i, j, k)] << "\t" << aux[ID(Aux::A, i, j, k)] << "\t" << aux[ID(Aux::Theta, i, j, k)] << "\t" << std::endl;

          // beta_n = -tau_q*(d->gamma - 1) - kappa*pow(aux[ID(Aux::T, i, j, k)],4)*(prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) 
          //   / (prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]);

          beta_n = -tau_q*(d->gamma - 1) 
            + (kappa*(prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)])/(pow(prims[ID(Prims::n, i, j, k)],3)*prims[ID(Prims::p, i, j, k)]))
            * ((d->gamma - 1)*pow(prims[ID(Prims::rho, i, j, k)],2) + pow(prims[ID(Prims::p, i, j, k)],2));

          aux[ID(Aux::q1, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v1, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxux + 
            prims[ID(Prims::v2, i, j, k)]*dyux + 
            prims[ID(Prims::v3, i, j, k)]*dzux )
            
            + beta_epsilon*( dxrho 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( aux[ID(Aux::drhodt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxrho 
            + prims[ID(Prims::v2, i, j, k)]*dyrho
            + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            + beta_n*( (prims[ID(Prims::n, i+1, j, k)] - prims[ID(Prims::n, i-1, j, k)])/(2*d->dx) 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( aux[ID(Aux::dndt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxn 
            + prims[ID(Prims::v2, i, j, k)]*dyn
            + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


          aux[ID(Aux::q2, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv2dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v2, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxuy + 
            prims[ID(Prims::v2, i, j, k)]*dyuy + 
            prims[ID(Prims::v3, i, j, k)]*dzuy )
            
            + beta_epsilon*( dyrho 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( aux[ID(Aux::drhodt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxrho 
            + prims[ID(Prims::v2, i, j, k)]*dyrho
            + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            + beta_n*( dyn
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( aux[ID(Aux::dndt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxn 
            + prims[ID(Prims::v2, i, j, k)]*dyn
            + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


          aux[ID(Aux::q3, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv3dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v3, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxuz + 
            prims[ID(Prims::v2, i, j, k)]*dyuz + 
            prims[ID(Prims::v3, i, j, k)]*dzuz )
            
            + beta_epsilon*( dzrho 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( aux[ID(Aux::drhodt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxrho 
            + prims[ID(Prims::v2, i, j, k)]*dyrho
            + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            + beta_n*( dzn 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( aux[ID(Aux::dndt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxn 
            + prims[ID(Prims::v2, i, j, k)]*dyn
            + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;

          // These expressions need expanding and improving from HP!
          aux[ID(Aux::pi11, i, j, k)] = -2*eta*( 2*dxux
            - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi12, i, j, k)] = -2*eta*( dxuy
            + dyux
            - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi13, i, j, k)] = -2*eta*( dxuz
            + dzux
            - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi22, i, j, k)] = -2*eta*( 2*dyuy
            - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi23, i, j, k)] = -2*eta*( dyuz
            + dzuy
            - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi33, i, j, k)] = -2*eta*( dzuz
            - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          // Again, repeating this here once the correct values for v1,v2,v3 have been set...
          aux[ID(Aux::qv, i, j, k)] = (aux[ID(Aux::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (aux[ID(Aux::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                                  + (aux[ID(Aux::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
          aux[ID(Aux::pi00, i, j, k)] = aux[ID(Aux::pi11, i, j, k)] + aux[ID(Aux::pi22, i, j, k)] + aux[ID(Aux::pi33, i, j, k)];
          aux[ID(Aux::pi01, i, j, k)] = aux[ID(Aux::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                    + aux[ID(Aux::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                    + aux[ID(Aux::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi02, i, j, k)] = aux[ID(Aux::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                    + aux[ID(Aux::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                    + aux[ID(Aux::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi03, i, j, k)] = aux[ID(Aux::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                    + aux[ID(Aux::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                    + aux[ID(Aux::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation

          // Spatially-calculated time derivatives!
          aux[ID(Aux::dtp, i, j, k)]  =  -d->gamma*(prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))*pow(aux[ID(Aux::W, i, j, k)],2)/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
          aux[ID(Aux::dtrho, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)))/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
          aux[ID(Aux::dtv1, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v1, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
          aux[ID(Aux::dtv2, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v2, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
          aux[ID(Aux::dtv3, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v3, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) + d->gamma*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));

          // Copy into the variables actually used to do the advancement
          aux[ID(Aux::dpdt, i, j, k)]  = aux[ID(Aux::dtp, i, j, k)];
          aux[ID(Aux::drhodt, i, j, k)]  = aux[ID(Aux::dtrho, i, j, k)];
          aux[ID(Aux::dv1dt, i, j, k)]  = aux[ID(Aux::dtv1, i, j, k)];
          aux[ID(Aux::dv2dt, i, j, k)]  = aux[ID(Aux::dtv2, i, j, k)];
          aux[ID(Aux::dv3dt, i, j, k)]  = aux[ID(Aux::dtv3, i, j, k)];

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
        aux[ID(Aux::vsqrd, i, j, k)] = prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - aux[ID(Aux::vsqrd, i, j, k)] );
        prev_vars[ID(0, i, j, k)] = aux[ID(Aux::W, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        prev_vars[ID(1, i, j, k)] = prims[ID(Prims::v1, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        prev_vars[ID(2, i, j, k)] = prims[ID(Prims::v2, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        prev_vars[ID(3, i, j, k)] = prims[ID(Prims::v3, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        prev_vars[ID(4, i, j, k)] = prims[ID(Prims::p, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        prev_vars[ID(5, i, j, k)] = prims[ID(Prims::rho, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        prev_vars[ID(6, i, j, k)] = prims[ID(Prims::n, i, j, k)]; // Set here for time-differencing - necessary? what about the others?
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        prims[ID(Prims::rho, i, j, k)] = prims[ID(Prims::n, i, j, k)]*(1+aux[ID(Aux::e, i, j, k)]);
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
        aux[ID(Aux::h, i, j, k)] = 1 + aux[ID(Aux::e, i, j, k)] + prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
      }
    }
  }

  // double kappa = this->data->optionalSimArgs[0];
  // double tau_q = this->data->optionalSimArgs[1];
  // double zeta = this->data->optionalSimArgs[2];
  // double tau_Pi = this->data->optionalSimArgs[3];
  // double eta = this->data->optionalSimArgs[4];
  // double tau_epsilon = this->data->optionalSimArgs[5];

  double kappa { d->optionalSimArgs[0] };
  double tau_q;
  double zeta { d->optionalSimArgs[2] };
  double eta;
  double tau_epsilon;
  double tau_Pi;

  double beta_epsilon;
  double beta_n;  

  // Addition BDNK variables

  // for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
  //   for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
  //     for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

  double dxrho;
  double dyrho;
  double dzrho;

  double dxn;
  double dyn;
  double dzn;

  double dxux;
  double dyuy;
  double dzuz;

  double dxuy;
  double dxuz;
  double dyux;
  double dyuz;
  double dzux;
  double dzuy;

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
 
          dxrho = minmodGradFO(prims[ID(Prims::rho, i-1, j, k)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i+1, j, k)], d->dx);
          dyrho = minmodGradFO(prims[ID(Prims::rho, i, j-1, k)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i, j+1, k)], d->dy);
          dzrho = minmodGradFO(prims[ID(Prims::rho, i, j, k-1)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i, j, k+1)], d->dz);
          
          dxn = minmodGradFO(prims[ID(Prims::n, i-1, j, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::n, i+1, j, k)], d->dx);
          dyn = minmodGradFO(prims[ID(Prims::n, i, j-1, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::n, i, j+1, k)], d->dy);
          dzn = minmodGradFO(prims[ID(Prims::n, i, j, k-1)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::n, i, j, k+1)], d->dz);

          dxux = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                              aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)], d->dx);
          dyuy = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v2, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                              aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)], d->dy);
          dzuz = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                              aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);                          

          dxuy = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                              aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)], d->dx);
          dxuz = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                              aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)], d->dx);
        
          dyux = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                              aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)], d->dy);
          dyuz = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                              aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)], d->dy);

          dzux = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                              aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)], d->dz);  
          dzuy = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                              aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);  

          tau_q = (3.0/4.0)*d->optionalSimArgs[1]*pow(prims[ID(Prims::rho, i, j, k)],0.25);
          eta = d->optionalSimArgs[3]*pow(prims[ID(Prims::rho, i, j, k)],0.25);
          tau_epsilon = (3.0/4.0)*d->optionalSimArgs[4]*pow(prims[ID(Prims::rho, i, j, k)],0.25);
          tau_Pi = tau_epsilon/3.0;
          beta_epsilon = tau_q*(d->gamma -1) - (d->gamma -1)*(kappa*prims[ID(Prims::rho, i, j, k)]*
            (prims[ID(Prims::rho, i, j, k)]+prims[ID(Prims::p, i, j, k)]))/(prims[ID(Prims::n, i, j, k)]**2 * prims[ID(Prims::p, i, j, k)]);

          aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + dxux + dyuy + dzuz;

          aux[ID(Aux::A, i, j, k)] = -tau_epsilon * ( aux[ID(Aux::W, i, j, k)]*(aux[ID(Aux::drhodt, i, j, k)] + 
          prims[ID(Prims::v1, i, j, k)]*dxrho + 
          prims[ID(Prims::v2, i, j, k)]*dyrho + 
          prims[ID(Prims::v3, i, j, k)]*dzrho ) +
          (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*aux[ID(Aux::Theta, i, j, k)]  );

          aux[ID(Aux::Pi, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)] + (tau_Pi/tau_epsilon)*aux[ID(Aux::A, i, j, k)];

          // if (i==4 && j ==0 && k==0)
          //   std::cout << zeta << "\t" << aux[ID(Aux::Pi, i, j, k)] << "\t" << aux[ID(Aux::A, i, j, k)] << "\t" << aux[ID(Aux::Theta, i, j, k)] << "\t" << std::endl;

          // beta_n = -tau_q*(d->gamma - 1) - kappa*pow(aux[ID(Aux::T, i, j, k)],4)*(prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) 
          //   / (prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]);

          beta_n = -tau_q*(d->gamma - 1) 
            + (kappa*(prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)])/(pow(prims[ID(Prims::n, i, j, k)],3)*prims[ID(Prims::p, i, j, k)]))
            * ((d->gamma - 1)*pow(prims[ID(Prims::rho, i, j, k)],2) + pow(prims[ID(Prims::p, i, j, k)],2));

          aux[ID(Aux::q1, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v1, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxux + 
            prims[ID(Prims::v2, i, j, k)]*dyux + 
            prims[ID(Prims::v3, i, j, k)]*dzux )
            
            + beta_epsilon*( dxrho 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( aux[ID(Aux::drhodt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxrho 
            + prims[ID(Prims::v2, i, j, k)]*dyrho
            + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            + beta_n*( (prims[ID(Prims::n, i+1, j, k)] - prims[ID(Prims::n, i-1, j, k)])/(2*d->dx) 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( aux[ID(Aux::dndt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxn 
            + prims[ID(Prims::v2, i, j, k)]*dyn
            + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


          aux[ID(Aux::q2, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv2dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v2, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxuy + 
            prims[ID(Prims::v2, i, j, k)]*dyuy + 
            prims[ID(Prims::v3, i, j, k)]*dzuy )
            
            + beta_epsilon*( dyrho 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( aux[ID(Aux::drhodt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxrho 
            + prims[ID(Prims::v2, i, j, k)]*dyrho
            + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            + beta_n*( dyn
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( aux[ID(Aux::dndt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxn 
            + prims[ID(Prims::v2, i, j, k)]*dyn
            + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


          aux[ID(Aux::q3, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv3dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v3, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxuz + 
            prims[ID(Prims::v2, i, j, k)]*dyuz + 
            prims[ID(Prims::v3, i, j, k)]*dzuz )
            
            + beta_epsilon*( dzrho 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( aux[ID(Aux::drhodt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxrho 
            + prims[ID(Prims::v2, i, j, k)]*dyrho
            + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            + beta_n*( dzn 
            + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( aux[ID(Aux::dndt, i, j, k)]
            + prims[ID(Prims::v1, i, j, k)]*dxn 
            + prims[ID(Prims::v2, i, j, k)]*dyn
            + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;

          // These expressions need expanding and improving from HP!
          aux[ID(Aux::pi11, i, j, k)] = -2*eta*( 2*dxux
            - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi12, i, j, k)] = -2*eta*( dxuy
            + dyux
            - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi13, i, j, k)] = -2*eta*( dxuz
            + dzux
            - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi22, i, j, k)] = -2*eta*( 2*dyuy
            - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi23, i, j, k)] = -2*eta*( dyuz
            + dzuy
            - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::pi33, i, j, k)] = -2*eta*( dzuz
            - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

          aux[ID(Aux::qv, i, j, k)] = (aux[ID(Aux::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (aux[ID(Aux::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                                  + (aux[ID(Aux::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
          aux[ID(Aux::pi00, i, j, k)] = aux[ID(Aux::pi11, i, j, k)] + aux[ID(Aux::pi22, i, j, k)] + aux[ID(Aux::pi33, i, j, k)];
          aux[ID(Aux::pi01, i, j, k)] = aux[ID(Aux::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                    + aux[ID(Aux::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                    + aux[ID(Aux::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi02, i, j, k)] = aux[ID(Aux::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                    + aux[ID(Aux::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                    + aux[ID(Aux::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi03, i, j, k)] = aux[ID(Aux::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                    + aux[ID(Aux::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                    + aux[ID(Aux::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation

      }
    }
  }  

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(Cons::D, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)];
        // S1,2,3
        cons[ID(Cons::S1, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)] 
          + (aux[ID(Aux::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];

        // std::cout << cons[ID(Cons::S1, i, j, k)] << std::endl;

        cons[ID(Cons::S2, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)] 
          + (aux[ID(Aux::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
        cons[ID(Cons::S3, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)] 
          + (aux[ID(Aux::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v3, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];
        // Tau
        cons[ID(Cons::Tau, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] 
        - (prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)] + prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]) 
        + 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi00, i, j, k)];

        aux[ID(Aux::dtp, i, j, k)]  =  -d->gamma*(prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))*pow(aux[ID(Aux::W, i, j, k)],2)/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
        aux[ID(Aux::dtrho, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)))/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
        aux[ID(Aux::dtv1, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v1, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
        aux[ID(Aux::dtv2, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v2, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
        aux[ID(Aux::dtv3, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v3, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) + d->gamma*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));

        // Copy into the variables actually used to do the advancement
        aux[ID(Aux::dpdt, i, j, k)]  = aux[ID(Aux::dtp, i, j, k)];
        aux[ID(Aux::drhodt, i, j, k)]  = aux[ID(Aux::dtrho, i, j, k)];
        aux[ID(Aux::dv1dt, i, j, k)]  = aux[ID(Aux::dtv1, i, j, k)];
        aux[ID(Aux::dv2dt, i, j, k)]  = aux[ID(Aux::dtv2, i, j, k)];
        aux[ID(Aux::dv3dt, i, j, k)]  = aux[ID(Aux::dtv3, i, j, k)];

        // if (i==3 && j==3 && k==3)
          // std::cout << cons[ID(Cons::S1, i, j, k)] << "\t" << prims[ID(Prims::v1, i, j, k)] << "\t" << aux[ID(Aux::qv, i, j, k)] << "\t" << "\t" << prims[ID(Prims::rho, i, j, k)]
          // << "\t" << aux[ID(Aux::W, i, j, k)] <<  "\t" << aux[ID(Aux::Pi, i, j, k)] <<  "\t" << aux[ID(Aux::pi11, i, j, k)] <<  "\t" << aux[ID(Aux::pi12, i, j, k)] << "\t" 
          // << aux[ID(Aux::pi13, i, j, k)] << "\t" <<  prims[ID(Prims::p, i, j, k)] << "\t" << aux[ID(Aux::A, i, j, k)] << "\t" << aux[ID(Aux::Theta, i, j, k)]
          // << "\t" << tau_Pi << "\t" << tau_epsilon << std::endl;
     
      }  
    }
  }


}

void IS::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
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
          f[ID(1+nvar, i, j, k)] = cons[ID(Cons::S1+nvar, i, j, k)]*prims[ID(dir, i, j, k)] + ( aux[ID(Aux::q1+dir, i, j, k)] * prims[ID(Prims::v1+nvar, i, j, k)]  
            - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+nvar, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)] ) * aux[ID(Aux::W, i, j, k)];
          // (p+Pi)delta_ij
          if (dir == nvar) {
            f[ID(1+nvar, i, j, k)] += (prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)]);
          }
        }
        //  pi^i_j  
        if (dir == 0) {
          for (int nvar(0); nvar < 3; nvar++) {
            f[ID(1+nvar, i, j, k)] += aux[ID(Aux::pi11+nvar, i, j, k)];
          }
        } else if (dir == 1) {
          f[ID(1, i, j, k)] += aux[ID(Aux::pi12, i, j, k)];
          f[ID(2, i, j, k)] += aux[ID(Aux::pi22, i, j, k)];
          f[ID(3, i, j, k)] += aux[ID(Aux::pi23, i, j, k)];
        } else if (dir == 2) {
          f[ID(1, i, j, k)] += aux[ID(Aux::pi13, i, j, k)];
          f[ID(2, i, j, k)] += aux[ID(Aux::pi23, i, j, k)];
          f[ID(3, i, j, k)] += aux[ID(Aux::pi33, i, j, k)];
        } else {
          throw std::runtime_error("Flux direction is not 0, 1 or 2");
        }

        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Cons::Tau, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (aux[ID(Aux::q1+dir, i, j, k)] - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)])*aux[ID(Aux::W, i, j, k)]
          + aux[ID(Aux::pi01+dir, i, j, k)];
      
        // if (i==4 && j==0 && k==0)
        //   std::cout << f[ID(1, i, j, k)] << "\t" << cons[ID(Cons::S1, i, j, k)] << "\t" << prims[ID(v1, i, j, k)] << "\t" << aux[ID(Aux::qv, i, j, k)] << "\t" 
        //   << aux[ID(Aux::W, i, j, k)] <<  "\t" << aux[ID(Aux::Pi, i, j, k)] <<  "\t" << aux[ID(Aux::pi11, i, j, k)] <<  "\t" << aux[ID(Aux::pi12, i, j, k)] << "\t" 
        //   << aux[ID(Aux::pi13, i, j, k)] << "\t" <<  prims[ID(Prims::p, i, j, k)] << "\t" << aux[ID(Aux::A, i, j, k)] << "\t" << aux[ID(Aux::Theta, i, j, k)] << std::endl;
      
      } // End k loop
    } // End j loop
  } // End i loop
}

