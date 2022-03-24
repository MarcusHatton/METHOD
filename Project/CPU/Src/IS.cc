#include "IS.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

IS::IS() : Model()
{
  this->Ncons = 10;
  this->Nprims = 10;
  this->Naux = 34;
}

IS::IS(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 10;
  this->Nprims = (this->data)->Nprims = 10;
  this->Naux = (this->data)->Naux = 34;

  // Solutions for C2P all cells
  solution = (double *) malloc(sizeof(double)*4*data->Nx*data->Ny*data->Nz);

  smartGuesses = 0;
  
  // 0  
  this->data->consLabels.push_back("D");
  this->data->consLabels.push_back("S1");  this->data->consLabels.push_back("S2");  
  this->data->consLabels.push_back("S3");  this->data->consLabels.push_back("Tau");
  // 5
  this->data->consLabels.push_back("v1_C");  this->data->consLabels.push_back("v2_C");
  this->data->consLabels.push_back("v3_C");  this->data->consLabels.push_back("p_C");
  this->data->consLabels.push_back("rho_C");
  
  // 0
  this->data->primsLabels.push_back("v1");   this->data->primsLabels.push_back("v2");
  this->data->primsLabels.push_back("v3");
  // 3
  this->data->primsLabels.push_back("p");   this->data->primsLabels.push_back("rho");
  // 5
  this->data->primsLabels.push_back("dv1dt");
  this->data->primsLabels.push_back("dv2dt");  this->data->primsLabels.push_back("dv3dt");
  this->data->primsLabels.push_back("dpdt");   this->data->primsLabels.push_back("drhodt");
  // 10

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
  this->data->auxLabels.push_back("n");
  // 19
  this->data->auxLabels.push_back("h");      this->data->auxLabels.push_back("T");
  this->data->auxLabels.push_back("e");      this->data->auxLabels.push_back("W");
  // 23
  this->data->auxLabels.push_back("dndt"); this->data->auxLabels.push_back("dWdt");     
  // 25
  this->data->auxLabels.push_back("S1_PF");  this->data->auxLabels.push_back("S2_PF");  
  this->data->auxLabels.push_back("S3_PF");  this->data->auxLabels.push_back("Tau_PF");
  // 29 - spatially calculated time derivatives for comparison purposes!
  this->data->auxLabels.push_back("dtv1");  this->data->auxLabels.push_back("dtv2");
  this->data->auxLabels.push_back("dtv3");
  this->data->auxLabels.push_back("dtp");  this->data->auxLabels.push_back("dtrho");  

  // 34

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
  for(int ncon(0); ncon < Ncons; ncon++) {
    source[ncon] = 0.0;
  }
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

  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int k(0); k < this->data->Nz; k++) {
        for(int ncon(0); ncon < 4; ncon++) {
          source[ID(ncon, i, j, k)] = 0;
        }
        for(int ncon(0); ncon < 5; ncon++) {
          source[ID(Cons::v1_C+ncon, i, j, k)] = prims[ID(Prims::dv1dt+ncon, i, j, k)];
        }
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
  
  // Chain rule on W to get dWdt
  double dWdt_rf = pow(args->W_rf,3)*args->v1_rf*x[0] + pow(args->W_rf,3)*args->v2_rf*x[1]
                + pow(args->W_rf,3)*args->v3_rf*x[2]; 
  double Theta_rf = dWdt_rf + args->dxux_rf + args->dyuy_rf + args->dzuz_rf;
  double dndt_rf = -( args->v1_rf*args->dndx_rf + args->v2_rf*args->dndy_rf + args->v3_rf*args->dndz_rf 
                   + args->n_rf*Theta_rf/args->W_rf );
  double drhodt_rf = dndt_rf + x[3]/(args->gamma-1);

  double E_rf = args->Tau_NI_rf + args->n_rf*args->W_rf; // Replacement for D
  double A_rf = -args->tau_epsilon_rf * ( args->W_rf*(drhodt_rf + 
                args->v1_rf*args->dxrho_rf + args->v2_rf*args->dyrho_rf + args->v3_rf*args->dzrho_rf ) +
                (args->p_rf + args->rho_rf)*Theta_rf  );
  double Pi_rf = -args->zeta_rf*Theta_rf + (args->tau_Pi_rf/args->tau_epsilon_rf)*A_rf;
  double q1_rf = -args->tau_q_rf * (args->rho_rf + args->p_rf) * 
                  args->W_rf*( (args->W_rf*x[0] + dWdt_rf*args->v1_rf) + // chain rule
                  args->v1_rf*args->dxux_rf + args->v2_rf*args->dyux_rf + args->v3_rf*args->dzux_rf );
  double q2_rf = -args->tau_q_rf * (args->rho_rf + args->p_rf) * 
                  args->W_rf*( (args->W_rf*x[1] + dWdt_rf*args->v2_rf) + // chain rule
                  args->v1_rf*args->dxuy_rf + args->v2_rf*args->dyuy_rf + args->v3_rf*args->dzuy_rf );
  double q3_rf = -args->tau_q_rf * (args->rho_rf + args->p_rf) * 
                  args->W_rf*( (args->W_rf*x[2] + dWdt_rf*args->v3_rf) + // chain rule
                  args->v1_rf*args->dxuz_rf + args->v2_rf*args->dyuz_rf + args->v3_rf*args->dzuz_rf );
  
  // These expressions still need improvement!
  double pi11_rf = -2*args->eta_rf*( 2*args->dxux_rf - (2/3)*(1 + (args->W_rf*args->v1_rf)*(args->W_rf*args->v1_rf))*Theta_rf );
  double pi12_rf = -2*args->eta_rf*( args->dxuy_rf + args->dyux_rf - (2/3)*(args->W_rf*args->v1_rf)*(args->W_rf*args->v2_rf)*Theta_rf );
  double pi13_rf = -2*args->eta_rf*( args->dxuz_rf + args->dzux_rf - (2/3)*(args->W_rf*args->v1_rf)*(args->W_rf*args->v3_rf)*Theta_rf );
  double pi22_rf = -2*args->eta_rf*( 2*args->dyuy_rf - (2/3)*(1 + (args->W_rf*args->v2_rf)*(args->W_rf*args->v2_rf))*Theta_rf );
  double pi23_rf = -2*args->eta_rf*( args->dyuz_rf + args->dzuy_rf - (2/3)*(args->W_rf*args->v2_rf)*(args->W_rf*args->v3_rf)*Theta_rf );
  double pi33_rf = -2*args->eta_rf*( 2*args->dzuz_rf - (2/3)*(1 + (args->W_rf*args->v3_rf)*(args->W_rf*args->v3_rf))*Theta_rf );

  // Values should be sensible    
  if (args->p_rf < 0 || args->rho_rf < 0 || args->W_rf < 1 || args->n_rf < 0 || abs(args->v1_rf) >= 1 || abs(args->v2_rf) >= 1 || abs(args->v3_rf) >= 1 ) {
    printf("EEK");
    fvec[0] = fvec[1] = fvec[2] = fvec[3] = 1e6;
    return 0;
  }
  
  double pi00_rf = pi11_rf + pi22_rf + pi33_rf;
  double qv_rf = q1_rf*args->v1_rf + q2_rf*args->v2_rf + q3_rf*args->v3_rf;
  double pi01_rf = pi11_rf*args->v1_rf + pi12_rf*args->v2_rf + pi13_rf*args->v3_rf; // dbl check sign on orthogonality relation
  double pi02_rf = pi12_rf*args->v1_rf + pi22_rf*args->v2_rf + pi23_rf*args->v3_rf;
  double pi03_rf = pi13_rf*args->v1_rf + pi23_rf*args->v2_rf + pi33_rf*args->v3_rf;

  fvec[0] = args->S1_NI_rf - ( (Pi_rf + A_rf)*args->W_rf*args->W_rf*args->v1_rf 
            + (q1_rf + qv_rf*args->v1_rf)*args->W_rf + pi01_rf );
  fvec[1] = args->S2_NI_rf - ( (Pi_rf + A_rf)*args->W_rf*args->W_rf*args->v2_rf 
            + (q2_rf + qv_rf*args->v2_rf)*args->W_rf + pi02_rf );
  fvec[2] = args->S3_NI_rf - ( (Pi_rf + A_rf)*args->W_rf*args->W_rf*args->v3_rf 
            + (q3_rf + qv_rf*args->v3_rf)*args->W_rf + pi03_rf );
  fvec[3] = E_rf - ( (Pi_rf + A_rf)*args->W_rf*args->W_rf - (Pi_rf + A_rf) 
            + 2*qv_rf*args->W_rf + pi00_rf );

  return 0;
}



void IS::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  printf("Warning: Not implemented (or indeed implementable...) for BDNKvPP");
  return;

  /*

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
  
  sol[0] = prims[Prims::p] + aux[Aux::Pi] - 2*aux[Aux::qv]*aux[Aux::W] - aux[Aux::pi00];
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
  */

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
  
  // for (int i(d->is); i < d->ie; i++) {
  //   for (int j(d->js); j < d->je; j++) {
  //     for (int k(d->ks); k < d->ke; k++) {

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

  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

        // Trivially recover basic variables
        for(int nbasic(0); nbasic < 5; nbasic++) {
          prims[ID(Prims::v1+nbasic, i, j, k)] = cons[ID(Cons::v1_C+nbasic, i, j, k)];
        }

        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - (prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
          + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)]
          + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)]) );
        aux[ID(Aux::n, i, j, k)] = cons[ID(Cons::D, i, j, k)]/aux[ID(Aux::W, i, j, k)];
        // aux[ID(Aux::n, i, j, k)] = prims[ID(Prims::rho, i, j, k)] - prims[ID(Prims::p, i, j, k)]/(d->gamma-1);
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (aux[ID(Aux::n, i, j, k)]*(d->gamma-1));
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / aux[ID(Aux::n, i, j, k)];
        aux[ID(Aux::h, i, j, k)] = 1 + aux[ID(Aux::e, i, j, k)] + prims[ID(Prims::p, i, j, k)] / aux[ID(Aux::n, i, j, k)];

      }
    }
  }          
           
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

          dxrho = minmodGradFO(prims[ID(Prims::rho, i-1, j, k)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i+1, j, k)], d->dx);
          dyrho = minmodGradFO(prims[ID(Prims::rho, i, j-1, k)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i, j+1, k)], d->dy);
          dzrho = minmodGradFO(prims[ID(Prims::rho, i, j, k-1)], prims[ID(Prims::rho, i, j, k)], prims[ID(Prims::rho, i, j, k+1)], d->dz);
          
          dxn = minmodGradFO(aux[ID(Aux::n, i-1, j, k)], aux[ID(Aux::n, i, j, k)], aux[ID(Aux::n, i+1, j, k)], d->dx);
          dyn = minmodGradFO(aux[ID(Aux::n, i, j-1, k)], aux[ID(Aux::n, i, j, k)], aux[ID(Aux::n, i, j+1, k)], d->dy);
          dzn = minmodGradFO(aux[ID(Aux::n, i, j, k-1)], aux[ID(Aux::n, i, j, k)], aux[ID(Aux::n, i, j, k+1)], d->dz);

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
          beta_epsilon = tau_q*(d->gamma -1);

          // Initial guesses - could replace these with spatial-calc'd
          sol[0] = prims[ID(Prims::dv1dt, i, j, k)];
          sol[1] = prims[ID(Prims::dv2dt, i, j, k)];
          sol[2] = prims[ID(Prims::dv3dt, i, j, k)];
          sol[3] = prims[ID(Prims::dpdt, i, j, k)];

          // Set additional args for rootfind
          args.v1_rf = prims[ID(Prims::v1, i, j, k)];
          args.v2_rf = prims[ID(Prims::v2, i, j, k)];
          args.v3_rf = prims[ID(Prims::v3, i, j, k)];
          args.W_rf = aux[ID(Aux::W, i, j, k)];
          args.p_rf = prims[ID(Prims::p, i, j, k)];
          args.rho_rf = prims[ID(Prims::rho, i, j, k)];
          args.S1_NI_rf = cons[ID(Cons::S1, i, j, k)]  - (args.rho_rf + args.p_rf)*args.W_rf*args.W_rf*args.v1_rf;
          args.S2_NI_rf = cons[ID(Cons::S2, i, j, k)]  - (args.rho_rf + args.p_rf)*args.W_rf*args.W_rf*args.v2_rf;
          args.S3_NI_rf = cons[ID(Cons::S3, i, j, k)]  - (args.rho_rf + args.p_rf)*args.W_rf*args.W_rf*args.v3_rf;
          args.Tau_NI_rf = cons[ID(Cons::Tau, i, j, k)]  - ((args.rho_rf + args.p_rf)*args.W_rf*args.W_rf - args.p_rf);
          args.n_rf = aux[ID(Aux::n, i, j, k)];
          args.dndx_rf = dxn;
          args.dndy_rf = dyn;
          args.dndz_rf = dzn;
          args.dxrho_rf = dxrho;
          args.dyrho_rf = dyrho;
          args.dzrho_rf = dzrho;
          args.dxux_rf = dxux;
          args.dxuy_rf = dxuy;
          args.dxuz_rf = dxuz;
          args.dyux_rf = dyux;
          args.dyuy_rf = dyuy;
          args.dyuz_rf = dyuz;
          args.dzux_rf = dzux;
          args.dzuy_rf = dzuy;
          args.dzuz_rf = dzuz;
          args.tau_epsilon_rf = tau_epsilon;
          args.tau_Pi_rf = tau_Pi;
          args.tau_q_rf = tau_q;
          args.eta_rf = eta;
          args.zeta_rf = zeta;

          args.gamma = d->gamma;
  
          // Solve residual = 0
          info = __cminpack_func__(hybrd1) (&ISresidual, &args, sys_size, sol, res,
                                            tol, wa, lwa);
                                                                    
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
            //printf("(%f, %f, %f, %f, %f) prims\n",  prims[ID(Prims::p, i, j, k)], aux[ID(Aux::Pi, i, j, k)], aux[ID(Aux::n, i, j, k)], prims[ID(Prims::v1, i, j, k)], aux[ID(Aux::q1, i, j, k)]);
            //printf("(%f, %f, %f, %f) aux  \n",  aux[ID(Aux::W, i, j, k)], aux[ID(Aux::qv, i, j, k)], aux[ID(Aux::pi00, i, j, k)], aux[ID(Aux::pi01, i, j, k)]);
            //printf("(%f, %f, %f, %f) cons \n",  cons[ID(Cons::Y1, i, j, k)], cons[ID(Cons::D, i, j, k)], cons[ID(Cons::S1, i, j, k)], cons[ID(Cons::Tau, i, j, k)]);
            exit(0);
            Failed fail = {i, j, k};
            fails.push_back(fail);
          } else {
            // Now have the correct values for time derivs
            // solution[ID(0, i, j, k)] = sol[0];
            // solution[ID(1, i, j, k)] = sol[1];
            // solution[ID(2, i, j, k)] = sol[2];
            // solution[ID(3, i, j, k)] = sol[3];
            prims[ID(Prims::dv1dt, i, j, k)] = sol[0];
            prims[ID(Prims::dv2dt, i, j, k)] = sol[1];
            prims[ID(Prims::dv3dt, i, j, k)] = sol[2];
            prims[ID(Prims::dpdt, i, j, k)] = sol[3];
          }      
        
  //       } // End k-loop
  //     } // End j-loop
  //   } // End i-loop

  // // double kappa = d->optionalSimArgs[0];
  // // double tau_q = d->optionalSimArgs[1];
  // // double zeta = d->optionalSimArgs[2];
  // // double tau_Pi = d->optionalSimArgs[3];
  // // double eta = d->optionalSimArgs[4];
  // // double tau_epsilon = d->optionalSimArgs[5];


  // for (int i(d->is); i < d->ie; i++) {
  //   for (int j(d->js); j < d->je; j++) {
  //     for (int k(d->ks); k < d->ke; k++) {
        
        aux[ID(Aux::dWdt, i, j, k)] = pow(aux[ID(Aux::W, i, j, k)],3)*(prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::dv1dt, i, j, k)] + 
          prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::dv2dt, i, j, k)] + 
          prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::dv3dt, i, j, k)]);        

        aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + dxux + dyuy + dzuz;

        aux[ID(Aux::dndt, i, j, k)] = prims[ID(Prims::v1, i, j, k)]*dxn + prims[ID(Prims::v2, i, j, k)]*dyn
          + prims[ID(Prims::v3, i, j, k)]*dzn + aux[ID(Aux::n, i, j, k)]*aux[ID(Aux::Theta, i, j, k)]/aux[ID(Aux::W, i, j, k)];
      
        prims[ID(Prims::drhodt, i, j, k)] = aux[ID(Aux::dndt, i, j, k)] + prims[ID(Prims::dpdt, i, j, k)]/(d->gamma-1);

        aux[ID(Aux::A, i, j, k)] = -tau_epsilon * ( aux[ID(Aux::W, i, j, k)]*(prims[ID(Prims::drhodt, i, j, k)] + 
        prims[ID(Prims::v1, i, j, k)]*dxrho + 
        prims[ID(Prims::v2, i, j, k)]*dyrho + 
        prims[ID(Prims::v3, i, j, k)]*dzrho ) +
        (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*aux[ID(Aux::Theta, i, j, k)]  );

        aux[ID(Aux::Pi, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)] + (tau_Pi/tau_epsilon)*aux[ID(Aux::A, i, j, k)];

        beta_n = -tau_q*(d->gamma - 1) - kappa*aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::T, i, j, k)]
          *(prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) / (prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]);

        aux[ID(Aux::q1, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
          aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::dv1dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v1, i, j, k)]) + // chain rule
          prims[ID(Prims::v1, i, j, k)]*dxux + 
          prims[ID(Prims::v2, i, j, k)]*dyux + 
          prims[ID(Prims::v3, i, j, k)]*dzux );
          
          // + beta_epsilon*( dxrho 
          // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( -prims[ID(Prims::drhodt, i, j, k)]
          // + prims[ID(Prims::v1, i, j, k)]*dxrho 
          // + prims[ID(Prims::v2, i, j, k)]*dyrho
          // + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
          
          // + beta_n*( (aux[ID(Aux::n, i+1, j, k)] - aux[ID(Aux::n, i-1, j, k)])/(2*d->dx) 
          // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( -aux[ID(Aux::dndt, i, j, k)]
          // + prims[ID(Prims::v1, i, j, k)]*dxn 
          // + prims[ID(Prims::v2, i, j, k)]*dyn
          // + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


        aux[ID(Aux::q2, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
          aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::dv2dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v2, i, j, k)]) + // chain rule
          prims[ID(Prims::v1, i, j, k)]*dxuy + 
          prims[ID(Prims::v2, i, j, k)]*dyuy + 
          prims[ID(Prims::v3, i, j, k)]*dzuy );
          
          // + beta_epsilon*( dyrho 
          // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( -prims[ID(Prims::drhodt, i, j, k)]
          // + prims[ID(Prims::v1, i, j, k)]*dxrho 
          // + prims[ID(Prims::v2, i, j, k)]*dyrho
          // + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
          
          // + beta_n*( dyn
          // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( -aux[ID(Aux::dndt, i, j, k)]
          // + prims[ID(Prims::v1, i, j, k)]*dxn 
          // + prims[ID(Prims::v2, i, j, k)]*dyn
          // + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


        aux[ID(Aux::q3, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
          aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::dv3dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v3, i, j, k)]) + // chain rule
          prims[ID(Prims::v1, i, j, k)]*dxuz + 
          prims[ID(Prims::v2, i, j, k)]*dyuz + 
          prims[ID(Prims::v3, i, j, k)]*dzuz );
          
          // + beta_epsilon*( dzrho 
          // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( -prims[ID(Prims::drhodt, i, j, k)]
          // + prims[ID(Prims::v1, i, j, k)]*dxrho 
          // + prims[ID(Prims::v2, i, j, k)]*dyrho
          // + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
          
          // + beta_n*( dzn 
          // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( -aux[ID(Aux::dndt, i, j, k)]
          // + prims[ID(Prims::v1, i, j, k)]*dxn 
          // + prims[ID(Prims::v2, i, j, k)]*dyn
          // + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;

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

        aux[ID(Aux::S1_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)];  
        aux[ID(Aux::S2_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)];  
        aux[ID(Aux::S3_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)];  
        aux[ID(Aux::Tau_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] - (prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]);

        // Just tracking here...
        aux[ID(Aux::dtp, i, j, k)]  =  -d->gamma*(prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))*pow(aux[ID(Aux::W, i, j, k)],2)/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
        aux[ID(Aux::dtrho, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)))/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
        aux[ID(Aux::dtv1, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v1, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
        aux[ID(Aux::dtv2, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v2, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
        aux[ID(Aux::dtv3, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v3, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) + d->gamma*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));

        } // End k-loop
      } // End j-loop
    } // End i-loop


  //if (i==4 && j ==0 && k==0)
  // std::cout << zeta << "\t" << aux[ID(Aux::Pi, 4, 0, 0)] << "\t" << aux[ID(Aux::A, 4, 0, 0)] << "\t" << aux[ID(Aux::Theta, 4, 0, 0)] << "\t" << std::endl;
  
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
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)] );
        aux[ID(Aux::n, i, j, k)] = prims[ID(Prims::rho, i, j, k)] - prims[ID(Prims::p, i, j, k)]/(d->gamma-1);                        
        // aux[ID(Aux::n, i, j, k)] = cons[ID(Cons::D, i, j, k)]/aux[ID(Aux::W, i, j, k)];
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (aux[ID(Aux::n, i, j, k)]*(d->gamma-1));
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / aux[ID(Aux::n, i, j, k)];
        aux[ID(Aux::h, i, j, k)] = 1 + aux[ID(Aux::e, i, j, k)] + prims[ID(Prims::p, i, j, k)] / aux[ID(Aux::n, i, j, k)];
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
          
          dxn = minmodGradFO(aux[ID(Aux::n, i-1, j, k)], aux[ID(Aux::n, i, j, k)], aux[ID(Aux::n, i+1, j, k)], d->dx);
          dyn = minmodGradFO(aux[ID(Aux::n, i, j-1, k)], aux[ID(Aux::n, i, j, k)], aux[ID(Aux::n, i, j+1, k)], d->dy);
          dzn = minmodGradFO(aux[ID(Aux::n, i, j, k-1)], aux[ID(Aux::n, i, j, k)], aux[ID(Aux::n, i, j, k+1)], d->dz);

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
          beta_epsilon = tau_q*(d->gamma -1);

          aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + dxux + dyuy + dzuz;

          aux[ID(Aux::A, i, j, k)] = -tau_epsilon * ( aux[ID(Aux::W, i, j, k)]*(prims[ID(Prims::drhodt, i, j, k)] + 
          prims[ID(Prims::v1, i, j, k)]*dxrho + 
          prims[ID(Prims::v2, i, j, k)]*dyrho + 
          prims[ID(Prims::v3, i, j, k)]*dzrho ) +
          (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*aux[ID(Aux::Theta, i, j, k)]  );

          aux[ID(Aux::Pi, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)] + (tau_Pi/tau_epsilon)*aux[ID(Aux::A, i, j, k)];

          // if (i==4 && j ==0 && k==0)
          //   std::cout << zeta << "\t" << aux[ID(Aux::Pi, i, j, k)] << "\t" << aux[ID(Aux::A, i, j, k)] << "\t" << aux[ID(Aux::Theta, i, j, k)] << "\t" << std::endl;

          beta_n = -tau_q*(d->gamma - 1) - kappa*aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::T, i, j, k)]
            *(prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) / (prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::p, i, j, k)]);

          aux[ID(Aux::q1, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::dv1dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v1, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxux + 
            prims[ID(Prims::v2, i, j, k)]*dyux + 
            prims[ID(Prims::v3, i, j, k)]*dzux );
            
            // + beta_epsilon*( dxrho 
            // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( -prims[ID(Prims::drhodt, i, j, k)]
            // + prims[ID(Prims::v1, i, j, k)]*dxrho 
            // + prims[ID(Prims::v2, i, j, k)]*dyrho
            // + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            // + beta_n*( (aux[ID(Aux::n, i+1, j, k)] - aux[ID(Aux::n, i-1, j, k)])/(2*d->dx) 
            // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*( -aux[ID(Aux::dndt, i, j, k)]
            // + prims[ID(Prims::v1, i, j, k)]*dxn 
            // + prims[ID(Prims::v2, i, j, k)]*dyn
            // + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


          aux[ID(Aux::q2, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::dv2dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v2, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxuy + 
            prims[ID(Prims::v2, i, j, k)]*dyuy + 
            prims[ID(Prims::v3, i, j, k)]*dzuy );
            
            // + beta_epsilon*( dyrho 
            // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( -prims[ID(Prims::drhodt, i, j, k)]
            // + prims[ID(Prims::v1, i, j, k)]*dxrho 
            // + prims[ID(Prims::v2, i, j, k)]*dyrho
            // + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            // + beta_n*( dyn
            // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*( -aux[ID(Aux::dndt, i, j, k)]
            // + prims[ID(Prims::v1, i, j, k)]*dxn 
            // + prims[ID(Prims::v2, i, j, k)]*dyn
            // + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;


          aux[ID(Aux::q3, i, j, k)] = -tau_q * (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * 
            aux[ID(Aux::W, i, j, k)]*( (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::dv3dt, i, j, k)] + aux[ID(Aux::dWdt, i, j, k)]*prims[ID(Prims::v3, i, j, k)]) + // chain rule
            prims[ID(Prims::v1, i, j, k)]*dxuz + 
            prims[ID(Prims::v2, i, j, k)]*dyuz + 
            prims[ID(Prims::v3, i, j, k)]*dzuz );
            
            // + beta_epsilon*( dzrho 
            // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( -prims[ID(Prims::drhodt, i, j, k)]
            // + prims[ID(Prims::v1, i, j, k)]*dxrho 
            // + prims[ID(Prims::v2, i, j, k)]*dyrho
            // + prims[ID(Prims::v3, i, j, k)]*dzrho ) ) 
            
            // + beta_n*( dzn 
            // + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*( -aux[ID(Aux::dndt, i, j, k)]
            // + prims[ID(Prims::v1, i, j, k)]*dxn 
            // + prims[ID(Prims::v2, i, j, k)]*dyn
            // + prims[ID(Prims::v3, i, j, k)]*dzn ) ) ;

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
        cons[ID(Cons::D, i, j, k)] = aux[ID(Aux::n, i, j, k)]*aux[ID(Aux::W, i, j, k)];
        // S1,2,3
        cons[ID(Cons::S1, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)] 
          + (aux[ID(Aux::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
        cons[ID(Cons::S2, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)] 
          + (aux[ID(Aux::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
        cons[ID(Cons::S3, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)] 
          + (aux[ID(Aux::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v3, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];
        // Tau
        cons[ID(Cons::Tau, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] 
        - (prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::Pi, i, j, k)] + aux[ID(Aux::A, i, j, k)] + aux[ID(Aux::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]) 
        + 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi00, i, j, k)];
     
        cons[ID(Cons::v1_C, i, j, k)] = prims[ID(Prims::v1, i, j, k)];
        cons[ID(Cons::v2_C, i, j, k)] = prims[ID(Prims::v2, i, j, k)];
        cons[ID(Cons::v3_C, i, j, k)] = prims[ID(Prims::v3, i, j, k)];
        cons[ID(Cons::p_C, i, j, k)] = prims[ID(Prims::p, i, j, k)];
        cons[ID(Cons::rho_C, i, j, k)] = prims[ID(Prims::rho, i, j, k)];

        aux[ID(Aux::S1_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)];  
        aux[ID(Aux::S2_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)];  
        aux[ID(Aux::S3_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)];  
        aux[ID(Aux::Tau_PF, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] - (prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]);

        // Spatially-calc'd time derivs
        aux[ID(Aux::dtp, i, j, k)]  =  -d->gamma*(prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))*pow(aux[ID(Aux::W, i, j, k)],2)/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
        aux[ID(Aux::dtrho, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) + d->gamma*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) + d->gamma*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::p, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::rho, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::rho, i+1, j, k)] - prims[ID(Prims::rho, i-1, j, k)])/(2*d->dx)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::rho, i, j+1, k)] - prims[ID(Prims::rho, i, j-1, k)])/(2*d->dy)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::rho, i, j, k+1)] - prims[ID(Prims::rho, i, j, k-1)])/(2*d->dz)))/(d->gamma*pow(aux[ID(Aux::W, i, j, k)],2) - d->gamma + 1);
        aux[ID(Aux::dtv1, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v1, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
        aux[ID(Aux::dtv2, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v2, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));
        aux[ID(Aux::dtv3, i, j, k)]  =  (-d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],4)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)) + 2*d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i+1, j, k)] - prims[ID(Prims::p, i-1, j, k)])/(2*d->dx)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::p, i, j+1, k)] - prims[ID(Prims::p, i, j-1, k)])/(2*d->dy)) + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*pow(prims[ID(Prims::v3, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) + d->gamma*((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)) - pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz)) - ((prims[ID(Prims::p, i, j, k+1)] - prims[ID(Prims::p, i, j, k-1)])/(2*d->dz)))/((d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::p, i, j, k)] + d->gamma*pow(aux[ID(Aux::W, i, j, k)],2)*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*pow(aux[ID(Aux::W, i, j, k)],2));

        // Copy them to do the first step's advancement (initial data)
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
        f[ID(0, i, j, k)] = cons[ID(Cons::D, i, j, k)]*prims[ID(dir, i, j, k)];;
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
          f[ID(0, i, j, k)] += aux[ID(Aux::pi12, i, j, k)];
          f[ID(1, i, j, k)] += aux[ID(Aux::pi22, i, j, k)];
          f[ID(2, i, j, k)] += aux[ID(Aux::pi23, i, j, k)];
        } else if (dir == 2) {
          f[ID(0, i, j, k)] += aux[ID(Aux::pi13, i, j, k)];
          f[ID(1, i, j, k)] += aux[ID(Aux::pi23, i, j, k)];
          f[ID(2, i, j, k)] += aux[ID(Aux::pi33, i, j, k)];
        } else {
          throw std::runtime_error("Flux direction is not 0, 1 or 2");
        }

        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Cons::Tau, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (aux[ID(Aux::q1+dir, i, j, k)] - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)])*aux[ID(Aux::W, i, j, k)]
          + aux[ID(Aux::pi01+dir, i, j, k)];
        
        // v1, v2, v3, p, rho
        for (int nvar(0); nvar < 5; nvar++) {
          f[ID(5+nvar, i, j, k)] = 0;
        }

      } // End k loop
    } // End j loop
  } // End i loop
}

