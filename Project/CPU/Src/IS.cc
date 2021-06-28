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

  prev_vars = (double *) malloc(sizeof(double)*1*data->Nx*data->Ny*data->Nz);

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
  source[5] = (prims[Prims::n] / tau_q) * (aux[Aux::q1NS] - prims[Prims::q1]);
  source[6] = (prims[Prims::n] / tau_q) * (aux[Aux::q2NS] - prims[Prims::q2]);
  source[7] = (prims[Prims::n] / tau_q) * (aux[Aux::q3NS] - prims[Prims::q3]);
  // U
  source[8] = (prims[Prims::n] / tau_Pi) * (aux[Aux::PiNS] - prims[Prims::Pi]);
  // Z11,12,13,22,23,33
  source[9] = (prims[Prims::n] / tau_pi) * (aux[Aux::pi11NS] - prims[Prims::pi11]);
  source[10] = (prims[Prims::n] / tau_pi) * (aux[Aux::pi12NS] - prims[Prims::pi12]);
  source[11] = (prims[Prims::n] / tau_pi) * (aux[Aux::pi13NS] - prims[Prims::pi13]);
  source[12] = (prims[Prims::n] / tau_pi) * (aux[Aux::pi22NS] - prims[Prims::pi22]);
  source[13] = (prims[Prims::n] / tau_pi) * (aux[Aux::pi23NS] - prims[Prims::pi23]);
  source[14] = (prims[Prims::n] / tau_pi) * (aux[Aux::pi33NS] - prims[Prims::pi33]);
  
}

void IS::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  double kappa = this->data->optionalSimArgs[0];
  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
  double tau_pi = this->data->optionalSimArgs[5];
        
  // q_j,NS 10
  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

        aux[ID(Aux::q1NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx) );
        aux[ID(Aux::q2NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy) );
        aux[ID(Aux::q3NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz) );
        
        //printf("(%i, %i, %i) ijk\n", i, j, k);
        //printf("(%f, %f, %f) \n", aux[ID(Aux::T, i+1, j, k)], aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::q1NS, i, j, k)]);
        //printf("(%f, %f, %f, %f) \n", prims[ID(Prims::p, i+1, j, k)], prims[ID(Prims::n, i+1, j, k)], prims[ID(Prims::p, i-1, j, k)], prims[ID(Prims::n, i-1, j, k)]);
        
//          aux[ID(Aux::q1NS, i, j, k)] = -kappa*aux[ID(Aux::T, i, j, k)] * ( (log(aux[ID(Aux::T, i+1, j, k)]) - log(aux[ID(Aux::T, i-1, j, k)]))/(2*d->dx) );
//          aux[ID(Aux::q2NS, i, j, k)] = -kappa*aux[ID(Aux::T, i, j, k)] * ( (log(aux[ID(Aux::T, i, j+1, k)]) - log(aux[ID(Aux::T, i, j-1, k)]))/(2*d->dy) );
//          aux[ID(Aux::q3NS, i, j, k)] = -kappa*aux[ID(Aux::T, i, j, k)] * ( (log(aux[ID(Aux::T, i, j, k+1)]) - log(aux[ID(Aux::T, i, j, k-1)]))/(2*d->dz) );
 
        // Theta 20 then Pi,NS 13 
        aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)] - aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v3, i, j, k+1)] - aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz);
        // Pi,NS = -zeta*Theta
        aux[ID(Aux::PiNS, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)];
  
        // pi^l_j,NS 14
        // 11
        aux[ID(Aux::pi11NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 12
        aux[ID(Aux::pi12NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 13
        aux[ID(Aux::pi13NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 22
        aux[ID(Aux::pi22NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 23
        aux[ID(Aux::pi23NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 33
        aux[ID(Aux::pi33NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
      }
    }
  }   

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
        source[ID(Y1, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_q) * (aux[ID(Aux::q1NS, i, j, k)] - prims[ID(Prims::q1, i, j, k)]);
        source[ID(Y2, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_q) * (aux[ID(Aux::q1NS, i, j, k)] - prims[ID(Prims::q2, i, j, k)]);
        source[ID(Y3, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_q) * (aux[ID(Aux::q1NS, i, j, k)] - prims[ID(Prims::q3, i, j, k)]);        
        // U
        source[ID(U, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_Pi) * (aux[ID(Aux::PiNS, i, j, k)] - prims[ID(Prims::Pi, i, j, k)]);
        // Z11,12,13,22,23,33
        source[ID(Z11, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_pi) * (aux[ID(Aux::pi11NS, i, j, k)] - prims[ID(Prims::pi11, i, j, k)]);
        source[ID(Z12, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_pi) * (aux[ID(Aux::pi12NS, i, j, k)] - prims[ID(Prims::pi12, i, j, k)]);
        source[ID(Z13, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_pi) * (aux[ID(Aux::pi13NS, i, j, k)] - prims[ID(Prims::pi13, i, j, k)]);
        source[ID(Z22, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_pi) * (aux[ID(Aux::pi22NS, i, j, k)] - prims[ID(Prims::pi22, i, j, k)]); 
        source[ID(Z23, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_pi) * (aux[ID(Aux::pi23NS, i, j, k)] - prims[ID(Prims::pi23, i, j, k)]);
        source[ID(Z33, i, j, k)] = (prims[ID(Prims::n, i, j, k)] / tau_pi) * (aux[ID(Aux::pi33NS, i, j, k)] - prims[ID(Prims::pi33, i, j, k)]);                       
        
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
  double rho_plus_p_rf = ((E_rf + x[0])/(W_rf*W_rf)) - args->Pi_rf;
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

  fvec[0] = p_rf + args->Pi_rf - 2*qv_rf*W_rf - pi00_rf - x[0];
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
  double p_rf = args->D_rf*(1/x[0] -1) - args->Pi_rf + 2*qv_rf*W_rf + pi00_rf - args->Tau_rf;
  double rho_rf = n_rf + p_rf/(args->gamma-1);
  double H_rf = 1 + (p_rf*(args->gamma/(args->gamma-1)) + args->Pi_rf)/n_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 0 || v1_rf >= 1 || v2_rf >= 1 || v3_rf >= 1) {
    printf("EEK");
    fvec[0] = fvec[1] = 1e6;
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

  // Y1-3,U,Z11-33
  for (int ndissvar(0); ndissvar < 10; ndissvar++) {
    prims[Prims::q1+ndissvar] = cons[Cons::Y1+ndissvar] / cons[Cons::D];
  }

  aux[Aux::pi00] = prims[Prims::pi11] + prims[Prims::pi22] + prims[Prims::pi33]; // this one can be done here fine
  // what about these? need them in the guesses...
  aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
  aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation

  // Hybrd1 set-up
  Args args;                      // Additional arguments structure
  const int sys_size(4);                     // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.4e-7;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Set additional args for rootfind
  args.D_rf = cons[Cons::D];
  args.S1_rf = cons[Cons::S1];
  args.S2_rf = cons[Cons::S2];
  args.S3_rf = cons[Cons::S3];
  args.Tau_rf = cons[Cons::Tau];
  args.q1_rf = prims[Prims::q1];
  args.q2_rf = prims[Prims::q2];
  args.q3_rf = prims[Prims::q3];
  args.Pi_rf = prims[Prims::Pi];
  args.pi11_rf = prims[Prims::pi11];
  args.pi12_rf = prims[Prims::pi12];
  args.pi13_rf = prims[Prims::pi13];
  args.pi22_rf = prims[Prims::pi22];
  args.pi23_rf = prims[Prims::pi23];
  args.pi33_rf = prims[Prims::pi33];
  args.gamma = d->gamma;
  
  bool alternative_C2P = true;
  
  if (alternative_C2P) {
  
  sol[0] = 1/(aux[Aux::W]*(1 + (prims[Prims::p]*(d->gamma/(d->gamma-1)) + prims[Prims::Pi])/prims[Prims::n]));
  sol[1] = (prims[Prims::q1] + aux[Aux::qv]*prims[Prims::v1])*aux[Aux::W] + aux[Aux::pi01];
  sol[2] = (prims[Prims::q2] + aux[Aux::qv]*prims[Prims::v2])*aux[Aux::W] + aux[Aux::pi02];
  sol[3] = (prims[Prims::q3] + aux[Aux::qv]*prims[Prims::v3])*aux[Aux::W] + aux[Aux::pi03];

  // Solve residual = 0
  info = __cminpack_func__(hybrd1) (&ISAlternativeResidual, &args, sys_size, sol, res,
                                    tol, wa, lwa);
  // If root find fails, add failed cell to the list
  if (info!=1) {
    //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
    throw std::runtime_error("C2P could not converge.\n");
  }
  aux[Aux::vsqrd] = sol[0]*sol[0]*((args->S1_rf - sol[1])*(args->S1_rf - sol[1]) + (args->S2_rf - sol[2])*(args->S2_rf - sol[2]) 
                    + (args->S3_rf - sol[3])*(args->S3_rf - sol[3]))/(cons[Cons::D])*(cons[Cons::D]);
  aux[Aux::W] = (1 / sqrt(1 - aux[Aux::vsqrd]));
  prims[Prims::n] = cons[Cons::D] / aux[Aux::W];
  prims[Prims::v1] = sol[0]*(cons[Cons::S1] - x[1])/cons[Cons::D];
  prims[Prims::v2] = sol[0]*(cons[Cons::S2] - x[2])/cons[Cons::D];
  prims[Prims::v3] = sol[0]*(cons[Cons::S3] - x[3])/cons[Cons::D];
  aux[Aux::pi00] = prims[Prims::pi11] + prims[Prims::pi22] + prims[Prims::pi33]; // not sure we need this again here tbh
  aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
  prims[Prims::p] = cons[Cons::D]*(1/sol[0] -1) - prims[Prims::Pi] + 2*aux[Aux::qv]*aux[Aux::W] + aux[Aux::pi00] - cons[Cons::Tau];
  prims[Prims::rho] = prims[Prims::n] + prims[Prims::p]/(d->gamma-1);

  
  // Repeating the ones here that depend on v1,v2,v3...
  aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
        
  aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
  aux[Aux::T] = prims[Prims::p] / prims[Prims::n]; 
    
  
  } else {
  
  sol[0] = prims[Prims::p] + prims[Prims::Pi] - 2*aux[Aux::qv]*aux[Aux::W] - aux[Aux::pi00];
  sol[1] = (prims[Prims::q1] + aux[Aux::qv]*prims[Prims::v1])*aux[Aux::W] + aux[Aux::pi01];
  sol[2] = (prims[Prims::q2] + aux[Aux::qv]*prims[Prims::v2])*aux[Aux::W] + aux[Aux::pi02];
  sol[3] = (prims[Prims::q3] + aux[Aux::qv]*prims[Prims::v3])*aux[Aux::W] + aux[Aux::pi03];

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
  aux[Aux::rho_plus_p] = (cons[Cons::Tau] + cons[Cons::D] + sol[0])/(aux[Aux::W]*aux[Aux::W]) - prims[Prims::Pi];
  prims[Prims::v1] = (cons[Cons::S1] - sol[1])/((aux[Aux::rho_plus_p] + prims[Prims::Pi])*aux[Aux::W]*aux[Aux::W]);
  prims[Prims::v2] = (cons[Cons::S2] - sol[2])/((aux[Aux::rho_plus_p] + prims[Prims::Pi])*aux[Aux::W]*aux[Aux::W]);  
  prims[Prims::v3] = (cons[Cons::S3] - sol[3])/((aux[Aux::rho_plus_p] + prims[Prims::Pi])*aux[Aux::W]*aux[Aux::W]);  
  prims[Prims::p] = (aux[Aux::rho_plus_p] - prims[Prims::n])*((d->gamma-1)/d->gamma);
  prims[Prims::rho] = aux[Aux::rho_plus_p] - prims[Prims::p];
  
  // Repeating the ones here that depend on v1,v2,v3...
  aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
  aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
        
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
  const double tol = 1.4e-7;          // Tolerance of rootfinder
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
          prims[ID(Prims::q1+ndissvar, i, j, k)] = cons[ID(Cons::Y1+ndissvar, i, j, k)] / cons[ID(Cons::D, i, j, k)];
        }
        //aux[ID(Aux::W, i, j, k)] = 1 / ( 1 - (prims[ID(Prims::v1, i, j, k)]**2 + prims[ID(Prims::v2, i, j, k)]**2 + prims[ID(Prims::v3, i, j, k)]**2) )**0.5; // Point in re-calcing this here?
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) 
                               + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
        aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)];
        aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
            
        // Set additional args for rootfind
        args.D_rf = cons[ID(Cons::D, i, j, k)];
        args.S1_rf = cons[ID(Cons::S1, i, j, k)];
        args.S2_rf = cons[ID(Cons::S2, i, j, k)];
        args.S3_rf = cons[ID(Cons::S3, i, j, k)];
        args.Tau_rf = cons[ID(Cons::Tau, i, j, k)];
        args.q1_rf = prims[ID(Prims::q1, i, j, k)];
        args.q2_rf = prims[ID(Prims::q2, i, j, k)];
        args.q3_rf = prims[ID(Prims::q3, i, j, k)];
        args.Pi_rf = prims[ID(Prims::Pi, i, j, k)];
        args.pi11_rf = prims[ID(Prims::pi11, i, j, k)];
        args.pi12_rf = prims[ID(Prims::pi12, i, j, k)];
        args.pi13_rf = prims[ID(Prims::pi13, i, j, k)];
        args.pi22_rf = prims[ID(Prims::pi22, i, j, k)];
        args.pi23_rf = prims[ID(Prims::pi23, i, j, k)];
        args.pi33_rf = prims[ID(Prims::pi33, i, j, k)];
        args.gamma = d->gamma;
        
        sol[0] = prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)] - 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] - aux[ID(Aux::pi00, i, j, k)];
        sol[1] = (prims[ID(Prims::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
        sol[2] = (prims[ID(Prims::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
        sol[3] = (prims[ID(Prims::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];

        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&ISresidual, &args, sys_size, sol, res,
                                          tol, wa, lwa);
//        info = __cminpack_func__(hybrd) (&ISresidual, &args, sys_size, sol, res,
//                                          tol, maxfev, ml, mu, epsfcn, &diag[0], mode, factor, nprint, &nfev, &fjac[0][0], ldfjac, &r[0], lr, &qtf[0], &wa1[0], &wa2[0], &wa3[0], &wa4[0]);        
                                                                   
        // If root find fails, add failed cell to the list
        if (info!=1) {
          printf("(%f, %f) Y1,q1\n",  cons[ID(Cons::Y1, i, j, k)], prims[ID(Prims::q1, i, j, k)]);
          printf("(%f, %f, %f, %f) sol\n", sol[0],sol[1],sol[2],sol[3]);
          printf("%i info\n",info);
          printf("(%i, %i, %i) failed\n", i, j, k);
          printf("(%f, %f, %f, %f) res\n", res[0], res[1], res[2], res[3]);
          printf("(%f, %f, %f, %f) sol\n", sol[0], sol[1], sol[2], sol[3]);
          printf("(%f, %f, %f, %f, %f) prims\n",  prims[ID(Prims::p, i, j, k)], prims[ID(Prims::Pi, i, j, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::v1, i, j, k)], prims[ID(Prims::q1, i, j, k)]);
          printf("(%f, %f, %f, %f) aux\n",  aux[ID(Aux::W, i, j, k)], aux[ID(Aux::qv, i, j, k)], aux[ID(Aux::pi00, i, j, k)], aux[ID(Aux::pi01, i, j, k)]);
          printf("(%f, %f, %f, %f) cons\n",  cons[ID(Cons::Y1, i, j, k)], cons[ID(Cons::D, i, j, k)], cons[ID(Cons::S1, i, j, k)], cons[ID(Cons::Tau, i, j, k)]);
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
        aux[ID(Aux::vsqrd, i, j, k)] = ((cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])*(cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)]) 
                                  + (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])*(cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])
                                  + (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])*(cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)]))
                                  /((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])*(cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)]));
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt((1-aux[ID(Aux::vsqrd, i, j, k)]));
        prims[ID(Prims::n, i, j, k)] = cons[ID(Cons::D, i, j, k)] / aux[ID(Aux::W, i, j, k)];
        aux[ID(Aux::rho_plus_p, i, j, k)] = ((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])/(aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)])) - prims[ID(Prims::Pi, i, j, k)];
        prims[ID(Prims::v1, i, j, k)] = (cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                 + prims[ID(Prims::Pi, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);
        prims[ID(Prims::v2, i, j, k)] = (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                 + prims[ID(Prims::Pi, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);  
        prims[ID(Prims::v3, i, j, k)] = (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                 + prims[ID(Prims::Pi, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);  
        prims[ID(Prims::p, i, j, k)] = (aux[ID(Aux::rho_plus_p, i, j, k)] - prims[ID(Prims::n, i, j, k)])*((d->gamma-1)/d->gamma);
        prims[ID(Prims::rho, i, j, k)] = aux[ID(Aux::rho_plus_p, i, j, k)] - prims[ID(Prims::p, i, j, k)];

        // Again, repeating this here once the correct values for v1,v2,v3 have been set...
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
        // aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)]; // Should be unncessary here
        aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];

      } // End k-loop
    } // End j-loop
  } // End i-loop  

        
  // Recompute NS variables here? Let's try it

/*

  double kappa = this->data->optionalSimArgs[0];
//  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
//  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
//  double tau_pi = this->data->optionalSimArgs[5];
        
  // q_j,NS 10
  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {

        aux[ID(Aux::q1NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx) );
        aux[ID(Aux::q2NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy) );
        aux[ID(Aux::q3NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz) );
        
        //printf("(%i, %i, %i) ijk\n", i, j, k);
        //printf("(%f, %f, %f) \n", aux[ID(Aux::T, i+1, j, k)], aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::q1NS, i, j, k)]);
        //printf("(%f, %f, %f, %f) \n", prims[ID(Prims::p, i+1, j, k)], prims[ID(Prims::n, i+1, j, k)], prims[ID(Prims::p, i-1, j, k)], prims[ID(Prims::n, i-1, j, k)]);
        
//          aux[ID(Aux::q1NS, i, j, k)] = -kappa*aux[ID(Aux::T, i, j, k)] * ( (log(aux[ID(Aux::T, i+1, j, k)]) - log(aux[ID(Aux::T, i-1, j, k)]))/(2*d->dx) );
//          aux[ID(Aux::q2NS, i, j, k)] = -kappa*aux[ID(Aux::T, i, j, k)] * ( (log(aux[ID(Aux::T, i, j+1, k)]) - log(aux[ID(Aux::T, i, j-1, k)]))/(2*d->dy) );
//          aux[ID(Aux::q3NS, i, j, k)] = -kappa*aux[ID(Aux::T, i, j, k)] * ( (log(aux[ID(Aux::T, i, j, k+1)]) - log(aux[ID(Aux::T, i, j, k-1)]))/(2*d->dz) );
 
        // Theta 20 then Pi,NS 13 
        aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)] - aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v3, i, j, k+1)] - aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz);
        // Pi,NS = -zeta*Theta
        aux[ID(Aux::PiNS, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)];
  
        // pi^l_j,NS 14
        // 11
        aux[ID(Aux::pi11NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 12
        aux[ID(Aux::pi12NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 13
        aux[ID(Aux::pi13NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 22
        aux[ID(Aux::pi22NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 23
        aux[ID(Aux::pi23NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 33
        aux[ID(Aux::pi33NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
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

  // W, q_kv^k, pi^0_0
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(Aux::vsqrd, i, j, k)] = prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - aux[ID(Aux::vsqrd, i, j, k)] );
        prev_vars[ID(0, i, j, k)] = aux[ID(Aux::W, i, j, k)]; // Set here for time-differencing
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
        aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)];

        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        prims[ID(Prims::rho, i, j, k)] = prims[ID(Prims::n, i, j, k)]*(1+aux[ID(Aux::e, i, j, k)]);
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
        aux[ID(Aux::h, i, j, k)] = 1 + aux[ID(Aux::e, i, j, k)] + prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
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
          aux[ID(Aux::q1NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx) );
          aux[ID(Aux::q2NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy) );
          aux[ID(Aux::q3NS, i, j, k)] = -kappa* ( (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz) );
      }
    }
  }  
  // Theta 20 then Pi,NS 13 
  for (int i(1); i < d->Nx-1; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        aux[ID(Aux::Theta, i, j, k)] = aux[ID(Aux::dWdt, i, j, k)] + (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)] - aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v3, i, j, k+1)] - aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k-1)])/(2*d->dz);
        // Pi,NS = -zeta*Theta
        aux[ID(Aux::PiNS, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)];
      }
    }
  }  
  // pi^l_j,NS 14
  for (int i(1); i < d->Nx-1; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        // 11
        aux[ID(Aux::pi11NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 12
        aux[ID(Aux::pi12NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 13
        aux[ID(Aux::pi13NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 22
        aux[ID(Aux::pi22NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 23
        aux[ID(Aux::pi23NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy)
          + (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 33
        aux[ID(Aux::pi33NS, i, j, k)] = -2*eta*( (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
      }
    }
  }  

  // pi^0_j
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) { // Minus signs here (?)
         aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
         aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
         aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(Cons::D, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)];
        // S1,2,3
        cons[ID(Cons::S1, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)] 
          + (prims[ID(Prims::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
        cons[ID(Cons::S2, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)] 
          + (prims[ID(Prims::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
        cons[ID(Cons::S3, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)] 
          + (prims[ID(Prims::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v3, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];
        // Tau
        cons[ID(Cons::Tau, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] 
        - (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)] + prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]) 
        + 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi00, i, j, k)];
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

/*

  for (int i(1); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(Aux::0, i, j, k)] = (prims[ID(Prims::v1, i+1, j, k)]-prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(Aux::1, i, j, k)] = (prims[ID(Prims::v1, i, j+1, k)]-prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        aux[ID(Aux::2, i, j, k)] = (prims[ID(Prims::v1, i, j, k+1)]-prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);
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
        f[ID(0, i, j, k)] = cons[ID(Cons::D, i, j, k)]*prims[ID(dir, i, j, k)];
        // Sv + ..
        for (int nvar(0); nvar < 3; nvar++) {
          f[ID(1+nvar, i, j, k)] = cons[ID(Cons::S1+nvar, i, j, k)]*prims[ID(dir, i, j, k)] + ( prims[ID(Prims::q1+dir, i, j, k)] * prims[ID(Prims::v1+nvar, i, j, k)]  
            - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+nvar, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)] ) * aux[ID(Aux::W, i, j, k)];
          // (p+Pi)delta_ij
          if (dir == nvar) {
            f[ID(1+nvar, i, j, k)] += (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]);
          }
        }
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

        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Cons::Tau, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (prims[ID(Prims::q1+dir, i, j, k)] - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)])*aux[ID(Aux::W, i, j, k)]
          + aux[ID(Aux::pi01+dir, i, j, k)];
        // Y1-3,U,Z11-33 *v
        for (int ndissvar(0); ndissvar < 10; ndissvar++) {                
          f[ID(Y1+ndissvar, i, j, k)] = cons[ID(Cons::Y1+ndissvar, i, j, k)]*prims[ID(dir, i, j, k)];
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
          source[ID(1+dir, i, j, k)] = -(kappa_of_T(cons[ID(Cons::0, i, j, k)], kappa_0) * aux[ID(Aux::dir, i, j, k)] +
                                         prims[ID(Prims::v2+dir, i, j, k)]) / tau_q_of_T(cons[ID(Cons::0, i, j, k)], tau_q_0);
        }
      }
    }
  }
}

*/

