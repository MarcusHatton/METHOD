#include "IS.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

IS::IS() : Model()
{
  this->Ncons = 4;
  this->Nprims = 4;
  this->Naux = 3;
}

IS::IS(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 4;
  this->Nprims = (this->data)->Nprims = 4;
  this->Naux = (this->data)->Naux = 3;

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
  this->data->auxLabels.push_back("h");     this->data->auxLabels.push_back("T")     
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
  this->data->auxLabels.push_back("Theta");  this->data->auxLabels.push_back("a");
       
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
  for (int dir(0); dir < 3; dir++) {
    source[5+dir] = (n / tau_q) * (aux[dir+4] - prims[dir+8]) 
  }
  // U
  source[8] = (n / tau_Pi) * (aux[7] - prims[11])
  // Z11,12,13,22,23,33
  for (int dir(0); dir < 6; dir++) {
    source[9+dir] = (n / tau_pi) * (aux[dir+8] - prims[dir+12])
  }

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
        source[ID(0, i, j, k)] = 0.0;
        // S1,2,3
        source[ID(1, i, j, k)] = 0.0;
        source[ID(2, i, j, k)] = 0.0;
        source[ID(3, i, j, k)] = 0.0;
        // Tau
        source[ID(4, i, j, k)] = 0.0;
        
        // Y1,2,3
        for (int dir(0); dir < 3; dir++) {
          source[ID(5+dir, i, j, k)] = (n / tau_q) * (aux[ID(dir+4, i, j, k)] - prims[ID(dir+8, i, j, k)]) 
        }
        // U
        source[8] = (n / tau_Pi) * (aux[ID(7, i, j, k)] - prims[ID(11, i, j, k)])
        // Z11,12,13,22,23,33
        for (int dir(0); dir < 6; dir++) {
          source[ID(9+dir] = (n / tau_pi) * (aux[ID(dir+8, i, j, k)] - prims[ID(dir+12, i, j, k)])
        }
        
        
      }
    }
  }
}

void IS::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  // printf("ToyQ model does not implement getPrimitiveVarsSingleCell\n");
  // exit(1);
  // Y1-3,U,Z11-33
  for (int ndissvar(0); ndissvar < 10; nvar++) {
    prims[ndissvar+8] = cons[ndissvar+5] / cons[0];
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
        aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
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
        aux[ID(3, i, j, k)] = 1 / ( 1 - (prims[ID(0, i, j, k)]**2 + prims[ID(1, i, j, k)]**2 + prims[ID(2, i, j, k)]**2) )**0.5;
      }
    }
  }

  // q_k v^k
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(6, i, j, k)] = (prims[ID(6, i, j, k)] * prims[ID(0, i, j, k)]) + (prims[ID(7, i, j, k)] * prims[ID(0, i, j, k)]) + (prims[ID(8, i, j, k)] * prims[ID(0, i, j, k)]);
      }
    }
  }
  
  // pi^0_0
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(6, i, j, k)] = prims[ID(10, i, j, k)] + prims[ID(13, i, j, k)] + prims[ID(15, i, j, k)];
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
          aux[ID(10, i, j, k)] = -kappa*aux[ID(1, i, j, k)] * ( (ln(aux[ID(1, i+1, j, k)]) - ln(aux[ID(1, i-1, j, k)]))/(2*d->dx) );
          aux[ID(11, i, j, k)] = -kappa*aux[ID(1, i, j, k)] * ( (ln(aux[ID(1, i, j+1, k)]) - ln(aux[ID(1, i, j-1, k)]))/(2*d->dy) );
          aux[ID(12, i, j, k)] = -kappa*aux[ID(1, i, j, k)] * ( (ln(aux[ID(1, i, j, k+1)]) - ln(aux[ID(1, i, j, k-1)]))/(2*d->dz) );
      }
    }
  }  
  // Theta 20 then Pi,NS 13 
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(20, i, j, k)] = ( dtu0 + (aux[ID(3, i+1, j, k)]*prims[ID(0, i+1, j, k)] - aux[ID(3, i-1, j, k)]*prims[ID(0, i-1, j, k)])/(2*d->dx) 
          + (aux[ID(3, i, j+1, k)]*prims[ID(1, i, j+1, k)] - aux[ID(3, i, j, k)]*prims[ID(1, i, j-1, k)])/(2*d->dy)
          + (aux[ID(3, i, j, k+1)]*prims[ID(2, i, j, k+1)] - aux[ID(3, i, j, k)]*prims[ID(2, i, j, k-1)])/(2*d->dz) );
        // Pi,NS = -zeta*Theta
        aux[ID(13, i, j, k)] = -zeta * aux[ID(20, i, j, k)];
      }
    }
  }  
  // pi^l_j,NS 14
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // 11
        aux[ID(14, i, j, k)] = -2*eta*( (aux[ID(3, i+1, j, k)]*prims[ID(0, i+1, j, k)] - aux[ID(3, i-1, j, k)]*prims[ID(0, i-1, j, k)])/(d->dx)
          - (2/3)*(1 + (aux[ID(3, i, j, k)]*prims[ID(0, i, j, k)])**2)*aux[ID(20, i, j, k)] );
        // 12
        aux[ID(15, i, j, k)] = -2*eta*( (aux[ID(3, i+1, j, k)]*prims[ID(1, i+1, j, k)] - aux[ID(3, i-1, j, k)]*prims[ID(1, i-1, j, k)])/(2*d->dx)
          + (aux[ID(3, i, j+1, k)]*prims[ID(0, i, j+1, k)] - aux[ID(3, i, j-1, k)]*prims[ID(0, i, j-1, k)])/(2*d->dy)
          - (2/3)*((aux[ID(3, i, j, k)]*prims[ID(0, i, j, k)])*(aux[ID(3, i, j, k)]*prims[ID(1, i, j, k)]))*aux[ID(20, i, j, k)] );
        // 13
        aux[ID(16, i, j, k)] = -2*eta*( (aux[ID(3, i+1, j, k)]*prims[ID(2, i+1, j, k)] - aux[ID(3, i-1, j, k)]*prims[ID(2, i-1, j, k)])/(2*d->dx)
          + (aux[ID(3, i, j, k+1)]*prims[ID(0, i, j, k+1)] - aux[ID(3, i, j, k-1)]*prims[ID(0, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(3, i, j, k)]*prims[ID(0, i, j, k)])*(aux[ID(3, i, j, k)]*prims[ID(2, i, j, k)]))*aux[ID(20, i, j, k)] );
        // 22
        aux[ID(17, i, j, k)] = -2*eta*( (aux[ID(3, i, j+1, k)]*prims[ID(0, i, j+1, k)] - aux[ID(3, i, j-1, k)]*prims[ID(0, i, j-1, k)])/(d->dy)
          - (2/3)*(1 + (aux[ID(3, i, j, k)]*prims[ID(1, i, j, k)])**2)*aux[ID(20, i, j, k)] );
        // 23
        aux[ID(18, i, j, k)] = -2*eta*( (aux[ID(3, i, j+1, k)]*prims[ID(2, i, j+1, k)] - aux[ID(3, i, j-1, k)]*prims[ID(2, i, j-1, k)])/(2*d->dy)
          + (aux[ID(3, i, j, k+1)]*prims[ID(1, i, j, k+1)] - aux[ID(3, i, j, k-1)]*prims[ID(1, i, j, k-1)])/(2*d->dz)
          - (2/3)*((aux[ID(3, i, j, k)]*prims[ID(1, i, j, k)])*(aux[ID(3, i, j, k)]*prims[ID(2, i, j, k)]))*aux[ID(20, i, j, k)] );
        // 33
        aux[ID(19, i, j, k)] = -2*eta*( (aux[ID(3, i, j, k+1)]*prims[ID(0, i, j, k+1)] - aux[ID(3, i, j, k-1)]*prims[ID(0, i, j, k-1)])/(d->dz)
          - (2/3)*(1 + (aux[ID(3, i, j, k)]*prims[ID(2, i, j, k)])**2)*aux[ID(20, i, j, k)] );
      }
    }
  }  
  
  

  // pi^0_j
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
         aux[ID(7, i, j, k)] = prims[ID(10, i, j, k)]*prims[ID(0, i, j, k)] + prims[ID(11, i, j, k)]*prims[ID(1, i, j, k)] + prims[ID(12, i, j, k)]*prims[ID(2, i, j, k)];
         aux[ID(8, i, j, k)] = prims[ID(11, i, j, k)]*prims[ID(0, i, j, k)] + prims[ID(13, i, j, k)]*prims[ID(1, i, j, k)] + prims[ID(14, i, j, k)]*prims[ID(2, i, j, k)];
         aux[ID(9, i, j, k)] = prims[ID(12, i, j, k)]*prims[ID(0, i, j, k)] + prims[ID(14, i, j, k)]*prims[ID(1, i, j, k)] + prims[ID(15, i, j, k)]*prims[ID(2, i, j, k)];
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(0, i, j, k)] = prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)];
        // S1,2,3
        for (int nvar(0); nvar < 3; nvar++) {
          cons[ID(nvar+1, i, j, k)] = (prims[ID(4, i, j, k)] + prims[ID(3, i, j, k)] + prims[ID(?, i, j, k)])Pi * prims[ID(6, i, j, k)]**2 * prims[ID(nvar, i, j, k)] 
            + (prims[ID(nvar+7, i, j, k)] + aux[ID(6, i, j, k)] * prims[ID(nvar, i, j, k)]) * prims[ID(6, i, j, k)] + aux[ID(nvar+7, i, j, k)];
        }
        // Tau
        cons[ID(4, i, j, k)] = (prims[ID(4, i, j, k)] + prims[ID(3, i, j, k)] + prims[ID(?, i, j, k)]) * prims[ID(6, i, j, k)]**2 
        - (prims[ID(3, i, j, k)] + prims[ID(?, i, j, k)]Pi + prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)]) 
        + 2*aux[ID(6, i, j, k)]*prims[ID(6, i, j, k)] + aux[ID(6, i, j, k)];
        // Y1-3
        for (int nvar(0); nvar < 3; nvar++) {
           cons[ID(5+nvar, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)] * prims[ID(nvar+7, i, j, k)];
        }   
        // U
        cons[ID(8, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)] * prims[ID(10, i, j, k)];
        // Z11-33
        for (int nvar(0); nvar < 6; nvar++) {
           cons[ID(9+nvar, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)] * prims[ID(nvar+11, i, j, k)];
        }      
      }
    }
  }

/*

  for (int i(1); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
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
        f[ID(0, i, j, k)] = cons[ID(0, i, j, k)]*prims[ID(dir, i, j, k)];
        // Sv + ..
        for (int nvar(1); nvar < 4; nvar++) {
          f[ID(nvar, i, j, k)] = cons[ID(nvar, i, j, k)]*prims[ID(dir, i, j, k)] + ( prims[ID(dir+6, i, j, k)] * prims[ID(nvar-1, i, j, k)]  
            - aux[ID(6, i, j, k)]*prims[ID(nvar-1, i, j, k)]*prims[ID(dir, i, j, k)] ) * aux[ID(3, i, j, k)] + prims[ID(nvar+dir+10, i, j, k)]pi_ij ; // need to add pi21,pi31,pi32 for this to work 
          // (p+Pi)delta_ij
          if (dir == nvar-1) {
            f[ID(nvar, i, j, k)] += (prims[ID(6, i, j, k)] + prims[ID(9, i, j, k)]);
          }
        }
        // Tau*v + ...
        f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(3, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (prims[ID(dir+6, i, j, k)] - aux[ID(6, i, j, k)]*prims[ID(dir, i, j, k)]) + aux[ID(7+dir, i, j, k)];
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
                                         prims[ID(1+dir, i, j, k)]) / tau_q_of_T(cons[ID(0, i, j, k)], tau_q_0);
        }
      }
    }
  }
}

*/

