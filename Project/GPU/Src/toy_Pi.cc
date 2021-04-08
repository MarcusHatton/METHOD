#include "toy_Pi.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>

ToyPi::ToyPi() : Model()
{
  this->Ncons = 9;
  this->Nprims = 9;
  this->Naux = 7;
}

ToyPi::ToyPi(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 9;
  this->Nprims = (this->data)->Nprims = 9;
  this->Naux = (this->data)->Naux = 7;

  this->data->consLabels.push_back("vx"); this->data->consLabels.push_back("vy");
  this->data->consLabels.push_back("vz");  this->data->consLabels.push_back("pixx");
  this->data->consLabels.push_back("piyy"); this->data->consLabels.push_back("pizz");
  this->data->consLabels.push_back("pixy"); this->data->consLabels.push_back("pixz");
  this->data->consLabels.push_back("piyz");

  this->data->primsLabels.push_back("vx"); this->data->primsLabels.push_back("vy");
  this->data->primsLabels.push_back("vz");  this->data->primsLabels.push_back("pixx");
  this->data->primsLabels.push_back("piyy"); this->data->primsLabels.push_back("pizz");
  this->data->primsLabels.push_back("pixy"); this->data->primsLabels.push_back("pixz");
  this->data->primsLabels.push_back("piyz");

  this->data->auxLabels.push_back("Theta"); this->data->auxLabels.push_back("sigmaxx");  
  this->data->auxLabels.push_back("sigmayy"); this->data->auxLabels.push_back("sigmazz"); 
  this->data->auxLabels.push_back("sigmaxy"); this->data->auxLabels.push_back("sigmaxz"); 
  this->data->auxLabels.push_back("sigmayz");
}

ToyPi::~ToyPi()
{
}


void ToyPi::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  // printf("ToyPi model does not implement sourceTermSingleCell\n");
  // exit(1);

  float eta = this->data->optionalSimArgs[0];
  float tau_Pi = this->data->optionalSimArgs[1];

  for (int nvar(0); nvar < 9; nvar++) {
    if (nvar <= 2) {
      source[nvar] = 0.0;
      } else {
        source[nvar] = ( 1 / tau_Pi ) * ( -2 * eta * aux[nvar-2] - prims[nvar] );
    }
  }
}

void ToyPi::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  float eta = this->data->optionalSimArgs[0];
  float tau_Pi = this->data->optionalSimArgs[1];

//  printf("Calling source\n");

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int nvar(0); nvar < 9; nvar++) {
          if (nvar <= 2) {
            source[ID(nvar, i, j, k)] = 0.0;
          } else {
            source[ID(nvar, i, j, k)] = ( 1 / tau_Pi ) * ( -2 * eta * aux[ID(nvar-2, i, j, k)] - prims[ID(nvar, i, j, k)] );
          }
        }
      }
    }
  }
}

void ToyPi::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // printf("ToyPi model does not implement getPrimitiveVarsSingleCell\n");
  // exit(1);
  for (int nvar(0); nvar < 9; nvar++) {
    prims[nvar] = cons[nvar];
  }
}

void ToyPi::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

//  printf("Calling getPrimVars %i %i %i %i %i %i\n",
//         d->is, d->ie, d->js, d->je, d->ks, d->ke);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int nvar(0); nvar < 9; nvar++) {
          prims[ID(nvar, i, j, k)] = cons[ID(nvar, i, j, k)];
        }
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
          // Theta
          aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
          // sigmaxx,yy,zz                      
          aux[ID(1, i, j, k)] = 2*((prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx));

          // sigmaxy,xz,yz
          aux[ID(4, i, j, k)] = ((prims[ID(1, i+1, j, k)]-prims[ID(1, i-1, j, k)])/(2*d->dx));   
          aux[ID(5, i, j, k)] = ((prims[ID(2, i+1, j, k)]-prims[ID(2, i-1, j, k)])/(2*d->dx));  
        }
      }
    }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
            // Theta
            aux[ID(0, i, j, k)] += (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
            // sigmaxx,yy,zz                      
            aux[ID(2, i, j, k)] = 2*((prims[ID(1, i, j+1, k)]-prims[ID(1, i, j-1, k)])/(2*d->dy));
  
            // sigmaxy,xz,yz
            aux[ID(4, i, j, k)] += ((prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy));   
            aux[ID(6, i, j, k)] = ((prims[ID(2, i, j+1, k)]-prims[ID(2, i, j-1, k)])/(2*d->dy));  
          }
        }
      }    
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
              // Theta
              aux[ID(0, i, j, k)] += (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
              // sigmazz                      
              aux[ID(3, i, j, k)] = 2*((prims[ID(2, i, j, k+1)]-prims[ID(2, i, j, k-1)])/(2*d->dz));
    
              // sigmaxz,yz
              aux[ID(5, i, j, k)] += ((prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz) );   
              aux[ID(6, i, j, k)] += ((prims[ID(1, i, j, k+1)]-prims[ID(1, i, j, k-1)])/(2*d->dz));   
            }
          }
        }
      }
    } 
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // sigmaxx,yy,zz                      
        aux[ID(1, i, j, k)] -= (2/3)*aux[ID(0, i, j, k)];
        aux[ID(2, i, j, k)] -= (2/3)*aux[ID(0, i, j, k)];    
        aux[ID(3, i, j, k)] -= (2/3)*aux[ID(0, i, j, k)];    
      }
    }
  }      

/*  
    for (int i(1); i < d->Nx-1; i++) {
      for (int j(1); j < d->Ny-1; j++) {
        for (int k(1); k < d->Nz-1; k++) {
          // Theta
          aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx)
                                + (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy)
                                + (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
          // sigmaxx,yy,zz                      
          aux[ID(1, i, j, k)] = (2*((prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx))
                                - (2/3)*aux[ID(0, i, j, k)]);
          aux[ID(2, i, j, k)] = (2*((prims[ID(1, i, j+1, k)]-prims[ID(1, i, j-1, k)])/(2*d->dy))
                                - (2/3)*aux[ID(0, i, j, k)]);                                
          aux[ID(3, i, j, k)] = (2*((prims[ID(2, i, j, k+1)]-prims[ID(2, i, j, k-1)])/(2*d->dz))
                                - (2/3)*aux[ID(0, i, j, k)]);
          // sigmaxy,xz,yz
          aux[ID(4, i, j, k)] = ( (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy) 
                                + (prims[ID(1, i+1, j, k)]-prims[ID(1, i-1, j, k)])/(2*d->dx) );   
          aux[ID(5, i, j, k)] = ( (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz) 
                                + (prims[ID(2, i+1, j, k)]-prims[ID(2, i-1, j, k)])/(2*d->dx) );  
          aux[ID(6, i, j, k)] = ( (prims[ID(2, i, j+1, k)]-prims[ID(2, i, j-1, k)])/(2*d->dy) 
                                + (prims[ID(1, i, j, k+1)]-prims[ID(1, i, j, k-1)])/(2*d->dz) );  
        }
      }
    }
    
  for (int i(d->is+1); i < d->ie-1; i++) {
    if (d->js+1 > d->je-1) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
        }
      }
    }
    else {
      for (int j(d->js+1); j < d->je-1; j++) {
        if (d->ks+1 > d->ke-1) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx)
                                + (prims[ID(1, i, j+1, k)]-prims[ID(1, i, j-1, k)])/(2*d->dy);
          }
        }
        else {
          for (int k(d->ks+1); k < d->ke-1; k++) {
            aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx)
                                + (prims[ID(1, i, j+1, k)]-prims[ID(1, i, j-1, k)])/(2*d->dy)
                                + (prims[ID(2, i, j, k+1)]-prims[ID(2, i, j, k-1)])/(2*d->dz);
          }
        }
      }
    }
  }
*/
}

void ToyPi::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  printf("Calling primsToAll\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int nvar(0); nvar < 9; nvar++) {
          cons[ID(nvar, i, j, k)] = prims[ID(nvar, i, j, k)];
        }
      }
    }
  }

    for (int i(1); i < d->Nx-1; i++) {
      for (int j(1); j < d->Ny-1; j++) {
        for (int k(1); k < d->Nz-1; k++) {
          // Theta
          aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx)
                                + (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy)
                                + (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
          // sigmaxx,yy,zz                      
          aux[ID(1, i, j, k)] = (2*((prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx))
                                - (2/3)*aux[ID(0, i, j, k)]);
          aux[ID(2, i, j, k)] = (2*((prims[ID(1, i, j+1, k)]-prims[ID(1, i, j-1, k)])/(2*d->dy))
                                - (2/3)*aux[ID(0, i, j, k)]);                                
          aux[ID(3, i, j, k)] = (2*((prims[ID(2, i, j, k+1)]-prims[ID(2, i, j, k-1)])/(2*d->dz))
                                - (2/3)*aux[ID(0, i, j, k)]);
          // sigmaxy,xz,yz
          aux[ID(4, i, j, k)] = ( (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy) 
                                + (prims[ID(1, i+1, j, k)]-prims[ID(1, i-1, j, k)])/(2*d->dx) );   
          aux[ID(5, i, j, k)] = ( (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz) 
                                + (prims[ID(2, i+1, j, k)]-prims[ID(2, i-1, j, k)])/(2*d->dx) );  
          aux[ID(6, i, j, k)] = ( (prims[ID(2, i, j+1, k)]-prims[ID(2, i, j-1, k)])/(2*d->dy) 
                                + (prims[ID(1, i, j, k+1)]-prims[ID(1, i, j, k-1)])/(2*d->dz) );  
        }
      }
    }

}

void ToyPi::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
  // v
  f[ID(0, i, j, k)] = cons[ID(dir, i, j, k)] * cons[ID(0, i, j, k)];	
  f[ID(1, i, j, k)] = cons[ID(dir, i, j, k)] * cons[ID(1, i, j, k)];	
  f[ID(2, i, j, k)] = cons[ID(dir, i, j, k)] * cons[ID(2, i, j, k)];	
  // pi
  f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * cons[ID(dir, i, j, k)];
  f[ID(4, i, j, k)] = cons[ID(4, i, j, k)] * cons[ID(dir, i, j, k)];
  f[ID(5, i, j, k)] = cons[ID(5, i, j, k)] * cons[ID(dir, i, j, k)];
  f[ID(6, i, j, k)] = cons[ID(6, i, j, k)] * cons[ID(dir, i, j, k)];
  f[ID(7, i, j, k)] = cons[ID(7, i, j, k)] * cons[ID(dir, i, j, k)];
  f[ID(8, i, j, k)] = cons[ID(8, i, j, k)] * cons[ID(dir, i, j, k)];
  
  if (dir == 0) {
    f[ID(0, i, j, k)] += cons[ID(3, i, j, k)]; // pi_xx
    f[ID(1, i, j, k)] += cons[ID(6, i, j, k)]; // pi_xy
    f[ID(2, i, j, k)] += cons[ID(7, i, j, k)]; // pi_xz
  } else if (dir == 1) {
    f[ID(0, i, j, k)] += cons[ID(6, i, j, k)]; // pi_xy
    f[ID(1, i, j, k)] += cons[ID(4, i, j, k)]; // pi_yy
    f[ID(2, i, j, k)] += cons[ID(8, i, j, k)]; // pi_yz
  } else if (dir == 2) {
    f[ID(0, i, j, k)] += cons[ID(7, i, j, k)]; // pi_xz
    f[ID(1, i, j, k)] += cons[ID(8, i, j, k)]; // pi_yz
    f[ID(2, i, j, k)] += cons[ID(5, i, j, k)]; // pi_zz   
  } else {
    printf("eek");
  }
      } // End k loop
    } // End j loop
  } // End i loop
}
