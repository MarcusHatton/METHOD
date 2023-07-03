#include "DEIFY.h"
#include <cstdio>
#include <cmath>
#include <ctime>


// Dissipative Extension for Ideal Fluid dYnamics
DEIFY::DEIFY()
{
  sourceExists = true;
}

DEIFY::DEIFY(Data * data, FluxMethod * fluxMethod) : ModelExtension(data), fluxMethod(fluxMethod)
{
  //  Syntax
  Data * d(this->data);

  sourceExists = true;
//  fluxExists = true;

  // Allocate arrays
  // Additional Flux vectors from NS Dissipation Terms
  Fx = new double[d->Nx*d->Ny*d->Nz*5] ();
  Fy = new double[d->Nx*d->Ny*d->Nz*5] ();
  Fz = new double[d->Nx*d->Ny*d->Nz*5] ();
  dtH = new double[d->Nx*d->Ny*d->Nz*5] ();
  d->sourceExtension = new double[d->Nx*d->Ny*d->Nz*5] ();
}

DEIFY::~DEIFY()
{
  //  Syntax
  Data * d(this->data);

  delete[] Fx;
  delete[] Fy;
  delete[] Fz;
  delete[] dtH;
  delete[] d->sourceExtension;
}

void DEIFY::sourceExtension(double * cons, double * prims, double * aux, double * source)
{
  // Syntax
  Data * d(this->data);

  // Set vars - dissipative NS forms
  this->set_vars(cons, prims, aux);
 
  // Determine the diffusion vectors
  this->set_Fx(cons, prims, aux);
  for (int var(0); var<5; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
            source[ID(var, i, j, k)] = -(Fx[ID(var, i+1, j, k)] - Fx[ID(var, i-1, j, k)]) / (2*d->dx);
        }
      }
    }
  }

  if (d->dims>1) {
    this->set_Fy(cons, prims, aux);
    for (int var(0); var<5; var++) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
              source[ID(var, i, j, k)] += -(Fy[ID(var, i, j+1, k)] - Fy[ID(var, i, j-1, k)]) / (2*d->dy);
          }
        }
      }
    }
  }

  if (d->dims==3) {
    this->set_Fz(cons, prims, aux);
    for (int var(0); var<5; var++) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            source[ID(var, i, j, k)] += -(Fz[ID(var, i, j, k+1)] - Fz[ID(var, i, j, k-1)]) / (2*d->dz);
          }
        }
      }
    }
  }

  this->set_dtH(cons, prims, aux);
  for (int var(0); var<5; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          source[ID(var, i, j, k)] += -dtH[ID(var, i, j, k)];
        }
      }
    }
  }

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

double minmidGradSOGeneral(int enum_number, int dir, int c_p_or_a, double * cons, double * prims, double * aux, int i, int j, int k, Data * d) {
  double FDSum = 0.0;
  double BDSum = 0.0;
  float Coeff_stencil [3] = {-1.5, 2.0, -0.5};
  int stencil [3] = {0,0,0};
  double dX = 0.0;
  switch (dir) {
    case 0:
      dX = d->dx;
      stencil[0] = 1;
    case 1:
      dX = d->dy;
      stencil[1] = 1;
    case 2:
      dX = d->dz;
      stencil[2] = 1;
    }

  if (c_p_or_a == 0) {
    for (int step=0; step<3; step++) {
      FDSum += Coeff_stencil[step]*cons[ID(enum_number, i+stencil[0], j+stencil[1], k+stencil[2])];
      BDSum -= Coeff_stencil[step]*cons[ID(enum_number, i-stencil[0], j-stencil[1], k-stencil[2])];
    }
  //   FDGrad = ( -1.5*cons[ID(enum_number, i, j, k)] + 2*cons[ID(enum_number, i+stencil[0], j+stencil[1], k+stencil[2])] 
  //            - 0.5*cons[ID(enum_number, i+2*stencil[0], j+2*stencil[1], k+2*stencil[2])] )/ dx;
  //   BDGrad = 1.5*cons[ID(enum_number, i, j, k)] - 2*cons[ID(enum_number, i-stencil[0], j-stencil[1], k-stencil[2])] 
  //            + 0.5*cons[ID(enum_number, i-2*stencil[0], j-2*stencil[1], k-2*stencil[2])];
  } else if (c_p_or_a == 1) {
    for (int step=0; step<3; step++) {
      FDSum += Coeff_stencil[step]*prims[ID(enum_number, i+stencil[0], j+stencil[1], k+stencil[2])];
      BDSum -= Coeff_stencil[step]*prims[ID(enum_number, i-stencil[0], j-stencil[1], k-stencil[2])];
    }
  } else if (c_p_or_a == 2) {
    for (int step=0; step<3; step++) {
      FDSum += Coeff_stencil[step]*aux[ID(enum_number, i+stencil[0], j+stencil[1], k+stencil[2])];
      BDSum -= Coeff_stencil[step]*aux[ID(enum_number, i-stencil[0], j-stencil[1], k-stencil[2])];
    }
  } else {
      throw std::invalid_argument("Need to pass {0,1,2} == {C,P,A}\n");
  }
  double FDGrad = FDSum/dX;
  double BDGrad = BDSum/dX;
  if ( (FDGrad < 0 && BDGrad > 0) || (FDGrad > 0 && BDGrad < 0) ) {
    return 0;
  } else {
    return abs(FDGrad) < abs(BDGrad) ? FDGrad : BDGrad;
  }
}

void DEIFY::set_vars(double * cons, double * prims, double * aux)
{
  Data * d(this->data);

  double kappa = this->data->optionalSimArgs[0];
//  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
//  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
  printf("eta: %5.2f\n",eta);
//  double tau_pi = this->data->optionalSimArgs[5];

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

  double u_x;
  double u_y;
  double u_z;
 
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        // dxT = (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx);
        // dyT = (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy);
        // dzT = (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz);

        // dxux = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx);
        // dyuy = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        // dzuz = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);

        // dxuy = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx);
        // dxuz = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx);
        // dyux = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        // dyuz = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy);
        // dzux = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);
        // dzuy = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz);

        u_x = aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)];
        u_y = aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)];
        u_z = aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)];

        // dxT = minmodGradSO(aux[ID(Aux::T, i-2, j, k)], aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i+1, j, k)], aux[ID(Aux::T, i+2, j, k)], d->dx);
        // dyT = minmodGradSO(aux[ID(Aux::T, i, j-2, k)], aux[ID(Aux::T, i, j-1, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j+1, k)], aux[ID(Aux::T, i, j+2, k)], d->dy);
        // dzT = minmodGradSO(aux[ID(Aux::T, i, j, k-2)], aux[ID(Aux::T, i, j, k-1)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j, k+1)], aux[ID(Aux::T, i, j, k+2)], d->dz);
        
        // dxux = minmodGradSO(aux[ID(Aux::W, i-2, j, k)]*prims[ID(Prims::v1, i-2, j, k)], aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
        //                     aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)], aux[ID(Aux::W, i+2, j, k)]*prims[ID(Prims::v1, i+2, j, k)], d->dx);
        // dyuy = minmodGradSO(aux[ID(Aux::W, i, j-2, k)]*prims[ID(Prims::v2, i, j-2, k)], aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v2, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
        //                     aux[ID(Aux::W, i, j+2, k)]*prims[ID(Prims::v2, i, j+2, k)], d->dy);
        // dzuz = minmodGradSO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
        //                     aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);                          

        // dxuy = minmodGradSO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
        //                     aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)], d->dx);
        // dxuz = minmodGradSO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
        //                     aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)], d->dx);
       
        // dyux = minmodGradSO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
        //                     aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)], d->dy);
        // dyuz = minmodGradSO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
        //                     aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)], d->dy);

        // dzux = minmodGradSO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
        //                     aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)], d->dz);  
        // dzuy = minmodGradSO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
        //                     aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);  

        // dxT = minmodGradFO(aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i+1, j, k)], d->dx);
        // dyT = minmodGradFO(aux[ID(Aux::T, i, j-1, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j+1, k)], d->dy);
        // dzT = minmodGradFO(aux[ID(Aux::T, i, j, k-1)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j, k+1)], d->dz);
        
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

        dxT = minmidGradSOGeneral(Aux::T, 0, 2, cons, prims, aux, i, j, k, d);
        dyT = minmidGradSOGeneral(Aux::T, 1, 2, cons, prims, aux, i, j, k, d);
        dzT = minmidGradSOGeneral(Aux::T, 2, 2, cons, prims, aux, i, j, k, d);

        // dxvx = minmidGradSOGeneral(Prims::v1, 0, 1, cons, prims, aux, i, j, k, d);
        // dyvx = minmidGradSOGeneral(Prims::v1, 1, 1, cons, prims, aux, i, j, k, d);
        // dzvx = minmidGradSOGeneral(Prims::v1, 2, 1, cons, prims, aux, i, j, k, d);
        // dxvy = minmidGradSOGeneral(Prims::v2, 0, 1, cons, prims, aux, i, j, k, d);
        // dyvy = minmidGradSOGeneral(Prims::v2, 1, 1, cons, prims, aux, i, j, k, d);
        // dzvy = minmidGradSOGeneral(Prims::v2, 2, 1, cons, prims, aux, i, j, k, d);
        // dxvz = minmidGradSOGeneral(Prims::v3, 0, 1, cons, prims, aux, i, j, k, d);
        // dyvz = minmidGradSOGeneral(Prims::v3, 1, 1, cons, prims, aux, i, j, k, d);
        // dzvz = minmidGradSOGeneral(Prims::v3, 2, 1, cons, prims, aux, i, j, k, d);


        aux[ID(Aux::a1, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( aux[ID(Aux::W, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] 
          + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*dxux
          + prims[ID(Prims::v2, i, j, k)]*dyux + prims[ID(Prims::v3, i, j, k)]*dzux );
        
        aux[ID(Aux::a2, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( aux[ID(Aux::W, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] 
          + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*dxuy
          + prims[ID(Prims::v2, i, j, k)]*dyuy + prims[ID(Prims::v3, i, j, k)]*dzuy );
        
        aux[ID(Aux::a3, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( aux[ID(Aux::W, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] 
          + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*dxuz
          + prims[ID(Prims::v2, i, j, k)]*dyuz + prims[ID(Prims::v3, i, j, k)]*dzuz );

        aux[ID(Aux::q1NS, i, j, k)] = -kappa* ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*dxT 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*dyT
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*dzT );
        aux[ID(Aux::q2NS, i, j, k)] = -kappa* ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*dyT 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*dxT
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*dzT );
        aux[ID(Aux::q3NS, i, j, k)] = -kappa* ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*dzT 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*dxT
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*dyT );

        // Theta 20 then Pi,NS 13 
        aux[ID(Aux::Theta, i, j, k)] = aux[ID(TDerivs::dtW, i, j, k)] + dxux + dyuy + dzuz;
        // Pi,NS = -zeta*Theta
        aux[ID(Aux::PiNS, i, j, k)] = -zeta * aux[ID(Aux::Theta, i, j, k)];
  
        // pi^l_j,NS 14 - STILL NOT FULLY CORRECT (but much better, I think only missing O(u^4) terms now!)
        // 11
        aux[ID(Aux::pi11NS, i, j, k)] = -eta*( 2*dxux - (2/3)*(1 + u_x*u_x)*aux[ID(Aux::Theta, i, j, k)] 
          + 2*u_x * ( 1 + u_x*u_x + u_y*u_y + u_z*u_z ) * ( 2*u_x*dxux + u_y*dxuy + u_z*dxuz + u_y*dyux + u_z*dzux ) );
        // 12
        aux[ID(Aux::pi12NS, i, j, k)] = -eta*( dxuy + dyux - (2/3)*(2 + u_x*u_x + u_y*u_y + u_z*u_z )*aux[ID(Aux::Theta, i, j, k)] );
        // 13
        aux[ID(Aux::pi13NS, i, j, k)] = -eta*( dxuz + dzux - (2/3)*(2 + u_x*u_x + u_y*u_y + u_z*u_z )*aux[ID(Aux::Theta, i, j, k)] );
        // 22
        aux[ID(Aux::pi22NS, i, j, k)] = -eta*( 2*dyuy - (2/3)*(1 + u_y*u_y)*aux[ID(Aux::Theta, i, j, k)] 
          + 2*u_y * ( 1 + u_x*u_x + u_y*u_y + u_z*u_z ) * ( 2*u_y*dyuy + u_z*dyuz + u_x*dyux + u_z*dzuy + u_x*dxuy ) );
        // 23
        aux[ID(Aux::pi23NS, i, j, k)] = -eta*( dyuz + dzuy - (2/3)*(2 + u_x*u_x + u_y*u_y + u_z*u_z )*aux[ID(Aux::Theta, i, j, k)] );
        // 33
        aux[ID(Aux::pi33NS, i, j, k)] = -eta*( 2*dzuz - (2/3)*(1 + u_z*u_z)*aux[ID(Aux::Theta, i, j, k)] 
          + 2*u_z * ( 1 + u_x*u_x + u_y*u_y + u_z*u_z ) * ( 2*u_z*dzuz + u_x*dzux + u_y*dzuy + u_x*dxuz + u_y*dyuz ) );

        aux[ID(TDerivs::dtE, i, j, k)] = aux[ID(TDerivs::dtTau, i, j, k)] + aux[ID(TDerivs::dtD, i, j, k)];

        // Begin Jacobian inversion solution results:
        aux[ID(TDerivs::dtp, i, j, k)] = (d->gamma - 1)*aux[ID(TDerivs::dtE, i, j, k)]/(d->gamma*sqr(aux[ID(Aux::W, i, j, k)]) - d->gamma + 1) + (-d->gamma*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::W, i, j, k)])*aux[ID(TDerivs::dtD, i, j, k)]/(d->gamma*sqr(aux[ID(Aux::W, i, j, k)]) - d->gamma + 1);
        aux[ID(TDerivs::dtrho, i, j, k)] = (d->gamma*sqr(aux[ID(Aux::W, i, j, k)]) - d->gamma - sqr(aux[ID(Aux::W, i, j, k)]) + 1)*aux[ID(TDerivs::dtD, i, j, k)]/(d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)]) - d->gamma*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::W, i, j, k)]) + aux[ID(TDerivs::dtE, i, j, k)]/(d->gamma*sqr(aux[ID(Aux::W, i, j, k)]) - d->gamma + 1);
        aux[ID(TDerivs::dtv1, i, j, k)] = -d->gamma*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtE, i, j, k)]/(d->gamma*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + d->gamma*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)]) + (d->gamma*prims[ID(Prims::v1, i, j, k)] - prims[ID(Prims::v1, i, j, k)])*aux[ID(TDerivs::dtD, i, j, k)]/(d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] - d->gamma*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::p, i, j, k)] - d->gamma*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::rho, i, j, k)] + aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::rho, i, j, k)]) + aux[ID(TDerivs::dtS1, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)]);
        aux[ID(TDerivs::dtv2, i, j, k)] = -d->gamma*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtE, i, j, k)]/(d->gamma*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + d->gamma*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)]) + (d->gamma*prims[ID(Prims::v2, i, j, k)] - prims[ID(Prims::v2, i, j, k)])*aux[ID(TDerivs::dtD, i, j, k)]/(d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] - d->gamma*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::p, i, j, k)] - d->gamma*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::rho, i, j, k)] + aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::rho, i, j, k)]) + aux[ID(TDerivs::dtS2, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)]);
        aux[ID(TDerivs::dtv3, i, j, k)] = -d->gamma*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtE, i, j, k)]/(d->gamma*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + d->gamma*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] - d->gamma*prims[ID(Prims::p, i, j, k)] - d->gamma*prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)]) + (d->gamma*prims[ID(Prims::v3, i, j, k)] - prims[ID(Prims::v3, i, j, k)])*aux[ID(TDerivs::dtD, i, j, k)]/(d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + d->gamma*aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)] - d->gamma*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::p, i, j, k)] - d->gamma*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::rho, i, j, k)] + aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::p, i, j, k)] + aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::rho, i, j, k)]) + aux[ID(TDerivs::dtS3, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::p, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::rho, i, j, k)]);

        // Use EoS to write dtn
        aux[ID(TDerivs::dtn, i, j, k)] = aux[ID(TDerivs::dtrho, i, j, k)] - (1/(d->gamma-1))*aux[ID(TDerivs::dtp, i, j, k)];
        // Now use chain rule to write dtW(dtv1,dtv2,dtv3...)
        aux[ID(TDerivs::dtW, i, j, k)] = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*(prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]
                                             + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)]);

        // Old formulation!
        // Using chain rule for n rather than matrix-inversion result - NOT GOOD!
        //aux[ID(TDerivs::dtn, i, j, k)] = aux[ID(TDerivs::dtD, i, j, k)]/aux[ID(Aux::W, i, j, k)]
        //                                     - (prims[ID(Prims::n, i, j, k)]/aux[ID(Aux::W, i, j, k)])*aux[ID(TDerivs::dtW, i, j, k)];
        // Heavily simplified result
        //aux[ID(TDerivs::dtrho, i, j, k)] = aux[ID(TDerivs::dtE, i, j, k)];
        // Using EoS
        //aux[ID(TDerivs::dtp, i, j, k)] = (d->gamma-1)*(aux[ID(TDerivs::dtrho, i, j, k)] + aux[ID(TDerivs::dtn, i, j, k)]);
        //aux[ID(TDerivs::dtrho, i, j, k)] = aux[ID(TDerivs::dtn, i, j, k)] + (1/(d->gamma-1))*aux[ID(TDerivs::dtp, i, j, k)];

        aux[ID(TDerivs::dtT, i, j, k)] = (1/prims[ID(Prims::n, i, j, k)])*aux[ID(TDerivs::dtp, i, j, k)] 
                                             - (prims[ID(Prims::p, i, j, k)]/sqr(prims[ID(Prims::n, i, j, k)]))*aux[ID(TDerivs::dtn, i, j, k)];


        aux[ID(TDerivs::dtq1NS, i, j, k)] = -kappa*(sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v1, i, j, k)])*((aux[ID(TDerivs::dtT, i+1, j, k)] - aux[ID(TDerivs::dtT, i-1, j, k)])/(d->dx)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtT, i, j+1, k)] - aux[ID(TDerivs::dtT, i, j-1, k)])/(d->dy)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtT, i, j, k+1)] - aux[ID(TDerivs::dtT, i, j, k-1)])/(d->dz)) + 2*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv1, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv1, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)])*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtW, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtW, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtW, i, j, k)] + ((aux[ID(TDerivs::dtT, i+1, j, k)] - aux[ID(TDerivs::dtT, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtq2NS, i, j, k)] = -kappa*(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtT, i+1, j, k)] - aux[ID(TDerivs::dtT, i-1, j, k)])/(d->dx)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v2, i, j, k)])*((aux[ID(TDerivs::dtT, i, j+1, k)] - aux[ID(TDerivs::dtT, i, j-1, k)])/(d->dy)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtT, i, j, k+1)] - aux[ID(TDerivs::dtT, i, j, k-1)])/(d->dz)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + 2*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv2, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtW, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)])*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtW, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtW, i, j, k)] + ((aux[ID(TDerivs::dtT, i, j+1, k)] - aux[ID(TDerivs::dtT, i, j-1, k)])/(d->dy)));
        aux[ID(TDerivs::dtq3NS, i, j, k)] = -kappa*(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtT, i+1, j, k)] - aux[ID(TDerivs::dtT, i-1, j, k)])/(d->dx)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtT, i, j+1, k)] - aux[ID(TDerivs::dtT, i, j-1, k)])/(d->dy)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v3, i, j, k)])*((aux[ID(TDerivs::dtT, i, j, k+1)] - aux[ID(TDerivs::dtT, i, j, k-1)])/(d->dz)) + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + 2*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtW, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtW, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)])*((aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtW, i, j, k)] + ((aux[ID(TDerivs::dtT, i, j, k+1)] - aux[ID(TDerivs::dtT, i, j, k-1)])/(d->dz)));

        aux[ID(TDerivs::dtPiNS, i, j, k)] = -zeta*(aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + 0 + ((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + ((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + ((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)]);

        double C = 2.0/3.0;
        aux[ID(TDerivs::dtpi11NS, i, j, k)] = 2*eta*(C*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) - 2*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtpi12NS, i, j, k)] = -2*eta*(((aux[ID(TDerivs::dtv1, i, j+1, k)] - aux[ID(TDerivs::dtv1, i, j-1, k)])/(d->dy)) + ((aux[ID(TDerivs::dtv2, i+1, j, k)] - aux[ID(TDerivs::dtv2, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtpi13NS, i, j, k)] = -2*eta*(((aux[ID(TDerivs::dtv1, i, j, k+1)] - aux[ID(TDerivs::dtv1, i, j, k-1)])/(d->dz)) + ((aux[ID(TDerivs::dtv3, i+1, j, k)] - aux[ID(TDerivs::dtv3, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtpi22NS, i, j, k)] = 2*eta*(C*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) - 2*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)));
        aux[ID(TDerivs::dtpi23NS, i, j, k)] = -2*eta*(((aux[ID(TDerivs::dtv2, i, j, k+1)] - aux[ID(TDerivs::dtv2, i, j, k-1)])/(d->dz)) + ((aux[ID(TDerivs::dtv3, i, j+1, k)] - aux[ID(TDerivs::dtv3, i, j-1, k)])/(d->dy)));
        aux[ID(TDerivs::dtpi33NS, i, j, k)] = 2*eta*(C*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) - 2*((aux[ID(TDerivs::dtv3, i+1, j, k)] - aux[ID(TDerivs::dtv3, i-1, j, k)])/(d->dx)));
        
        // Expression with full W(t) dependence - ones above are for simplification of W=1
        // dtpi11NS = 2*eta*(C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**3*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*0 + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**3*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]**2 + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + C*prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*0 + C*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + C*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + C*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] - 2*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) - 2*prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) - 2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) - 2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)])
        // dtpi12NS = 2*eta*(C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*0 + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]**2 - aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i, j+1, k)] - aux[ID(TDerivs::dtv1, i, j-1, k)])/(d->dy)) - aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i+1, j, k)] - aux[ID(TDerivs::dtv2, i-1, j, k)])/(d->dx)) - prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) - prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) - aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i, j+1, k)] - prims[ID(Prims::v1, i, j-1, k)])/(d->dy)) - aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i+1, j, k)] - prims[ID(Prims::v2, i-1, j, k)])/(d->dx)) - ((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] - ((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv1, i, j, k)])
        // dtpi13NS = 2*eta*(C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]**2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*0 + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]**2 - aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i, j, k+1)] - aux[ID(TDerivs::dtv1, i, j, k-1)])/(d->dz)) - aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i+1, j, k)] - aux[ID(TDerivs::dtv3, i-1, j, k)])/(d->dx)) - prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) - prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) - aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i, j, k+1)] - prims[ID(Prims::v1, i, j, k-1)])/(d->dz)) - aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(d->dx)) - ((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] - ((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv1, i, j, k)])
        // dtpi22NS = 2*eta*(C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**3*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*0 + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv2, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**3*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]**2 + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + C*prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*0 + C*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + C*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + C*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] - 2*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) - 2*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) - 2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) - 2*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)])
        // dtpi23NS = 2*eta*(C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]**2*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*0 + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv2, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]**2*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]**2 - aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j, k+1)] - aux[ID(TDerivs::dtv2, i, j, k-1)])/(d->dz)) - aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j+1, k)] - aux[ID(TDerivs::dtv3, i, j-1, k)])/(d->dy)) - prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) - prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) - aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j, k+1)] - prims[ID(Prims::v2, i, j, k-1)])/(d->dz)) - aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j+1, k)] - prims[ID(Prims::v3, i, j-1, k)])/(d->dy)) - ((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] - ((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv2, i, j, k)])
        // dtpi33NS = 2*eta*(C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]**3*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv3, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**3*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*0 + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + 3*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]**2*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] + 2*C*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**3*aux[ID(TDerivs::dtW, i, j, k)]*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz)) + 2*C*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]**2*aux[ID(TDerivs::dtW, i, j, k)]**2 + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) + C*prims[ID(Prims::v1, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) + C*prims[ID(Prims::v2, i, j, k)]*((aux[ID(TDerivs::dtW, i, j+1, k)] - aux[ID(TDerivs::dtW, i, j-1, k)])/(d->dy)) + C*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i, j, k+1)] - aux[ID(TDerivs::dtW, i, j, k-1)])/(d->dz)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v1, i+1, j, k)] - prims[ID(Prims::v1, i-1, j, k)])/(d->dx)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v2, i, j+1, k)] - prims[ID(Prims::v2, i, j-1, k)])/(d->dy)) + C*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i, j, k+1)] - prims[ID(Prims::v3, i, j, k-1)])/(d->dz)) + C*0 + C*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv1, i, j, k)] + C*((aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(d->dy))*aux[ID(TDerivs::dtv2, i, j, k)] + C*((aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(d->dz))*aux[ID(TDerivs::dtv3, i, j, k)] - 2*aux[ID(Aux::W, i, j, k)]*((aux[ID(TDerivs::dtv3, i+1, j, k)] - aux[ID(TDerivs::dtv3, i-1, j, k)])/(d->dx)) - 2*prims[ID(Prims::v3, i, j, k)]*((aux[ID(TDerivs::dtW, i+1, j, k)] - aux[ID(TDerivs::dtW, i-1, j, k)])/(d->dx)) - 2*aux[ID(TDerivs::dtW, i, j, k)]*((prims[ID(Prims::v3, i+1, j, k)] - prims[ID(Prims::v3, i-1, j, k)])/(d->dx)) - 2*((aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(d->dx))*aux[ID(TDerivs::dtv3, i, j, k)])

      }
    }
  }   

}

void DEIFY::set_dtH(double * cons, double * prims, double * aux)
{
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // Now can get the state vector NS correction
        dtH[ID(0, i, j, k)] = 0;
        dtH[ID(1, i, j, k)] = ((aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q1NS, i, j, k)])*aux[ID(TDerivs::dtW, i, j, k)] + ((aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(TDerivs::dtv1, i, j, k)] + (aux[ID(Aux::q1NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtq1NS, i, j, k)] + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtq2NS, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtq3NS, i, j, k)])*prims[ID(Prims::v1, i, j, k)] + aux[ID(TDerivs::dtq1NS, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*aux[ID(TDerivs::dtv1, i, j, k)] + 2*aux[ID(Aux::PiNS, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtPiNS, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtpi11NS, i, j, k)] + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtpi12NS, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtpi13NS, i, j, k)];
        dtH[ID(2, i, j, k)] = ((aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)])*aux[ID(TDerivs::dtW, i, j, k)] + ((aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(TDerivs::dtv2, i, j, k)] + (aux[ID(Aux::q1NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtq1NS, i, j, k)] + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtq2NS, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtq3NS, i, j, k)])*prims[ID(Prims::v2, i, j, k)] + aux[ID(TDerivs::dtq2NS, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*aux[ID(TDerivs::dtv2, i, j, k)] + 2*aux[ID(Aux::PiNS, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtPiNS, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtpi12NS, i, j, k)] + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtpi22NS, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtpi23NS, i, j, k)];
        dtH[ID(3, i, j, k)] = ((aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)])*aux[ID(TDerivs::dtW, i, j, k)] + ((aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(TDerivs::dtv3, i, j, k)] + (aux[ID(Aux::q1NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + aux[ID(Aux::q2NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + aux[ID(Aux::q3NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtq1NS, i, j, k)] + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtq2NS, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtq3NS, i, j, k)])*prims[ID(Prims::v3, i, j, k)] + aux[ID(TDerivs::dtq3NS, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*aux[ID(TDerivs::dtv3, i, j, k)] + 2*aux[ID(Aux::PiNS, i, j, k)]*aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtPiNS, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtpi13NS, i, j, k)] + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtpi23NS, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtpi33NS, i, j, k)];
        dtH[ID(4, i, j, k)] = (sqr(aux[ID(Aux::W, i, j, k)]) - 1)*aux[ID(TDerivs::dtPiNS, i, j, k)] + (2*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + 2*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + 2*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(TDerivs::dtW, i, j, k)] + (2*aux[ID(Aux::q1NS, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)] + 2*aux[ID(Aux::q2NS, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + 2*aux[ID(Aux::q3NS, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)] + 2*prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtq1NS, i, j, k)] + 2*prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtq2NS, i, j, k)] + 2*prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtq3NS, i, j, k)])*aux[ID(Aux::W, i, j, k)] + 2*aux[ID(Aux::PiNS, i, j, k)]*aux[ID(Aux::W, i, j, k)]*aux[ID(TDerivs::dtW, i, j, k)] + aux[ID(TDerivs::dtpi11NS, i, j, k)] + aux[ID(TDerivs::dtpi22NS, i, j, k)] + aux[ID(TDerivs::dtpi33NS, i, j, k)];
      }
    }
  }

}


void DEIFY::set_Fx(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // Now can set the diffusion vector
        Fx[ID(0, i, j, k)] = 0;
        Fx[ID(1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::PiNS, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::pi11NS, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(2, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::pi12NS, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(4, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)] - aux[ID(Aux::PiNS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + 2*aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)];
      }
    }
  }

}

void DEIFY::set_Fy(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // Now can get the diffusion vector
        Fx[ID(0, i, j, k)] = 0;
        Fx[ID(1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::pi12NS, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(2, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::PiNS, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::pi22NS, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::pi23NS, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(4, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)] - aux[ID(Aux::PiNS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + 2*aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)];
      }
    }
  }

}

void DEIFY::set_Fz(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // Now can get the diffusion vector
        Fx[ID(0, i, j, k)] = 0;
        Fx[ID(1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::pi13NS, i, j, k)];
        Fx[ID(2, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::pi23NS, i, j, k)];
        Fx[ID(3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::PiNS, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::pi33NS, i, j, k)];
        Fx[ID(4, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)] - aux[ID(Aux::PiNS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + 2*aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
      }
    }
  }

}
