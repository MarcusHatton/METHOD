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

void DEIFY::set_vars(double * cons, double * prims, double * aux)
{
  Data * d(this->data);

  double kappa = this->data->optionalSimArgs[0];
//  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
//  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
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

/*
  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
*/

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        /*

        dxT = (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx);
        dyT = (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy);
        dzT = (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz);

        dxux = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx);
        dyuy = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        dzuz = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);

        dxuy = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx);
        dxuz = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx);
        dyux = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        dyuz = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy);
        dzux = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);
        dzuy = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz);

        */

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
  
        // pi^l_j,NS 14 - STILL NOT FULLY CORRECT - NEED MORE h_munu factors in front of big bracket!
        // 11
        aux[ID(Aux::pi11NS, i, j, k)] = -2*eta*( 2*dxux 
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 12
        aux[ID(Aux::pi12NS, i, j, k)] = -2*eta*( dxuy + dyux
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 13
        aux[ID(Aux::pi13NS, i, j, k)] = -2*eta*( dxuz + dzux
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 22
        aux[ID(Aux::pi22NS, i, j, k)] = -2*eta*( 2*dyuy
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 23
        aux[ID(Aux::pi23NS, i, j, k)] = -2*eta*( dyuz + dzuy
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );
        // 33
        aux[ID(Aux::pi33NS, i, j, k)] = -2*eta*( 2*dzuz
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::Theta, i, j, k)] );

        aux[ID(TDerivs::dtE, i, j, k)] = aux[ID(TDerivs::dtTau, i, j, k)] + aux[ID(TDerivs::dtD, i, j, k)];

        // Write out time-derivative forms for primitives and dissNS variables in terms of spatial calcs
        aux[ID(TDerivs::dtv1, i, j, k)] = -prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtE, i, j, k)]/(prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)]) +  aux[ID(TDerivs::dtS1, i, j, k)]/((prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*sqr(aux[ID(Aux::W, i, j, k)]));
        aux[ID(TDerivs::dtv2, i, j, k)] = -prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtE, i, j, k)]/(prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)]) + aux[ID(TDerivs::dtS2, i, j, k)]/((prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*sqr(aux[ID(Aux::W, i, j, k)]));
        aux[ID(TDerivs::dtv3, i, j, k)] = -prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtE, i, j, k)]/(prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)]) + aux[ID(TDerivs::dtS3, i, j, k)]/((prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::rho, i, j, k)])*sqr(aux[ID(Aux::W, i, j, k)]));

        aux[ID(TDerivs::dtW, i, j, k)] = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]*(prims[ID(Prims::v1, i, j, k)]*aux[ID(TDerivs::dtv1, i, j, k)]
                                             + prims[ID(Prims::v2, i, j, k)]*aux[ID(TDerivs::dtv2, i, j, k)] + prims[ID(Prims::v3, i, j, k)]*aux[ID(TDerivs::dtv3, i, j, k)]);

        // Using chain rule rather than matrix-inversion result
        aux[ID(TDerivs::dtn, i, j, k)] = aux[ID(TDerivs::dtD, i, j, k)]/aux[ID(Aux::W, i, j, k)]
                                             - (prims[ID(Prims::n, i, j, k)]/aux[ID(Aux::W, i, j, k)])*aux[ID(TDerivs::dtW, i, j, k)];
        
        aux[ID(TDerivs::dtrho, i, j, k)] = aux[ID(TDerivs::dtE, i, j, k)];
        aux[ID(TDerivs::dtp, i, j, k)] = (d->gamma-1)*(aux[ID(TDerivs::dtrho, i, j, k)] + aux[ID(TDerivs::dtn, i, j, k)]);

        // aux[ID(TDerivs::dtp, i, j, k)] = aux[ID(TDerivs::dtS3, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]) + aux[ID(TDerivs::dtS2, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]) + aux[ID(TDerivs::dtS1, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]) + aux[ID(TDerivs::dtE, i, j, k)]/(sqr(aux[ID(Aux::W, i, j, k)]) - 1);
        // Using EoS
        // aux[ID(TDerivs::dtrho, i, j, k)] = aux[ID(TDerivs::dtn, i, j, k)] + (1/(d->gamma-1))*aux[ID(TDerivs::dtp, i, j, k)];
        
        
        aux[ID(TDerivs::dtT, i, j, k)] = (1/prims[ID(Prims::n, i, j, k)])*aux[ID(TDerivs::dtp, i, j, k)] 
                                             - (prims[ID(Prims::p, i, j, k)]/sqr(prims[ID(Prims::n, i, j, k)]))*aux[ID(TDerivs::dtn, i, j, k)];

        // Heavily simplified expressions for now (actually, just low-velocity limit... v**2 = 0, W = 1)      
        aux[ID(TDerivs::dtq1NS, i, j, k)] = -kappa*((aux[ID(TDerivs::dtT, i+1, j, k)] - aux[ID(TDerivs::dtT, i-1, j, k)])/(d->dx));
        aux[ID(TDerivs::dtq2NS, i, j, k)] = -kappa*((aux[ID(TDerivs::dtT, i, j+1, k)] - aux[ID(TDerivs::dtT, i, j-1, k)])/(d->dy));
        aux[ID(TDerivs::dtq3NS, i, j, k)] = -kappa*((aux[ID(TDerivs::dtT, i, j, k+1)] - aux[ID(TDerivs::dtT, i, j, k-1)])/(d->dz));

        aux[ID(TDerivs::dtPiNS, i, j, k)] = -zeta*(((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + ((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + ((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)));

        double C = 2.0/3.0;

        aux[ID(TDerivs::dtpi11NS, i, j, k)] = 2*eta*(C*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) - 2*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtpi12NS, i, j, k)] = -2*eta*(((aux[ID(TDerivs::dtv1, i, j+1, k)] - aux[ID(TDerivs::dtv1, i, j-1, k)])/(d->dy)) + ((aux[ID(TDerivs::dtv2, i+1, j, k)] - aux[ID(TDerivs::dtv2, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtpi13NS, i, j, k)] = -2*eta*(((aux[ID(TDerivs::dtv1, i, j, k+1)] - aux[ID(TDerivs::dtv1, i, j, k-1)])/(d->dz)) + ((aux[ID(TDerivs::dtv3, i+1, j, k)] - aux[ID(TDerivs::dtv3, i-1, j, k)])/(d->dx)));
        aux[ID(TDerivs::dtpi22NS, i, j, k)] = 2*eta*(C*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) - 2*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)));
        aux[ID(TDerivs::dtpi23NS, i, j, k)] = -2*eta*(((aux[ID(TDerivs::dtv2, i, j, k+1)] - aux[ID(TDerivs::dtv2, i, j, k-1)])/(d->dz)) + ((aux[ID(TDerivs::dtv3, i, j+1, k)] - aux[ID(TDerivs::dtv3, i, j-1, k)])/(d->dy)));
        aux[ID(TDerivs::dtpi33NS, i, j, k)] = 2*eta*(C*((aux[ID(TDerivs::dtv1, i+1, j, k)] - aux[ID(TDerivs::dtv1, i-1, j, k)])/(d->dx)) + C*((aux[ID(TDerivs::dtv2, i, j+1, k)] - aux[ID(TDerivs::dtv2, i, j-1, k)])/(d->dy)) + C*((aux[ID(TDerivs::dtv3, i, j, k+1)] - aux[ID(TDerivs::dtv3, i, j, k-1)])/(d->dz)) - 2*((aux[ID(TDerivs::dtv3, i+1, j, k)] - aux[ID(TDerivs::dtv3, i-1, j, k)])/(d->dx)));

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
        Fx[ID(1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::PiNS, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(2, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(4, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)] - aux[ID(Aux::PiNS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*sqr(prims[ID(Prims::v1, i, j, k)]) + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)];
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
        Fx[ID(1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(2, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::PiNS, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        Fx[ID(4, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)] - aux[ID(Aux::PiNS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*sqr(prims[ID(Prims::v2, i, j, k)]) + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)];
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
        Fx[ID(1, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]);
        Fx[ID(2, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)] + aux[ID(Aux::pi12NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]);
        Fx[ID(3, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::PiNS, i, j, k)] + 2*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi13NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi23NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]);
        Fx[ID(4, i, j, k)] = aux[ID(Aux::PiNS, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)] - aux[ID(Aux::PiNS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q1NS, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q2NS, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)]*sqr(prims[ID(Prims::v3, i, j, k)]) + aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::q3NS, i, j, k)] + aux[ID(Aux::pi11NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi22NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)] + aux[ID(Aux::pi33NS, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
      }
    }
  }

}
