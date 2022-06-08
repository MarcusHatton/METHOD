//! Hybrid SRRMHD/SRMHD model
/*!
    This script contains the function definitions for the hybrid model.
    For detailed documentation about the methods contained herein, see hybrid.h
  and model.h.
*/

#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "hybrid.h"

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

Hybrid::Hybrid() : Model()
{
  this->Ncons = 15;
  this->Nprims = 16;
  this->Naux = 17;
}

Hybrid::Hybrid(Data * data, double tauCrossOver, double tauSpan, bool useDEIFY) : Model(data), tauCrossOver(tauCrossOver), tauSpan(tauSpan), useDEIFY(useDEIFY)
{
  // The hybrid model is basically a dissipative model in disguise, i.e. its cons
  // prims and aux are the same as SRRMHD. This model contains pointers to both
  // SRMHD and SRRMHD models, and REGIME if requested.

  dissipativeModel = new IS(data);
  idealModel = new ISCE(data);

  this->Ncons = (this->data)->Ncons = 15;
  this->Nprims = (this->data)->Nprims = 16;
  this->Naux = (this->data)->Naux = 30;

  // Allocate ideal arrays
  icons    = new double[data->Nx*data->Ny*data->Nz*idealModel->Ncons];
  iprims   = new double[data->Nx*data->Ny*data->Nz*idealModel->Nprims];
  iaux     = new double[data->Nx*data->Ny*data->Nz*idealModel->Naux];
  sicons   = new double[idealModel->Ncons];
  siprims  = new double[idealModel->Nprims];
  siaux    = new double[idealModel->Naux];
  iflux    = new double[data->Nx*data->Ny*data->Nz*idealModel->Ncons];
  dflux    = new double[data->Nx*data->Ny*data->Nz*dissipativeModel->Ncons];
  isource  = new double[idealModel->Ncons];
  dsource  = new double[dissipativeModel->Ncons];

  // 0
  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("S1");
  this->data->consLabels.push_back("S2");  this->data->consLabels.push_back("S3");
  this->data->consLabels.push_back("Tau");  this->data->consLabels.push_back("Y1");
  this->data->consLabels.push_back("Y2");  this->data->consLabels.push_back("Y3");
  this->data->consLabels.push_back("U");  this->data->consLabels.push_back("Z11");
  this->data->consLabels.push_back("Z12");  this->data->consLabels.push_back("Z13");
  this->data->consLabels.push_back("Z22");  this->data->consLabels.push_back("Z23");
  this->data->consLabels.push_back("Z33");
  // 15

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
  // 16

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
  // 30

}

Hybrid::~Hybrid()
{
  delete dissipativeModel;
  delete idealModel;
  if (useDEIFY)
  {
    delete subgridModel;
    delete[] deifySource;
    delete[] mask;
  }
  delete[] icons;
  delete[] iprims;
  delete[] iaux;
  delete[] sicons;
  delete[] siprims;
  delete[] siaux;
  delete[] iflux;
  delete[] dflux;
  delete[] isource;
  delete[] dsource;
}

void Hybrid::setupDEIFY(FluxMethod * fluxMethod)
{
  // Syntax
  Data * d(this->data);

  // Store pointer DEIFY and allocate work arrays
  if (useDEIFY)
  {
    subgridModel = new DEIFY(d, fluxMethod);
    deifySource = new double[d->Nx*d->Ny*d->Nz*idealModel->Ncons]; // 5
    mask         = new int[d->Nx*d->Ny*d->Nz];
  }
}

double Hybrid::idealWeight(double * cons, double * prims, double * aux)
{
  // Penalty function for a given cell
  return data->tauFunc(cons, prims, aux) < tauCrossOver-tauSpan ? 0 :
         data->tauFunc(cons, prims, aux) < tauCrossOver+tauSpan?
        (tanh((data->tauFunc(cons, prims, aux) - tauCrossOver) / (tauSpan/3))+1)/2 :
        1;
}

double Hybrid::idealWeightID(double * cons, double * prims, double * aux, int i, int j, int k)
{
  // Penalty function for a given cell, given global cons prims and aux
  return data->tauFunc(cons, prims, aux, i, j, k) < tauCrossOver-tauSpan ? 0 :
         data->tauFunc(cons, prims, aux, i, j, k) < tauCrossOver+tauSpan?
         (tanh((data->tauFunc(cons, prims, aux, i, j, k) - tauCrossOver) / (tauSpan/3))+1)/2 :
         1;
}

bool Hybrid::useDissipative(double * cons, double * prims, double * aux)
{
  // Should we use the Resistive C2P?
  return data->tauFunc(cons, prims, aux) < tauCrossOver;
}

void Hybrid::setIdealCPAs(double * dcons, double * dprims, double * daux)
{
  // Set the ideal cons prims and aux from the dissipative versions (single cell)
  for(int ncon(0); ncon < 5; ncon++) {
    sicons[ncon] = dcons[ncon];
  }
  // sicons[0] = dcons[0]; sicons[1] = dcons[1]; sicons[2] = dcons[2]; sicons[3] = dcons[3];
  // sicons[4] = dcons[4];

  for(int nprim(0); nprim < 6; nprim++) {
    siprims[nprim] = dprims[nprim];
  }
  // siprims[0] = dprims[0]; siprims[1] = dprims[1]; siprims[2] = dprims[2]; 
  // siprims[3] = dprims[3]; siprims[4] = dprims[4]; siprims[5] = dprims[5]; 

  for(int naux(0); naux < 4; naux++) {
    siaux[naux] = daux[naux];
  }
  for(int naux(10); naux < 20; naux++) {
    siaux[naux] = daux[naux];
  }
  // Need to be more careful here
  // siaux[0] = daux[0]; siaux[1] = daux[1]; siaux[2] = daux[2]; siaux[3] = daux[3];

  // double b0(daux[1]*(dprims[5]*dprims[1] + dprims[6]*dprims[2] + dprims[7]*dprims[3]));
  // double bx(dprims[5] / daux[1] + b0*dprims[1]);
  // double by(dprims[6] / daux[1] + b0*dprims[2]);
  // double bz(dprims[7] / daux[1] + b0*dprims[3]);
  // double bsq((dprims[5]*dprims[5] + dprims[6]*dprims[6] + dprims[7]*dprims[7] + b0*b0)/(daux[1]*daux[1]));
  // double BS(daux[1]*(dprims[5]*dcons[1] + dprims[6]*dcons[2] + dprims[7]*dcons[3]));
  // double Ssq(dcons[1]*dcons[1] + dcons[2]*dcons[2] + dcons[3]*dcons[3]);

  // siaux[4]  = b0;
  // siaux[5]  = bx;
  // siaux[6]  = by;
  // siaux[7]  = bz;
  // siaux[8]  = bsq;
  // siaux[9]  = daux[9];
  // siaux[10] = BS;
  // siaux[11] = daux[7];
  // siaux[12] =  Ssq;
}

void Hybrid::setIdealCPAsAll(double * dcons, double * dprims, double * daux)
{
  // Syntax
  Data * d(this->data);

  // Set the ideal cons prims and aux from the dissipative versions (all cells)
  for (int i(0); i < data->Nx; i++) {
    for (int j(0); j < data->Ny; j++) {
      for (int k(0); k < data->Nz; k++) {
        for(int ncon(0); ncon < 5; ncon++) {
          icons[ID(ncon, i, j, k)] = dcons[ID(ncon, i, j, k)];
        } 
        // icons[ID(0, i, j, k)] = dcons[ID(0, i, j, k)]; icons[ID(1, i, j, k)] = dcons[ID(1, i, j, k)]; icons[ID(2, i, j, k)] = dcons[ID(2, i, j, k)]; icons[ID(3, i, j, k)] = dcons[ID(3, i, j, k)];
        // icons[ID(4, i, j, k)] = dcons[ID(4, i, j, k)]; icons[ID(5, i, j, k)] = dcons[ID(5, i, j, k)]; icons[ID(6, i, j, k)] = dcons[ID(6, i, j, k)]; icons[ID(7, i, j, k)] = dcons[ID(7, i, j, k)];
        // icons[ID(8, i, j, k)] = dcons[ID(12, i, j, k)];
        for(int nprim(0); nprim < 6; nprim++) {
          iprims[ID(nprim, i, j, k)] = dprims[ID(nprim, i, j, k)];
        } 
        // iprims[ID(0, i, j, k)] = dprims[ID(0, i, j, k)]; iprims[ID(1, i, j, k)] = dprims[ID(1, i, j, k)]; iprims[ID(2, i, j, k)] = dprims[ID(2, i, j, k)]; iprims[ID(3, i, j, k)] = dprims[ID(3, i, j, k)];
        // iprims[ID(4, i, j, k)] = dprims[ID(4, i, j, k)]; iprims[ID(5, i, j, k)] = dprims[ID(5, i, j, k)]; iprims[ID(6, i, j, k)] = dprims[ID(6, i, j, k)]; iprims[ID(7, i, j, k)] = dprims[ID(7, i, j, k)];
        for(int naux(0); naux < 4; naux++) {
          iaux[ID(naux, i, j, k)] = daux[ID(naux, i, j, k)];
        }
        for(int naux(10); naux < 20; naux++) {
          iaux[ID(naux, i, j, k)] = daux[ID(naux, i, j, k)];
        }
        // iaux[ID(0, i, j, k)] = daux[ID(0, i, j, k)]; iaux[ID(1, i, j, k)] = daux[ID(1, i, j, k)]; iaux[ID(2, i, j, k)] = daux[ID(2, i, j, k)]; iaux[ID(3, i, j, k)] = daux[ID(3, i, j, k)];

        // double b0(daux[ID(1, i, j, k)]*(dprims[ID(5, i, j, k)]*dprims[ID(1, i, j, k)] + dprims[ID(6, i, j, k)]*dprims[ID(2, i, j, k)] + dprims[ID(7, i, j, k)]*dprims[ID(3, i, j, k)]));
        // double bx(dprims[ID(5, i, j, k)] / daux[ID(1, i, j, k)] + b0*dprims[ID(1, i, j, k)]);
        // double by(dprims[ID(6, i, j, k)] / daux[ID(1, i, j, k)] + b0*dprims[ID(2, i, j, k)]);
        // double bz(dprims[ID(7, i, j, k)] / daux[ID(1, i, j, k)] + b0*dprims[ID(3, i, j, k)]);
        // double bsq((dprims[ID(5, i, j, k)]*dprims[ID(5, i, j, k)] + dprims[ID(6, i, j, k)]*dprims[ID(6, i, j, k)] + dprims[ID(7, i, j, k)]*dprims[ID(7, i, j, k)] + b0*b0)/(daux[ID(1, i, j, k)]*daux[ID(1, i, j, k)]));
        // double BS(daux[ID(1, i, j, k)]*(dprims[ID(5, i, j, k)]*dcons[ID(1, i, j, k)] + dprims[ID(6, i, j, k)]*dcons[ID(2, i, j, k)] + dprims[ID(7, i, j, k)]*dcons[ID(3, i, j, k)]));
        // double Ssq(dcons[ID(1, i, j, k)]*dcons[ID(1, i, j, k)] + dcons[ID(2, i, j, k)]*dcons[ID(2, i, j, k)] + dcons[ID(3, i, j, k)]*dcons[ID(3, i, j, k)]);

        // iaux[ID(4, i, j, k)]  = b0;
        // iaux[ID(5, i, j, k)]  = bx;
        // iaux[ID(6, i, j, k)]  = by;
        // iaux[ID(7, i, j, k)]  = bz;
        // iaux[ID(8, i, j, k)]  = bsq;
        // iaux[ID(9, i, j, k)]  = daux[ID(9, i, j, k)];
        // iaux[ID(10, i, j, k)] = BS;
        // iaux[ID(11, i, j, k)] = daux[ID(7, i, j, k)];
        // iaux[ID(12, i, j, k)] =  Ssq;
      }
    }
  }
}

void Hybrid::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  // Set ideal cons prims and aux
  setIdealCPAsAll(cons, prims, aux);

  // Calculate the ideal and dissipative flux vectors
  idealModel->fluxVector(icons, iprims, iaux, iflux, dir);
  dissipativeModel->fluxVector(cons, prims, aux, dflux, dir);

  // Add dissipative contribution to hybrid->f
  for (int var(0); var < dissipativeModel->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {

          double iW = idealWeightID(cons, prims, aux, i, j, k);
          f[ID(var, i, j, k)] = (1-iW)*dflux[ID(var, i, j, k)];

        }
      }
    }
  }
  // Add ideal contribution to hybrid->f
  for (int var(0); var < idealModel->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {

          double iW = idealWeightID(cons, prims, aux, i, j, k);
          f[ID(var, i, j, k)] += iW*iflux[ID(var, i, j, k)];

        }
      }
    }
  }

}

void Hybrid::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  double iW = idealWeight(cons, prims, aux);
  // Set ideal cons prims and aux
  setIdealCPAs(cons, prims, aux);

  // Calculate the ideal and dissipative source vectors
  dissipativeModel->sourceTermSingleCell(cons, prims, aux, dsource, i, j, k);
  idealModel->sourceTermSingleCell(sicons, siprims, siaux, isource, i, j, k);

  // Do IS dissipative variables' source
  for (int var(5); var < 15; var++) {
    source[var] = (1-iW)*dsource[var];
  }

  // Do ISCE source for fluid part
  for (int var(0); var < 5; var++) {
    source[var] = iW*isource[var];
  }

  // // Add ideal contribution
  // for (int var(0); var < 8; var++) {
  //   source[var] = iW*isource[var] + (1-iW)*dsource[var];
  // }
  // // Add dissipative contribution
  // for (int var(8); var < dissipativeModel->Ncons; var++)
  //   source[var] = (1-iW)*dsource[var];
}

void Hybrid::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  // Work arrays
  double * singleCons;
  double * singlePrims;
  double * singleAux;
  double * singleSource;

  singleCons = (double *) malloc(sizeof(double) * d->Ncons);
  singlePrims = (double *) malloc(sizeof(double) * d->Nprims);
  singleAux = (double *) malloc(sizeof(double) * d->Naux);
  singleSource = (double *) malloc(sizeof(double) * d->Ncons);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Copy data to work arrays
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          singlePrims[var] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          singleAux[var] = aux[ID(var, i, j, k)];
        }

        // Get source for this cell
        this->sourceTermSingleCell(singleCons, singlePrims, singleAux, singleSource, i, j, k);
        // Copy result back
        for (int var(0); var < d->Ncons; var++) {
          source[ID(var, i, j, k)] = singleSource[var];
        }

      }
    }
  }

  // Now add DEIFY source
  if (useDEIFY)
  {
    // Calculate the ideal cons prims and aux vectors
    setIdealCPAsAll(cons, prims, aux);

    // Set the DEIFY mask
    setMasks(cons, prims, aux);
    // Calculate the REGIEME source
    subgridModel->sourceExtension(icons, iprims, iaux, deifySource);

    // Add DEIFY contribution
    for (int var(0); var < idealModel->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            double iW = idealWeightID(icons, iprims, iaux, i, j, k);
            source[ID(var, i, j, k)] += deifySource[ID(var, i, j, k)] * iW * mask[ID(0, i, j, k)];
          }
        }
      }
    }
  }

  // Free up
  free(singleCons);
  free(singlePrims);
  free(singleAux);
  free(singleSource);
}

void Hybrid::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  if (useDissipative(cons, prims, aux))
  {
    // Resistive cons2prims
    dissipativeModel->getPrimitiveVarsSingleCell(cons, prims, aux, i, j, k);
  }
  else
  {
    // Ideal cons2prims
    setIdealCPAs(cons, prims, aux);
    idealModel->getPrimitiveVarsSingleCell(sicons, siprims, siaux, i, j, k);

    for(int nprim(0); nprim < 6; nprim++) {
      prims[nprim] = siprims[nprim];
    } 
    // And convert from ideal to dissipative prims and aux
    // prims[0] = siprims[0]; prims[1] = siprims[1]; prims[2] = siprims[2];
    // prims[3] = siprims[3]; prims[4] = siprims[4]; prims[5] = siprims[5];
    // prims[6] = siprims[6]; prims[7] = siprims[7];

    for(int naux(0); naux < 4; naux++) {
      aux[naux] = siaux[naux];
    } 
    for(int naux(10); naux < 20; naux++) {
      aux[naux] = siaux[naux];
    }
    // aux[0] = siaux[0]; aux[1] = siaux[1]; aux[2] = siaux[2]; aux[3] = siaux[3];
    // aux[4] = 0; aux[5] = 0; aux[6] = 0; aux[7] = siaux[11];
    // aux[8] = cons[8]*cons[8] + cons[9]*cons[9] + cons[10]*cons[10];
    // aux[9] = siaux[9]; aux[10] = prims[0]*aux[0]*aux[1]*aux[1];
    // aux[11] = prims[1]*cons[8] + prims[2]*cons[9] + prims[3]*cons[10];
    // aux[12] = aux[10] * prims[1];
    // aux[13] = aux[10] * prims[2];
    // aux[14] = aux[10] * prims[3];
    // aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
    // aux[16] = aux[10] - prims[4] - cons[0];

  }
}

void Hybrid::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Work arrays
  double * singleCons;
  double * singlePrims;
  double * singleAux;
  singleCons = (double *) malloc(sizeof(double) * d->Ncons);
  singlePrims = (double *) malloc(sizeof(double) * d->Nprims);
  singleAux = (double *) malloc(sizeof(double) * d->Naux);


  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        // Store this cell's cons data and rhohWsq from last step
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          singlePrims[var] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          singleAux[var] = aux[ID(var, i, j, k)];
        }

        singlePrims[5] = singleCons[5]; singlePrims[6] = singleCons[6]; singlePrims[7] = singleCons[7];
        singlePrims[8] = singleCons[8]; singlePrims[9] = singleCons[9]; singlePrims[10] = singleCons[10];

        this->getPrimitiveVarsSingleCell(singleCons, singlePrims, singleAux, i, j, k);

        // Copy cell's prim and aux back to data class
        // Store this cell's cons data
        for (int var(0); var < dissipativeModel->Nprims; var++) {
          prims[ID(var, i, j, k)] = singlePrims[var];
        }
        for (int var(0); var < dissipativeModel->Naux; var++) {
          aux[ID(var, i, j, k)] = singleAux[var];
        }
      }
    }
  }

  // Free up
  free(singleCons);
  free(singlePrims);
  free(singleAux);
}

void Hybrid::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  setIdealCPAsAll(cons, prims, aux);

  // Calculate ideal and dissipative variables
  dissipativeModel->primsToAll(cons, prims, aux);
  idealModel->primsToAll(icons, iprims, iaux);

  // Compute the hybrid variables using the penalty function
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double iW = idealWeightID(cons, prims, aux, i, j, k);
        for (int var(0); var < dissipativeModel->Ncons; var++) {
          cons[ID(var, i, j, k)] *= (1-iW);
          if (var < idealModel->Ncons)
            cons[ID(var, i, j, k)] += iW*icons[ID(var, i, j, k)];
        }
        for (int var(0); var < dissipativeModel->Nprims; var++) {
          prims[ID(var, i, j, k)] *= (1-iW);
          if (var < idealModel->Nprims)
            prims[ID(var, i, j, k)] += iW*iprims[ID(var, i, j, k)];
        }
        for (int var(0); var < dissipativeModel->Naux; var++) {
          aux[ID(var, i, j, k)] *= (1-iW);
          if (var < idealModel->Naux)
            aux[ID(var, i, j, k)] += iW*iaux[ID(var, i, j, k)];
        }
      }
    }
  }
}

void Hybrid::finalise(double *cons, double *prims, double *aux)
{

  // Time derivs?

  //printf("final_step: %d", final_step);
  if (!final_step) return;

  Data * d(this->data);

  // Get timestep
  double dt=d->dt;
  //printf("dt: %.17g", dt);
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        // dW/dt \equiv du0/dt
        aux[ID(Aux::dWdt, i, j, k)] = (aux[ID(Aux::W, i, j, k)] - prev_vars[ID(0, i, j, k)])/dt;
        aux[ID(Aux::dv1dt, i, j, k)] = (prims[ID(Prims::v1, i, j, k)] - prev_vars[ID(1, i, j, k)])/dt;
        aux[ID(Aux::dv2dt, i, j, k)] = (prims[ID(Prims::v2, i, j, k)] - prev_vars[ID(2, i, j, k)])/dt;
        aux[ID(Aux::dv3dt, i, j, k)] = (prims[ID(Prims::v3, i, j, k)] - prev_vars[ID(3, i, j, k)])/dt;
        // Update previous values
        prev_vars[ID(0, i, j, k)] = aux[ID(Aux::W, i, j, k)]; 
        prev_vars[ID(1, i, j, k)] = prims[ID(Prims::v1, i, j, k)]; 
        prev_vars[ID(2, i, j, k)] = prims[ID(Prims::v2, i, j, k)];
        prev_vars[ID(3, i, j, k)] = prims[ID(Prims::v3, i, j, k)]; 
      } // End k-loop
    } // End j-loop
  } // End i-loop

}

void Hybrid::setMasks(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Assume DEIFY is not valid...
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        mask[ID(0, i, j, k)] = 0;
      }
    }
  }

  // Why were these hard-coded to 3 instead of d->Ng
  // Only loop over interior points
  int is(d->Ng); int ie(d->Nx-d->Ng);
  int js(d->Ng); int je(d->Ny-d->Ng);
  int ks(d->Ng); int ke(d->Nz-d->Ng);
  if (d->dims<3) {
    ks = 0; ke = 1;
  }
  if (d->dims<2) {
    js = 0; je = 1;
  }


  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {

        // Assume this cell's DEIFY source is valid
        bool termsPossible(true);

        // If this cell needs the DEIFY source....
        if (d->tauFunc(icons, iprims, iaux, i, j, k) > tauCrossOver-tauSpan)
        {
          // Can we compute all of the terms too? I.e. and neighbours' terms be calculated
          int nn_req {2}; // MOVE THIS SOMEWHERE BETTER - also, value?
          for (int l(-nn_req); l < nn_req; l++) {
            for (int m(-nn_req); m < nn_req; m++) {
              for (int n(-nn_req); n < nn_req; n++) {
                if (d->tauFunc(icons, iprims, iaux, i+l, j+m, k+n) < tauCrossOver-tauSpan)
                  // If this neighbour is too dissipative then we cannot calculate DEIFY for the original cell
                  termsPossible = false;
              }
            }
          }
          // And we can calculate all the terms too... set the masks
          if (termsPossible)
          {
            mask[ID(0, i, j, k)] = 1;
          }
        }
      }
    }
  }
}
