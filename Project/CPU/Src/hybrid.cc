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
  this->Ncons = 14;
  this->Nprims = 11;
  this->Naux = 17;
}

Hybrid::Hybrid(Data * data, double sigmaCrossOver, double sigmaSpan, bool useDEIFY) : Model(data), sigmaCrossOver(sigmaCrossOver), sigmaSpan(sigmaSpan), useDEIFY(useDEIFY)
{
  // The hybrid model is basically a dissipative model in disguise, i.e. its cons
  // prims and aux are the same as SRRMHD. This model contains pointers to both
  // SRMHD and SRRMHD models, and REGIME if requested.

  dissipativeModel = new IS(data);
  idealModel = new ISCE(data);

  this->Ncons = (this->data)->Ncons = 15;
  this->Nprims = (this->data)->Nprims = 16;
  this->Naux = (this->data)->Naux = 17;

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

  // Store pointer REGIME and allocate work arrays
  if (useDEIFY)
  {
    subgridModel = new DEIFY(d, fluxMethod);
    deifySource = new double[d->Nx*d->Ny*d->Nz*idealModel->Ncons];
    mask         = new int[d->Nx*d->Ny*d->Nz];
  }
}

double Hybrid::idealWeight(double * cons, double * prims, double * aux)
{
  // Penalty function for a given cell
  return data->sigmaFunc(cons, prims, aux) < sigmaCrossOver-sigmaSpan ? 0 :
         data->sigmaFunc(cons, prims, aux) < sigmaCrossOver+sigmaSpan?
        (tanh((data->sigmaFunc(cons, prims, aux) - sigmaCrossOver) / (sigmaSpan/3))+1)/2 :
        1;
}

double Hybrid::idealWeightID(double * cons, double * prims, double * aux, int i, int j, int k)
{
  // Penalty function for a given cell, given global cons prims and aux
  return data->sigmaFunc(cons, prims, aux, i, j, k) < sigmaCrossOver-sigmaSpan ? 0 :
         data->sigmaFunc(cons, prims, aux, i, j, k) < sigmaCrossOver+sigmaSpan?
         (tanh((data->sigmaFunc(cons, prims, aux, i, j, k) - sigmaCrossOver) / (sigmaSpan/3))+1)/2 :
         1;
}

bool Hybrid::useDissipative(double * cons, double * prims, double * aux)
{
  // Should we use the Resistive C2P?
  return data->sigmaFunc(cons, prims, aux) < sigmaCrossOver;
}

void Hybrid::setIdealCPAs(double * rcons, double * rprims, double * raux)
{
  // Set the ideal cons prims and aux from the dissipative versions (single cell)
  sicons[0] = rcons[0]; sicons[1] = rcons[1]; sicons[2] = rcons[2]; sicons[3] = rcons[3];
  sicons[4] = rcons[4]; sicons[5] = rcons[5]; sicons[6] = rcons[6]; sicons[7] = rcons[7];
  sicons[8] = rcons[12];

  siprims[0] = rprims[0]; siprims[1] = rprims[1]; siprims[2] = rprims[2]; siprims[3] = rprims[3];
  siprims[4] = rprims[4]; siprims[5] = rprims[5]; siprims[6] = rprims[6]; siprims[7] = rprims[7];

  siaux[0] = raux[0]; siaux[1] = raux[1]; siaux[2] = raux[2]; siaux[3] = raux[3];

  double b0(raux[1]*(rprims[5]*rprims[1] + rprims[6]*rprims[2] + rprims[7]*rprims[3]));
  double bx(rprims[5] / raux[1] + b0*rprims[1]);
  double by(rprims[6] / raux[1] + b0*rprims[2]);
  double bz(rprims[7] / raux[1] + b0*rprims[3]);
  double bsq((rprims[5]*rprims[5] + rprims[6]*rprims[6] + rprims[7]*rprims[7] + b0*b0)/(raux[1]*raux[1]));
  double BS(raux[1]*(rprims[5]*rcons[1] + rprims[6]*rcons[2] + rprims[7]*rcons[3]));
  double Ssq(rcons[1]*rcons[1] + rcons[2]*rcons[2] + rcons[3]*rcons[3]);

  siaux[4]  = b0;
  siaux[5]  = bx;
  siaux[6]  = by;
  siaux[7]  = bz;
  siaux[8]  = bsq;
  siaux[9]  = raux[9];
  siaux[10] = BS;
  siaux[11] = raux[7];
  siaux[12] =  Ssq;
}

void Hybrid::setIdealCPAsAll(double * rcons, double * rprims, double * raux)
{
  // Syntax
  Data * d(this->data);

  // Set the ideal cons prims and aux from the dissipative versions (all cells)
  for (int i(0); i < data->Nx; i++) {
    for (int j(0); j < data->Ny; j++) {
      for (int k(0); k < data->Nz; k++) {
        icons[ID(0, i, j, k)] = rcons[ID(0, i, j, k)]; icons[ID(1, i, j, k)] = rcons[ID(1, i, j, k)]; icons[ID(2, i, j, k)] = rcons[ID(2, i, j, k)]; icons[ID(3, i, j, k)] = rcons[ID(3, i, j, k)];
        icons[ID(4, i, j, k)] = rcons[ID(4, i, j, k)]; icons[ID(5, i, j, k)] = rcons[ID(5, i, j, k)]; icons[ID(6, i, j, k)] = rcons[ID(6, i, j, k)]; icons[ID(7, i, j, k)] = rcons[ID(7, i, j, k)];
        icons[ID(8, i, j, k)] = rcons[ID(12, i, j, k)];

        iprims[ID(0, i, j, k)] = rprims[ID(0, i, j, k)]; iprims[ID(1, i, j, k)] = rprims[ID(1, i, j, k)]; iprims[ID(2, i, j, k)] = rprims[ID(2, i, j, k)]; iprims[ID(3, i, j, k)] = rprims[ID(3, i, j, k)];
        iprims[ID(4, i, j, k)] = rprims[ID(4, i, j, k)]; iprims[ID(5, i, j, k)] = rprims[ID(5, i, j, k)]; iprims[ID(6, i, j, k)] = rprims[ID(6, i, j, k)]; iprims[ID(7, i, j, k)] = rprims[ID(7, i, j, k)];

        iaux[ID(0, i, j, k)] = raux[ID(0, i, j, k)]; iaux[ID(1, i, j, k)] = raux[ID(1, i, j, k)]; iaux[ID(2, i, j, k)] = raux[ID(2, i, j, k)]; iaux[ID(3, i, j, k)] = raux[ID(3, i, j, k)];

        double b0(raux[ID(1, i, j, k)]*(rprims[ID(5, i, j, k)]*rprims[ID(1, i, j, k)] + rprims[ID(6, i, j, k)]*rprims[ID(2, i, j, k)] + rprims[ID(7, i, j, k)]*rprims[ID(3, i, j, k)]));
        double bx(rprims[ID(5, i, j, k)] / raux[ID(1, i, j, k)] + b0*rprims[ID(1, i, j, k)]);
        double by(rprims[ID(6, i, j, k)] / raux[ID(1, i, j, k)] + b0*rprims[ID(2, i, j, k)]);
        double bz(rprims[ID(7, i, j, k)] / raux[ID(1, i, j, k)] + b0*rprims[ID(3, i, j, k)]);
        double bsq((rprims[ID(5, i, j, k)]*rprims[ID(5, i, j, k)] + rprims[ID(6, i, j, k)]*rprims[ID(6, i, j, k)] + rprims[ID(7, i, j, k)]*rprims[ID(7, i, j, k)] + b0*b0)/(raux[ID(1, i, j, k)]*raux[ID(1, i, j, k)]));
        double BS(raux[ID(1, i, j, k)]*(rprims[ID(5, i, j, k)]*rcons[ID(1, i, j, k)] + rprims[ID(6, i, j, k)]*rcons[ID(2, i, j, k)] + rprims[ID(7, i, j, k)]*rcons[ID(3, i, j, k)]));
        double Ssq(rcons[ID(1, i, j, k)]*rcons[ID(1, i, j, k)] + rcons[ID(2, i, j, k)]*rcons[ID(2, i, j, k)] + rcons[ID(3, i, j, k)]*rcons[ID(3, i, j, k)]);

        iaux[ID(4, i, j, k)]  = b0;
        iaux[ID(5, i, j, k)]  = bx;
        iaux[ID(6, i, j, k)]  = by;
        iaux[ID(7, i, j, k)]  = bz;
        iaux[ID(8, i, j, k)]  = bsq;
        iaux[ID(9, i, j, k)]  = raux[ID(9, i, j, k)];
        iaux[ID(10, i, j, k)] = BS;
        iaux[ID(11, i, j, k)] = raux[ID(7, i, j, k)];
        iaux[ID(12, i, j, k)] =  Ssq;
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
  for (int var(0); var < 8; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {

          double iW = idealWeightID(cons, prims, aux, i, j, k);
          f[ID(var, i, j, k)] += iW*iflux[ID(var, i, j, k)];

        }
      }
    }
  }
  // And the divergence cleaning part separately
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        double iW = idealWeightID(cons, prims, aux, i, j, k);
        f[ID(12, i, j, k)] += iW*iflux[ID(8, i, j, k)];

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

  // Add ideal contribution
  for (int var(0); var < 8; var++) {
    source[var] = iW*isource[var] + (1-iW)*dsource[var];
  }
  // Add dissipative contribution
  for (int var(8); var < dissipativeModel->Ncons; var++)
    source[var] = (1-iW)*dsource[var];
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

  // Now add REGIME source
  if (useDEIFY)
  {
    // Calculate the ideal cons prims and aux vectors
    setIdealCPAsAll(cons, prims, aux);

    // Set the REGIME mask
    setMasks(cons, prims, aux);
    // Calculate the REGIEME source
    subgridModel->sourceExtension(icons, iprims, iaux, deifySource);

    // Add REGIME contribution
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

    // And convert from ideal to dissipative prims and aux
    prims[0] = siprims[0]; prims[1] = siprims[1]; prims[2] = siprims[2];
    prims[3] = siprims[3]; prims[4] = siprims[4]; prims[5] = siprims[5];
    prims[6] = siprims[6]; prims[7] = siprims[7];

    aux[0] = siaux[0]; aux[1] = siaux[1]; aux[2] = siaux[2]; aux[3] = siaux[3];
    aux[4] = 0; aux[5] = 0; aux[6] = 0; aux[7] = siaux[11];
    aux[8] = cons[8]*cons[8] + cons[9]*cons[9] + cons[10]*cons[10];
    aux[9] = siaux[9]; aux[10] = prims[0]*aux[0]*aux[1]*aux[1];
    aux[11] = prims[1]*cons[8] + prims[2]*cons[9] + prims[3]*cons[10];
    aux[12] = aux[10] * prims[1];
    aux[13] = aux[10] * prims[2];
    aux[14] = aux[10] * prims[3];
    aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
    aux[16] = aux[10] - prims[4] - cons[0];

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
  // Syntax
  Data * d(this->data);

  // Set ideal electric fields
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double iW = idealWeightID(cons, prims, aux, i, j, k);
        double iEx = -(prims[ID(2, i, j, k)]*prims[ID(7, i, j, k)] - prims[ID(3, i, j, k)]*prims[ID(6, i, j, k)]);
        double iEy = -(prims[ID(3, i, j, k)]*prims[ID(5, i, j, k)] - prims[ID(1, i, j, k)]*prims[ID(7, i, j, k)]);
        double iEz = -(prims[ID(1, i, j, k)]*prims[ID(6, i, j, k)] - prims[ID(2, i, j, k)]*prims[ID(5, i, j, k)]);

        cons[ID(8, i, j, k)]  *= (1-iW);
        cons[ID(9, i, j, k)]  *= (1-iW);
        cons[ID(10, i, j, k)] *= (1-iW);

        cons[ID(8, i, j, k)]  += iW*iEx;
        cons[ID(9, i, j, k)]  += iW*iEy;
        cons[ID(10, i, j, k)] += iW*iEz;

        prims[ID(8, i, j, k)]  = cons[ID(8, i, j, k)];
        prims[ID(9, i, j, k)]  = cons[ID(9, i, j, k)];
        prims[ID(10, i, j, k)] = cons[ID(10, i, j, k)];

      }
    }
  }
}

void Hybrid::setMasks(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Assume REGIME is not valid...
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        mask[ID(0, i, j, k)] = 0;
      }
    }
  }

  // Only loop over interior points
  int is(3); int ie(d->Nx-3);
  int js(3); int je(d->Ny-3);
  int ks(3); int ke(d->Nz-3);
  if (d->dims<3) {
    ks = 0; ke = 1;
  }
  if (d->dims<2) {
    js = 0; je = 1;
  }


  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {

        // Assume this cell's REGIME source is valid
        bool termsPossible(true);

        // If this cell needs the REGIME source....
        if (d->sigmaFunc(icons, iprims, iaux, i, j, k) > sigmaCrossOver-sigmaSpan)
        {
          // Can we compute all of the terms too? I.e. and neighbours' terms be calculated
          for (int l(-3); l < 3; l++) {
            for (int m(-3); m < 3; m++) {
              for (int n(-3); n < 3; n++) {
                if (d->sigmaFunc(icons, iprims, iaux, i+l, j+m, k+n) < sigmaCrossOver-sigmaSpan)
                  // If this neighbour is too dissipative then we cannot calculate REGIME for the original cell
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
