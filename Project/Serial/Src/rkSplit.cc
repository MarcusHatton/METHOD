#include "rkSplit.h"
#include <cstdio>
#include <iostream>


void RKSplit::setSource(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Set source contribution
  this->model->sourceTerm(cons, prims, aux, d->source);

  // If there is a subgrid model, set that contribution
  if (modelExtension != NULL && modelExtension->sourceExists) {
    modelExtension->sourceExtension(cons, prims, aux, d->sourceExtension);

    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            d->source[ID(var, i, j, k)] += d->sourceExtension[ID(var, i, j, k)];
          }
        }
      }
    }
  }

}

void RKSplit::step(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);
  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Predictor + source
  RK2::predictorStep(cons, prims, aux, dt);
  RK2::finalise(p1cons, p1prims, p1aux);

  // Corrector + source
  RK2::correctorStep(cons, prims, aux, dt);
  RK2::finalise(cons, prims, aux);

  // Set and add source
  this->setSource(p1cons, p1prims, p1aux);
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          cons[ID(var, i, j, k)] +=  dt * d->source[ID(var, i, j, k)];
        }
      }
    }
  }
  // RK2::finalise(cons, prims, aux);
  model->finalise(cons, prims, aux);
  RK2::finalise(cons, prims, aux);
}
