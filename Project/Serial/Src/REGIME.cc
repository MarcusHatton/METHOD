#include "REGIME.h"
#include <cstdio>
#include <cmath>
#include <ctime>

/*
    This script defines the functions for REGIME. More information can be found
  in the documentation for the REGIME class. Alternatively, you can look in the
  REGIME.h header file.
    Please see Wright & Hawke 2019 - A resistive extension for ideal MHD for
  more.
*/

#define ID(variable, idx, jdx, kdx)  ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

// dwdsb
#define IDWS(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
// dfxdw, dfydw, dfzdw
#define IDFW(ldx, mdx, idx, jdx, kdx)  ((ldx)*(12)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
// Mx, My, and Mz matrix
#define IDM(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

REGIME::REGIME(Data * data, FluxMethod * fluxMethod) : ModelExtension(data), fluxMethod(fluxMethod)
{
  //  Syntax
  Data * d(this->data);

  sourceExists = true;

  // Allocate arrays
  dfxdw = new double[d->Nx*d->Ny*d->Nz*8*12] ();
  dfydw = new double[d->Nx*d->Ny*d->Nz*8*12] ();
  dfzdw = new double[d->Nx*d->Ny*d->Nz*8*12] ();
  dwdsb = new double[d->Nx*d->Ny*d->Nz*12*3] ();
  E = new double[d->Nx*d->Ny*d->Nz*3] ();
  q = new double[d->Nx*d->Ny*d->Nz] ();
  K = new double[d->Nx*d->Ny*d->Nz*3] ();
  Mx = new double[d->Nx*d->Ny*d->Nz*8*3] ();
  My = new double[d->Nx*d->Ny*d->Nz*8*3] ();
  Mz = new double[d->Nx*d->Ny*d->Nz*8*3] ();
  fbx = new double[d->Nx*d->Ny*d->Nz*3] ();
  fby = new double[d->Nx*d->Ny*d->Nz*3] ();
  fbz = new double[d->Nx*d->Ny*d->Nz*3] ();
  diffuX = new double[d->Nx*d->Ny*d->Nz*8] ();
  diffuY = new double[d->Nx*d->Ny*d->Nz*8] ();
  diffuZ = new double[d->Nx*d->Ny*d->Nz*8] ();
  alpha = new double[d->Nx*d->Ny*d->Nz] ();
  d->sourceExtension = new double[d->Nx*d->Ny*d->Nz*d->Ncons] ();
  // vortX = new double[d->Nx*d->Ny*d->Nz] ();
  // vortY = new double[d->Nx*d->Ny*d->Nz] ();

}

REGIME::~REGIME()
{
  //  Syntax
  Data * d(this->data);

  delete[] dfxdw;
  delete[] dfydw;
  delete[] dfzdw;
  delete[] dwdsb;
  delete[] E;
  delete[] q;
  delete[] K;
  delete[] Mx;
  delete[] My;
  delete[] Mz;
  delete[] fbx;
  delete[] fby;
  delete[] fbz;
  delete[] diffuX;
  delete[] diffuY;
  delete[] diffuZ;
  delete[] alpha;
  delete[] d->sourceExtension;
  // delete[] vortX;
  // delete[] vortY;
}

// void REGIME::init(double * prims)
// {
//   // Syntax
//   Data * d(this->data);
//
//   // X-direction vorticity
//   maxCurVortX = -1;
//   for (int i(1); i < d->Nx-1; i++) {
//     for (int j(0); j < d->Ny; j++) {
//       for (int k(0); k < d->Nz; k++) {
//         vortX[ID(0, i, j, k)] = (prims[ID(2, i+1, j, k)] - prims[ID(2, i-1, j, k)]) / (2 * d->dx);
//         if (fabs(vortX[ID(0, i, j, k)]) > maxCurVortX)
//           maxCurVortX = fabs(vortX[ID(0, i, j, k)]);
//       }
//     }
//   }
//   // Y-direction vorticity
//   maxCurVortY = -1;
//   for (int i(0); i < d->Nx; i++) {
//     for (int j(1); j < d->Ny-1; j++) {
//       for (int k(0); k < d->Nz; k++) {
//         vortY[ID(0, i, j, k)] = (prims[ID(1, i, j+1, k)] - prims[ID(1, i, j-1, k)]) / (2 * d->dy);
//         if (fabs(vortY[ID(0, i, j, k)]) > maxCurVortY)
//           maxCurVortY = fabs(vortY[ID(0, i, j, k)]);
//       }
//     }
//   }
// }

// void REGIME::set_vorts(double * prims)
// {
//   // Syntax
//   Data * d(this->data);
//
//   // Update previous maximums
//   maxPrevVortX = maxCurVortX;
//   maxPrevVortY = maxCurVortY;
//
//   // X-direction vorticity
//   maxCurVortX = -1;
//   for (int i(1); i < d->Nx-1; i++) {
//     for (int j(0); j < d->Ny; j++) {
//       for (int k(0); k < d->Nz; k++) {
//         vortX[ID(0, i, j, k)] = (prims[ID(2, i+1, j, k)] - prims[ID(2, i-1, j, k)]) / (2 * d->dx);
//         if (fabs(vortX[ID(0, i, j, k)]) > maxCurVortX)
//           maxCurVortX = fabs(vortX[ID(0, i, j, k)]);
//       }
//     }
//   }
//   // Y-direction vorticity
//   maxCurVortY = -1;
//   for (int i(0); i < d->Nx; i++) {
//     for (int j(1); j < d->Ny-1; j++) {
//       for (int k(0); k < d->Nz; k++) {
//         vortY[ID(0, i, j, k)] = (prims[ID(1, i, j+1, k)] - prims[ID(1, i, j-1, k)]) / (2 * d->dy);
//         if (fabs(vortY[ID(0, i, j, k)]) > maxCurVortY)
//           maxCurVortY = fabs(vortY[ID(0, i, j, k)]);
//       }
//     }
//   }
// }

void REGIME::sourceExtension(double * cons, double * prims, double * aux, double * source)
{
  // Syntax
  Data * d(this->data);

  // Ahead and behind difference for MINMOD
  double ahead(0);
  double behind(0);
  // Tolerance for MINMOD switch
  // double tol(0.001);

  // Set vars
  this->set_vars(cons, prims, aux);
  // Ensure K and dwdsb are set
  this->set_K(cons, prims, aux);
  this->set_dwdsb(cons, prims, aux);

  // if (d->dims > 1)
  //   this->set_vorts(prims);

  // Get Da, determine its gradient and add to source
  {
    // Determine the diffusion vectors
    this->set_Dx(cons, prims, aux);
    for (int var(0); var<8; var++) {
      for (int i(1); i<d->Nx-1; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0);k<d->Nz; k++) {

            // // Mindmod
            // ahead =  (diffuX[ID(var, i+1, j, k)] - diffuX[ID(var, i, j, k)]) / (d->dx);
            // behind = (diffuX[ID(var, i, j, k)] - diffuX[ID(var, i-1, j, k)]) / (d->dx);
            //
            // if (ahead * behind > 0){
            //   if (fabs(ahead) < fabs(behind))
            //     source[ID(var, i, j, k)] = ahead;
            //   else
            //     source[ID(var, i, j, k)] = behind;
            // }
            // else
            //   source[ID(var, i, j, k)] = 0;//

            // Hybrid
            // ahead =  (diffuX[ID(var, i+1, j, k)] - diffuX[ID(var, i, j, k)]) / (d->dx);
            // behind = (diffuX[ID(var, i, j, k)] - diffuX[ID(var, i-1, j, k)]) / (d->dx);
            //
            // if (ahead * behind > 0){
            //   if (fabs(ahead) < fabs(behind))
            //     source[ID(var, i, j, k)] = ahead;
            //   else
            //     source[ID(var, i, j, k)] = behind;
            // }
            // else
            //   source[ID(var, i, j, k)] = (diffuX[ID(var, i+1, j, k)] - diffuX[ID(var, i-1, j, k)]) / (2 * d->dx);


            // // Hybrid v2
            ahead =  (diffuX[ID(var, i+1, j, k)] - diffuX[ID(var, i-1, j, k)]) / (2 * d->dx);
            behind = (diffuX[ID(var, i, j, k)] - diffuX[ID(var, i-2, j, k)]) / (2 * d->dx);

            if (ahead * behind > 0){
              if (fabs(ahead) < fabs(behind))
                source[ID(var, i, j, k)] = ahead;
              else
                source[ID(var, i, j, k)] = behind;
            }
            else
              source[ID(var, i, j, k)] = ahead;

              // // Hybrid v2
              // ahead =  (diffuX[ID(var, i+1, j, k)] - diffuX[ID(var, i-1, j, k)]) / (2 * d->dx);
              // behind = (diffuX[ID(var, i, j, k)] - diffuX[ID(var, i-2, j, k)]) / (2 * d->dx);
              //
              // if (ahead * behind > 0){
              //   if (fabs(ahead) < fabs(behind))
              //     source[ID(var, i, j, k)] = ahead;
              //   else
              //     source[ID(var, i, j, k)] = behind;
              // }
              // else
              //   source[ID(var, i, j, k)] = 0;

          }
        }
      }
    }
    if (d->dims>1)
    {
      this->set_Dy(cons, prims, aux);
      for (int var(0); var<8; var++) {
        for (int i(0); i<d->Nx; i++) {
          for (int j(1); j<d->Ny-1; j++) {
            for (int k(0);k<d->Nz; k++) {

              // // Mindmod
              // ahead =  (diffuY[ID(var, i, j+1, k)] - diffuY[ID(var, i, j, k)]) / (d->dy);
              // behind = (diffuY[ID(var, i, j, k)] - diffuY[ID(var, i, j-1, k)]) / (d->dy);
              //
              // if (ahead * behind > 0){
              //   if (fabs(ahead) < fabs(behind))
              //     source[ID(var, i, j, k)] += ahead;
              //   else
              //     source[ID(var, i, j, k)] += behind;
              // }
              // else
              //   source[ID(var, i, j, k)] += 0;

              // // Hybrid
              // ahead =  (diffuY[ID(var, i, j+1, k)] - diffuY[ID(var, i, j, k)]) / (d->dy);
              // behind = (diffuY[ID(var, i, j, k)] - diffuY[ID(var, i, j-1, k)]) / (d->dy);
              //
              // if (ahead * behind > 0){
              //   if (fabs(ahead) < fabs(behind))
              //     source[ID(var, i, j, k)] += ahead;
              //   else
              //     source[ID(var, i, j, k)] += behind;
              // }
              // else
              //   source[ID(var, i, j, k)] += (diffuY[ID(var, i, j+1, k)] - diffuY[ID(var, i, j-1, k)]) / (2 * d->dy);

              // // Hybrid v2
              ahead =  (diffuY[ID(var, i, j+1, k)] - diffuY[ID(var, i, j-1, k)]) / (2 * d->dy);
              behind = (diffuY[ID(var, i, j, k)] - diffuY[ID(var, i, j-2, k)]) / (2 * d->dy);

              if (ahead * behind > 0){
                if (fabs(ahead) < fabs(behind))
                  source[ID(var, i, j, k)] += ahead;
                else
                  source[ID(var, i, j, k)] += behind;
              }
              else
                source[ID(var, i, j, k)] += ahead;

                // // Hybrid v3
                // ahead =  (diffuY[ID(var, i, j+1, k)] - diffuY[ID(var, i, j-1, k)]) / (2 * d->dy);
                // behind = (diffuY[ID(var, i, j, k)] - diffuY[ID(var, i, j-2, k)]) / (2 * d->dy);
                //
                // if (ahead * behind > 0){
                //   if (fabs(ahead) < fabs(behind))
                //     source[ID(var, i, j, k)] += ahead;
                //   else
                //     source[ID(var, i, j, k)] += behind;
                // }
                // else
                //   source[ID(var, i, j, k)] += 0;
            }
          }
        }
      }
    }
    if (d->dims==3)
    {
      this->set_Dz(cons, prims, aux);
      for (int var(0); var<8; var++) {
        for (int i(0); i<d->Nx; i++) {
          for (int j(0); j<d->Ny; j++) {
            for (int k(1);k<d->Nz-1; k++) {
              source[ID(var, i, j, k)] += (diffuZ[ID(var, i, j, k+1)] - diffuZ[ID(var, i, j, k-1)]) / (2*d->dz);
            }
          }
        }
      }
    }
  }


}

void REGIME::set_vars(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // First set E = - v cross B, zero the K and diffuXYZ vectors, the Mxyz matrices, and set the prefactor
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        E[ID(0, i, j, k)] = - (prims[ID(2, i, j, k)] * prims[ID(7, i, j, k)] - prims[ID(3, i, j, k)] * prims[ID(6, i, j, k)]);
        E[ID(1, i, j, k)] = - (prims[ID(3, i, j, k)] * prims[ID(5, i, j, k)] - prims[ID(1, i, j, k)] * prims[ID(7, i, j, k)]);
        E[ID(2, i, j, k)] = - (prims[ID(1, i, j, k)] * prims[ID(6, i, j, k)] - prims[ID(2, i, j, k)] * prims[ID(5, i, j, k)]);

        K[ID(0, i, j, k)] = 0.0;
        K[ID(1, i, j, k)] = 0.0;
        K[ID(2, i, j, k)] = 0.0;

        diffuX[ID(0, i, j, k)] = diffuX[ID(1, i, j, k)] = diffuX[ID(2, i, j, k)] = diffuX[ID(3, i, j, k)] = diffuX[ID(4, i, j, k)] = diffuX[ID(5, i, j, k)] =
        diffuX[ID(6, i, j, k)] = diffuX[ID(7, i, j, k)] = 0.0;
        diffuY[ID(0, i, j, k)] = diffuY[ID(1, i, j, k)] = diffuY[ID(2, i, j, k)] = diffuY[ID(3, i, j, k)] = diffuY[ID(4, i, j, k)] = diffuY[ID(5, i, j, k)] =
        diffuY[ID(6, i, j, k)] = diffuY[ID(7, i, j, k)] = 0.0;
        diffuZ[ID(0, i, j, k)] = diffuZ[ID(1, i, j, k)] = diffuZ[ID(2, i, j, k)] = diffuZ[ID(3, i, j, k)] = diffuZ[ID(4, i, j, k)] = diffuZ[ID(5, i, j, k)] =
        diffuZ[ID(6, i, j, k)] = diffuZ[ID(7, i, j, k)] = 0.0;

        for (int l(0); l<8; l++) {
          for (int m(0); m<3; m++) {
            Mx[IDM(l, m, i, j, k)] = 0.0;
            My[IDM(l, m, i, j, k)] = 0.0;
            Mz[IDM(l, m, i, j, k)] = 0.0;
          }
        }

        alpha[ID(0, i, j, k)] = 1 / ((q[ID(0, i, j, k)]*q[ID(0, i, j, k)] + d->sigma*d->sigma)*((q[ID(0, i, j, k)]*q[ID(0, i, j, k)] + d->sigma*d->sigma) + d->sigma*d->sigma*(prims[ID(5, i, j, k)]*prims[ID(5, i, j, k)] +
                   prims[ID(6, i, j, k)]*prims[ID(6, i, j, k)] + prims[ID(7, i, j, k)]*prims[ID(7, i, j, k)])));
      }
    }
  }

  // Charge density, q, is the gradient of the electric field. Central differencing
  if (d->dims == 3) {
    for (int i(2); i<d->Nx-2; i++) {
      for (int j(2); j<d->Ny-2; j++) {
        for (int k(2); k<d->Nz-2; k++) {
          q[ID(0, i, j, k)] = (-E[ID(0, i+2, j, k)] + 8*E[ID(0, i+1, j, k)] - 8*E[ID(0, i-1, j, k)] + E[ID(0, i-2, j, k)]) / (12*d->dx) +
                              (-E[ID(1, i, j+2, k)] + 8*E[ID(1, i, j+1, k)] - 8*E[ID(1, i, j-1, k)] + E[ID(1, i, j-2, k)]) / (12*d->dy) +
                              (-E[ID(2, i, j, k+2)] + 8*E[ID(2, i, j, k+1)] - 8*E[ID(2, i, j, k-1)] + E[ID(2, i, j, k-2)]) / (12*d->dz);
        }
      }
    }
  }
  else if (d->dims == 2) {
    for (int i(2); i<d->Nx-2; i++) {
      for (int j(2); j<d->Ny-2; j++) {
        q[ID(0, i, j, 0)] = (-E[ID(0, i+2, j, 0)] + 8*E[ID(0, i+1, j, 0)] - 8*E[ID(0, i-1, j, 0)] + E[ID(0, i-2, j, 0)]) / (12*d->dx) +
                            (-E[ID(1, i, j+2, 0)] + 8*E[ID(1, i, j+1, 0)] - 8*E[ID(1, i, j-1, 0)] + E[ID(1, i, j-2, 0)]) / (12*d->dy);
      }
    }
  }
  else {
    for (int i(2); i<d->Nx-2; i++) {
      q[ID(0, i, 0, 0)] = (-E[ID(0, i+2, 0, 0)] + 8*E[ID(0, i+1, 0, 0)] - 8*E[ID(0, i-1, 0, 0)] + E[ID(0, i-2, 0, 0)]) / (12*d->dx);
    }
  }

}


void REGIME::set_dwdsb(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Save some typing
        double qsigsq(q[ID(0, i, j, k)]*q[ID(0, i, j, k)] + d->sigma*d->sigma);
        double sigcu(d->sigma*d->sigma*d->sigma);
        double sig(d->sigma);
        double qch(q[ID(0, i, j, k)]);


        // First, do A
        {
          dwdsb[IDWS(1, 0, i, j, k)] = -qch * (qch*qch + sig*sig * (1 + prims[ID(5, i, j, k)]*prims[ID(5, i, j, k)]));
          dwdsb[IDWS(1, 1, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)]*qch*sig - prims[ID(7, i, j, k)]*qsigsq);
          dwdsb[IDWS(1, 2, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig + prims[ID(6, i, j, k)]*qsigsq);
          dwdsb[IDWS(2, 0, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)]*qch*sig + prims[ID(7, i, j, k)]*qsigsq);
          dwdsb[IDWS(2, 1, i, j, k)] = -qch * (qch*qch + sig*sig * (1 + prims[ID(6, i, j, k)]*prims[ID(6, i, j, k)]));
          dwdsb[IDWS(2, 2, i, j, k)] = -sig * (prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig - prims[ID(5, i, j, k)]*qsigsq);
          dwdsb[IDWS(3, 0, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig - prims[ID(6, i, j, k)]*qsigsq);
          dwdsb[IDWS(3, 1, i, j, k)] = -sig * (prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig + prims[ID(5, i, j, k)]*qsigsq);
          dwdsb[IDWS(3, 2, i, j, k)] = -qch * (qch*qch + sig*sig * (1 + prims[ID(7, i, j, k)]*prims[ID(7, i, j, k)]));
        }

        // Now do B
        {
          dwdsb[IDWS(5, 0, i, j, k)] = -prims[ID(5, i, j, k)] * sigcu * E[ID(0, i, j, k)];
          dwdsb[IDWS(5, 1, i, j, k)] = -prims[ID(6, i, j, k)] * sigcu * E[ID(0, i, j, k)] - sig*qsigsq*prims[ID(3, i, j, k)];
          dwdsb[IDWS(5, 2, i, j, k)] = -prims[ID(7, i, j, k)] * sigcu * E[ID(0, i, j, k)] + sig*qsigsq*prims[ID(2, i, j, k)];
          dwdsb[IDWS(6, 0, i, j, k)] = -prims[ID(5, i, j, k)] * sigcu * E[ID(1, i, j, k)] + sig*qsigsq*prims[ID(3, i, j, k)];
          dwdsb[IDWS(6, 1, i, j, k)] = -prims[ID(6, i, j, k)] * sigcu * E[ID(1, i, j, k)];
          dwdsb[IDWS(6, 2, i, j, k)] = -prims[ID(7, i, j, k)] * sigcu * E[ID(1, i, j, k)] - sig*qsigsq*prims[ID(1, i, j, k)];
          dwdsb[IDWS(7, 0, i, j, k)] = -prims[ID(5, i, j, k)] * sigcu * E[ID(2, i, j, k)] - sig*qsigsq*prims[ID(2, i, j, k)];
          dwdsb[IDWS(7, 1, i, j, k)] = -prims[ID(6, i, j, k)] * sigcu * E[ID(2, i, j, k)] + sig*qsigsq*prims[ID(1, i, j, k)];
          dwdsb[IDWS(7, 2, i, j, k)] = -prims[ID(7, i, j, k)] * sigcu * E[ID(2, i, j, k)];
        }

        // Now, C
        {
          dwdsb[IDWS(8, 0, i, j, k)] = -prims[ID(5, i, j, k)] * prims[ID(5, i, j, k)] * sigcu - sig*qsigsq;
          dwdsb[IDWS(8, 1, i, j, k)] = -prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)] * sigcu;
          dwdsb[IDWS(8, 2, i, j, k)] = -prims[ID(5, i, j, k)] * prims[ID(7, i, j, k)] * sigcu;
          dwdsb[IDWS(9, 0, i, j, k)] = -prims[ID(6, i, j, k)] * prims[ID(5, i, j, k)] * sigcu;
          dwdsb[IDWS(9, 1, i, j, k)] = -prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] * sigcu - sig*qsigsq;
          dwdsb[IDWS(9, 2, i, j, k)] = -prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)] * sigcu;
          dwdsb[IDWS(10, 0, i, j, k)] = -prims[ID(7, i, j, k)] * prims[ID(5, i, j, k)] * sigcu;
          dwdsb[IDWS(10, 1, i, j, k)] = -prims[ID(7, i, j, k)] * prims[ID(6, i, j, k)] * sigcu;
          dwdsb[IDWS(10, 2, i, j, k)] = -prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)] * sigcu - sig*qsigsq;
        }

        // Dont bother with D, it is always multiplied by zero.

        // Dont forget to multiply by the prefactor!
        for (int l(0); l<12; l++) {
          for (int m(0); m<3; m++) {
            dwdsb[IDWS(l, m, i, j, k)] *= alpha[ID(0, i, j, k)];
          }
        }
      }
    }
  }
}

void REGIME::set_K(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

// First do x-direction
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        fbx[ID(0, i, j, k)] = 0;
        fbx[ID(1, i, j, k)] = prims[ID(7, i, j, k)];
        fbx[ID(2, i, j, k)] = -prims[ID(6, i, j, k)];
      }
    }
  }
  // Reconstruct stiff fluxes
  this->fluxMethod->fluxReconstruction(E, NULL, NULL, fbx, d->fnet, 0, 3);
  // Add flux differencing to K
  for (int var(0); var<3; var++) {
    for (int i(1); i<d->Nx-1; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          K[ID(var, i, j, k)] += d->fnet[ID(var, i+1, j, k)]/(2*d->dx) - d->fnet[ID(var, i-1, j, k)]/(2*d->dx);
        }
      }
    }
  }

  // Now add y-contribution
  if (d->dims>1)
  {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          fby[ID(0, i, j, k)] = -prims[ID(7, i, j, k)];
          fby[ID(1, i, j, k)] = 0;
          fby[ID(2, i, j, k)] = prims[ID(5, i, j, k)];
        }
      }
    }
    // Reconstruct stiff fluxes
    this->fluxMethod->fluxReconstruction(E, NULL, NULL, fby, d->fnet, 1, 3);
    // Add flux differencing to K
    for (int var(0); var<3; var++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(1); j<d->Ny-1; j++) {
          for (int k(0); k<d->Nz; k++) {
            K[ID(var, i, j, k)] += d->fnet[ID(var, i, j+1, k)]/(2*d->dy) - d->fnet[ID(var, i, j-1, k)]/(2*d->dy);
          }
        }
      }
    }
  }

  // Finally, add z-contribution
  if (d->dims==3)
  {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          fbz[ID(0, i, j, k)] = prims[ID(6, i, j, k)];
          fbz[ID(1, i, j, k)] = -prims[ID(5, i, j, k)];
          fbz[ID(2, i, j, k)] = 0;
        }
      }
    }
    // Reconstruct stiff fluxes
    this->fluxMethod->fluxReconstruction(E, NULL, NULL, fbz, d->fnet, 2, 3);
    // Add flux differencing to K
    for (int var(0); var<3; var++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(1); k<d->Nz-1; k++) {
            K[ID(var, i, j, k)] += d->fnet[ID(var, i, j, k+1)]/(2*d->dz) - d->fnet[ID(var, i, j, k-1)]/(2*d->dz);
          }
        }
      }
    }
  }
}


void REGIME::set_dfxdw(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Row 0
        dfxdw[IDFW(0, 0, i, j, k)] = prims[ID(1, i, j, k)];
        dfxdw[IDFW(0, 1, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfxdw[IDFW(1, 4, i, j, k)] = 1;
        dfxdw[IDFW(1, 5, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(1, 6, i, j, k)] = prims[ID(6, i, j, k)];
        dfxdw[IDFW(1, 7, i, j, k)] = prims[ID(7, i, j, k)];
        dfxdw[IDFW(1, 8, i, j, k)] = -E[ID(0, i, j, k)];
        dfxdw[IDFW(1, 9, i, j, k)] = E[ID(1, i, j, k)];
        dfxdw[IDFW(1, 10, i, j, k)] = E[ID(2, i, j, k)];
        // Row 2
        dfxdw[IDFW(2, 5, i, j, k)] = -prims[ID(6, i, j, k)];
        dfxdw[IDFW(2, 6, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(2, 8, i, j, k)] = -E[ID(1, i, j, k)];
        dfxdw[IDFW(2, 9, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 3
        dfxdw[IDFW(3, 5, i, j, k)] = -prims[ID(7, i, j, k)];
        dfxdw[IDFW(3, 7, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(3, 8, i, j, k)] = -E[ID(2, i, j, k)];
        dfxdw[IDFW(3, 10, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 4
        dfxdw[IDFW(4, 1, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfxdw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(1, i, j, k)] / (d->gamma - 1);
        dfxdw[IDFW(4, 6, i, j, k)] = -E[ID(2, i, j, k)];
        dfxdw[IDFW(4, 7, i, j, k)] = E[ID(1, i, j, k)];
        dfxdw[IDFW(4, 9, i, j, k)] = prims[ID(7, i, j, k)];
        dfxdw[IDFW(4, 10, i, j, k)] = -prims[ID(6, i, j, k)];
        // Row 6
        dfxdw[IDFW(6, 10, i, j, k)] = -1;
        // Row 7
        dfxdw[IDFW(7, 9, i, j, k)] = 1;
      }
    }
  }
}


void REGIME::set_dfydw(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Row 0
        dfydw[IDFW(0, 0, i, j, k)] = prims[ID(2, i, j, k)];
        dfydw[IDFW(0, 2, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfydw[IDFW(1, 5, i, j, k)] = -prims[ID(6, i, j, k)];
        dfydw[IDFW(1, 6, i, j, k)] = -prims[ID(5, i, j, k)];
        dfydw[IDFW(1, 8, i, j, k)] = -E[ID(1, i, j, k)];
        dfydw[IDFW(1, 9, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 2
        dfydw[IDFW(2, 4, i, j, k)] = 1;
        dfydw[IDFW(2, 5, i, j, k)] = prims[ID(5, i, j, k)];
        dfydw[IDFW(2, 6, i, j, k)] = -prims[ID(6, i, j, k)];
        dfydw[IDFW(2, 7, i, j, k)] = prims[ID(7, i, j, k)];
        dfydw[IDFW(2, 8, i, j, k)] = E[ID(0, i, j, k)];
        dfydw[IDFW(2, 9, i, j, k)] = -E[ID(1, i, j, k)];
        dfydw[IDFW(2, 10, i, j, k)] = E[ID(2, i, j, k)];
        // Row 3
        dfydw[IDFW(3, 6, i, j, k)] = -prims[ID(7, i, j, k)];
        dfydw[IDFW(3, 7, i, j, k)] = -prims[ID(6, i, j, k)];
        dfydw[IDFW(3, 9, i, j, k)] = -E[ID(2, i, j, k)];
        dfydw[IDFW(3, 10, i, j, k)] = -E[ID(1, i, j, k)];
        // Row 4
        dfydw[IDFW(4, 2, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfydw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(2, i, j, k)] / (d->gamma - 1);
        dfydw[IDFW(4, 5, i, j, k)] = E[ID(2, i, j, k)];
        dfydw[IDFW(4, 7, i, j, k)] = -E[ID(0, i, j, k)];
        dfydw[IDFW(4, 8, i, j, k)] = -prims[ID(7, i, j, k)];
        dfydw[IDFW(4, 10, i, j, k)] = prims[ID(5, i, j, k)];
        // Row 6
        dfydw[IDFW(5, 10, i, j, k)] = 1;
        // Row 7
        dfydw[IDFW(7, 8, i, j, k)] = -1;
      }
    }
  }
}


void REGIME::set_dfzdw(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Row 0
        dfzdw[IDFW(0, 0, i, j, k)] = prims[ID(3, i, j, k)];
        dfzdw[IDFW(0, 3, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfzdw[IDFW(1, 5, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(1, 7, i, j, k)] = -prims[ID(5, i, j, k)];
        dfzdw[IDFW(1, 8, i, j, k)] = -E[ID(2, i, j, k)];
        dfzdw[IDFW(1, 10, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 2
        dfzdw[IDFW(2, 6, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(2, 7, i, j, k)] = -prims[ID(6, i, j, k)];
        dfzdw[IDFW(2, 9, i, j, k)] = -E[ID(2, i, j, k)];
        dfzdw[IDFW(2, 10, i, j, k)] = -E[ID(1, i, j, k)];
        // Row 3
        dfzdw[IDFW(3, 4, i, j, k)] = 1;
        dfzdw[IDFW(3, 5, i, j, k)] = prims[ID(5, i, j, k)];
        dfzdw[IDFW(3, 6, i, j, k)] = prims[ID(6, i, j, k)];
        dfzdw[IDFW(3, 7, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(3, 8, i, j, k)] = E[ID(0, i, j, k)];
        dfzdw[IDFW(3, 9, i, j, k)] = E[ID(1, i, j, k)];
        dfzdw[IDFW(3, 10, i, j, k)] = -E[ID(2, i, j, k)];
        // Row 4
        dfzdw[IDFW(4, 3, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfzdw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(3, i, j, k)] / (d->gamma - 1);
        dfzdw[IDFW(4, 5, i, j, k)] = -E[ID(1, i, j, k)];
        dfzdw[IDFW(4, 6, i, j, k)] = E[ID(0, i, j, k)];
        dfzdw[IDFW(4, 8, i, j, k)] = prims[ID(6, i, j, k)];
        dfzdw[IDFW(4, 9, i, j, k)] = -prims[ID(5, i, j, k)];
        // Row 6
        dfzdw[IDFW(5, 9, i, j, k)] = -1;
        // Row 7
        dfzdw[IDFW(6, 8, i, j, k)] = 1;
      }
    }
  }
}



void REGIME::set_Dx(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // We have moved the set_dfxdw function in to this loop ---> saves ~few percent on execution
  // this->set_dfxdw(cons, prims, aux);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // First set dfxdw
        // Row 0
        dfxdw[IDFW(0, 0, i, j, k)] = prims[ID(1, i, j, k)];
        dfxdw[IDFW(0, 1, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfxdw[IDFW(1, 4, i, j, k)] = 1;
        dfxdw[IDFW(1, 5, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(1, 6, i, j, k)] = prims[ID(6, i, j, k)];
        dfxdw[IDFW(1, 7, i, j, k)] = prims[ID(7, i, j, k)];
        dfxdw[IDFW(1, 8, i, j, k)] = -E[ID(0, i, j, k)];
        dfxdw[IDFW(1, 9, i, j, k)] = E[ID(1, i, j, k)];
        dfxdw[IDFW(1, 10, i, j, k)] = E[ID(2, i, j, k)];
        // Row 2
        dfxdw[IDFW(2, 5, i, j, k)] = -prims[ID(6, i, j, k)];
        dfxdw[IDFW(2, 6, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(2, 8, i, j, k)] = -E[ID(1, i, j, k)];
        dfxdw[IDFW(2, 9, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 3
        dfxdw[IDFW(3, 5, i, j, k)] = -prims[ID(7, i, j, k)];
        dfxdw[IDFW(3, 7, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(3, 8, i, j, k)] = -E[ID(2, i, j, k)];
        dfxdw[IDFW(3, 10, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 4
        dfxdw[IDFW(4, 1, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfxdw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(1, i, j, k)] / (d->gamma - 1);
        dfxdw[IDFW(4, 6, i, j, k)] = -E[ID(2, i, j, k)];
        dfxdw[IDFW(4, 7, i, j, k)] = E[ID(1, i, j, k)];
        dfxdw[IDFW(4, 9, i, j, k)] = prims[ID(7, i, j, k)];
        dfxdw[IDFW(4, 10, i, j, k)] = -prims[ID(6, i, j, k)];
        // Row 6
        dfxdw[IDFW(6, 10, i, j, k)] = -1;
        // Row 7
        dfxdw[IDFW(7, 9, i, j, k)] = 1;

        // Now compute Mx: Mx = -1 * DOT(dfxdw, dwdsb)
        // Optimised version: many entries of dfxdw are zero so explicitly do these
        // Note: element (0, 0) of dfxdw multiplies a zero in dwdsb so ignore, same for [:, 4]
        Mx[IDM(0, 0, i, j, k)] -= (dfxdw[IDFW(0, 1, i, j, k)] * dwdsb[IDWS(1, 0, i, j, k)]);
        Mx[IDM(0, 1, i, j, k)] -= (dfxdw[IDFW(0, 1, i, j, k)] * dwdsb[IDWS(1, 1, i, j, k)]);
        Mx[IDM(0, 2, i, j, k)] -= (dfxdw[IDFW(0, 1, i, j, k)] * dwdsb[IDWS(1, 2, i, j, k)]);

        Mx[IDM(4, 0, i, j, k)] -= (dfxdw[IDFW(4, 1, i, j, k)] * dwdsb[IDWS(1, 0, i, j, k)]);
        Mx[IDM(4, 1, i, j, k)] -= (dfxdw[IDFW(4, 1, i, j, k)] * dwdsb[IDWS(1, 1, i, j, k)]);
        Mx[IDM(4, 2, i, j, k)] -= (dfxdw[IDFW(4, 1, i, j, k)] * dwdsb[IDWS(1, 2, i, j, k)]);

        Mx[IDM(6, 0, i, j, k)] -= (- dwdsb[IDWS(10, 0, i, j, k)]);
        Mx[IDM(6, 1, i, j, k)] -= (- dwdsb[IDWS(10, 1, i, j, k)]);
        Mx[IDM(6, 2, i, j, k)] -= (- dwdsb[IDWS(10, 2, i, j, k)]);

        Mx[IDM(7, 0, i, j, k)] -= (dwdsb[IDWS(9, 0, i, j, k)]);
        Mx[IDM(7, 1, i, j, k)] -= (dwdsb[IDWS(9, 1, i, j, k)]);
        Mx[IDM(7, 2, i, j, k)] -= (dwdsb[IDWS(9, 2, i, j, k)]);

        for (int l(1); l<5; l++) {
          for (int m(0); m<3; m++) {
            for (int n(5); n<11; n++) {
              Mx[IDM(l, m, i, j, k)] -= dfxdw[IDFW(l, n, i, j, k)] * dwdsb[IDWS(n, m, i, j, k)];
            }
          }
        }

        // Now can get the diffusion vector
        // Dx = DOT(Mx, K)
        for (int l(0); l<8; l++) {
          for (int m(0); m<3; m++) {
            diffuX[ID(l, i, j, k)] += Mx[IDM(l, m, i, j, k)] * K[ID(m, i, j, k)];
          }
        }
      }
    }
  }
}




void REGIME::set_Dy(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // We have moved the set_dfydw function in to this loop ---> saves ~few percent on execution
  // this->set_dfydw(cons, prims, aux);


  // My = -1 * DOT(dfydw, dwdsb)
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {

          // First set dfydw
          // Row 0
          dfydw[IDFW(0, 0, i, j, k)] = prims[ID(2, i, j, k)];
          dfydw[IDFW(0, 2, i, j, k)] = prims[ID(0, i, j, k)];
          // Row 1
          dfydw[IDFW(1, 5, i, j, k)] = -prims[ID(6, i, j, k)];
          dfydw[IDFW(1, 6, i, j, k)] = -prims[ID(5, i, j, k)];
          dfydw[IDFW(1, 8, i, j, k)] = -E[ID(1, i, j, k)];
          dfydw[IDFW(1, 9, i, j, k)] = -E[ID(0, i, j, k)];
          // Row 2
          dfydw[IDFW(2, 4, i, j, k)] = 1;
          dfydw[IDFW(2, 5, i, j, k)] = prims[ID(5, i, j, k)];
          dfydw[IDFW(2, 6, i, j, k)] = -prims[ID(6, i, j, k)];
          dfydw[IDFW(2, 7, i, j, k)] = prims[ID(7, i, j, k)];
          dfydw[IDFW(2, 8, i, j, k)] = E[ID(0, i, j, k)];
          dfydw[IDFW(2, 9, i, j, k)] = -E[ID(1, i, j, k)];
          dfydw[IDFW(2, 10, i, j, k)] = E[ID(2, i, j, k)];
          // Row 3
          dfydw[IDFW(3, 6, i, j, k)] = -prims[ID(7, i, j, k)];
          dfydw[IDFW(3, 7, i, j, k)] = -prims[ID(6, i, j, k)];
          dfydw[IDFW(3, 9, i, j, k)] = -E[ID(2, i, j, k)];
          dfydw[IDFW(3, 10, i, j, k)] = -E[ID(1, i, j, k)];
          // Row 4
          dfydw[IDFW(4, 2, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
          dfydw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(2, i, j, k)] / (d->gamma - 1);
          dfydw[IDFW(4, 5, i, j, k)] = E[ID(2, i, j, k)];
          dfydw[IDFW(4, 7, i, j, k)] = -E[ID(0, i, j, k)];
          dfydw[IDFW(4, 8, i, j, k)] = -prims[ID(7, i, j, k)];
          dfydw[IDFW(4, 10, i, j, k)] = prims[ID(5, i, j, k)];
          // Row 6
          dfydw[IDFW(5, 10, i, j, k)] = 1;
          // Row 7
          dfydw[IDFW(7, 8, i, j, k)] = -1;


          // Now compute My: My = -1 * DOT(dfydw, dwdsb)
          // Optimised version: many entries of dfxdw are zero so explicitly do these
          // Note: element (0, 0) of dfydw multiplies a zero in dwdsb so ignore, same for [:, 4]
          My[IDM(0, 0, i, j, k)] -= (dfydw[IDFW(0, 2, i, j, k)] * dwdsb[IDWS(2, 0, i, j, k)]);
          My[IDM(0, 1, i, j, k)] -= (dfydw[IDFW(0, 2, i, j, k)] * dwdsb[IDWS(2, 1, i, j, k)]);
          My[IDM(0, 2, i, j, k)] -= (dfydw[IDFW(0, 2, i, j, k)] * dwdsb[IDWS(2, 2, i, j, k)]);

          My[IDM(4, 0, i, j, k)] -= (dfydw[IDFW(4, 2, i, j, k)] * dwdsb[IDWS(2, 0, i, j, k)]);
          My[IDM(4, 1, i, j, k)] -= (dfydw[IDFW(4, 2, i, j, k)] * dwdsb[IDWS(2, 1, i, j, k)]);
          My[IDM(4, 2, i, j, k)] -= (dfydw[IDFW(4, 2, i, j, k)] * dwdsb[IDWS(2, 2, i, j, k)]);

          My[IDM(5, 0, i, j, k)] -= (dwdsb[IDWS(10, 0, i, j, k)]);
          My[IDM(5, 1, i, j, k)] -= (dwdsb[IDWS(10, 1, i, j, k)]);
          My[IDM(5, 2, i, j, k)] -= (dwdsb[IDWS(10, 2, i, j, k)]);

          My[IDM(7, 0, i, j, k)] -= (- dwdsb[IDWS(8, 0, i, j, k)]);
          My[IDM(7, 1, i, j, k)] -= (- dwdsb[IDWS(8, 1, i, j, k)]);
          My[IDM(7, 2, i, j, k)] -= (- dwdsb[IDWS(8, 2, i, j, k)]);

          for (int l(1); l<5; l++) {
            for (int m(0); m<3; m++) {
              for (int n(5); n<11; n++) {
                My[IDM(l, m, i, j, k)] -= dfydw[IDFW(l, n, i, j, k)] * dwdsb[IDWS(n, m, i, j, k)];
              }
            }
          }

        // Now can get the diffusion vector
        // Dy = DOT(My, K)
        for (int l(0); l<8; l++) {
          for (int m(0); m<3; m++) {
            diffuY[ID(l, i, j, k)] += My[IDM(l, m, i, j, k)] * K[ID(m, i, j, k)];
          }
        }
      }
    }
  }
}


void REGIME::set_Dz(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // We have moved the set_dfzdw function in to this loop ---> saves ~few percent on execution
  // this->set_dfzdw(cons, prims, aux);


  // Mz = -1 * DOT(dfzdw, dwdsb)
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {

        // First set dfzdw

        // Row 0
        dfzdw[IDFW(0, 0, i, j, k)] = prims[ID(3, i, j, k)];
        dfzdw[IDFW(0, 3, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfzdw[IDFW(1, 5, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(1, 7, i, j, k)] = -prims[ID(5, i, j, k)];
        dfzdw[IDFW(1, 8, i, j, k)] = -E[ID(2, i, j, k)];
        dfzdw[IDFW(1, 10, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 2
        dfzdw[IDFW(2, 6, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(2, 7, i, j, k)] = -prims[ID(6, i, j, k)];
        dfzdw[IDFW(2, 9, i, j, k)] = -E[ID(2, i, j, k)];
        dfzdw[IDFW(2, 10, i, j, k)] = -E[ID(1, i, j, k)];
        // Row 3
        dfzdw[IDFW(3, 4, i, j, k)] = 1;
        dfzdw[IDFW(3, 5, i, j, k)] = prims[ID(5, i, j, k)];
        dfzdw[IDFW(3, 6, i, j, k)] = prims[ID(6, i, j, k)];
        dfzdw[IDFW(3, 7, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(3, 8, i, j, k)] = E[ID(0, i, j, k)];
        dfzdw[IDFW(3, 9, i, j, k)] = E[ID(1, i, j, k)];
        dfzdw[IDFW(3, 10, i, j, k)] = -E[ID(2, i, j, k)];
        // Row 4
        dfzdw[IDFW(4, 3, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfzdw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(3, i, j, k)] / (d->gamma - 1);
        dfzdw[IDFW(4, 5, i, j, k)] = -E[ID(1, i, j, k)];
        dfzdw[IDFW(4, 6, i, j, k)] = E[ID(0, i, j, k)];
        dfzdw[IDFW(4, 8, i, j, k)] = prims[ID(6, i, j, k)];
        dfzdw[IDFW(4, 9, i, j, k)] = -prims[ID(5, i, j, k)];
        // Row 6
        dfzdw[IDFW(5, 9, i, j, k)] = -1;
        // Row 7
        dfzdw[IDFW(6, 8, i, j, k)] = 1;

        // Now compute My: My = -1 * DOT(dfzdw, dwdsb)
        // Optimised version: many entries of dfxdw are zero so explicitly do these
        // Note: element (0, 0) of dfzdw multiplies a zero in dwdsb so ignore, same for [:, 4]
        Mz[IDM(0, 0, i, j, k)] -= (dfzdw[IDFW(0, 3, i, j, k)] * dwdsb[IDWS(3, 0, i, j, k)]);
        Mz[IDM(0, 1, i, j, k)] -= (dfzdw[IDFW(0, 3, i, j, k)] * dwdsb[IDWS(3, 1, i, j, k)]);
        Mz[IDM(0, 2, i, j, k)] -= (dfzdw[IDFW(0, 3, i, j, k)] * dwdsb[IDWS(3, 2, i, j, k)]);

        Mz[IDM(4, 0, i, j, k)] -= (dfzdw[IDFW(4, 3, i, j, k)] * dwdsb[IDWS(3, 0, i, j, k)]);
        Mz[IDM(4, 1, i, j, k)] -= (dfzdw[IDFW(4, 3, i, j, k)] * dwdsb[IDWS(3, 1, i, j, k)]);
        Mz[IDM(4, 2, i, j, k)] -= (dfzdw[IDFW(4, 3, i, j, k)] * dwdsb[IDWS(3, 2, i, j, k)]);

        Mz[IDM(5, 0, i, j, k)] -= (- dwdsb[IDWS(9, 0, i, j, k)]);
        Mz[IDM(5, 1, i, j, k)] -= (- dwdsb[IDWS(9, 1, i, j, k)]);
        Mz[IDM(5, 2, i, j, k)] -= (- dwdsb[IDWS(9, 2, i, j, k)]);

        Mz[IDM(6, 0, i, j, k)] -= (dwdsb[IDWS(8, 0, i, j, k)]);
        Mz[IDM(6, 1, i, j, k)] -= (dwdsb[IDWS(8, 1, i, j, k)]);
        Mz[IDM(6, 2, i, j, k)] -= (dwdsb[IDWS(8, 2, i, j, k)]);
        for (int l(1); l<5; l++) {
          for (int m(0); m<3; m++) {
            for (int n(5); n<11; n++) {
              Mz[IDM(l, m, i, j, k)] -= dfzdw[IDFW(l, n, i, j, k)] * dwdsb[IDWS(n, m, i, j, k)];
            }
          }
        }

        // Dz = DOT(Mz, K)
        for (int l(0); l<8; l++) {
          for (int m(0); m<3; m++) {
            diffuZ[ID(l, i, j, k)] += Mz[IDM(l, m, i, j, k)] * K[ID(m, i, j, k)];
          }
        }
      }
    }
  }
}
