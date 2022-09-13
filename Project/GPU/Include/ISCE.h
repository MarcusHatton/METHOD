#ifndef ISCE_H
#define ISCE_H

#include "model.h"
#include "C2PArgs.h"
#include "deviceArguments.h"

#include <stdio.h>

/*
This is the human readable description of this models variables.

  ISCE has 5 conserved variables:
[0-5]  D, Sx, Sy, Sz, tau
*/

//! <b> Special Relativistic Dissipative HydroDynamics </b>
/*!
    @par
      @todo COMPLETE DESCRIPTION
*/
class ISCE : public Model
{

  public:

    // enums to save looking up numbering of C/P/As when using ID accessor.
    enum Cons { D, S1, S2, S3, Tau };
    enum Prims { v1, v2, v3, p, rho, n, q1, q2, q3, Pi, pi11, pi12, pi13, pi22, pi23, pi33 };
    enum Aux { h, T, e, W, q0, qv, pi00, pi01, pi02, pi03, Theta, vsqrd, rhohWsq, S_sqrd,
               q1NS, q2NS, q3NS, PiNS, pi11NS, pi12NS, pi13NS, pi22NS, pi23NS, pi33NS,
               q1LO, q2LO, q3LO, PiLO, pi11LO, pi12LO, pi13LO, pi22LO, pi23LO, pi33LO,  
               a1, a2, a3 };
    enum TDerivs { dtp = 37, dtrho, dtn, dtv1, dtv2, dtv3, dtW, dtT, dtq1NS, dtq2NS, dtq3NS, dtPiNS,
               dtpi11NS, dtpi12NS, dtpi13NS, dtpi22NS, dtpi23NS, dtpi33NS, dtD, dtS1, dtS2, dtS3,
               dtTau, dtE};

    // Work arrays
    double * singleCons;
    double * singlePrims;
    double * singleAux;
    double * singleSource;
    C2PArgs * c2pArgs;

    ISCE(); //!< Default constructor

    //! Parameterized constructor
    /*!
      @par
        Stores a pointer to the Data class for reference in its members

      @param *data pointer to Data class containing global simulation data
    */
    ISCE(Data * data);

    ~ISCE();  //!< Destructor

    //! Single cell source term contribution
    /*!
        Determines the source vector due the the cond prims and aux vector
      of a single compute cell.
      Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxiliary vector work array. Size is Naux
      @param *source pointer to source vector work array. Size is Ncons
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::sourceTermSingleCell
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
    Source terms arise from the evolution of the electric fields
    and from implementing the divergence cleaning method. This function
    determines the source contribution to the change in the conserved vector
    for all cells. This function calls sourceTermSingleCell for every compute
    cell.

    @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
    @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
    @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
    @param *source pointer to source vector work array. Size is Ncons*Nx*Ny*Nz
    @sa Model::sourceTerm
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell spectral decomposition
    /*!
        Generates the values for aux and prims for the given cons vector for only
      a single cell (i, j, k).
        Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.
        Single celled version required for the inside of the residual functions
      for the (semi) implicit time integrators.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxiliary vector work array. Size is Naux
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);

    //! Spectral decomposition
    /*!
    Generates the values of the primitive and auxiliary variables consistent
    with the conservative vector given. Method first subtracts the EM fields
    from the conserved quantities, reducing the problem to the hydrodynamic
    properties only before determining their values via a newton method.
    Function calls single celled version (below) for each compute cell.

    @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
    @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
    @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
    @sa Model::getPrimitiveVars
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Dionysopoulou 2016.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
        Method determines the value of the conserved vector in the specified direction.
      For the form of the fluxes see Dionysopoulou 2016 with the inclusion of
      divergence cleaning from Muddle 2015.

      @note We are assuming that all primitive and auxiliary variables are up-to-date
      at the time of this function execution.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *f pointer to flux vector work array. Size is Ncons*Nx*Ny*Nz
      @param dir direction in which to generate flux vector. (x, y, z) = (0, 1, 2)
      @sa Model::fluxVector
    */
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);
};

//! <b> ISCE class on the device </b>
/*!
  @par
    Device class for ISCE
*/
class ISCE_D : public Model_D
{
  public:
    __device__ ISCE_D(TimeIntAndModelArgs * args) : Model_D(args) { }

    //!< @sa Model::sourceTermSingleCell
    __device__ void sourceTermSingleCell(double * cons, double * prims, double * aux, double * source);

    //!< @sa Model::getPrimitiveVarsSingleCell
    __device__ void getPrimitiveVarsSingleCell(double * cons, double * prims, double * aux);
};

#endif
