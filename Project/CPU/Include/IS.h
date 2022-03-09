#ifndef IS_H
#define IS_H

#include "model.h"

class IS : public Model
{

  public:

    // enums to save looking up numbering of C/P/As when using ID accessor.
    enum Cons { S1, S2, S3, Tau, v1_C, v2_C, v3_C, p_C, rho_C };
    enum Prims { v1, v2, v3, p, rho, dv1dt, dv2dt, dv3dt, dpdt, drhodt };
    enum Aux { A, Theta, Pi, q0, q1, q2, q3, qv, pi11, pi12, pi13, pi22, pi23, pi33, 
               pi00, pi01, pi02, pi03, n, h, T, e, W, dndt, dWdt };


    int smartGuesses;     //!< Number of smart guess required

    double * solution;    //!< Pointer to array to hold solution of C2P for every cell. Size is 2*Nx*Ny*Nz
    
    IS();     //!< Default constructor

    //! Parameterized constructor
    /*!
      @parm
      @param *data Pointer to Data class containing global simulation data
    */
    IS(Data * data, bool alternative_C2P);

    virtual ~IS();     //!< Destructor


    //! Single cell source term contribution
    /*!
      @par
        Models that can posess a stiff source term and hence (semi-)implicit time
      integrators will require a source contribution (and cons2prims method) that
      applies to a single cell.
      @par
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}\f$
      @param *source pointer to source vector work array. Size is \f$N_{cons}\f$
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::sourceTermSingleCell
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
      @par
        Non-zero flux for cons[8], phi, as a result of implementing divergence
      cleaning. For details see Muddle.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param *source pointer to source vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @sa Model::sourceTerm
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell cons2prims conversion
    /*!
      @par
        For the same reason as outlined in sourceTermSingleCell, some models will
      require a single celled primitive conversion method.
      @par
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}\f$
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);

    //! Spectral decomposition
    /*!
      @par
        Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
      Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
      an N=2 rootfind using cminpack library.
      @par
        Initial inputs will be the current values of the conserved vector and the
      <b> old </b> values for the prims and aux vectors.
      Output will be the current values of cons, prims and aux.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @sa Model::getPrimitiveVars
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
      @par
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Anton 2010.
      Riemann Solver`

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
      @par
        Method determines the value for the flux vector.
      @par
        For the form of the fluxes see Anton 2010 with the inclusion of divergence
      cleaning from Dedner et al. 2002.
      interfaces, John Muddle.

      @note We are assuming that all primitive and auxiliary variables are up-to-date
      at the time of this function execution.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param *f pointer to flux vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param dir direction in which to generate flux vector. \f$(x, y, z) = (0, 1, 2)\f$
      @sa Model::fluxVector
    */
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);

    //! Finalise the simulation variables
    /*!
      @par
        Mostly, this probably wont be needed, but if there is any final steps to finish
      off a timestep, this can be done here.
    */
    void finalise(double *cons, double *prims, double *aux, bool final_step=false) {};

    //! <b> Additional arguments for the IS residual function </b>
    /*!
      @par
        The conservative to primitive transformation for the IS class requires an
      N=2 dimensional nonlinear rootfind and thus requires the multi-dimensional Newton-
      Secant solver of the Cminpack library, i.e. the function @e hybrd1. This function
      can take additional arguments in the form of a void pointer to some data object,
      and thus for the data required in the cons2prims solver, the additional data is
      located in this Args data structure.
    
    */
    typedef struct
    {
      double
      S1_rf,
      S2_rf,
      S3_rf,
      Tau_rf,
      v1_rf,
      v2_rf,
      v3_rf,
      W_rf,
      p_rf,
      rho_rf,
      n_rf,
      dndx_rf,
      dndy_rf,
      dndz_rf,
      dxrho_rf,
      dyrho_rf,
      dzrho_rf,
      dxux_rf,
      dxuy_rf,
      dxuz_rf,
      dyux_rf,
      dyuy_rf,
      dyuz_rf,
      dzux_rf,
      dzuy_rf,
      dzuz_rf,
      tau_epsilon_rf,
      tau_Pi_rf,
      tau_q_rf,
      eta_rf,
      zeta_rf,
      gamma;
      int i;
    } Args;

    //! <b> Stores data of the failed cons2prims rootfinder </b>
    /*!
      @par
        When the cons2prims rootfinder fails, we can take note of the cell, continue
      throughout the domain and come back to that failed cell, using the average of
      the successfully completed surrounding cells as an initial estimate for the
      solution of the failed cell. This struct holds the failed cells data, and thus
      we can use a vector type to hold instances of this structure when an unknown
      number of cells inevitably fail.
    */
    typedef struct
    {
      // Store coordinates of the failed cell
      int
      //@{
      x, y, z;  //!< Cell number of failed C2P conversion
      //@}
    } Failed;

};

//! <b> Residual function for spectral analysis </b>
/*!
  @par
    IS requires N=4 rootfind, therefore need to implement the hybrd cminpack
  multiD Newton solver. Things may get ugly.
  @par
    Cant do anything about the arguments of this function, cminpack demands you
  not enjoy anything it offers...

  @param *p void pointer to the additional arguments struct, Args
  @param n size of system (n=2 for IS)
  @param *x pointer to array containing initial estimate of solution, will also hold solution
  @param *fvec pointer to array to hold residual values. These should be 0 +- tol
  @param iflag Cminpack error flag

  @note For more information regarding the form of this function and its parameters see the URL below
  @sa [Original source](https://github.com/devernay/cminpack)
*/
int residual(void *p, int n, const double *x, double *fvec, int iflag);




  



#endif
