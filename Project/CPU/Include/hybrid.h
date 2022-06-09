#ifndef HYBRID_H
#define HYBRID_H

#include "model.h"
#include "IS.h"
#include "ISCE.h"
#include "DEIFY.h"
#include "flux.h"

//! <b> Hybrid model </b>
/*!
  @par
    This hybrid model combines the ideal and dissipative models, dynamically
    switching between the two dependant upon the local (and nearest neighbours'
    if using DEIFY) conductivity. This model is SRRMHD in disguise, that is,
    for any cons, prims, aux arguments, these refer to the resistive versions.
    Ideal vectors are prefixed with `i'.
  @par
    Switching between the models uses the tanh function, centred around some
    crossover conductivity, \f$\sigma_{crossover}\f$, and with in the span of
    \f$\pm\sigma_{span}\f$. For conductivities less than \f$\sigma_c-\sigma_s\f$
    the resistive model is used, and for \f$\simga_c+\sigma_s\f$ ideal is used.
    If DEIFY is used, it is calculated across the entire domain and added to
    the ideal MHD contribution.
  @par
    The source and flux vectors of the hybrid model are a combination of the
    SRMHD and SRRMHD source and flux vectors, weighted according to the tanh
    penalty function.
  @par
    The C2P of this model utilises the resistive C2P for \f$\sigma<\sigma_c\f$
    and the ideal C2P otherwise.
  @par
    If specified in the constructor, he DEIFY source is calculated for the
    whole domain, but only added to the source vector if the cell, and its
    nearest 3 neighbours have conductivites \f$\sigma\ge\sigma_c-\sigma_s\f$.
*/
class Hybrid : public Model
{

  public:

    // enums to save looking up numbering of C/P/As when using ID accessor.
    enum Cons { D, S1, S2, S3, Tau, Y1, Y2, Y3, U, Z11, Z12, Z13, Z22, Z23, Z33 };
    enum Prims { v1, v2, v3, p, rho, n, q1, q2, q3, Pi, pi11, pi12, pi13, pi22, pi23, pi33 };
    enum Aux { h, T, e, W, q0, qv, pi00, pi01, pi02, pi03, q1NS, q2NS, q3NS, PiNS, 
               pi11NS, pi12NS, pi13NS, pi22NS, pi23NS, pi33NS, Theta, dv1dt, 
               dv2dt, dv3dt, a1, a2, a3, vsqrd, dWdt, rho_plus_p };
               
    double
    *icons, *iprims, *iaux,                     // Ideal cons, prims and aux. Size is \f$N_{var}*N_x*N_y*N_z\f$.
    *sicons, *siprims, *siaux,                  // Ideal cons, prims and aux. Size is \f$N_{var}\f$.
    *iflux, *dflux,                             // Flux vectors for ideal, resistive and DEIFY. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
    *isource, *dsource, *deifySource,          // Source vectors for ideal, resistive and DEIFY. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
    tauCrossOver,                             // Centre conductivity of penalty function
    tauSpan,                                  // Span of conductivity of penalty function


    bool useDEIFY;                             // Should we use DEIFY? (Default to true)

    IS * dissipativeModel;                    // Pointer to SRRMHD model

    ISCE * idealModel;                         // Pointer to SRMHD model

    DEIFY * subgridModel = NULL;               // Pointer to DEIFY model
    
    int *mask;                                  // Flag: can we set DEIFY source?

    Hybrid();                                   //!< Default constructor


    //! Parameterized constructor
    /*!
      @par
        Stores a pointer to the Data class for reference in its methods. Sets
        up the parameters for the hybrid model.

        @param[in *data Pointer to Data class containing global simulation data
        @param[in] tauCrossOver double centre conductivity of penalty function (defaults to 150)
        @param[in] tauSpan double conductivity span of penalty function (defaults to 50)
        @param[in] useDEIFY bool should we use DEIFY? (defaults to true)
    */
    Hybrid(Data * data, double tauCrossOver=150, double tauSpan=50, bool useDEIFY=true);


    virtual ~Hybrid();  //!< Destructor

    //! Setup the DEIFY model
    /*!
      @par
        Initialise the DEIFY model and store a pointer it. Also allocate memory
        for the DEIFY source vector and the mask.

        @param[in] fluxMethod pointer to simulations FluxMethod
    */
    void setupDEIFY(FluxMethod * fluxMethod);

  private:
    //! Penalty function: ideal contribution
    /*!
      @par
        Compute the contribution from idealMHD for a single cell. The arguments
        are the cons prims and aux of only one cell.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}\f$.
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}\f$.
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}\f$.
        @return iW double fractional contribution of the ideal model. (1-iW) gives the resistive contribution.
    */
    double idealWeight(double * cons, double * prims, double * aux);

    //! Penalty function: ideal contribution
    /*!
      @par
        Compute the contribution from idealMHD for a single cell. The arguments
        are the cons prims and aux of the whole domain, and a cell position to
        return the weighting for.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$.
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$.
        @param i int x-coordinate of cell to get weighting for.
        @param j int y-coordinate of cell to get weighting for.
        @param k int z-coordinate of cell to get weighting for.
        @return iW double fractional contribution of the ideal model. (1-iW) gives the resistive contribution.
    */
    double idealWeightID(double * cons, double * prims, double * aux, int i, int j, int k);

    //! Use resistive C2P?
    /*!
      @par
        Should we use the resistive C2P? If the local conductivity is less than
        the crossover conductivity, we use the resistive C2P.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$.
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$.
        @return use bool do we use the resistive C2P?
    */
    bool useDissipative(double * cons, double * prims, double * aux);

    //! Set ideal cons, prims and aux vectors for a single cell
    /*!
      @par
        Set the ideal conserved, primitive and auxiliary vectors from the
        (given) resistive conserved, primitive and auxiliary vectors, for a
        single computational cell.

        @param[in] cons pointer to resistive conserved vector. Size is \f$N_{cons}\f$.
        @param[in] prims pointer to resistive primitive vector. Size is \f$N_{prims}\f$.
        @param[in] aux pointer to resistive auxiliary vector. Size is \f$N_{aux}\f$.
        @param[out] sicons pointer to ideal conserved vector. Size is \f$N_{cons}\f$.
        @param[out] siprims pointer to ideal primitive vector. Size is \f$N_{prims}\f$.
        @param[out] siaux pointer to ideal auxiliary vector. Size is \f$N_{aux}\f$.
    */
    void setIdealCPAs(double *cons, double * prims, double * aux);

    //! Set ideal cons, prims and aux vectors for all cells
    /*!
      @par
        Set the ideal conserved, primitive and auxiliary vectors from the
        (given) resistive conserved, primitive and auxiliary vectors, for all
        computational cells.

        @param[in] cons pointer to resistive conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
        @param[in] prims pointer to resistive primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$.
        @param[in] aux pointer to resistive auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$.
        @param[out] icons pointer to ideal conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
        @param[out] iprims pointer to ideal primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$.
        @param[out] iaux pointer to ideal auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$.
    */
    void setIdealCPAsAll(double *cons, double * prims, double * aux);

    //! Set the DEIFY source mask
    /*!
      @par
        Set the DEIFY mask to indicate whether the DEIFY source should be
        added to the hybrid->source vector.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$.
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$.
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$.
        @param[out] mask pointer to mask array.
    */
    void setMasks(double * cons, double * prims, double * aux);

  public:
    //! Single cell source term contribution
    /*!
      @par
        Source term calculation for a single computational cell. Weighted sum
        of ideal and resistive contributions

      @param[in] cons pointer to conserved vector. Size is \f$N_{cons}\f$
      @param[in] prims pointer to primitive vector. Size is \f$N_{prims}\f$
      @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}\f$
      @param[in, out] source pointer to source vector. Size is \f$N_{cons}\f$
      @param i int cell number in x-direction (optional)
      @param j int cell number in y-direction (optional)
      @param k int cell number in z-direction (optional)
      @sa Model::sourceTermSingleCell
      @sa SRMHD::sourceTermSingleCell
      @sa SRRMHD::sourceTermSingleCell
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
      @par
        Source term calculation for all cells. Weighted sum of ideal and
        resistive contributions.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$
        @param[in, out] source pointer to source vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$
        @sa Model::sourceTermSingleCell
        @sa SRMHD::sourceTermSingleCell
        @sa SRRMHD::sourceTermSingleCell
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell conservative to primitive transformation
    /*!
      @par
        Conservative to primitive transformation for a single computational
        cell. If local conductivity is less than the crossover conductivity
        we use the resistive C2P, otherwise we use ideal.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}\f$
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}\f$
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}\f$
        @param i int cell number in x-direction (optional)
        @param j int cell number in y-direction (optional)
        @param k int cell number in z-direction (optional)
        @sa Model::getPrimitiveVarsSingleCell
        @sa SRMHD::getPrimitiveVarsSingleCell
        @sa SRRMHD::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);

    //! Conservative to primitive transformation for all cells
    /*!
      @par
        Conservative to primitive transformation for all cells. If local
        conductivity is less than the crossover conductivity we use the
        resistive C2P, otherwise we use ideal.

      @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @sa Model::getPrimitiveVarsSingleCell
      @sa SRMHD::getPrimitiveVarsSingleCell
      @sa SRRMHD::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
      @par
        Primitive to everything transformation. We weight the contributions
        from ideal and resistive MHD for this using the penalty function.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$
        @sa Model::primsToAll
        @sa SRMHD::primsToAll
        @sa SRRMHD::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
      @par
        Compute the flux vector for all cells. We weight the contributions
        from ideal and resistive MHD for this using the penalty function.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$
        @param[in, out] f pointer to flux vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
        @param dir direction in which to generate flux vector. \f$(x, y, z) = (0, 1, 2)\f$
        @sa Model::fluxVector
        @sa SRMHD::fluxVector
        @sa SRRMHD::fluxVector
    */
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);

    //! Finalise the simulation variables
    /*!
      @par
        At the end of each step, we need to calculate the ideal electric fields
        and add their contribution to hybrids electric fields.

        @param[in] cons pointer to conserved vector. Size is \f$N_{cons}*N_x*N_y*N_z\f$
        @param[in] prims pointer to primitive vector. Size is \f$N_{prims}*N_x*N_y*N_z\f$
        @param[in] aux pointer to auxiliary vector. Size is \f$N_{aux}*N_x*N_y*N_z\f$
    */
   
    // void finalise(double *cons, double *prims, double *aux);
    void finalise(double *cons, double *prims, double *aux, bool final_step=false);

};

#endif
