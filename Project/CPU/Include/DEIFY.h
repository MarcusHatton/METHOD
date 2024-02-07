#ifndef DEIFY_H
#define DEIFY_H

#include "modelExtension.h"
#include "flux.h"

//! <b> DEIFY: // Dissipative Extension for Ideal Fluid dYnamics </b>
/*!
    This class represents the implementation of DEIFY, a resistive extension
  to the special relativistic, ideal MHD equations. Details can be found in
  Wright & Hawke 2019 `A resistive extension to ideal MHD`.

    DEIFY extends the equations of ideal MHD by means of an additional, diffusive
  source term. The new system has the following form:
  \f{align}{
    \frac{\partial q}{\partial t} + \frac{\partial f^i(q)}{\partial x^i} = \frac{\partial D^i}{\partial x^i}
  \f}
  where the LHS corresponds to the standard, special relativistic ideal MHD (SRMHD)
  equations, and the RHD is the source term extension of DEIFY.

  The diffusion vector is defined by the following:
  \f{align}{
    D^i = -\frac{\partial f^i}{\partial \overline{q}} \bigg( \frac{\partial \overline{s}_0}{\partial \overline{q}} \bigg)^{-1} \frac{\partial \overline{f}^j_0}{\partial x^j}
  \f}
  where \f$ f^i \f$ is the \f$i^{th}\f$ direction, ideal MHD flux, \f$ \overline{q} \f$ is the stiff variables (electric fields of resistive MHD),
  \f$ \overline{s} \f$ is the source of the stiff variables and \f$\overline{f}^j\f$ is the \f$j^{th}\f$ flux of the stiff variables.

    Within this code, we use the following naming conventions:
    <ul>
     <li> `dfxdw`, `dfydw` and `dfzdw` are the derivatives of the \f$x\f$, \f$y\f$ and \f$z\f$ direction fluxes
    with respect to the primitive variables. </li>
     <li> `dwdsb` is the inverse of the stiff source vector with respect to the primitive variables. </li>
     <li> `Mx`, `My`, and `Mz` are the matrices \f$ -\frac{\partial f^i}{\partial w} \bigg( \frac{\partial \overline{s}_0}{\partial w} \bigg)^{-1}  \f$. </li>
     <li> `diffuX`, `diffuY`, and `diffuZ` are the \f$D^i\f$ vectors.
    </ul>

  To understand the elements of this extension, please view the paper. Including this
  extension is as simple as declaring `DEIFY modelExtension(&data, &fluxMethod);` in
  `main`, and including a pointer to `modelExtension` in the time integrator's
  constructor. With this, the SRMHD model will be able to capture resistive
  effects. You should play with the model to get a feel for how it behaves with
  \f$ \sigma \f$, and make sure that for high resolution simulations you check
  that the resolution criterion is met (available in the paper). You may need to
  reduce the CFL factor if simulations look unstable. Although simulations may
  converge with conductivities as low as \f$ \sigma = 10 \f$, I would be very
  careful for values less than 50. Conductivities of 100 or more are reliably
  converging, but remember for smaller conductivities, higher order moments
  may become important and solutions may differ from resistive MHD.

  Optimisations
  -------------
    To improve the performance of DEIFY, we have implemented a number of optimisations.
  These will make the source faster to compute, but may have led to a slightly
  less readable code (and less modular than I would have liked). We list some of
  these below:
  <ul>
    <li> First note that the vector elements D1, D2, and D3 inside dwdsb (Note:
  not the diffusion vector!) are not calculated. This is because in all cases
  the elements are multiplied by zero. </li>
    <li> When dotting `dfdw` with `dwdsb`, many of the elements are zero and so
  we manually multiply many of the non-zero entries to save on computation. </li>
    <li> A number of loops have been fused together to improve the memory
  management. </li>
  </ul>

  Examples
  --------
    The best way to understand this extension is to view the examples we have
  provided. Run the main programs with `make run` and observe the results using
  the `interactivePlot` script.
*/
class DEIFY : public ModelExtension
{
  public:

    // Bool controlling whether to add NLO contributions from the
    // CE-expanded fluxes (state vector contribution not yet implemented anyway)
    bool NLO = false;

    // Whether to use simple central differencing for gradients.
    // When false, used a MinMod slope-limiting scheme.
    bool SimpleCentreDifference = true;

    template<typename T>
    T sqr(T x) { return ((x) * (x)); }

    template<typename T>
    T cube(T x) { return ((x) * (x) * (x)); }

    // enums to save looking up numbering of C/P/As when using ID accessor.
    enum Cons { D, S1, S2, S3, Tau };
    enum Prims { v1, v2, v3, p, rho, n, q1, q2, q3, Pi, pi11, pi12, pi13, pi22, pi23, pi33 };
    enum Aux { h, T, e, W, q0, qv, pi00, pi01, pi02, pi03,
               q1NS, q2NS, q3NS, PiNS, pi11NS, pi12NS, pi13NS, pi22NS, pi23NS, pi33NS,
               q1LO, q2LO, q3LO, PiLO, pi11LO, pi12LO, pi13LO, pi22LO, pi23LO, pi33LO,  
               Theta, vsqrd, a1, a2, a3 };
    enum TDerivs { dtp = 35, dtrho, dtn, dtv1, dtv2, dtv3, dtW, dtT, dtq1NS, dtq2NS, dtq3NS, dtPiNS,
               dtpi11NS, dtpi12NS, dtpi13NS, dtpi22NS, dtpi23NS, dtpi33NS, dtD, dtS1, dtS2, dtS3,
               dtTau, dtE};

    double Grad(int enum_n, int direction, int C_P_or_A, double * cons, double * prims, double * aux, int i, int j, int k) {

      Data * d(this->data);

      int stencil[3] = {0, 0, 0};
      int stencil2[3] = {0, 0, 0};
      stencil[direction] += 1;
      stencil2[direction] += 2;
      double dX[3] = {data->dx, data->dy, data->dz};
      double var_cen, var_fw, var_fw2, var_bw, var_bw2;

      if (C_P_or_A == 0) {
        var_cen = cons[ID(enum_n, i, j, k)];
        var_fw = cons[ID(enum_n, i+stencil[0], j+stencil[1], k+stencil[2])];
        var_fw2 = cons[ID(enum_n, i+stencil2[0], j+stencil2[1], k+stencil2[2])];
        var_bw = cons[ID(enum_n, i-stencil[0], j-stencil[1], k-stencil[2])];
        var_bw2 = cons[ID(enum_n, i-stencil2[0], j-stencil2[1], k-stencil2[2])];
      }
      else if (C_P_or_A == 1) {
        var_cen = prims[ID(enum_n, i, j, k)];
        var_fw = prims[ID(enum_n, i+stencil[0], j+stencil[1], k+stencil[2])];
        var_fw2 = prims[ID(enum_n, i+stencil2[0], j+stencil2[1], k+stencil2[2])];
        var_bw = prims[ID(enum_n, i-stencil[0], j-stencil[1], k-stencil[2])];
        var_bw2 = prims[ID(enum_n, i-stencil2[0], j-stencil2[1], k-stencil2[2])];
      }
      else if (C_P_or_A == 2) {
        var_cen = aux[ID(enum_n, i, j, k)];
        var_fw = aux[ID(enum_n, i+stencil[0], j+stencil[1], k+stencil[2])];
        var_fw2 = aux[ID(enum_n, i+stencil2[0], j+stencil2[1], k+stencil2[2])];
        var_bw = aux[ID(enum_n, i-stencil[0], j-stencil[1], k-stencil[2])];
        var_bw2 = aux[ID(enum_n, i-stencil2[0], j-stencil2[1], k-stencil2[2])];
      }
      
      //double CDGrad = (var_fw - var_bw)/(2*dX[direction]);
      double CDGrad = (-0.08333*var_fw2 + 0.6667*var_fw - 0.6667*var_bw + 0.08333*var_bw2)/(dX[direction]);

      // Min-Mod First-Order
      double FDGrad = (-1.0*var_cen + 1.0*var_fw)/dX[direction];
      double BDGrad = (1.0*var_cen - 1.0*var_bw)/dX[direction];
      //printf("CDGrad: %f", CDGrad);
      //printf("FDGrad: %f", FDGrad);
      //printf("BDGrad: %f", BDGrad);

      if (SimpleCentreDifference) {
        return CDGrad;
      }

      printf("EEK EEK EEK");

      if ( (FDGrad < 0 && BDGrad > 0) || (FDGrad > 0 && BDGrad < 0) ) {
        return 0;
      } else {
        return abs(FDGrad) < abs(BDGrad) ? FDGrad : BDGrad;
      }
    }

    double *Fx, *Fy, *Fz;   //!< Diffusion vectors

    double *dtH;
    
    FluxMethod * fluxMethod;     //!< Pointer to the flux method class

    DEIFY();

    //! Constructor
    DEIFY(Data * data, FluxMethod * fluxMethod);

    //! Destructor
    virtual ~DEIFY();

    //! Main user function for the modified source
    /*!
        This method implements the modified source for DEIFY. Given the current
      state of the system in cons, prims and aux, it writes into source the
      contribution from the derivative of the diffusion vector.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *source pointer to source vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @sa ModelExtension
    */
    void sourceExtension(double * cons, double * prims, double * aux, double * source);

    //! Sets up variables including the electric field and charge density
    void set_vars(double * cons, double * prims, double * aux);

    //{
    //! Set the diffusion vectors
    void set_Fx(double * cons, double * prims, double * aux);
    void set_Fy(double * cons, double * prims, double * aux);
    void set_Fz(double * cons, double * prims, double * aux);
    //}

    //{
    //! Set the time derivative of the state vector's NS contribution,
    //! purely in terms of spatial derivatives
    void set_dtH(double * cons, double * prims, double * aux);
    //}

};

#endif
