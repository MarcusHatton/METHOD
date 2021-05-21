#ifndef IS_H
#define IS_H

#include "model.h"
#include "simData.h"

//! Full IS model for (relativistic) non-ideal fluids
/*!
  @par
    Full IS model with no second-order terms in source vector
*/
class IS : public Model
{
  public:
    IS();
    IS(Data * data);
    virtual ~IS();

    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);
    void sourceTerm(double *cons, double *prims, double *aux, double *source);
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);
    void getPrimitiveVars(double *cons, double *prims, double *aux);
    void primsToAll(double *cons, double *prims, double *aux);
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);
    void finalise(double *cons, double *prims, double *aux) { };

};

#endif
