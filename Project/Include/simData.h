#ifndef SIMDATA_H
#define SIMDATA_H

class Data
/*
  Class contains all the data of the simulation relevant to any of the other
  modules. Containing it in this way prevents issues of cyclic includes, also
  results in Simulation as more of an interface than a class that needs to be
  known to lower objects---good practice.
*/
{
  public:
    int
    nx, ny;               // Number of physical cells in x- and y-direction
    double
    xmin, xmax,           // x-axis limits
    ymin, ymax,           // y-axis limits
    endTime,              // End time of simulation
    cfl;                  // Courant factor
    int Ng;               // Number of ghost cells
    double
    gamma,                // Adiabatic index
    sigma;                // Resistivity
    int
    memSet,              // Indicator that memory has been allocated for state vectors
    Ncons, Nprims, Naux;  // Number of conserved, primitive and auxilliary variables
    double
    cp,                   // Constant divergence cleaning term
    *cons, *prims, *aux,  // State vectors of conserved, primitive and auxilliary vars
    *f, *fnet,            // Flux vector and net numerical flux vector
    *source,              // Source vector
    *x, *y;               // x- and y-coordinate positions (incl ghost cells)
    double alphaX, alphaY,// Max wave speed in x and y direction
    t, dt,                // Current time and timestep
    dx, dy;               // Spatial x and y step
    int
    iters,                // Number of interations that have been completed
    Nx, Ny;               // Total number of compute cells in domain


    //! Element ID function
    /*!
        To access the 2nd conserved variable at (x, y) = (12, 4) for example,
      we call elem=data.id(2, 12, 4) and use this in d.cons[elem].
    */
    int id(int var, int i, int j) {
      return var * this->Nx * this->Ny + i * this->Ny + j;
    }

    //! Constructor
    /*!
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object.
    */
    Data(int nx, int ny,
         double xmin, double xmax,
         double ymin, double ymax,
         double endTime, double cfl=0.5, int Ng=4,
         double gamma=5.0/3.0, double sigma=0,
         int dataSet=0,
         int Ncons=0, int Nprims=0, int Naux=0,
         double cp=1.0);

};

#endif
