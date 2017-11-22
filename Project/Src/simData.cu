#include "simData.h"
#include <stdexcept>

Data::Data(int nx, int ny, int nz,
           double xmin, double xmax,
           double ymin, double ymax,
           double zmin, double zmax,
           double endTime, double cfl, int Ng,
           double gamma, double sigma,
           double cp,
           double mu1, double mu2)
           :
           nx(nx), ny(ny), nz(nz),
           xmin(xmin), xmax(xmax),
           ymin(ymin), ymax(ymax),
           zmin(zmin), zmax(zmax),
           endTime(endTime), cfl(cfl), Ng(Ng),
           gamma(gamma), sigma(sigma),
           memSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(cp),
           mu1(mu1), mu2(mu2)
{
  // Ensure even number of cells to prevent zero in initial data at centre of domain
  if (nx%2 || ny%2 || nz%2) {
    throw std::invalid_argument("Please ensure even number of cells in each direction.\n");
  }


  this->Nx = nx + 2 * Ng;
  this->Ny = ny + 2 * Ng;
  this->Nz = nz + 2 * Ng;

  // Catch 2D case
  if (nz == 0) {
    this->Nz = 1;
    zmin = -1e20;
    zmax = 1e20;
  }
  // Catch 1D case
  if (ny == 0) {
    this->Nz = this->Ny = 1;
    zmin = ymin = -1e20;
    zmax = ymax = 1e20;
  }
  // Ensure there is some Resistivity
  if (this->sigma <= 0.0) {
    throw std::invalid_argument("Conductivity must be positive, sigma > 0.\n");
  }
}
