#include "SubgridNS.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

template<typename T>
T sqr(T x) { return ((x) * (x)); }

NS::NS() : Model()
{
  this->Ncons = 5;
  this->Nprims = 16;
  this->Naux = 52;
}

NS::NS(Data * data, bool alt_C2P=false) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 5;
  this->Nprims = (this->data)->Nprims = 16;
  this->Naux = (this->data)->Naux = 52;

  // Solutions for C2P all cells
  solution = (double *) malloc(sizeof(double)*4*data->Nx*data->Ny*data->Nz);

  // Vector for storing variable at previous time-step...
  // the 4 here is for the 4 time-deriv variables currently needed... should be automated really not hard-set
  prev_vars = (double *) malloc(sizeof(double)*5*data->Nx*data->Ny*data->Nz); 

  smartGuesses = 4;
  
  alternative_C2P = alt_C2P;

  scale_ratio = data->optionalSimArgs[6];
  
  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("S1");
  this->data->consLabels.push_back("S2");  this->data->consLabels.push_back("S3");
  this->data->consLabels.push_back("Tau");  

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
  this->data->auxLabels.push_back("theta");  this->data->auxLabels.push_back("dv1dt");
  this->data->auxLabels.push_back("dv2dt");  this->data->auxLabels.push_back("dv3dt");
  this->data->auxLabels.push_back("a0");     this->data->auxLabels.push_back("a1");   
  this->data->auxLabels.push_back("a2");     this->data->auxLabels.push_back("a3");     
  this->data->auxLabels.push_back("vsqrd"); this->data->auxLabels.push_back("dWdt");  
  this->data->auxLabels.push_back("dTdt");  this->data->auxLabels.push_back("rho_plus_p");
  // 22 - Shear
  this->data->auxLabels.push_back("sigma00"); this->data->auxLabels.push_back("sigma01");
  this->data->auxLabels.push_back("sigma02"); this->data->auxLabels.push_back("sigma03");
  this->data->auxLabels.push_back("sigma11");  this->data->auxLabels.push_back("sigma12");
  this->data->auxLabels.push_back("sigma13");  this->data->auxLabels.push_back("sigma22");
  this->data->auxLabels.push_back("sigma23");  this->data->auxLabels.push_back("sigma33");  
  this->data->auxLabels.push_back("sigmasqrd");  this->data->auxLabels.push_back("detsigma");
  // 34 - Vorticity
  this->data->auxLabels.push_back("omega00"); this->data->auxLabels.push_back("omega01");
  this->data->auxLabels.push_back("omega02"); this->data->auxLabels.push_back("omega03");
  this->data->auxLabels.push_back("omega11");  this->data->auxLabels.push_back("omega12");
  this->data->auxLabels.push_back("omega13");  this->data->auxLabels.push_back("omega22");
  this->data->auxLabels.push_back("omega23");  this->data->auxLabels.push_back("omega33");
  this->data->auxLabels.push_back("omegasqrd");
  // 45 - Heat/Acc. Combo
  this->data->auxLabels.push_back("Theta0");   this->data->auxLabels.push_back("Theta1");
  this->data->auxLabels.push_back("Theta2");   this->data->auxLabels.push_back("Theta3");
  // 49 - Coefficients
  this->data->auxLabels.push_back("zeta");  this->data->auxLabels.push_back("kappa");
  this->data->auxLabels.push_back("eta");
 
}

NS::~NS()
{
  free(solution);
}

// First and second order "minmod" functions for slope-limiting
double minmodGradFO(double im1, double i, double ip1, double dX) {

  double FDGrad = (-1.0*i + 1*ip1)/dX;
  double BDGrad = (1.0*i - 1*im1)/dX;
  if ( (FDGrad < 0 && BDGrad > 0) || (FDGrad > 0 && BDGrad < 0) ) {
    return 0;
  } else {
    return abs(FDGrad) < abs(BDGrad) ? FDGrad : BDGrad;
  }
}

double minmodGradSO(double im2, double im1, double i, double ip1, double ip2, double dX) {
  
  double FDGrad = (-1.5*i + 2*ip1 - 0.5*ip2)/dX;
  double BDGrad = (1.5*i - 2*im1 + 0.5*im2)/dX;
  if ( (FDGrad < 0 && BDGrad > 0) || (FDGrad > 0 && BDGrad < 0) ) {
    return 0;
  } else {
    return abs(FDGrad) < abs(BDGrad) ? FDGrad : BDGrad;
  }
}

double CalcDeterminant(std::vector<std::vector<double>> Matrix)
{
  //this function is written in c++ to calculate the determinant of matrix
  // it's a recursive function that can handle matrix of any dimension
  double det = 0; // the determinant value will be stored here
  if (Matrix.size() == 1)
  {
      return Matrix[0][0]; // no calculation needed
  }
  else if (Matrix.size() == 2)
  {
      //in this case we calculate the determinant of a 2-dimensional matrix in a 
      //default procedure
      det = (Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0]);
      return det;
  }
  else
  {
      //in this case we calculate the determinant of a squared matrix that have 
      // for example 3x3 order greater than 2
      for (int p = 0; p < Matrix[0].size(); p++)
      {
          //this loop iterate on each elements of the first row in the matrix.
          //at each element we cancel the row and column it exist in
          //and form a matrix from the rest of the elements in the matrix
          std::vector<std::vector<double>> TempMatrix; // to hold the shaped matrix;
          for (int i = 1; i < Matrix.size(); i++)
          {
              // iteration will start from row one cancelling the first row values
              std::vector<double> TempRow;
              for (int j = 0; j < Matrix[i].size(); j++)
              {
                  // iteration will pass all cells of the i row excluding the j 
                  //value that match p column
                  if (j != p)
                  {
                     TempRow.push_back(Matrix[i][j]);//add current cell to TempRow 
                  }
              }
              if (TempRow.size() > 0)
                  TempMatrix.push_back(TempRow);
              //after adding each row of the new matrix to the vector tempx
              //we add it to the vector temp which is the vector where the new 
              //matrix will be formed
          }
          det = det + Matrix[0][p] * pow(-1, p) * CalcDeterminant(TempMatrix);
          //then we calculate the value of determinant by using a recursive way
          //where we re-call the function by passing to it the new formed matrix
          //we keep doing this until we get our determinant
      }
//      for (int i = 0; i < Matrix.size(); i++) {
//        for (int j = 0; j < Matrix.size(); j++) {
//          if (Matrix[i][j] != 0.0) {
//          printf("(%i, %i, %g) i, j, k: \n", i, j, Matrix[i][j]);
//          printf("(%g) det:\n", det);
//          exit(0);
//          }
//        }
//      }
      return det;
  }
};

double calculate4Determinant(double** matrix) {

  double a = matrix[0][0];
  double b = matrix[0][1];
  double c = matrix[0][2];
  double d = matrix[0][3];
  double e = matrix[1][1];
  double f = matrix[1][2];
  double g = matrix[1][3];
  double h = matrix[2][2];
  double i = matrix[2][3];
  double j = matrix[3][3];

  double det = a*e*h*j - (a*e*i*i + a*h*g*g + a*j*f*f + e*h*d*d + e*j*c*c + h*j*b*b)
  + 2*(a*f*g*i * e*c*d*i + h*b*d*g + j*b*c*f) - 2*(b*f*i*d + b*c*g*i + c*g*f*d)
  + (b*i)*(b*i) + (c*g)*(c*g) + (d*f)*(d*f);
  
  return det;

}


// Function to calculate the determinant
// of a matrix
double calculateDeterminant(double** matrix, int size)
{
    double det = 0;
    int sign = 1;
 
    // Base Case
    if (size == 1) {
        det = matrix[0][0];
    }
    else if (size == 2) {
        det = (matrix[0][0] * matrix[1][1])
              - (matrix[0][1] * matrix[1][0]);
    }
 
    // Perform the Laplace Expansion
    else {
        for (int i = 0; i < size; i++) {
 
            // Stores the cofactor matrix
            double** cofactor = new double*[size - 1];
            for (int j = 0; j < size - 1; j++) {
                cofactor[j] = new double[size - 1];
            }
            int sub_i = 0, sub_j = 0;
            for (int j = 1; j < size; j++) {
                for (int k = 0; k < size; k++) {
                    if (k == i) {
                        continue;
                    }
                    cofactor[sub_i][sub_j] = matrix[j][k];
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }
 
            // Update the determinant value
            det += sign * matrix[0][i]
                   * calculateDeterminant(cofactor, size - 1);
            sign = -sign;
            for (int j = 0; j < size - 1; j++) {
                delete[] cofactor[j];
            }
            delete[] cofactor;
        }
    }
 
    // Return the final determinant value
    return det;
}


void NS::calculateDissipativeCoefficients(double *cons, double *prims, double *aux)
{
  Data * d(this->data);

  // scale_ratio = this->scale_ratio;
  
  scale_ratio = 800 / this->data->nx; // should be using this as calibrated on 800x800

  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int k(0); k < this->data->Nz; k++) {
            aux[ID(Aux::zeta, i, j, k)] = pow(this->scale_ratio, 2) * pow(10,-5.7) * pow(abs(aux[ID(Aux::omegasqrd, i, j, k)]), 0.1) * pow(aux[ID(Aux::T, i, j, k)], 0.4) * pow(prims[ID(Prims::n, i, j, k)], 0.5)
                                          *  pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] - aux[ID(Aux::omegasqrd, i, j, k)]), 0.46) * pow(abs(aux[ID(Aux::theta, i, j, k)]), -0.85);
            aux[ID(Aux::kappa, i, j, k)] = pow(this->scale_ratio, 2) * pow(10,-6.3) * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)]), 0.15) * pow(prims[ID(Prims::n, i, j, k)], 0.3)
                                          * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] - aux[ID(Aux::omegasqrd, i, j, k)]), 0.23); // * ...
            aux[ID(Aux::eta, i, j, k)] = pow(this->scale_ratio, 2) * pow(10,-3.4) * pow(abs(aux[ID(Aux::detsigma, i, j, k)]), 0.15)
                                          * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] - aux[ID(Aux::omegasqrd, i, j, k)]), 0.045)
                                          * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] / aux[ID(Aux::omegasqrd, i, j, k)]), -0.13);

            //printf("(%g, %g, %g) zeta, kappa, eta\n", aux[ID(Aux::zeta, i, j, k)], aux[ID(Aux::kappa, i, j, k)], aux[ID(Aux::eta, i, j, k)]);

            if( isnan(aux[ID(Aux::zeta, i, j, k)]) ) {
              aux[ID(Aux::zeta, i, j, k)] = 0.0;
            } if( isnan(aux[ID(Aux::kappa, i, j, k)]) ) {
              aux[ID(Aux::kappa, i, j, k)] = 0.0;
            } if( isnan(aux[ID(Aux::eta, i, j, k)]) ) {
              aux[ID(Aux::eta, i, j, k)] = 0.0;
            } 

            // Check for infinities - from negative powers of terms that are 0
            if( isinf(aux[ID(Aux::zeta, i, j, k)]) ) {
              aux[ID(Aux::zeta, i, j, k)] = 0.0; // Should probably have a better solution than setting them to zero...
            } if( isinf(aux[ID(Aux::kappa, i, j, k)]) ) {
              aux[ID(Aux::kappa, i, j, k)] = 0.0;
            } if( isinf(aux[ID(Aux::eta, i, j, k)]) ) {
              aux[ID(Aux::eta, i, j, k)] = 0.0;
            } 

      }
    }
  }

}

void NS::calculateDissipativeCoefficientsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  Data * d(this->data);

  // scale_ratio = this->scale_ratio;
  
  scale_ratio = 800 / this->data->nx; // should be using this as calibrated on 800x800

  aux[ID(Aux::zeta, i, j, k)] = pow(this->scale_ratio, 2) * pow(10,-5.7) * pow(abs(aux[ID(Aux::omegasqrd, i, j, k)]), 0.1) * pow(aux[ID(Aux::T, i, j, k)], 0.4) * pow(prims[ID(Prims::n, i, j, k)], 0.5)
                                *  pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] - aux[ID(Aux::omegasqrd, i, j, k)]), 0.46) * pow(abs(aux[ID(Aux::theta, i, j, k)]), -0.85);
  aux[ID(Aux::kappa, i, j, k)] = pow(this->scale_ratio, 2) * pow(10,-6.3) * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)]), 0.15) * pow(prims[ID(Prims::n, i, j, k)], 0.3)
                                * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] - aux[ID(Aux::omegasqrd, i, j, k)]), 0.23); // * ...
  aux[ID(Aux::eta, i, j, k)] = pow(this->scale_ratio, 2) * pow(10,-3.4) * pow(abs(aux[ID(Aux::detsigma, i, j, k)]), 0.15)
                                * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] - aux[ID(Aux::omegasqrd, i, j, k)]), 0.045)
                                * pow(abs(aux[ID(Aux::sigmasqrd, i, j, k)] / aux[ID(Aux::omegasqrd, i, j, k)]), -0.13);

  // printf("(%g, %g, %g) zeta, kappa, eta\n", aux[ID(Aux::zeta, i, j, k)], aux[ID(Aux::kappa, i, j, k)], aux[ID(Aux::eta, i, j, k)]);
    //printf("(%d, %d, %d) i, j, k\n", i, j, k);
    //printf("(%g, %g, %g) zeta, kappa, eta\n", aux[ID(Aux::zeta, i, j, k)], aux[ID(Aux::kappa, i, j, k)], aux[ID(Aux::eta, i, j, k)]);
    //printf("(%g, %g, %g, %g) omegasqrd, sigmasqrd, theta, detsigma \n", aux[ID(Aux::omegasqrd, i, j, k)], aux[ID(Aux::sigmasqrd, i, j, k)], aux[ID(Aux::theta, i, j, k)], aux[ID(Aux::detsigma, i, j, k)]);
    //printf("(%g, %g, %g) zeta, kappa, eta\n", aux[ID(Aux::detsigma, i, j, k)], aux[ID(Aux::kappa, i, j, k)], aux[ID(Aux::eta, i, j, k)]);

    //exit(0);

  // Check for nan values - from fractional powers of terms that are negative
  if( isnan(aux[ID(Aux::zeta, i, j, k)]) ) {
    aux[ID(Aux::zeta, i, j, k)] = 0.0;
  } if( isnan(aux[ID(Aux::kappa, i, j, k)]) ) {
    aux[ID(Aux::kappa, i, j, k)] = 0.0;
  } if( isnan(aux[ID(Aux::eta, i, j, k)]) ) {
    aux[ID(Aux::eta, i, j, k)] = 0.0;
  } 

  // Check for infinities - from negative powers of terms that are 0
  if( isinf(aux[ID(Aux::zeta, i, j, k)]) ) {
    aux[ID(Aux::zeta, i, j, k)] = 0.0; // Should probably have a better solution than setting them to zero...
  } if( isinf(aux[ID(Aux::kappa, i, j, k)]) ) {
    aux[ID(Aux::kappa, i, j, k)] = 0.0;
  } if( isinf(aux[ID(Aux::eta, i, j, k)]) ) {
    aux[ID(Aux::eta, i, j, k)] = 0.0;
  } 

}

void NS::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{

  // D
  source[0] = 0.0;
  // S1,2,3
  source[1] = 0.0; 
  source[2] = 0.0;
  source[3] = 0.0; 
  // Tau
  source[4] = 0.0;
 
}

void NS::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int k(0); k < this->data->Nz; k++) {
        // D
        source[ID(D, i, j, k)] = 0.0;
        // S1,2,3
        source[ID(S1, i, j, k)] = 0.0;
        source[ID(S2, i, j, k)] = 0.0;
        source[ID(S3, i, j, k)] = 0.0;
        // Tau
        source[ID(Tau, i, j, k)] = 0.0;
        
      }
    }
  }

}

//! Residual function to minimize in the format required by cminpack
/*!
    I know, its a horrible layout, alas we're stuck with it.

    Parameters
    ----------
    p : pointer to struct
      Struct contains additional arguments that are required (if any)
    n : int
      Size of system
    x : pointer to double
      The array containing the initial guess
    fvec : pointer to double
      The array containing the solution
    iflag : int
      Error flag
*/
int NSresidual(void *ptr, int n, const double *x, double *fvec, int iflag)
{

//  Data * d(this->data);
  
  // Retrieve additional arguments
  NS::Args * args = (NS::Args*) ptr;

  // Values must make sense
  // Think this should be taken out for now - need new sensible values
  /*
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  */
  
  double E_rf = args->Tau_rf + args->D_rf; 
  double vsqrd_rf = ((args->S1_rf - x[1])*(args->S1_rf - x[1]) + (args->S2_rf - x[2])*(args->S2_rf - x[2]) + (args->S3_rf - x[3])*(args->S3_rf - x[3]))/((E_rf + x[0])*(E_rf + x[0]));
  double W_rf(1 / sqrt(1 - vsqrd_rf));
  double n_rf(args->D_rf / W_rf);
  double rho_plus_p_rf = ((E_rf + x[0])/(W_rf*W_rf)) - args->Pi_rf;
  double v1_rf = (args->S1_rf - x[1])/((rho_plus_p_rf + args->Pi_rf)*W_rf*W_rf);
  double v2_rf = (args->S2_rf - x[2])/((rho_plus_p_rf + args->Pi_rf)*W_rf*W_rf);
  double v3_rf = (args->S3_rf - x[3])/((rho_plus_p_rf + args->Pi_rf)*W_rf*W_rf);
  double p_rf = (rho_plus_p_rf - n_rf)*((args->gamma-1)/args->gamma);
  double rho_rf = rho_plus_p_rf - p_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 1 || n_rf < 0 || abs(v1_rf) >= 1 || abs(v2_rf) >= 1 || abs(v3_rf) >= 1 || vsqrd_rf >= 1) {
    printf("EEK");
    fvec[0] = fvec[1] = fvec[2] = fvec[3] = 1e6;
    return 0;
  }
  
  double pi00_rf = args->pi11_rf + args->pi22_rf + args->pi33_rf;
  double qv_rf = args->q1_rf*v1_rf + args->q2_rf*v2_rf + args->q3_rf*v3_rf;
  double pi01_rf = args->pi11_rf*v1_rf + args->pi12_rf*v2_rf + args->pi13_rf*v3_rf; // dbl check sign on orthogonality relation
  double pi02_rf = args->pi12_rf*v1_rf + args->pi22_rf*v2_rf + args->pi23_rf*v3_rf;
  double pi03_rf = args->pi13_rf*v1_rf + args->pi23_rf*v2_rf + args->pi33_rf*v3_rf;

  fvec[0] = p_rf + args->Pi_rf - 2*qv_rf*W_rf - pi00_rf - x[0];
  fvec[1] = (args->q1_rf + qv_rf*v1_rf)*W_rf + pi01_rf - x[1];
  fvec[2] = (args->q2_rf + qv_rf*v2_rf)*W_rf + pi02_rf - x[2];
  fvec[3] = (args->q3_rf + qv_rf*v3_rf)*W_rf + pi03_rf - x[3];

  return 0;
}

int NSAlternativeResidual(void *ptr, int n, const double *x, double *fvec, int iflag)
{

//  Data * d(this->data);
  
  // Retrieve additional arguments
  NS::Args * args = (NS::Args*) ptr;

  // Values must make sense
  // Think this should be taken out for now - need new sensible values
  /*
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  */
  
  double vsqrd_rf = x[0]*x[0]*((args->S1_rf - x[1])*(args->S1_rf - x[1]) + (args->S2_rf - x[2])*(args->S2_rf - x[2]) + (args->S3_rf - x[3])*(args->S3_rf - x[3]))/((args->D_rf)*(args->D_rf));
  double W_rf(1 / sqrt(1 - vsqrd_rf));
  double n_rf(args->D_rf / W_rf);
  double v1_rf = x[0]*(args->S1_rf - x[1])/args->D_rf;
  double v2_rf = x[0]*(args->S2_rf - x[2])/args->D_rf;
  double v3_rf = x[0]*(args->S3_rf - x[3])/args->D_rf;
  double pi00_rf = args->pi11_rf + args->pi22_rf + args->pi33_rf;
  double qv_rf = args->q1_rf*v1_rf + args->q2_rf*v2_rf + args->q3_rf*v3_rf;
  double p_rf = args->D_rf*((1/x[0]) -1) - args->Pi_rf + 2*qv_rf*W_rf + pi00_rf - args->Tau_rf;
  double rho_rf = n_rf + p_rf/(args->gamma-1);
  double H_rf = 1 + (p_rf*(args->gamma/(args->gamma-1)) + args->Pi_rf)/n_rf;

  // Values should be sensible    
  if (p_rf < 0 || rho_rf < 0 || W_rf < 1 || n_rf < 0 || abs(v1_rf) >= 1 || abs(v2_rf) >= 1 || abs(v3_rf) >= 1 || vsqrd_rf >= 1) {
    printf("EEK");
    fvec[0] = fvec[1] = fvec[2] = fvec[3] = 1e6;
    return 0;
  }
  
  double pi01_rf = args->pi11_rf*v1_rf + args->pi12_rf*v2_rf + args->pi13_rf*v3_rf; // dbl check sign on orthogonality relation
  double pi02_rf = args->pi12_rf*v1_rf + args->pi22_rf*v2_rf + args->pi23_rf*v3_rf;
  double pi03_rf = args->pi13_rf*v1_rf + args->pi23_rf*v2_rf + args->pi33_rf*v3_rf;

  fvec[0] = 1/(W_rf*H_rf) - x[0];
  fvec[1] = (args->q1_rf + qv_rf*v1_rf)*W_rf + pi01_rf - x[1];
  fvec[2] = (args->q2_rf + qv_rf*v2_rf)*W_rf + pi02_rf - x[2];
  fvec[3] = (args->q3_rf + qv_rf*v3_rf)*W_rf + pi03_rf - x[3];

  return 0;
}

void NS::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  Data * d(this->data);

  // Do what we can first before root-find
  // No calculation of diss prims being done for single-cell currently...

  aux[Aux::pi00] = prims[Prims::pi11] + prims[Prims::pi22] + prims[Prims::pi33]; // this one can be done here fine
  // what about these? need them in the guesses...
  aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
  aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
  aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation

  // Hybrd1 set-up
  Args args;                      // Additional arguments structure
  const int sys_size(4);                     // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1e-6;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Set additional args for rootfind
  args.D_rf = cons[Cons::D];
  args.S1_rf = cons[Cons::S1];
  args.S2_rf = cons[Cons::S2];
  args.S3_rf = cons[Cons::S3];
  args.Tau_rf = cons[Cons::Tau];
  args.q1_rf = prims[Prims::q1];
  args.q2_rf = prims[Prims::q2];
  args.q3_rf = prims[Prims::q3];
  args.Pi_rf = prims[Prims::Pi];
  args.pi11_rf = prims[Prims::pi11];
  args.pi12_rf = prims[Prims::pi12];
  args.pi13_rf = prims[Prims::pi13];
  args.pi22_rf = prims[Prims::pi22];
  args.pi23_rf = prims[Prims::pi23];
  args.pi33_rf = prims[Prims::pi33];
  args.gamma = d->gamma;
  
  if (alternative_C2P) {
  
    sol[0] = 1/(aux[Aux::W]*(1 + (prims[Prims::p]*(d->gamma/(d->gamma-1)) + prims[Prims::Pi])/prims[Prims::n]));
    sol[1] = (prims[Prims::q1] + aux[Aux::qv]*prims[Prims::v1])*aux[Aux::W] + aux[Aux::pi01];
    sol[2] = (prims[Prims::q2] + aux[Aux::qv]*prims[Prims::v2])*aux[Aux::W] + aux[Aux::pi02];
    sol[3] = (prims[Prims::q3] + aux[Aux::qv]*prims[Prims::v3])*aux[Aux::W] + aux[Aux::pi03];
  
    // Solve residual = 0
    info = __cminpack_func__(hybrd1) (&NSAlternativeResidual, &args, sys_size, sol, res,
                                      tol, wa, lwa);
    // If root find fails, add failed cell to the list
    if (info!=1) {
      //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
      throw std::runtime_error("C2P could not converge.\n");
    }
    aux[Aux::vsqrd] = sol[0]*sol[0]*((cons[Cons::S1] - sol[1])*(cons[Cons::S1] - sol[1]) + (cons[Cons::S2] - sol[2])*(cons[Cons::S2] - sol[2]) 
                      + (cons[Cons::S3] - sol[3])*(cons[Cons::S3] - sol[3]))/(cons[Cons::D]*cons[Cons::D]);
    aux[Aux::W] = (1 / sqrt(1 - aux[Aux::vsqrd]));
    prims[Prims::n] = cons[Cons::D] / aux[Aux::W];
    prims[Prims::v1] = sol[0]*(cons[Cons::S1] - sol[1])/cons[Cons::D];
    prims[Prims::v2] = sol[0]*(cons[Cons::S2] - sol[2])/cons[Cons::D];
    prims[Prims::v3] = sol[0]*(cons[Cons::S3] - sol[3])/cons[Cons::D];
    aux[Aux::pi00] = prims[Prims::pi11] + prims[Prims::pi22] + prims[Prims::pi33]; // not sure we need this again here tbh
    aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
    prims[Prims::p] = cons[Cons::D]*(1/sol[0] -1) - prims[Prims::Pi] + 2*aux[Aux::qv]*aux[Aux::W] + aux[Aux::pi00] - cons[Cons::Tau];
    prims[Prims::rho] = prims[Prims::n] + prims[Prims::p]/(d->gamma-1);
    
    // Repeating the ones here that depend on v1,v2,v3...
    aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
          
    aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
    aux[Aux::T] = prims[Prims::p] / prims[Prims::n]; 
    
  
  } else {
  
    sol[0] = prims[Prims::p] + prims[Prims::Pi] - 2*aux[Aux::qv]*aux[Aux::W] - aux[Aux::pi00];
    sol[1] = (prims[Prims::q1] + aux[Aux::qv]*prims[Prims::v1])*aux[Aux::W] + aux[Aux::pi01];
    sol[2] = (prims[Prims::q2] + aux[Aux::qv]*prims[Prims::v2])*aux[Aux::W] + aux[Aux::pi02];
    sol[3] = (prims[Prims::q3] + aux[Aux::qv]*prims[Prims::v3])*aux[Aux::W] + aux[Aux::pi03];
  
    // Solve residual = 0
    info = __cminpack_func__(hybrd1) (&NSresidual, &args, sys_size, sol, res,
                                      tol, wa, lwa);
    // If root find fails, add failed cell to the list
    if (info!=1) {
      //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
      throw std::runtime_error("C2P could not converge.\n");
    }
    aux[Aux::vsqrd] = ((cons[Cons::S1] - sol[1])*(cons[Cons::S1] - sol[1]) + (cons[Cons::S2] - sol[2])*(cons[Cons::S2] - sol[2]) 
                       + (cons[Cons::S3] - sol[3])*(cons[Cons::S3] - sol[3]))
                       /((cons[Cons::Tau] + cons[Cons::D] + sol[0])*(cons[Cons::Tau]  + cons[Cons::D] + sol[0]));
    aux[Aux::W] = 1 / sqrt((1-aux[Aux::vsqrd]));
    prims[Prims::n] = cons[Cons::D] / aux[Aux::W];
    aux[Aux::rho_plus_p] = (cons[Cons::Tau] + cons[Cons::D] + sol[0])/(aux[Aux::W]*aux[Aux::W]) - prims[Prims::Pi];
    prims[Prims::v1] = (cons[Cons::S1] - sol[1])/((aux[Aux::rho_plus_p] + prims[Prims::Pi])*aux[Aux::W]*aux[Aux::W]);
    prims[Prims::v2] = (cons[Cons::S2] - sol[2])/((aux[Aux::rho_plus_p] + prims[Prims::Pi])*aux[Aux::W]*aux[Aux::W]);  
    prims[Prims::v3] = (cons[Cons::S3] - sol[3])/((aux[Aux::rho_plus_p] + prims[Prims::Pi])*aux[Aux::W]*aux[Aux::W]);  
    prims[Prims::p] = (aux[Aux::rho_plus_p] - prims[Prims::n])*((d->gamma-1)/d->gamma);
    prims[Prims::rho] = aux[Aux::rho_plus_p] - prims[Prims::p];
    
    // Repeating the ones here that depend on v1,v2,v3...
    aux[Aux::qv] = prims[Prims::q1]*prims[Prims::v1] + prims[Prims::q2]*prims[Prims::v2] + prims[Prims::q3]*prims[Prims::v3];
    aux[Aux::pi01] = prims[Prims::pi11]*prims[Prims::v1] + prims[Prims::pi12]*prims[Prims::v2] + prims[Prims::pi13]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi02] = prims[Prims::pi12]*prims[Prims::v1] + prims[Prims::pi22]*prims[Prims::v2] + prims[Prims::pi23]*prims[Prims::v3]; // dbl check sign on orthogonality relation
    aux[Aux::pi03] = prims[Prims::pi13]*prims[Prims::v1] + prims[Prims::pi23]*prims[Prims::v2] + prims[Prims::pi33]*prims[Prims::v3]; // dbl check sign on orthogonality relation
          
    aux[Aux::e] = prims[Prims::p] / (prims[Prims::n]*(d->gamma-1));
    aux[Aux::T] = prims[Prims::p] / prims[Prims::n];     
 
  }
   
}

void NS::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int sys_size(4);              // Size of system
  double sol[sys_size];                      // Guess and solution vector
  double res[sys_size];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1e-4;          // Tolerance of rootfinder
  const int lwa = 50;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array
  std::vector<Failed> fails;          // Vector of failed structs. Stores location of failed cons2prims cells.
  
  /* - fix n->sys_size if ever used
  int maxfev(50);
  int ml(n);
  int mu(n);
  double epsfcn(0.1);
  double diag[n];
  int mode(1);
  double factor(100);
  int nprint(-1);
  int nfev(0);
  double fjac[n][n];
  int ldfjac(n);
  int lr(10);
  double r[lr];
  double qtf[n];
  double wa1[n];                     // Work array
  double wa2[n];                     // Work array
  double wa3[n];                     // Work array
  double wa4[n];                     // Work array
  */
  
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        //aux[ID(Aux::W, i, j, k)] = 1 / ( 1 - (prims[ID(Prims::v1, i, j, k)]**2 + prims[ID(Prims::v2, i, j, k)]**2 + prims[ID(Prims::v3, i, j, k)]**2) )**0.5; // Point in re-calcing this here?
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) 
                               + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
        aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)];
        aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
        aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                 + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                 + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
            
        // Set additional args for rootfind
        args.D_rf = cons[ID(Cons::D, i, j, k)];
        args.S1_rf = cons[ID(Cons::S1, i, j, k)];
        args.S2_rf = cons[ID(Cons::S2, i, j, k)];
        args.S3_rf = cons[ID(Cons::S3, i, j, k)];
        args.Tau_rf = cons[ID(Cons::Tau, i, j, k)];
        args.q1_rf = prims[ID(Prims::q1, i, j, k)];
        args.q2_rf = prims[ID(Prims::q2, i, j, k)];
        args.q3_rf = prims[ID(Prims::q3, i, j, k)];
        args.Pi_rf = prims[ID(Prims::Pi, i, j, k)];
        args.pi11_rf = prims[ID(Prims::pi11, i, j, k)];
        args.pi12_rf = prims[ID(Prims::pi12, i, j, k)];
        args.pi13_rf = prims[ID(Prims::pi13, i, j, k)];
        args.pi22_rf = prims[ID(Prims::pi22, i, j, k)];
        args.pi23_rf = prims[ID(Prims::pi23, i, j, k)];
        args.pi33_rf = prims[ID(Prims::pi33, i, j, k)];
        args.gamma = d->gamma;
        
        if (alternative_C2P) {
          
          sol[0] = 1/(aux[ID(Aux::W, i, j, k)]*(1 + (prims[ID(Prims::p, i, j, k)]*(d->gamma/(d->gamma-1)) + prims[ID(Prims::Pi, i, j, k)])/prims[ID(Prims::n, i, j, k)]));
          sol[1] = (prims[ID(Prims::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
          sol[2] = (prims[ID(Prims::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
          sol[3] = (prims[ID(Prims::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];
        
          // Solve residual = 0
          info = __cminpack_func__(hybrd1) (&NSAlternativeResidual, &args, sys_size, sol, res,
                                            tol, wa, lwa);        
        
        } else {
        
          sol[0] = prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)] - 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] - aux[ID(Aux::pi00, i, j, k)];
          sol[1] = (prims[ID(Prims::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
          sol[2] = (prims[ID(Prims::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
          sol[3] = (prims[ID(Prims::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];
  
          // Solve residual = 0
          info = __cminpack_func__(hybrd1) (&NSresidual, &args, sys_size, sol, res,
                                            tol, wa, lwa);
  //        info = __cminpack_func__(hybrd) (&NSresidual, &args, sys_size, sol, res,
  //                                          tol, maxfev, ml, mu, epsfcn, &diag[0], mode, factor, nprint, &nfev, &fjac[0][0], ldfjac, &r[0], lr, &qtf[0], &wa1[0], &wa2[0], &wa3[0], &wa4[0]);        
        }
                                                                   
        // If root find fails, add failed cell to the list
        if (info!=1) {
          printf("%i info\n",info);
          printf("(%i, %i, %i) failed\n", i, j, k);
          printf("(%g, %g, %g, %g) res\n", res[0], res[1], res[2], res[3]);
          printf("(%g, %g, %g, %g) sol\n", sol[0], sol[1], sol[2], sol[3]);
          std::cout << "Prims ";
          for (int vz(0); vz < d->Nprims; vz++) {
            std::cout << d->prims[ID(vz, i, j, k)] << " ";
          }
          std::cout << std::endl;
          std::cout << "Cons  ";
          for (int vz(0); vz < d->Ncons; vz++) {
            std::cout << d->cons[ID(vz, i, j, k)] << " ";
          }
          std::cout << std::endl;
          std::cout << "Aux   ";
          for (int vz(0); vz < d->Naux; vz++) {
            std::cout << d->aux[ID(vz, i, j, k)] << " ";
          }
          std::cout << std::endl;
          //printf("(%f, %f, %f, %f, %f) prims\n",  prims[ID(Prims::p, i, j, k)], prims[ID(Prims::Pi, i, j, k)], prims[ID(Prims::n, i, j, k)], prims[ID(Prims::v1, i, j, k)], prims[ID(Prims::q1, i, j, k)]);
          //printf("(%f, %f, %f, %f) aux  \n",  aux[ID(Aux::W, i, j, k)], aux[ID(Aux::qv, i, j, k)], aux[ID(Aux::pi00, i, j, k)], aux[ID(Aux::pi01, i, j, k)]);
          //printf("(%f, %f, %f, %f) cons \n",  cons[ID(Cons::Y1, i, j, k)], cons[ID(Cons::D, i, j, k)], cons[ID(Cons::S1, i, j, k)], cons[ID(Cons::Tau, i, j, k)]);
          exit(0);
          Failed fail = {i, j, k};
          fails.push_back(fail);
        }
        else {
          // Now have the correct values for Chi, Sigma123
          solution[ID(0, i, j, k)] = sol[0];
          solution[ID(1, i, j, k)] = sol[1];
          solution[ID(2, i, j, k)] = sol[2];
          solution[ID(3, i, j, k)] = sol[3];
        }      
      
      } // End k-loop
    } // End j-loop
  } // End i-loop
  
  /*
  
  // ################################## Smart guessing ########################### //
  // Are there any failures?
  if (fails.size() > 0) {
    int x, y, z;
    // Loop through any failed cells and try again, using the mean of successfull
    // surrounding cells solutions as an initial estimate
    for (Failed fail : fails) {
      x = fail.x;
      y = fail.y;
      z = fail.z;
      // Vector to contain successful neighbours
      std::vector<Failed> neighbours;
      if (x > 0) neighbours.push_back(Failed {x-1, y, z});
      if (y > 0) neighbours.push_back(Failed {x, y-1, z});
      if (z > 0) neighbours.push_back(Failed {x, y, z-1});
      if (x < d->Nx - 1) neighbours.push_back(Failed {x+1, y, z});
      if (y < d->Ny - 1) neighbours.push_back(Failed {x, y+1, z});
      if (z < d->Nz - 1) neighbours.push_back(Failed {x, y, z+1});

      sol[0] = 0;
      sol[1] = 0;
      sol[2] = 0;
      sol[3] = 0;
      for (Failed neighbour : neighbours) {
        sol[0] += solution[ID(0, neighbour.x, neighbour.y, neighbour.z)];
        sol[1] += solution[ID(1, neighbour.x, neighbour.y, neighbour.z)];
        sol[2] += solution[ID(2, neighbour.x, neighbour.y, neighbour.z)];
        sol[3] += solution[ID(3, neighbour.x, neighbour.y, neighbour.z)];
      }
      sol[0] /= neighbours.size();
      sol[1] /= neighbours.size();
      sol[2] /= neighbours.size();
      sol[3] /= neighbours.size();
      // Solve residual = 0
      info = __cminpack_func__(hybrd1) (&NSresidual, &args, n, sol, res,
                                        tol, wa, lwa);
      if (info != 1) {
        printf("Smart guessing did not work, exiting\n");
        printf("(%i) \n",info);
        printf("(%d, %d, %d) failed\n", fail.x, fail.y, fail.z);
        std::exit(1);
      } else {
        smartGuesses++;
        printf("Smart guessing worked!\n");
        solution[ID(0, x, y, z)] = sol[0];
        solution[ID(1, x, y, z)] = sol[1];
        solution[ID(2, x, y, z)] = sol[2];
        solution[ID(3, x, y, z)] = sol[3];
      } // else
    } // for

  } // if

  */

  if (alternative_C2P) {
  
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
   
          aux[ID(Aux::vsqrd, i, j, k)] = solution[ID(0, i, j, k)]*solution[ID(0, i, j, k)]*((cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])*(cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)]) 
                            + (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])*(cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)]) 
                            + (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])*(cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)]))/(cons[ID(Cons::D, i, j, k)]*cons[ID(Cons::D, i, j, k)]);
          aux[ID(Aux::W, i, j, k)] = (1 / sqrt(1 - aux[ID(Aux::vsqrd, i, j, k)]));
          prims[Prims::n] = cons[ID(Cons::D, i, j, k)] / aux[ID(Aux::W, i, j, k)];
          prims[ID(Prims::v1, i, j, k)] = solution[ID(0, i, j, k)]*(cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])/cons[ID(Cons::D, i, j, k)];
          prims[ID(Prims::v2, i, j, k)] = solution[ID(0, i, j, k)]*(cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])/cons[ID(Cons::D, i, j, k)];
          prims[ID(Prims::v3, i, j, k)] = solution[ID(0, i, j, k)]*(cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])/cons[ID(Cons::D, i, j, k)];
          prims[ID(Prims::p, i, j, k)] = cons[ID(Cons::D, i, j, k)]*(1/solution[ID(0, i, j, k)] -1) - prims[ID(Prims::Pi, i, j, k)] 
                                         + 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi00, i, j, k)] - cons[ID(Cons::Tau, i, j, k)];
          prims[ID(Prims::rho, i, j, k)] = prims[ID(Prims::n, i, j, k)] + prims[ID(Prims::p, i, j, k)]/(d->gamma-1);
          
          // Again, repeating this here once the correct values for v1,v2,v3 have been set...
          aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                                 + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
          aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)]; // Should be unncessary here
          aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                   + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                   + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                   + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                   + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                   + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                   + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          
          aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
          aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
  
        } // End k-loop
      } // End j-loop
    } // End i-loop    
  
  } else {
    
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {

          // C2P Scheme as outlined in HP/FYR
          aux[ID(Aux::vsqrd, i, j, k)] = ((cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])*(cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)]) 
                                    + (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])*(cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])
                                    + (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])*(cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)]))
                                    /((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])*(cons[ID(Cons::Tau, i, j, k)] 
                                    + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)]));
          aux[ID(Aux::W, i, j, k)] = 1 / sqrt((1-aux[ID(Aux::vsqrd, i, j, k)]));
          prims[ID(Prims::n, i, j, k)] = cons[ID(Cons::D, i, j, k)] / aux[ID(Aux::W, i, j, k)];
          aux[ID(Aux::rho_plus_p, i, j, k)] = ((cons[ID(Cons::Tau, i, j, k)] + cons[ID(Cons::D, i, j, k)] + solution[ID(0, i, j, k)])
                                              /(aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)])) - prims[ID(Prims::Pi, i, j, k)];
          prims[ID(Prims::v1, i, j, k)] = (cons[ID(Cons::S1, i, j, k)] - solution[ID(1, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                   + prims[ID(Prims::Pi, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);
          prims[ID(Prims::v2, i, j, k)] = (cons[ID(Cons::S2, i, j, k)] - solution[ID(2, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                   + prims[ID(Prims::Pi, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);  
          prims[ID(Prims::v3, i, j, k)] = (cons[ID(Cons::S3, i, j, k)] - solution[ID(3, i, j, k)])/((aux[ID(Aux::rho_plus_p, i, j, k)] 
                                   + prims[ID(Prims::Pi, i, j, k)])*aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)]);  
          prims[ID(Prims::p, i, j, k)] = (aux[ID(Aux::rho_plus_p, i, j, k)] - prims[ID(Prims::n, i, j, k)])*((d->gamma-1)/d->gamma);
          prims[ID(Prims::rho, i, j, k)] = aux[ID(Aux::rho_plus_p, i, j, k)] - prims[ID(Prims::p, i, j, k)];
  
          // Again, repeating this here once the correct values for v1,v2,v3 have been set...
          aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                                 + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
          aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)]; // Should be unncessary here
          aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                   + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                   + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                   + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                   + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                   + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                   + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)]; // dbl check sign on orthogonality relation
          
          aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
          aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
  
        } // End k-loop
      } // End j-loop
    } // End i-loop  
  } // else

  double dtT;
  double dxT;
  double dyT;
  double dzT;
  
  double dtW;
  double dxW;
  double dyW;
  double dzW;
  double dtux;
  double dtuy;
  double dtuz;
  double dxux;
  double dxuy;
  double dxuz;
  double dyux;
  double dyuy;
  double dyuz;
  double dzux;
  double dzuy;
  double dzuz;
  
  // double a,b,c,d,e,f,g,h,l,m;

  double kappa = this->data->optionalSimArgs[0];
  double zeta = this->data->optionalSimArgs[2];
  double eta = this->data->optionalSimArgs[4];

  double** sigmatensor = new double*[4];
  for (int i = 0; i < 4; i++) {
      sigmatensor[i] = new double[4];
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        dtT = aux[ID(Aux::dTdt, i, j, k)];
        dxT = (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx);
        dyT = (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy);
        dzT = (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz);

        dtW = aux[ID(Aux::dWdt, i, j, k)]; 
        dxW = (aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(2*d->dx);
        dyW = (aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(2*d->dy);
        dzW = (aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(2*d->dz);

       
        // Chain rule for time derivs
        dtux = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(Aux::dWdt, i, j, k)];
        dtuy = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(Aux::dWdt, i, j, k)];
        dtuz = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(Aux::dWdt, i, j, k)];
        
        dxux = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx);
        dyuy = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        dzuz = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);

        dxuy = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx);
        dxuz = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx);
        dyux = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        dyuz = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy);
        dzux = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);
        dzuy = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz);
        
        /*
        dxT = minmodGradFO(aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i+1, j, k)], d->dx);
        dyT = minmodGradFO(aux[ID(Aux::T, i, j-1, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j+1, k)], d->dy);
        dzT = minmodGradFO(aux[ID(Aux::T, i, j, k-1)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j, k+1)], d->dz);
        
        dxux = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)], d->dx);
        dyuy = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v2, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)], d->dy);
        dzuz = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);                          

        dxuy = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)], d->dx);
        dxuz = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)], d->dx);
       
        dyux = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)], d->dy);
        dyuz = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)], d->dy);

        dzux = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)], d->dz);  
        dzuy = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);  
        */

        aux[ID(Aux::omega00, i, j, k)] = 0.0;
        aux[ID(Aux::omega01, i, j, k)] = dtux - dxW;
        aux[ID(Aux::omega02, i, j, k)] = dtuy - dyW;
        aux[ID(Aux::omega03, i, j, k)] = dtuz - dzW;
        aux[ID(Aux::omega11, i, j, k)] = 0.0;
        aux[ID(Aux::omega12, i, j, k)] = dxuy - dyux;
        aux[ID(Aux::omega13, i, j, k)] = dxuz - dzux;
        aux[ID(Aux::omega22, i, j, k)] = 0.0;
        aux[ID(Aux::omega23, i, j, k)] = dyuz - dzuy;
        aux[ID(Aux::omega33, i, j, k)] = 0.0;
        aux[ID(Aux::omegasqrd, i, j, k)] = aux[ID(Aux::omega00, i, j, k)] * aux[ID(Aux::omega00, i, j, k)]
                                         + 2*aux[ID(Aux::omega01, i, j, k)] * aux[ID(Aux::omega01, i, j, k)]  
                                         + 2*aux[ID(Aux::omega02, i, j, k)] * aux[ID(Aux::omega02, i, j, k)]  
                                         + 2*aux[ID(Aux::omega03, i, j, k)] * aux[ID(Aux::omega03, i, j, k)]  
                                         + aux[ID(Aux::omega11, i, j, k)] * aux[ID(Aux::omega11, i, j, k)]  
                                         + 2*aux[ID(Aux::omega12, i, j, k)] * aux[ID(Aux::omega12, i, j, k)]
                                         + 2*aux[ID(Aux::omega13, i, j, k)] * aux[ID(Aux::omega13, i, j, k)]
                                         + aux[ID(Aux::omega22, i, j, k)] * aux[ID(Aux::omega22, i, j, k)]
                                         + 2*aux[ID(Aux::omega23, i, j, k)] * aux[ID(Aux::omega23, i, j, k)]
                                         + aux[ID(Aux::omega33, i, j, k)] * aux[ID(Aux::omega33, i, j, k)];

        // Shear
        aux[ID(Aux::sigma11, i, j, k)] = 2*dxux 
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma12, i, j, k)] = dxuy + dyux
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma13, i, j, k)] = dxuz + dzux
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma22, i, j, k)] = 2*dyuy
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma23, i, j, k)] = dyuz + dzuy
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma33, i, j, k)] = 2*dzuz
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        // Use orthogonality for sigma00, sigma01 etc.
        aux[ID(Aux::sigma00, i, j, k)] = aux[ID(Aux::sigma11, i, j, k)] + aux[ID(Aux::sigma22, i, j, k)] + aux[ID(Aux::sigma33, i, j, k)];
        aux[ID(Aux::sigma01, i, j, k)] = aux[ID(Aux::sigma11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                + aux[ID(Aux::sigma12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                + aux[ID(Aux::sigma13, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::sigma02, i, j, k)] = aux[ID(Aux::sigma12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                + aux[ID(Aux::sigma22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                + aux[ID(Aux::sigma23, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::sigma03, i, j, k)] = aux[ID(Aux::sigma13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                + aux[ID(Aux::sigma23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                + aux[ID(Aux::sigma33, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::sigmasqrd, i, j, k)] = aux[ID(Aux::sigma00, i, j, k)] * aux[ID(Aux::sigma00, i, j, k)] 
                                         + 2*aux[ID(Aux::sigma01, i, j, k)] * aux[ID(Aux::sigma01, i, j, k)]
                                         + 2*aux[ID(Aux::sigma02, i, j, k)] * aux[ID(Aux::sigma02, i, j, k)]
                                         + 2*aux[ID(Aux::sigma03, i, j, k)] * aux[ID(Aux::sigma03, i, j, k)]
                                         + aux[ID(Aux::sigma11, i, j, k)] * aux[ID(Aux::sigma11, i, j, k)]
                                         + 2*aux[ID(Aux::sigma12, i, j, k)] * aux[ID(Aux::sigma12, i, j, k)]
                                         + 2*aux[ID(Aux::sigma13, i, j, k)] * aux[ID(Aux::sigma13, i, j, k)]
                                         + aux[ID(Aux::sigma22, i, j, k)] * aux[ID(Aux::sigma22, i, j, k)]
                                         + 2*aux[ID(Aux::sigma23, i, j, k)] * aux[ID(Aux::sigma23, i, j, k)]
                                         + aux[ID(Aux::sigma33, i, j, k)] * aux[ID(Aux::sigma33, i, j, k)];

        /*
        std::vector<std::vector<double>> sigmatensor{ {aux[ID(Aux::sigma00, i, j, k)], aux[ID(Aux::sigma01, i, j, k)], aux[ID(Aux::sigma02, i, j, k)], aux[ID(Aux::sigma03, i, j, k)]},
                                                      {aux[ID(Aux::sigma01, i, j, k)], aux[ID(Aux::sigma11, i, j, k)], aux[ID(Aux::sigma12, i, j, k)], aux[ID(Aux::sigma13, i, j, k)]},
                                                      {aux[ID(Aux::sigma02, i, j, k)], aux[ID(Aux::sigma12, i, j, k)], aux[ID(Aux::sigma22, i, j, k)], aux[ID(Aux::sigma23, i, j, k)]},
                                                      {aux[ID(Aux::sigma03, i, j, k)], aux[ID(Aux::sigma13, i, j, k)], aux[ID(Aux::sigma23, i, j, k)], aux[ID(Aux::sigma33, i, j, k)]} };

        
        sigmatensor[0][0] = aux[ID(Aux::sigma00, i, j, k)];
        sigmatensor[0][1] = aux[ID(Aux::sigma01, i, j, k)];
        sigmatensor[0][2] = aux[ID(Aux::sigma02, i, j, k)];
        sigmatensor[0][3] = aux[ID(Aux::sigma03, i, j, k)];
        sigmatensor[1][0] = aux[ID(Aux::sigma01, i, j, k)];
        sigmatensor[1][1] = aux[ID(Aux::sigma11, i, j, k)];
        sigmatensor[1][2] = aux[ID(Aux::sigma12, i, j, k)];
        sigmatensor[1][3] = aux[ID(Aux::sigma13, i, j, k)];
        sigmatensor[2][0] = aux[ID(Aux::sigma02, i, j, k)];
        sigmatensor[2][1] = aux[ID(Aux::sigma12, i, j, k)];
        sigmatensor[2][2] = aux[ID(Aux::sigma22, i, j, k)];
        sigmatensor[2][3] = aux[ID(Aux::sigma23, i, j, k)];
        sigmatensor[3][0] = aux[ID(Aux::sigma03, i, j, k)];
        sigmatensor[3][1] = aux[ID(Aux::sigma13, i, j, k)];
        sigmatensor[3][2] = aux[ID(Aux::sigma23, i, j, k)];
        sigmatensor[3][3] = aux[ID(Aux::sigma33, i, j, k)];
        */
        
        // aux[ID(Aux::detsigma, i, j, k)] = calculateDeterminant(sigmatensor, 4);
        // aux[ID(Aux::detsigma, i, j, k)] = calculate4Determinant(sigmatensor);
        // aux[ID(Aux::detsigma, i, j, k)] = CalcDeterminant(sigmatensor);

        // Manual calculation of sigma 3-determinant
        aux[ID(Aux::detsigma, i, j, k)] = aux[ID(Aux::sigma00, i, j, k)]*(aux[ID(Aux::sigma11, i, j, k)]*aux[ID(Aux::sigma22, i, j, k)] - aux[ID(Aux::sigma12, i, j, k)]*aux[ID(Aux::sigma12, i, j, k)]) 
                                          - aux[ID(Aux::sigma01, i, j, k)]*(aux[ID(Aux::sigma01, i, j, k)]*aux[ID(Aux::sigma22, i, j, k)] - aux[ID(Aux::sigma12, i, j, k)]*aux[ID(Aux::sigma02, i, j, k)]) 
                                          + aux[ID(Aux::sigma02, i, j, k)]*(aux[ID(Aux::sigma01, i, j, k)]*aux[ID(Aux::sigma12, i, j, k)] - aux[ID(Aux::sigma11, i, j, k)]*aux[ID(Aux::sigma02, i, j, k)]);
        
        
        /*
        double a = aux[ID(Aux::sigma00, i, j, k)];
        double b = aux[ID(Aux::sigma01, i, j, k)];
        double c = aux[ID(Aux::sigma02, i, j, k)];
        double n = aux[ID(Aux::sigma03, i, j, k)];
        double e = aux[ID(Aux::sigma11, i, j, k)];
        double f = aux[ID(Aux::sigma12, i, j, k)];
        double g = aux[ID(Aux::sigma13, i, j, k)];
        double h = aux[ID(Aux::sigma22, i, j, k)];
        double l = aux[ID(Aux::sigma23, i, j, k)];
        double m = aux[ID(Aux::sigma33, i, j, k)];
      
        aux[ID(Aux::detsigma, i, j, k)] = a*e*h*m - (a*e*l*l + a*h*g*g + a*m*f*f + e*h*n*n + e*m*c*c + h*m*b*b)
        + 2*(a*f*g*l * e*c*n*l + h*b*n*g + m*b*c*f) - 2*(b*f*l*n + b*c*g*l + c*g*f*n)
        + (b*l)*(b*l) + (c*g)*(c*g) + (n*f)*(n*f);
        
        if (a != 0.0) {
        printf("(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) a,b,c,n,e,f,g,h,l,m\n", a,b,c,n,e,f,g,h,l,m);
        printf("(%g) detsigma\n", aux[ID(Aux::detsigma, i, j, k)]);
        exit(0);
        }
        */

        aux[ID(Aux::a0, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtW + prims[ID(Prims::v1, i, j, k)]*dxW
          + prims[ID(Prims::v2, i, j, k)]*dyW + prims[ID(Prims::v3, i, j, k)]*dzW );

        aux[ID(Aux::a1, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtux + prims[ID(Prims::v1, i, j, k)]*dxux
          + prims[ID(Prims::v2, i, j, k)]*dyux + prims[ID(Prims::v3, i, j, k)]*dzux );
        
        aux[ID(Aux::a2, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtuy + prims[ID(Prims::v1, i, j, k)]*dxuy
          + prims[ID(Prims::v2, i, j, k)]*dyuy + prims[ID(Prims::v3, i, j, k)]*dzuy );
        
        aux[ID(Aux::a3, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtuz + prims[ID(Prims::v1, i, j, k)]*dxuz
          + prims[ID(Prims::v2, i, j, k)]*dyuz + prims[ID(Prims::v3, i, j, k)]*dzuz );
        
        aux[ID(Aux::Theta0, i, j, k)] = ( (-1+ sqr(aux[ID(Aux::W, i, j, k)]))*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)]) );
        aux[ID(Aux::Theta1, i, j, k)] = ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)]) + 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)]) );
        aux[ID(Aux::Theta2, i, j, k)] = ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)]) + 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)]) );
        aux[ID(Aux::Theta3, i, j, k)] = ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)]) + 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)]) );


        calculateDissipativeCoefficientsSingleCell(cons, prims, aux, i, j, k);

        // q = -kappa*Theta
        prims[ID(Prims::q1, i, j, k)] = -kappa*aux[ID(Aux::kappa, i, j, k)]*aux[ID(Aux::Theta1, i, j, k)];
        prims[ID(Prims::q2, i, j, k)] = -kappa*aux[ID(Aux::kappa, i, j, k)]*aux[ID(Aux::Theta2, i, j, k)];
        prims[ID(Prims::q3, i, j, k)] = -kappa*aux[ID(Aux::kappa, i, j, k)]*aux[ID(Aux::Theta3, i, j, k)];
        // qv contraction
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);

        // theta 20 then Pi,NS 13 
        aux[ID(Aux::theta, i, j, k)] = dtW + dxux + dyuy + dzuz; // minus sign on dtW here??
        // Pi,NS = -zeta*theta
        prims[ID(Prims::Pi, i, j, k)] = -zeta*aux[ID(Aux::zeta, i, j, k)] * aux[ID(Aux::theta, i, j, k)];
  
        // pi^l_j,NS 14
        prims[ID(Prims::pi11, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma11, i, j, k)];
        prims[ID(Prims::pi12, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma12, i, j, k)];
        prims[ID(Prims::pi13, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma13, i, j, k)];
        prims[ID(Prims::pi22, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma22, i, j, k)];
        prims[ID(Prims::pi23, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma23, i, j, k)];
        prims[ID(Prims::pi33, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma33, i, j, k)];
        aux[ID(Aux::pi00, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma00, i, j, k)];
        aux[ID(Aux::pi01, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma01, i, j, k)];
        aux[ID(Aux::pi02, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma02, i, j, k)];
        aux[ID(Aux::pi03, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma03, i, j, k)];

      }
    }
  }

//  // pi^0_mu, qv
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) { // Minus signs here (?)
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
//        aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)];
//
//        aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
//                                + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
//                                + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
//        aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
//                                + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
//                                + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
//        aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
//                                + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
//                                + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
      }
    }
  }

//  int i = 100;
//  int j = 100;
//  int k = 0;
//
//  std::cout << "Prims ";
//  for (int vz(0); vz < d->Nprims; vz++) {
//    std::cout << d->prims[ID(vz, i, j, k)] << " ";
//  }
//  std::cout << std::endl;
//  std::cout << "Cons  ";
//  for (int vz(0); vz < d->Ncons; vz++) {
//    std::cout << d->cons[ID(vz, i, j, k)] << " ";
//  }
//  std::cout << std::endl;
//  std::cout << "Aux   ";
//  for (int vz(0); vz < d->Naux; vz++) {
//    std::cout << d->aux[ID(vz, i, j, k)] << " ";
//  }
//  std::cout << std::endl;


}

void NS::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  printf("Calling primsToAll\n");

  // W, q_kv^k, pi^0_0
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(Aux::vsqrd, i, j, k)] = prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                  + prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                  + prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::W, i, j, k)] = 1 / sqrt( 1 - aux[ID(Aux::vsqrd, i, j, k)] );
        prev_vars[ID(0, i, j, k)] = aux[ID(Aux::W, i, j, k)]; // Set here for time-differencing
        // set the others too??
        prims[ID(Prims::rho, i, j, k)] = prims[ID(Prims::p, i, j, k)]/(d->gamma -1) + prims[ID(Prims::n, i, j, k)];
        
        // Alternative expression I'd rather not use...
        //prims[ID(Prims::rho, i, j, k)] =  prims[ID(Prims::n, i, j, k)]*(1+aux[ID(Aux::e, i, j, k)]);
        aux[ID(Aux::T, i, j, k)] = prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
        aux[ID(Aux::e, i, j, k)] = prims[ID(Prims::p, i, j, k)] / (prims[ID(Prims::n, i, j, k)]*(d->gamma-1));
        aux[ID(Aux::h, i, j, k)] = 1 + aux[ID(Aux::e, i, j, k)] + prims[ID(Prims::p, i, j, k)] / prims[ID(Prims::n, i, j, k)];
       
      }
    }
  }

  double kappa = this->data->optionalSimArgs[0];
//  double tau_q = this->data->optionalSimArgs[1];
  double zeta = this->data->optionalSimArgs[2];
//  double tau_Pi = this->data->optionalSimArgs[3];
  double eta = this->data->optionalSimArgs[4];
//  double tau_pi = this->data->optionalSimArgs[5];

  double dtT;
  double dxT;
  double dyT;
  double dzT;
  
  double dtW;
  double dxW;
  double dyW;
  double dzW;
  double dtux;
  double dtuy;
  double dtuz;
  double dxux;
  double dxuy;
  double dxuz;
  double dyux;
  double dyuy;
  double dyuz;
  double dzux;
  double dzuy;
  double dzuz;
  
  double** sigmatensor = new double*[4];
  for (int i = 0; i < 4; i++) {
      sigmatensor[i] = new double[4];
  }


  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        dtT = aux[ID(Aux::dTdt, i, j, k)];
        dxT = (aux[ID(Aux::T, i+1, j, k)] - aux[ID(Aux::T, i-1, j, k)])/(2*d->dx);
        dyT = (aux[ID(Aux::T, i, j+1, k)] - aux[ID(Aux::T, i, j-1, k)])/(2*d->dy);
        dzT = (aux[ID(Aux::T, i, j, k+1)] - aux[ID(Aux::T, i, j, k-1)])/(2*d->dz);

        dtW = aux[ID(Aux::dWdt, i, j, k)];
        dxW = (aux[ID(Aux::W, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)])/(2*d->dx);
        dyW = (aux[ID(Aux::W, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)])/(2*d->dy);
        dzW = (aux[ID(Aux::W, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)])/(2*d->dz);

        // Chain rule for time derivs
        dtux = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(Aux::dWdt, i, j, k)];
        dtuy = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(Aux::dWdt, i, j, k)];
        dtuz = aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::dv1dt, i, j, k)] + prims[ID(Prims::v1, i, j, k)]*aux[ID(Aux::dWdt, i, j, k)];
        
        dxux = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)])/(2*d->dx);
        dyuy = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        dzuz = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);

        dxuy = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)])/(2*d->dx);
        dxuz = (aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)] - aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)])/(2*d->dx);
        dyux = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)])/(2*d->dy);
        dyuz = (aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)] - aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)])/(2*d->dy);
        dzux = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)])/(2*d->dz);
        dzuy = (aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)] - aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)])/(2*d->dz);
        
        /*
        dxT = minmodGradFO(aux[ID(Aux::T, i-1, j, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i+1, j, k)], d->dx);
        dyT = minmodGradFO(aux[ID(Aux::T, i, j-1, k)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j+1, k)], d->dy);
        dzT = minmodGradFO(aux[ID(Aux::T, i, j, k-1)], aux[ID(Aux::T, i, j, k)], aux[ID(Aux::T, i, j, k+1)], d->dz);
        
        dxux = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v1, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v1, i+1, j, k)], d->dx);
        dyuy = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v2, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v2, i, j+1, k)], d->dy);
        dzuz = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);                          

        dxuy = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v2, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v2, i+1, j, k)], d->dx);
        dxuz = minmodGradFO(aux[ID(Aux::W, i-1, j, k)]*prims[ID(Prims::v3, i-1, j, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i+1, j, k)]*prims[ID(Prims::v3, i+1, j, k)], d->dx);
       
        dyux = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v1, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v1, i, j+1, k)], d->dy);
        dyuz = minmodGradFO(aux[ID(Aux::W, i, j-1, k)]*prims[ID(Prims::v3, i, j-1, k)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)],
                            aux[ID(Aux::W, i, j+1, k)]*prims[ID(Prims::v3, i, j+1, k)], d->dy);

        dzux = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v1, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v1, i, j, k+1)], d->dz);  
        dzuy = minmodGradFO(aux[ID(Aux::W, i, j, k-1)]*prims[ID(Prims::v2, i, j, k-1)],  aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)],
                            aux[ID(Aux::W, i, j, k+1)]*prims[ID(Prims::v2, i, j, k+1)], d->dz);  
        */

        aux[ID(Aux::omega00, i, j, k)] = 0.0;
        aux[ID(Aux::omega01, i, j, k)] = dtux - dxW;
        aux[ID(Aux::omega02, i, j, k)] = dtuy - dyW;
        aux[ID(Aux::omega03, i, j, k)] = dtuz - dzW;
        aux[ID(Aux::omega11, i, j, k)] = 0.0;
        aux[ID(Aux::omega12, i, j, k)] = dxuy - dyux;
        aux[ID(Aux::omega13, i, j, k)] = dxuz - dzux;
        aux[ID(Aux::omega22, i, j, k)] = 0.0;
        aux[ID(Aux::omega23, i, j, k)] = dyuz - dzuy;
        aux[ID(Aux::omega33, i, j, k)] = 0.0;
        aux[ID(Aux::omegasqrd, i, j, k)] = aux[ID(Aux::omega00, i, j, k)] * aux[ID(Aux::omega00, i, j, k)]
                                         + 2*aux[ID(Aux::omega01, i, j, k)] * aux[ID(Aux::omega01, i, j, k)]  
                                         + 2*aux[ID(Aux::omega02, i, j, k)] * aux[ID(Aux::omega02, i, j, k)]  
                                         + 2*aux[ID(Aux::omega03, i, j, k)] * aux[ID(Aux::omega03, i, j, k)]  
                                         + aux[ID(Aux::omega11, i, j, k)] * aux[ID(Aux::omega11, i, j, k)]  
                                         + 2*aux[ID(Aux::omega12, i, j, k)] * aux[ID(Aux::omega12, i, j, k)]
                                         + 2*aux[ID(Aux::omega13, i, j, k)] * aux[ID(Aux::omega13, i, j, k)]
                                         + aux[ID(Aux::omega22, i, j, k)] * aux[ID(Aux::omega22, i, j, k)]
                                         + 2*aux[ID(Aux::omega23, i, j, k)] * aux[ID(Aux::omega23, i, j, k)]
                                         + aux[ID(Aux::omega33, i, j, k)] * aux[ID(Aux::omega33, i, j, k)];

        // Shear
        aux[ID(Aux::sigma11, i, j, k)] = 2*dxux 
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma12, i, j, k)] = dxuy + dyux
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma13, i, j, k)] = dxuz + dzux
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma22, i, j, k)] = 2*dyuy
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma23, i, j, k)] = dyuz + dzuy
          - (2/3)*((aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        aux[ID(Aux::sigma33, i, j, k)] = 2*dzuz
          - (2/3)*(1 + (aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)])*(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*aux[ID(Aux::theta, i, j, k)];
        // Use orthogonality for sigma00, sigma01 etc.
        aux[ID(Aux::sigma00, i, j, k)] = aux[ID(Aux::sigma11, i, j, k)] + aux[ID(Aux::sigma22, i, j, k)] + aux[ID(Aux::sigma33, i, j, k)];
        aux[ID(Aux::sigma01, i, j, k)] = aux[ID(Aux::sigma11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                + aux[ID(Aux::sigma12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                + aux[ID(Aux::sigma13, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::sigma02, i, j, k)] = aux[ID(Aux::sigma12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                + aux[ID(Aux::sigma22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                + aux[ID(Aux::sigma23, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::sigma03, i, j, k)] = aux[ID(Aux::sigma13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
                                + aux[ID(Aux::sigma23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
                                + aux[ID(Aux::sigma33, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
        aux[ID(Aux::sigmasqrd, i, j, k)] = aux[ID(Aux::sigma00, i, j, k)] * aux[ID(Aux::sigma00, i, j, k)] 
                                         + 2*aux[ID(Aux::sigma01, i, j, k)] * aux[ID(Aux::sigma01, i, j, k)]
                                         + 2*aux[ID(Aux::sigma02, i, j, k)] * aux[ID(Aux::sigma02, i, j, k)]
                                         + 2*aux[ID(Aux::sigma03, i, j, k)] * aux[ID(Aux::sigma03, i, j, k)]
                                         + aux[ID(Aux::sigma11, i, j, k)] * aux[ID(Aux::sigma11, i, j, k)]
                                         + 2*aux[ID(Aux::sigma12, i, j, k)] * aux[ID(Aux::sigma12, i, j, k)]
                                         + 2*aux[ID(Aux::sigma13, i, j, k)] * aux[ID(Aux::sigma13, i, j, k)]
                                         + aux[ID(Aux::sigma22, i, j, k)] * aux[ID(Aux::sigma22, i, j, k)]
                                         + 2*aux[ID(Aux::sigma23, i, j, k)] * aux[ID(Aux::sigma23, i, j, k)]
                                         + aux[ID(Aux::sigma33, i, j, k)] * aux[ID(Aux::sigma33, i, j, k)];
         /*
        std::vector<std::vector<double>> sigmatensor{ {aux[ID(Aux::sigma00, i, j, k)], aux[ID(Aux::sigma01, i, j, k)], aux[ID(Aux::sigma02, i, j, k)], aux[ID(Aux::sigma03, i, j, k)]},
                                                      {aux[ID(Aux::sigma01, i, j, k)], aux[ID(Aux::sigma11, i, j, k)], aux[ID(Aux::sigma12, i, j, k)], aux[ID(Aux::sigma13, i, j, k)]},
                                                      {aux[ID(Aux::sigma02, i, j, k)], aux[ID(Aux::sigma12, i, j, k)], aux[ID(Aux::sigma22, i, j, k)], aux[ID(Aux::sigma23, i, j, k)]},
                                                      {aux[ID(Aux::sigma03, i, j, k)], aux[ID(Aux::sigma13, i, j, k)], aux[ID(Aux::sigma23, i, j, k)], aux[ID(Aux::sigma33, i, j, k)]} };

       
        sigmatensor[0][0] = aux[ID(Aux::sigma00, i, j, k)];
        sigmatensor[0][1] = aux[ID(Aux::sigma01, i, j, k)];
        sigmatensor[0][2] = aux[ID(Aux::sigma02, i, j, k)];
        sigmatensor[0][3] = aux[ID(Aux::sigma03, i, j, k)];
        sigmatensor[1][0] = aux[ID(Aux::sigma01, i, j, k)];
        sigmatensor[1][1] = aux[ID(Aux::sigma11, i, j, k)];
        sigmatensor[1][2] = aux[ID(Aux::sigma12, i, j, k)];
        sigmatensor[1][3] = aux[ID(Aux::sigma13, i, j, k)];
        sigmatensor[2][0] = aux[ID(Aux::sigma02, i, j, k)];
        sigmatensor[2][1] = aux[ID(Aux::sigma12, i, j, k)];
        sigmatensor[2][2] = aux[ID(Aux::sigma22, i, j, k)];
        sigmatensor[2][3] = aux[ID(Aux::sigma23, i, j, k)];
        sigmatensor[3][0] = aux[ID(Aux::sigma03, i, j, k)];
        sigmatensor[3][1] = aux[ID(Aux::sigma13, i, j, k)];
        sigmatensor[3][2] = aux[ID(Aux::sigma23, i, j, k)];
        sigmatensor[3][3] = aux[ID(Aux::sigma33, i, j, k)];
        */
        
        // aux[ID(Aux::detsigma, i, j, k)] = calculateDeterminant(sigmatensor, 4);
        // aux[ID(Aux::detsigma, i, j, k)] = calculate4Determinant(sigmatensor);
        // aux[ID(Aux::detsigma, i, j, k)] = CalcDeterminant(sigmatensor);

        // Manual calculation of sigma 3-determinant
        aux[ID(Aux::detsigma, i, j, k)] = aux[ID(Aux::sigma00, i, j, k)]*(aux[ID(Aux::sigma11, i, j, k)]*aux[ID(Aux::sigma22, i, j, k)] - aux[ID(Aux::sigma12, i, j, k)]*aux[ID(Aux::sigma12, i, j, k)]) 
                                          - aux[ID(Aux::sigma01, i, j, k)]*(aux[ID(Aux::sigma01, i, j, k)]*aux[ID(Aux::sigma22, i, j, k)] - aux[ID(Aux::sigma12, i, j, k)]*aux[ID(Aux::sigma02, i, j, k)]) 
                                          + aux[ID(Aux::sigma02, i, j, k)]*(aux[ID(Aux::sigma01, i, j, k)]*aux[ID(Aux::sigma12, i, j, k)] - aux[ID(Aux::sigma11, i, j, k)]*aux[ID(Aux::sigma02, i, j, k)]);
        
        
        /*
        double a = aux[ID(Aux::sigma00, i, j, k)];
        double b = aux[ID(Aux::sigma01, i, j, k)];
        double c = aux[ID(Aux::sigma02, i, j, k)];
        double n = aux[ID(Aux::sigma03, i, j, k)];
        double e = aux[ID(Aux::sigma11, i, j, k)];
        double f = aux[ID(Aux::sigma12, i, j, k)];
        double g = aux[ID(Aux::sigma13, i, j, k)];
        double h = aux[ID(Aux::sigma22, i, j, k)];
        double l = aux[ID(Aux::sigma23, i, j, k)];
        double m = aux[ID(Aux::sigma33, i, j, k)];
      
        aux[ID(Aux::detsigma, i, j, k)] = a*e*h*m - (a*e*l*l + a*h*g*g + a*m*f*f + e*h*n*n + e*m*c*c + h*m*b*b)
        + 2*(a*f*g*l * e*c*n*l + h*b*n*g + m*b*c*f) - 2*(b*f*l*n + b*c*g*l + c*g*f*n)
        + (b*l)*(b*l) + (c*g)*(c*g) + (n*f)*(n*f);
        
        if (a != 0.0) {
        printf("(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) a,b,c,n,e,f,g,h,l,m\n", a,b,c,n,e,f,g,h,l,m);
        printf("(%g) detsigma\n", aux[ID(Aux::detsigma, i, j, k)]);
        exit(0);
        }
        */

        aux[ID(Aux::a0, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtW + prims[ID(Prims::v1, i, j, k)]*dxW
          + prims[ID(Prims::v2, i, j, k)]*dyW + prims[ID(Prims::v3, i, j, k)]*dzW );
        aux[ID(Aux::a1, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtux + prims[ID(Prims::v1, i, j, k)]*dxux
          + prims[ID(Prims::v2, i, j, k)]*dyux + prims[ID(Prims::v3, i, j, k)]*dzux );
        aux[ID(Aux::a2, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtuy + prims[ID(Prims::v1, i, j, k)]*dxuy
          + prims[ID(Prims::v2, i, j, k)]*dyuy + prims[ID(Prims::v3, i, j, k)]*dzuy );
        aux[ID(Aux::a3, i, j, k)] = aux[ID(Aux::W, i, j, k)] * ( dtuz + prims[ID(Prims::v1, i, j, k)]*dxuz
          + prims[ID(Prims::v2, i, j, k)]*dyuz + prims[ID(Prims::v3, i, j, k)]*dzuz );
        
        aux[ID(Aux::Theta0, i, j, k)] = ( (-1+ sqr(aux[ID(Aux::W, i, j, k)]))*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)]) );
        aux[ID(Aux::Theta1, i, j, k)] = ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v1, i, j, k)]))*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)]) + 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v1, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)]) );
        aux[ID(Aux::Theta2, i, j, k)] = ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v2, i, j, k)]))*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)]) + 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v2, i, j, k)]*prims[ID(Prims::v3, i, j, k)]*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)]) );
        aux[ID(Aux::Theta3, i, j, k)] = ( (1+ sqr(aux[ID(Aux::W, i, j, k)]*prims[ID(Prims::v3, i, j, k)]))*(dzT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a3, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*(dtT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a0, i, j, k)]) + 
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v1, i, j, k)]*(dxT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a1, i, j, k)])
          + sqr(aux[ID(Aux::W, i, j, k)])*prims[ID(Prims::v3, i, j, k)]*prims[ID(Prims::v2, i, j, k)]*(dyT + aux[ID(Aux::T, i, j, k)]*aux[ID(Aux::a2, i, j, k)]) );


        calculateDissipativeCoefficientsSingleCell(cons, prims, aux, i, j, k);

        // q = -kappa*Theta
        prims[ID(Prims::q1, i, j, k)] = -kappa*aux[ID(Aux::kappa, i, j, k)]*aux[ID(Aux::Theta1, i, j, k)];
        prims[ID(Prims::q2, i, j, k)] = -kappa*aux[ID(Aux::kappa, i, j, k)]*aux[ID(Aux::Theta2, i, j, k)];
        prims[ID(Prims::q3, i, j, k)] = -kappa*aux[ID(Aux::kappa, i, j, k)]*aux[ID(Aux::Theta3, i, j, k)];
        // qv contraction
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);

        // theta 20 then Pi,NS 13 
        aux[ID(Aux::theta, i, j, k)] = dtW + dxux + dyuy + dzuz; // minus sign here??
        // Pi,NS = -zeta*theta
        prims[ID(Prims::Pi, i, j, k)] = -zeta*aux[ID(Aux::zeta, i, j, k)] * aux[ID(Aux::theta, i, j, k)];
  
        // pi^l_j,NS 14
        prims[ID(Prims::pi11, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma11, i, j, k)];
        prims[ID(Prims::pi12, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma12, i, j, k)];
        prims[ID(Prims::pi13, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma13, i, j, k)];
        prims[ID(Prims::pi22, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma22, i, j, k)];
        prims[ID(Prims::pi23, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma23, i, j, k)];
        prims[ID(Prims::pi33, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma33, i, j, k)];
        aux[ID(Aux::pi00, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma00, i, j, k)];
        aux[ID(Aux::pi01, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma01, i, j, k)];
        aux[ID(Aux::pi02, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma02, i, j, k)];
        aux[ID(Aux::pi03, i, j, k)] = -2*eta*aux[ID(Aux::eta, i, j, k)]*aux[ID(Aux::sigma03, i, j, k)];

      }
    }
  }

//  // pi^0_mu, qv
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) { // Minus signs here (?)
        aux[ID(Aux::qv, i, j, k)] = (prims[ID(Prims::q1, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) + (prims[ID(Prims::q2, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) 
                               + (prims[ID(Prims::q3, i, j, k)] * prims[ID(Prims::v3, i, j, k)]);
//        aux[ID(Aux::pi00, i, j, k)] = prims[ID(Prims::pi11, i, j, k)] + prims[ID(Prims::pi22, i, j, k)] + prims[ID(Prims::pi33, i, j, k)];
//
//        aux[ID(Aux::pi01, i, j, k)] = prims[ID(Prims::pi11, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
//                                + prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
//                                + prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
//        aux[ID(Aux::pi02, i, j, k)] = prims[ID(Prims::pi12, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
//                                + prims[ID(Prims::pi22, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
//                                + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
//        aux[ID(Aux::pi03, i, j, k)] = prims[ID(Prims::pi13, i, j, k)]*prims[ID(Prims::v1, i, j, k)] 
//                                + prims[ID(Prims::pi23, i, j, k)]*prims[ID(Prims::v2, i, j, k)] 
//                                + prims[ID(Prims::pi33, i, j, k)]*prims[ID(Prims::v3, i, j, k)];
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // D
        cons[ID(Cons::D, i, j, k)] = prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)];
        // S1,2,3
        cons[ID(Cons::S1, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v1, i, j, k)] 
          + (prims[ID(Prims::q1, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v1, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi01, i, j, k)];
        cons[ID(Cons::S2, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v2, i, j, k)] 
          + (prims[ID(Prims::q2, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v2, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi02, i, j, k)];
        cons[ID(Cons::S3, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] * prims[ID(Prims::v3, i, j, k)] 
          + (prims[ID(Prims::q3, i, j, k)] + aux[ID(Aux::qv, i, j, k)] * prims[ID(Prims::v3, i, j, k)]) * aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi03, i, j, k)];
        // Tau
        cons[ID(Cons::Tau, i, j, k)] = (prims[ID(Prims::rho, i, j, k)] + prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]) * aux[ID(Aux::W, i, j, k)]*aux[ID(Aux::W, i, j, k)] 
        - (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)] + prims[ID(Prims::n, i, j, k)] * aux[ID(Aux::W, i, j, k)]) 
        + 2*aux[ID(Aux::qv, i, j, k)]*aux[ID(Aux::W, i, j, k)] + aux[ID(Aux::pi00, i, j, k)];
      }  
    }
  }

//  int i = 4;
//  int j = 4;
//  int k = 0;
//
//  std::cout << "Prims ";
//  for (int vz(0); vz < d->Nprims; vz++) {
//    std::cout << d->prims[ID(vz, i, j, k)] << " ";
//  }
//  std::cout << std::endl;
//  std::cout << "Cons  ";
//  for (int vz(0); vz < d->Ncons; vz++) {
//    std::cout << d->cons[ID(vz, i, j, k)] << " ";
//  }
//  std::cout << std::endl;
//  std::cout << "Aux   ";
//  for (int vz(0); vz < d->Naux; vz++) {
//    std::cout << d->aux[ID(vz, i, j, k)] << " ";
//  }
//  std::cout << std::endl;
//
//  exit(0);

}

void NS::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Dv
        f[ID(0, i, j, k)] = cons[ID(Cons::D, i, j, k)]*prims[ID(dir, i, j, k)];
        // Sv + ..
        for (int nvar(0); nvar < 3; nvar++) {
          f[ID(1+nvar, i, j, k)] = cons[ID(Cons::S1+nvar, i, j, k)]*prims[ID(dir, i, j, k)] + ( prims[ID(Prims::q1+dir, i, j, k)] * prims[ID(Prims::v1+nvar, i, j, k)]  
            - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+nvar, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)] ) * aux[ID(Aux::W, i, j, k)];
          // (p+Pi)delta_ij
          if (dir == nvar) {
            f[ID(1+nvar, i, j, k)] += (prims[ID(Prims::p, i, j, k)] + prims[ID(Prims::Pi, i, j, k)]);
          }
        }
        //  pi^i_j  
        if (dir == 0) {
          for (int nvar(0); nvar < 3; nvar++) {
            f[ID(1+nvar, i, j, k)] += prims[ID(Prims::pi11+nvar, i, j, k)];
          }
        } else if (dir == 1) {
          f[ID(1, i, j, k)] += prims[ID(Prims::pi12, i, j, k)];
          f[ID(2, i, j, k)] += prims[ID(Prims::pi22, i, j, k)];
          f[ID(3, i, j, k)] += prims[ID(Prims::pi23, i, j, k)];
        } else if (dir == 2) {
          f[ID(1, i, j, k)] += prims[ID(Prims::pi13, i, j, k)];
          f[ID(2, i, j, k)] += prims[ID(Prims::pi23, i, j, k)];
          f[ID(3, i, j, k)] += prims[ID(Prims::pi33, i, j, k)];
        } else {
          throw std::runtime_error("Flux direction is not 0, 1 or 2");
        }

        // (Tau+p)*v + ...
        f[ID(4, i, j, k)] = (cons[ID(Cons::Tau, i, j, k)] + prims[ID(Prims::p, i, j, k)]) * prims[ID(dir, i, j, k)] 
          + (prims[ID(Prims::q1+dir, i, j, k)] - aux[ID(Aux::qv, i, j, k)]*prims[ID(Prims::v1+dir, i, j, k)])*aux[ID(Aux::W, i, j, k)]
          + aux[ID(Aux::pi01+dir, i, j, k)];
      } // End k loop
    } // End j loop
  } // End i loop
}
