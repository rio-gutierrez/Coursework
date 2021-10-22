#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#ifndef M_PI
    #define M_PI 3.14159265358979323846;
#endif

// Global parameters
const double r0 {1.};
const double rn {10.};

// Function prototypes
VectorXd SOR_RES   (const MatrixXd &A, const VectorXd &b,
                    const VectorXd &x0, const int &dim,
                    const int max_it = 100,  const double omega = 1.);
VectorXd SOR       (const MatrixXd &A, const VectorXd &b,
                    const VectorXd &x0, const int &dim,
                    const int max_it = 1000, const double omega = 1.,
                    const double tol = 1.e-10);
VectorXd exact_sol (const VectorXd &x);
VectorXd rhs_vec   (const int &dim);
MatrixXd A_mat     (const int &dim);
MatrixXd Ah_mat    (const int &dim);
double compare     (const VectorXd &numerical, const VectorXd &analytical);

#endif