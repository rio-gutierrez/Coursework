#include "functions.hpp"
/* -------------------------------------------------------------------------------
                        FUNCTION IMPLEMENTATIONS
 -------------------------------------------------------------------------------*/

// Generate A matrix
MatrixXd A_mat(const int &dim){

    const double h {(rn-r0)/dim};
    MatrixXd A(dim-1,dim-1);

    double r_i {};
    double r_j {};
    double theta_p {};
    double theta_m {};

    // Generate (n-1)x(n-1) matrix A (from Eq.5)
    for (int i {0}; i < dim-1; i++) {

        r_i     = r0 + (i+1)*h;
        r_j     = r0 + (i+2)*h;
        theta_p = 1. + h/r_i;
        theta_m = 1. - h/r_j;

        if (i == 0)
            A(i,i) = -theta_p;
        else
            A(i,i) = -2.;

        if (i != dim-2){
            A(i,i+1) = theta_p;
            A(i+1,i) = theta_m;
        }
    }
    return A;
}

// Generate \hat{A} matrix
MatrixXd Ah_mat(const int &dim){

    const double h {(rn-r0)/dim};
    const double Xi_1 = 4./3. * (1. - h/(r0+h)) - 2.;           //\Xi_i from Eq. (8)
    const double Xi_2 = 1. + h/(r0+h) - 1./3. * (1. - h/(r0+h));
    double r_i {};
    double r_j {};
    double theta_p {};
    double theta_m {};

    MatrixXd Ah(dim-1,dim-1);

    // Generate (n-1)x(n-1) matrix \hat{A} (from Eq.9)
    for (int i {0}; i < dim-1; i++) {

        r_i     = r0 + (i+1)*h;
        r_j     = r0 + (i+2)*h;
        theta_p = 1. + h/r_i;
        theta_m = 1. - h/r_j;

        if (i == 0){
            Ah(i,i)   = Xi_1;
            Ah(i,i+1) = Xi_2;
            Ah(i+1,i) = theta_m;
        }
        else
            Ah(i,i) = -2.;

        if (i != dim-2 && i != 0){
            Ah(i,i+1) = theta_p;
            Ah(i+1,i) = theta_m;
        }
    }

    return Ah;
}

// Generate (\tilde{\varrho}}) vector
VectorXd rhs_vec(const int &dim){

    const double h {(rn-r0)/dim};
    VectorXd rhs(dim-1);
    double r_i {};
    double rho_i {};
    double theta_p {};

    // Generate (n-1)x1 rhs vector (\tilde{\varrho}}) (from Eq.5 or 9)
    for (int i {0}; i < dim-1; i++) {

        r_i     = r0 + (i+1)*h;
        rho_i   = - 4. * M_PI * h*h * 1./(pow(r_i,4));
        theta_p = 1. + h/r_i;

        if (i == dim-2)
            rhs(i) = rho_i - theta_p;
        else
            rhs(i) = rho_i;
    }
    return rhs;
}

// SOR_RES function implementation
VectorXd SOR_RES(const MatrixXd &A, const VectorXd &b,
                 const VectorXd &x0, const int &dim,
                 const int max_it, const double omega){
    /*
    ** SOR_RES routine that returns the residual b-Ax from a system Ax=b at the end of SOR algorithm**
    INPUTS:
        A: (dim-1)x(dim-1) matrix
        b: (dim-1)-vector
        x0: initial guess for iterative SOR solver
        dim: num of cells in the grid
        max_it: max number of itartions allowed
        omega: relaxation parameter (omega=1 => Gauss-Seidel)
    OUTPUT:
        diff: the residuals b - Ax
    */

    const double h {(rn-r0)/dim};
    VectorXd Psi_new = x0;
    size_t it {0};

    do{

        VectorXd Psi_old = Psi_new;

        for (int i {0}; i < dim-1; i++){

            VectorXd ai2  = A.row(i)(seq(i+1,dim-2));
            VectorXd psi2 = Psi_old(seq(i+1,dim-2));
            double sum2   = ai2.dot(psi2);

            if (i == 0)
                Psi_new(i) = omega/A(i,i) * (b(i) - sum2) + (1. - omega) * Psi_old(i);
            else{
                VectorXd ai1  = A.row(i)(seq(0,i-1));
                VectorXd psi1 = Psi_new(seq(0,i-1));
                double sum1   = ai1.dot(psi1);

                Psi_new(i) = omega/A(i,i) * (b(i) - sum1 - sum2) + (1. - omega) * Psi_old(i);
            }
        }

        it+=1;
    } while (it <= max_it);

    // residuals to be output
    VectorXd diff = b - A * Psi_new;
    diff = diff/(h*h);

    return diff;
}

// SOR function implementation
VectorXd SOR(const MatrixXd &A, const VectorXd &b,
             const VectorXd &x0, const int &dim,
             const int max_it, const double omega,
             const double tol){
    /*
    ** SOR routine that returns the solution x to the system Ax=b at the end of SOR algorithm **
    INPUTS:
        A: (dim-1)x(dim-1) matrix
        b: (dim-1)-vector
        x0: initial guess for iterative SOR solver
        dim: num of cells in the grid
        max_it: max number of itartions allowed
        omega: relaxation parameter (omega=1 => Gauss-Seidel)
        tol: user-defined error tolerance
    OUTPUT:
        Psi_full: the solution itself
    */
    VectorXd Psi_new = x0;
    size_t it {0};
    double residual_norm {};
    const double h {(rn-r0)/dim};

    do{

        VectorXd Psi_old = Psi_new;

        for (int i {0}; i < dim-1; i++){

            VectorXd ai2  = A.row(i)(seq(i+1,dim-2));
            VectorXd psi2 = Psi_old(seq(i+1,dim-2));
            double sum2   = ai2.dot(psi2);

            if (i == 0)
                Psi_new(i) = omega/A(i,i) * (b(i) - sum2) + (1. - omega) * Psi_old(i);
            else{
                VectorXd ai1  = A.row(i)(seq(0,i-1));
                VectorXd psi1 = Psi_new(seq(0,i-1));
                double sum1   = ai1.dot(psi1);

                Psi_new(i)    = omega/A(i,i) * (b(i) - sum1 - sum2) + (1. - omega) * Psi_old(i);
            }
        }

        VectorXd diff = b - A * Psi_new;
        diff = diff/(h*h);
        residual_norm = diff.lpNorm<Infinity>()/b.lpNorm<Infinity>();

        cout << "\nresidual = " << residual_norm << " at iteration " << it << endl;

        if (residual_norm <= tol){
            cout << "The relative residual norm " << residual_norm <<
                    " has reached the derired tolerance. \
                     Convergence successful after " << it <<
                     " iterations!\n" << endl;
            break;
        }

        it+=1;
    } while (it <= max_it);

    // append solution at the boundaries
    VectorXd Psi_full(dim+1);
    // double Psi_0 = Psi_new(0);                          // if using Neumann condition (4)
    double Psi_0 = (1./3.) * (4. * Psi_new(0) - Psi_new(1));   // if using Neumann condition (8)
    double Psi_n = 1.;
    Psi_full << Psi_0, Psi_new, Psi_n;

    return Psi_full;
}

// Evaluate exact solution (Eq. 17) on the grid
VectorXd exact_sol(const VectorXd &grid){
    VectorXd sol (grid.size());
    for (int i {0}; i < grid.size(); i++)
        sol(i) = - 2. * M_PI/(grid(i)*grid(i)) + 4. * M_PI/(grid(i)) + 1. - 19. * M_PI/50.;
    return sol;
}

// Compare L^{\infty} norm of the error between numerical and exact solutions
double compare(const VectorXd &numerical, const VectorXd &analytical){
    VectorXd diff = numerical - analytical;
    double error = diff.lpNorm<Infinity>();
    return error;
}