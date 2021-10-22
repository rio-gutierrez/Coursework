// main.cpp
// Successive Overrelaxation (SOR) applied to the Poisson equation
//    \nabla^2 \psi = - 4\pi\rho
// assuming spherical symmetry.
// Created by Mario L Gutierrez on 10/03/21.

#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "functions.hpp"

using namespace std;
using namespace Eigen;

/*  Set this bool to 'false' for parts a) and b) of the problem,
    or set to 'true' for part c)   */
const bool N_TEST {false};
// const bool N_TEST {true};

/* -------------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------------*/
// Start of main function
int main(int argc, const char * argv[]) {

    /* ---------------------------------------------------------
    Parts a), b) of the Problem. Set N_TEST=false in global pars.
    -----------------------------------------------------------*/
    if (!N_TEST){

        const int n {1024};
        const double h {(rn-r0)/n};

        // MatrixXd A   = A_mat(n);    // use either A or \hat{A}
        MatrixXd Ah  = Ah_mat(n);
        VectorXd rhs = rhs_vec(n);

        // set initial guess for the algorithm
        VectorXd guess = VectorXd::Ones(n-1);

        /* ------------------------------------------------------
        Save residuals after 100, 200, and 1000 iterations of GS
        (omega = 1)
        ------------------------------------------------------ */
        VectorXd residuals_100  = SOR_RES(Ah, rhs, guess, n, 100);
        VectorXd residuals_200  = SOR_RES(Ah, rhs, guess, n, 200);
        VectorXd residuals_1000 = SOR_RES(Ah, rhs, guess, n, 1000);

        ofstream res100file ("../Data/res100.csv");
            for (int j{0}; j < residuals_100.size(); j++)
                    res100file << abs(residuals_100(j)) << endl;
        res100file.close();

        ofstream res200file ("../Data/res200.csv");
            for (int j{0}; j < residuals_200.size(); j++)
                    res200file << abs(residuals_200(j)) << endl;
        res200file.close();

        ofstream res1000file ("../Data/res1000.csv");
            for (int j{0}; j < residuals_1000.size(); j++)
                    res1000file << abs(residuals_1000(j)) << endl;
        res1000file.close();

        /* -------------------------------------------------------
        Save residuals after 100, 200, and 1000 iterations of SOR_RES
        (omega = 1.5)
        ---------------------------------------------------------*/
        VectorXd residuals_SOR_100  = SOR_RES(Ah, rhs, guess, n, 100, 1.5);
        VectorXd residuals_SOR_200  = SOR_RES(Ah, rhs, guess, n, 200, 1.5);
        VectorXd residuals_SOR_1000 = SOR_RES(Ah, rhs, guess, n, 1000, 1.5);

        ofstream res100SORfile ("../Data/res_SOR_100.csv");
            for (int j{0}; j < residuals_SOR_100.size(); j++)
                    res100SORfile << abs(residuals_SOR_100(j)) << endl;
        res100SORfile.close();

        ofstream res200SORfile ("../Data/res_SOR_200.csv");
            for (int j{0}; j < residuals_SOR_200.size(); j++)
                    res200SORfile << abs(residuals_SOR_200(j)) << endl;
        res200SORfile.close();

        ofstream res1000SORfile ("../Data/res_SOR_1000.csv");
            for (int j{0}; j < residuals_SOR_1000.size(); j++)
                    res1000SORfile << abs(residuals_SOR_1000(j)) << endl;
        res1000SORfile.close();


        /* -------------------------------------------------------
        Save array with r-values
        ---------------------------------------------------------*/
        VectorXd r_vals = VectorXd::LinSpaced(n-1,r0+h,rn-h);

        ofstream rfile ("../Data/r_vals.csv");
            for (int j{0}; j < r_vals.size(); j++)
                    rfile << r_vals(j) << endl;
        rfile.close();
    }
    // end of parts a) and b)


   /* -------------------------------------------------------
    Part c) of the Problem. Set N_TEST=true in global pars.
    ---------------------------------------------------------*/
    if (N_TEST){

        // Do N_TEST for different values of n
        VectorXi n_array (13);
        n_array << 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800;

        ofstream nfile ("../Data/Partc/nvals.csv");
        ofstream n_errors_file ("../Data/Partc/n_errors.csv");

        for (int N : n_array){
            // MatrixXd A     = A_mat(N);    // use either A or \hat{A}
            MatrixXd Ah    = Ah_mat(N);
            VectorXd rhs   = rhs_vec(N);
            VectorXd guess = VectorXd::Ones(N-1);

            // Save the full grid, this time including the endpoints as well
            VectorXd r_grid       = VectorXd::LinSpaced(N+1,r0,rn);
            VectorXd solution     = exact_sol(r_grid);
            VectorXd num_solution = SOR(Ah, rhs, guess, N, 100000, 1.9);
            double   n_error      = compare(num_solution, solution);

            // output numerical solution for N=800
            if (N==800){
                ofstream Psifile ("../Data/Partc/psi.csv");
                for (int j{0}; j < num_solution.size(); j++)
                        Psifile << num_solution(j) << endl;
                Psifile.close();
            }

            nfile << N << endl;
            n_errors_file << n_error << endl;
        }
        nfile.close();
        n_errors_file.close();
    }
    // end of part c)

}  //end of main

/* -------------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------------*/