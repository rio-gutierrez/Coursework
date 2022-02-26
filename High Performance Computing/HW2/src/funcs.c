#include "funcs.h"

void update_potential(struct ngfs *gfs, double *psi_fun, double *phi_fun,
                      size_t potential) {
    /* The argument 'potential' is set to zero if we want a vanishing potential.
     * Otherwise we set it to 1 if we want the potentials given
     * by Eq. (5) from the Latex report.
     */
    DECLARE_VARS(gfs);

    for (int i = 0; i < nsize; i++) {
        switch (potential) {
            case 0: {
                psi_fun[i] = 0.;
                phi_fun[i] = 0.;
            } break;
            case 1: {
                psi_fun[i] = V[i] * V[i] * exp(-2. * U[i] * U[i] - V[i] * V[i]);
                phi_fun[i] = U[i] * U[i] * exp(-2. * U[i] * U[i] - V[i] * V[i]);
            } break;
            default: {
                psi_fun[i] = 0.;
                phi_fun[i] = 0.;
            }
        }
    }
}

/* Initialization step. The user provides the initialization */
void set_initial_data(struct ngfs *gfs, const double t) {
    DECLARE_VARS(gfs);

    for (int i = 0; i < nsize; i++) {
        // The first "ghost_zones" points are ghostzones, x=0 occurs at index
        // "ghost_zones"
        const double x = 0 + dx * (i - gfs->ghost_zones);
        const double theta_min = 2 * PI * (x - .25);
        const double theta_plus = 2 * PI * (x + .25);

        // Eqs 3a - 3f on Latex doc
        U[i] = exp(-10 * (1 - cos(theta_min)));
        Piu[i] = -20 * PI * sin(theta_min) * exp(-10 * (1 - cos(theta_min)));
        Ku[i] = 20 * PI * sin(theta_min) * exp(-10 * (1 - cos(theta_min)));
        V[i] = exp(-10 * (1 - cos(theta_plus)));
        Piv[i] = -20 * PI * sin(theta_plus) * exp(-10 * (1 - cos(theta_plus)));
        Kv[i] = -20 * PI * sin(theta_plus) * exp(-10 * (1 - cos(theta_plus)));
    }
}

// Calculate the time derivative of the evolved variables:
void wave_eq_time_deriv(struct ngfs *gfs, double *psi_fun, double *phi_fun,
                        const double t0) {
    DECLARE_VARS(gfs);

    const double h = dx;
    const double oo12h = 1. / (12. * h);

    for (int i = ghost_zones; i < nsize - ghost_zones; i++) {
        U_dot[i] = Ku[i];
        V_dot[i] = Kv[i];
    }

    for (int i = ghost_zones; i < nsize - ghost_zones; i++) {
        /* 4th-order centered stencil for first-derivatives */
        const double Piu_x = D4CEN(Piu, i, 1, h);
        const double Piv_x = D4CEN(Piv, i, 1, h);
        const double Ku_x = D4CEN(Ku, i, 1, h);
        const double Kv_x = D4CEN(Kv, i, 1, h);

        Piu_dot[i] = Ku_x;
        Piv_dot[i] = Kv_x;
        Ku_dot[i] = Piu_x + psi_fun[i];
        Kv_dot[i] = Piv_x + phi_fun[i];
    }

    /* dissipation */
    apply_dissipation(gfs);

    /* Boundary conditions  */
    periodic_sync(gfs);
}

// output gridfunctions U and V to .asc files
void output_gfs(struct ngfs *gfs) {
    DECLARE_VARS(gfs);
    static int counter = 0;
    char name_buff[1024];
    snprintf(name_buff, 1024, "%07d.asc", counter);
    FILE *ofile = fopen(name_buff, "w");

    // we don't output the ghostzones
    for (int i = ghost_zones; i < nsize - ghost_zones; ++i)
        fprintf(ofile, "%20.16e %20.16e %20.16e\n", (i - ghost_zones) * dx,
                U[i], V[i]);

    /* the output file contains the x_grid on the first column,
     * U values on the second, and V values on the third
     */
    ++counter;

    fclose(ofile);
    ofile = NULL;
}

/* The routines below generally don't have to be adjusted  */

/* periodic_sync
 * Applies periodic boundary conditions to all grid functions
 */
void periodic_sync(struct ngfs *gfs) {
    const int n = gfs->n;
    const int ghost_zones = gfs->ghost_zones;
    const int N = n - 2 * ghost_zones;  // number of subintervals
    const int nvars = gfs->nvars;

    for (int i = 0; i < ghost_zones; i++) {
        const int im = ((i - ghost_zones) + N) % N + ghost_zones;
        const int ip = ((n - i - 1 - ghost_zones)) % N + ghost_zones;

        for (int var = 0; var < nvars; var++) {
            gfs->vars[var]->dot[i] = gfs->vars[var]->dot[im];
            gfs->vars[var]->dot[n - i - 1] = gfs->vars[var]->dot[ip];
        }
    }
}

/* dissipation step. Damps high-frequency modes. Only needed in
 * nonlinear case. */
void apply_dissipation(struct ngfs *gfs) {
    const int n = gfs->n;
    const int ghost_zones = gfs->ghost_zones;
    const double dx = gfs->dx;
    const double coef = DISS_COEFF / dx;

    for (int var = 0; var < gfs->nvars; ++var) {
        double *var_dot = gfs->vars[var]->dot;
        const double *var_new = gfs->vars[var]->new;

        for (int i = ghost_zones; i < n - ghost_zones; i++) {
            var_dot[i] += coef * DISSIP_any_5(var_new, i, 1);
        }
    }
}
