#include "RK4.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "GF.h"

/* Generic RK4. Updates nvars different variables at once
 * Arguments : Pointer to struct ngfs
 * time value (generally unused)
 * timestep size
 * pointer to function that calculates time derivatives
 *
 * Return value is always zero
 *
 * Each step below consists if two nested loops
 * Loop over variables:
 *   Loop over individual points
 *     update variable at individual point
 */

int RK4_Step(struct ngfs* gfs, const double t0, const double dt,
             double* psi_fun, double* phi_fun, dotfunction time_deriv) {
    const int nvars = gfs->nvars;
    const int n = gfs->n;

    /* When entering this function, the 'old' values are actually stored
     * in the 'new' locations. So we copy them to the appropriate place.
     * We also make sure the dot pointer points to K1
     */
    for (int i = 0; i < nvars; i++) {
        gfs->vars[i]->dot = gfs->vars[i]->K1;

        for (int j = 0; j < n; j++) gfs->vars[i]->old[j] = gfs->vars[i]->new[j];
    }

    time_deriv(gfs, psi_fun, phi_fun, t0); /* this will fill in K1 */

    /* Step 1 of RK4. Also make sure dot points to
     *  K2 in preparation for step 2
     */

    for (int i = 0; i < nvars; i++) {
        gfs->vars[i]->dot = gfs->vars[i]->K2;

        for (int j = 0; j < n; j++) {
            gfs->vars[i]->new[j] =
                gfs->vars[i]->old[j] + .5 * dt * gfs->vars[i]->K1[j];
        }
    }

    time_deriv(gfs, psi_fun, phi_fun, t0 + .5 * dt); /* this will fill in K2 */

    /* Step 2 of RK4. Also make sure dot points
     *  to K3 in preparation for step 3
     */
    for (int i = 0; i < nvars; i++) {
        gfs->vars[i]->dot = gfs->vars[i]->K3;

        for (int j = 0; j < n; j++) {
            gfs->vars[i]->new[j] =
                gfs->vars[i]->old[j] + .5 * dt * gfs->vars[i]->K2[j];
        }
    }

    time_deriv(gfs, psi_fun, phi_fun, t0 + .5 * dt); /* this will fill in K3 */

    /* ditto */
    for (int i = 0; i < nvars; i++) {
        gfs->vars[i]->dot = gfs->vars[i]->K4;

        for (int j = 0; j < n; j++)
            gfs->vars[i]->new[j] =
                gfs->vars[i]->old[j] + dt * gfs->vars[i]->K3[j];
    }

    time_deriv(gfs, psi_fun, phi_fun, t0 + dt); /* this will fill in K4 */

    /* final step */

    for (int i = 0; i < nvars; i++) {
        for (int j = 0; j < n; j++) {
            gfs->vars[i]->new[j] =
                gfs->vars[i]->old[j] +
                1. / 6. *
                    dt*(gfs->vars[i]->K1[j] + 2. * gfs->vars[i]->K2[j] +
                        2. * gfs->vars[i]->K3[j] + gfs->vars[i]->K4[j]);
        }
    }

    /* 'new' now contains the new value */

    return 0;
}
