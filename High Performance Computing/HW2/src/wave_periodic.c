#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "GF.h"
#include "RK4.h"
#include "funcs.h"

int main(int argc, char** argv) {
    const int ghost_zones = 3;  // number of ghostzones
    struct ngfs gfs = {0, 0., 0.,
                       0, 0,  NULL};  // initialization; struct defined on GF.h
    const int nvars = 6;

    /* Read runtime gridsize */
    if (argc != 2) {
        fprintf(stderr,
                "Must supply exacly one argument, "
                "the number of points in the computational domain\n");
        return EXIT_FAILURE;
    }

    // n is the number of intervals in the computational domain [0,1)
    // The number of points will be n + 2 * ghost_zones

    char* char_ptr;
    const int n = strtoll(argv[1], &char_ptr, 10);

    if (errno == ERANGE) {
        printf("\n\terrno = %d\n", errno);
        fprintf(stderr,
                "We are out of range. Your n_max value is too large.\n");
        return EXIT_FAILURE;
    }

    if (n < 2 * ghost_zones + 1) {
        fprintf(stderr,
                "The number of points  must be larger twice the number of "
                "ghostzones\n");
        return EXIT_FAILURE;
    }

    /* physical extent of the grid is [0,1)  plus ghostzones */
    const double dx = 1. / (double)(n);
    const double x0 =
        -ghost_zones * dx;  // x0 is the coordinate of the first ghost point

    ngfs_allocate(nvars, n + 2 * ghost_zones, ghost_zones, x0, dx,
                  &gfs);  // the gridsize is n + 2 * ghost_zones

    double t = 0;

    /* Courant stability conditions requires dt < k * dx, where k is
     * usually of order 1
     */
    const double dt = .5 * dx;

    /* Initialization step. */
    set_initial_data(&gfs, t);

    // allocate memory for the potential
    double* Psi = calloc(n + 2 * ghost_zones, sizeof *Psi);
    double* Phi = calloc(n + 2 * ghost_zones, sizeof *Phi);

    // Make sure the memory has been successfully  allocated
    if ((Psi == NULL) || (Phi == NULL)) {
        printf("\tMemory not allocated. Perhaps n is too large...\n");
        return EXIT_FAILURE;
    }

    const double time_interval = .1;

    for (; t <= 30;) {
        /* first update values of potentials
         *  (last arg is 0 if we have a vanishing potential,
         *  otherwise we set it to 1...see function implementation)
         */

        update_potential(&gfs, Psi, Phi, 0);

        // then take a time-step with RK4
        RK4_Step(&gfs, t, dt, Psi, Phi, wave_eq_time_deriv);
        t += dt;

        // output every 0.1 time units
        if (t / time_interval == t * 10) output_gfs(&gfs);
    }

    /* We're done. Deallocate the gridfunctions and exit */

    ngfs_deallocate(&gfs);
    free(Psi);
    free(Phi);
    Psi = NULL;
    Phi = NULL;

    return 0;
}
