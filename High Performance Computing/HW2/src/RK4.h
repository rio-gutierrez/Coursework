#ifndef RK4_H
#define RK4_H

#include <stdlib.h>

#include "GF.h"

typedef void (*dotfunction)(struct ngfs* gfs, double* psi_fun, double* phi_fun,
                            const double t0);

int RK4_Step(struct ngfs* gfs, const double t0, const double dt,
             double* psi_fun, double* phi_fun, dotfunction time_deriv);

#endif
