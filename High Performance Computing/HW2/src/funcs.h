#ifndef FUNCS_H
#define FUNCS_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "GF.h"

#define PI (3.1415926535897932384626433832795)

/* Dissipation coefficient */
static const double DISS_COEFF = .01;

/* Fifth order dissipation */
#define DISSIP_any_5(f, i, di)                                             \
    (1.0 / 16.0 *                                                          \
     (-20 * f[i] + f[-3 * di + i] - 6 * f[-2 * di + i] + 15 * f[-di + i] + \
      15 * f[di + i] - 6 * f[2 * di + i] + f[3 * di + i]))

/* Fourth order accurate first-derivative operators */
#define D4CEN(f, i, di, h)                                               \
    ((f[-2 * di + i] - 8 * f[-di + i] + 8 * f[di + i] - f[2 * di + i]) * \
     oo12##h)

/* Here we consider an evolution system consisting of 6 evolved
 * variables:
 * U, Ku, Piu
 * V, Kv, Piv
 *
 * We provide a macro to identify which component of
 * the gridfunctions (gfs) corresponds to which variable.
 *
 * The 'unused' attribute tells the compiler not to issue a warning if
 * that variable is not used.
 *
 * With this macro, the user can directly read/write to arrays named
 * U, V, Ku, Kv, Piu, Piv, as well as U_dot, etc. The array size, ghost size,
 * and gridspacing is also given.
 */
#define DECLARE_VARS(_ptr_to_gfs)                                          \
    double __attribute__((unused)) *U = (_ptr_to_gfs)->vars[0]->new;       \
    double __attribute__((unused)) *Ku = (_ptr_to_gfs)->vars[1]->new;      \
    double __attribute__((unused)) *Piu = (_ptr_to_gfs)->vars[2]->new;     \
    double __attribute__((unused)) *V = (_ptr_to_gfs)->vars[3]->new;       \
    double __attribute__((unused)) *Kv = (_ptr_to_gfs)->vars[4]->new;      \
    double __attribute__((unused)) *Piv = (_ptr_to_gfs)->vars[5]->new;     \
    double __attribute__((unused)) *U_dot = (_ptr_to_gfs)->vars[0]->dot;   \
    double __attribute__((unused)) *Ku_dot = (_ptr_to_gfs)->vars[1]->dot;  \
    double __attribute__((unused)) *Piu_dot = (_ptr_to_gfs)->vars[2]->dot; \
    double __attribute__((unused)) *V_dot = (_ptr_to_gfs)->vars[3]->dot;   \
    double __attribute__((unused)) *Kv_dot = (_ptr_to_gfs)->vars[4]->dot;  \
    double __attribute__((unused)) *Piv_dot = (_ptr_to_gfs)->vars[5]->dot; \
    const int __attribute__((unused)) nsize = (_ptr_to_gfs)->n;            \
    const int __attribute__((unused)) ghost_zones =                        \
        (_ptr_to_gfs)->ghost_zones;                                        \
    const double __attribute__((unused)) dx = (_ptr_to_gfs)->dx;

void wave_eq_time_deriv(struct ngfs *gfs, double *psi_fun, double *phi_fun,
                        const double t0);
void periodic_sync(struct ngfs *gfs);
void apply_dissipation(struct ngfs *gfs);
void set_initial_data(struct ngfs *gfs, const double t);
void output_gfs(struct ngfs *gfs);
void update_potential(struct ngfs *gfs, double *psi_fun, double *phi_fun,
                      size_t potential);

#endif
