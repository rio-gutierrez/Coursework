#ifndef GF_H
#define GF_H

#include <stdlib.h>

/* struct gf holds the data for a single gridfunction */
struct gf {
    size_t n;        /* length of each array */
    int ghost_zones; /* Ghost size of algorithm */
    double *old;     /* for old values */
    double *new;     /* for updated values */
    double *dot;     /* time derivatives  */
    double *K1;      /* RK variables */
    double *K2;      /* RK variables */
    double *K3;      /* RK variables */
    double *K4;      /* RK variables */
};

/* struct ngfs holds data for all gridfunctions */
struct ngfs {
    int nvars;        /* How many variables */
    double dx;        /* x coordinate grid spacing */
    double x0;        /* x coordinate "origin" */
    size_t n;         /* length of arrays */
    int ghost_zones;  /* ghost size */
    struct gf **vars; /* pointer to pointer (see DECLARE_VARS on funcs.h) */
};

int gf_allocate(size_t n, int ghost_zones, struct gf *gptr);
int gf_deallocate(struct gf *gptr);

int ngfs_allocate(int nvars, size_t n, int ghost_zones, const double x0,
                  const double dx, struct ngfs *ptr);

int ngfs_deallocate(struct ngfs *ptr);
#endif
