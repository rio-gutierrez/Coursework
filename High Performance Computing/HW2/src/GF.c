#include "GF.h"

#include <assert.h>
#include <stdlib.h>

/* Allocate gridfunctions. External routines call this function */
int ngfs_allocate(int nvars, size_t n, int ghost_zones, const double x0,
                  const double dx, struct ngfs *ptr) {
    assert(ptr->vars == NULL);
    ptr->nvars = nvars;
    ptr->n = n;
    ptr->ghost_zones = ghost_zones;
    ptr->dx = dx;
    ptr->x0 = x0;

    ptr->vars = calloc(nvars, sizeof(struct gf *));
    assert(ptr->nvars);

    for (int i = 0; i < nvars; i++) {
        ptr->vars[i] = calloc(1, sizeof(struct gf));
        gf_allocate(n, ghost_zones, ptr->vars[i]);
    }

    return 0;
}

/* free gridfunctions */
int ngfs_deallocate(struct ngfs *ptr) {
    assert(ptr->vars);

    for (int i = 0; i < ptr->nvars; i++) {
        gf_deallocate(ptr->vars[i]);
        free(ptr->vars[i]);
    }
    free(ptr->vars);
    ptr->nvars = 0;
    ptr->n = 0;
    ptr->ghost_zones = 0;
    ptr->vars = NULL;
    return 0;
}

/* Allocate individual gridfunction. External routines generally don't call this
 * function */
int gf_allocate(size_t n, int ghost_zones, struct gf *gptr) {
    gptr->n = n;
    gptr->ghost_zones = ghost_zones;
    gptr->old = calloc(n, sizeof(double));
    gptr->new = calloc(n, sizeof(double));
    gptr->K1 = calloc(n, sizeof(double));
    gptr->K2 = calloc(n, sizeof(double));
    gptr->K3 = calloc(n, sizeof(double));
    gptr->K4 = calloc(n, sizeof(double));
    gptr->dot = NULL;

    assert(gptr->old);
    assert(gptr->new);
    assert(gptr->K1);
    assert(gptr->K2);
    assert(gptr->K3);
    assert(gptr->K4);

    return 0;
}

/* deallocate individual gridfunction. External routines generally don't call
 * this function */
int gf_deallocate(struct gf *gptr) {
    free(gptr->old);
    free(gptr->new);
    free(gptr->K1);
    free(gptr->K2);
    free(gptr->K3);
    free(gptr->K4);

    gptr->n = 0;
    gptr->ghost_zones = 0;
    gptr->old = NULL;
    gptr->new = NULL;
    gptr->dot = NULL;
    gptr->K1 = NULL;
    gptr->K2 = NULL;
    gptr->K3 = NULL;
    gptr->K4 = NULL;

    return 0;
}
