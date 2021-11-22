'''
    Define the Chebyshev polynomial and its first and second derivatives
'''

import numpy as np
from functools import lru_cache

@lru_cache
def Cheby_poly(n, x):
    assert -1. <= x <= 1.
    phi = np.arccos(x)
    return np.cos(n * phi)


@lru_cache
def Cheby_poly_x(n, x):
    assert -1. <= x <= 1.
    phi = np.arccos(x)
    eps = 1.e-10

    if np.fabs(1.-x) < eps:
        der = n*n
    elif np.fabs(-1.-x) < eps:
        der = (-1.)**(n+1) * n*n
    else:
        der = n * np.sin(n * phi) / np.sin(phi)

    return der


@lru_cache
def Cheby_poly_xx(n, x):
    assert -1. <= x <= 1.

    phi     = np.arccos(x)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    eps     = 1.e-10

    if np.fabs(1.-x) < eps:
        der = 1./3. * n*n * (n*n - 1.)
    elif np.fabs(-1.-x) < eps:
        der = (-1.)**(n+1) * 1./3. * n*n * (n*n - 1.)
    else:
        der1 = - n * n * np.cos(n * phi) / sin_phi**2
        der2 = n * cos_phi * np.sin(n * phi) / sin_phi**3
        der  = der1 + der2

    return der