import numpy as np
from scipy.optimize import minimize_scalar

'''
    Conjugate Gradient Method
'''

def conj_grad(p, func, grad, tol = 1.e-16, it_max = 10000):
    '''
    Initial parameters:
        p      (initial position)
        func   (function of interest)
        grad   (gradient of func)
        tol    (user-defined tolerance)
        it_max (max number of iterations allowed)

    Initial direction vector g0:
        g_0 = - grad(p)

    Set h_0 = g_0

    Then
        g_{i+1} = - grad(p_{i+1})
        h_{i+1} = g_{i+1} + \gamma_i h_i

    where
        gamma_i = dot(g_{i+1} - g_i, g_{i+1}) / dot(g_i, g_i)

    Output:
        p_f    (found local minimum of the function)
        f(p_f) (value of the function at p_f)
        it     (number of iterations needed to find p_f)
    '''

    # here we need to use '[:,0]' because of the way sympy's lambdify works
    # if not using lambdify, just remove '[0,:]'
    g = - grad(p[0], p[1])[:,0]
    h = g

    diff = 1.
    it   = 0

    while diff > tol:

        it += 1
        if it == it_max:
            print(f'Aborting after {it} iterations...Unable to find minimum.')
            break

        f_old = func(p[0],p[1])
        g_old = g

        # perform 1D optimization
        f_1dim = lambda s : func(p[0] + s * h[0], p[1] + s * h[1])
        brent  = minimize_scalar(f_1dim, method='brent')
        # update to new point p0 -> p0 + s*h
        s = brent.x
        p = p + s * h

        # evaluate function at the new location and compare with previous value
        f_new = func(p[0],p[1])
        diff  = np.abs(f_old - f_new)

        # if diff > tol not yet satisfied continue loop...
        g      = - grad(p[0], p[1])[:,0]
        gamma  = np.dot(g - g_old, g) / np.dot(g_old, g_old)
        h      = g + gamma * h

    return p, func(p[0], p[1]), it