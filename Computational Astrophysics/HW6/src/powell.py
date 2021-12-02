import numpy as np
from scipy.optimize import minimize_scalar

'''
    Powell's Set Method
'''

def powell(p, func, tol  = 1.e-16, it_max = 10000):
    '''
    Initial parameters:
        p      (initial position)
        func   (function of interest)
        tol    (user-defined tolerance)
        it_max (max number of iterations allowed)

    Choose initial direction vector sets u = [e0, e1],
    where e0, e1 are the unit coordinate vectors

    Then perform 1D optimization over the line
        func(p[0] + s * e0[0], p[1] + s * e0[1])

    Take the parameter s that yields the minimum of func
    along the above line and set a new point
       p1 = p + s * e0

    From p1 repeat the above process along the direction e1
    to arrive at new point p2. Then define a new direction
    p2 - p0, and optimize again along this direction to
    arrive at another point p3.

    Rinse and repeat until a minimum has been found within
    set tolerance.

    Output:
        p_f    (found local minimum of the function)
        f(p_f) (value of the function at p_f)
        it     (number of iterations needed to find p_f)
    '''

    # define unit vectors
    e0 = np.array((1., 0.))
    e1 = np.array((0., 1.))

    # initialize direction sets u_i
    u = [e0, e1]

    diff = 1.
    it   = 0

    while diff > tol:

        it += 1
        if it == it_max:
            print(f'Aborting after {it} iterations...Unable to find minimum.')
            break

        p_old = p
        f_old = func(p[0],p[1])

        for ui in u:
            # perform 1D optimization
            f_1dim = lambda s : func(p[0] + s * ui[0], p[1] + s * ui[1])
            brent  = minimize_scalar(f_1dim, method='brent')
            # update to new point p0 -> p0 + s*u
            s = brent.x
            p = p + s * ui

        # update direction sets u_i
        u[0] = u[1]
        u[1] = p - p_old

        # perform a further 1D optimization, along the new direction of u[1]
        f_1dim = lambda s : func(p[0] + s * u[1][0], p[1] + s * u[1][1])
        brent  = minimize_scalar(f_1dim, method='brent')
        s      = brent.x
        p      = p + s * u[1]

        # evaluate function at the new location and compare with previous value
        f_new = func(p[0],p[1])
        diff  = np.abs(f_old - f_new)

    return p, func(p[0], p[1]), it