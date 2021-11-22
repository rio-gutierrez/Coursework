import numpy as np
import chebyshev as ch

# establish some default parameters
N_MAX_DEFAULT = 100
M_MAX_DEFAULT = 1000
R_MIN_DEFAULT = 0.
R_MAX_DEFAULT = 20.

''' ----------------------------
    Create Spectral Solver Class
    ---------------------------- '''

class SpectralSystem():

    '''   Constructor  '''
    def __init__(self, rho,
                       N_max = N_MAX_DEFAULT,
                       M_max = M_MAX_DEFAULT,
                       r_min = R_MIN_DEFAULT,
                       r_max = R_MAX_DEFAULT
                ):
        self.rho         = rho
        self.N_max       = N_max
        self.M_max       = M_max
        self.r_min       = r_min
        self.r_max       = r_max
        self.dxdr        = 2. / (self.r_max - self.r_min)
        self.x           = self.gauss_lobatto_grid(self.N_max)
        self.r           = self.r_of_x(self.x, self.r_min, self.r_max)
        self.source      = self.source_func(self.r)
        self.coeffs      = self.solve(self.spectral_matrix(), self.source)
        self.r_test      = self.r_test_grid(self.r_min, self.r_max, self.M_max)
        self.x_test      = self.x_of_r(self.r_test, self.r_min, self.r_max)
        self.source_test = self.source_func(self.r_test)
        self.y           = self.y_func(self.x_test, self.N_max, self.coeffs)
        self.y_r         = self.y_r_func(self.x_test, self.N_max, self.dxdr, self.coeffs)
        self.y_rr        = self.y_rr_func(self.x_test, self.N_max, self.dxdr, self.coeffs)


    '''   Class Methods   '''
    # Change of coordinates mappings
    def x_of_r(self, r, rmin, rmax):
        return 2. * (r - rmin)/(rmax - rmin) - 1.

    def r_of_x(self, x, rmin, rmax):
        return rmin + (x + 1.) * (rmax - rmin) / 2.


    # generate grid with Guass-Lobatto roots
    def gauss_lobatto_grid(self, nmax):

        grid = np.zeros((nmax + 1))
        for i in range(nmax + 1):
            grid[i] = np.cos((nmax - i) * np.pi / nmax)

        return grid


    # build the source vector 4*pi*rho
    def source_func(self, r):
        source_vec = np.zeros_like(r)
        for i in range(1, len(source_vec)-1):
            source_vec[i] = 4. * np.pi * self.rho(r[i])
        return source_vec


    # build spectral matrix
    def spectral_matrix(self):

        T = np.zeros((self.N_max+1, self.N_max+1))

        # top and bottom rows
        for j in range(self.N_max + 1):
            T[0,j]          = self.dxdr * ch.Cheby_poly_x(j, self.x[0])
            T[self.N_max,j] = self.r_max * self.dxdr * \
                            ch.Cheby_poly_x(j, self.x[self.N_max]) + \
                            ch.Cheby_poly(j, self.x[self.N_max])
        # rest of the matrix
        for i in range(1,self.N_max):
            for j in range(self.N_max + 1):
                inner  = 0.5 * self.dxdr * ch.Cheby_poly_xx(j, self.x[i]) + \
                         1./self.r[i] * ch.Cheby_poly_x(j, self.x[i])
                T[i,j] = 2. * self.dxdr * inner

        return T


    # Compute the solution coefficients c_n
    def solve(self, mat, src):
        C = np.linalg.solve(mat, -src)
        return C



    '''
        The following methods are only used when testing
        the numerical solution on some test grid
    '''

    def r_test_grid(self, rmin, rmax, mmax):
        grid     = np.zeros(mmax+1)
        grid[0]  = rmin
        grid[-1] = rmax

        for i in range(1, mmax):
            grid[i] = i  * rmax / mmax

        return grid



    '''
        Reconstruct the solution y(r) and its first and second derivatives
        evaluated on test grid
    '''

    def y_func(self, xgrid, nmax, c_n):

        temp  = np.zeros(nmax+1)
        y_vec = np.zeros_like(xgrid)

        for i in range(len(y_vec)):
            for j in range(len(temp)):
                temp[j] = ch.Cheby_poly(j, xgrid[i])
            y_vec[i] = np.dot(c_n,temp)

        return y_vec


    def y_r_func(self, xgrid, nmax, der, c_n):

        temp  = np.zeros(nmax+1)
        y_r_vec = np.zeros_like(xgrid)

        for i in range(len(y_r_vec)):
            for j in range(len(temp)):
                temp[j] = der * ch.Cheby_poly_x(j, xgrid[i])
            y_r_vec[i] = np.dot(c_n,temp)

        return y_r_vec


    def y_rr_func(self, xgrid, nmax, der, c_n):

        temp  = np.zeros(nmax+1)
        y_rr_vec = np.zeros_like(xgrid)

        for i in range(len(y_rr_vec)):
            for j in range(len(temp)):
                temp[j] = der * der * ch.Cheby_poly_xx(j, xgrid[i])
            y_rr_vec[i] = np.dot(c_n,temp)

        return y_rr_vec