import numpy as np
import chebyshev as ch
import spectral_solver as spec

''' -----------------------------------------------
         Create Spectral Solver Child Class
    (splits the system into two adjacent intervals)
    ------------------------------------------------ '''

# some further default parameters
N_MAX_RIGHT = 20
M_MAX_RIGHT = 500
R_MIN_LEFT  = 0.
R_MAX_LEFT  = 1.
R_MIN_RIGHT = 1.
R_MAX_RIGHT = 5.


class SpectralSystemTwoDomains(spec.SpectralSystem):

    '''   Constructor  '''
    def __init__(self, rho,
                       N_max = spec.N_MAX_DEFAULT,
                       M_max = spec.M_MAX_DEFAULT,
                       r_min = R_MIN_LEFT,
                       r_max = R_MAX_LEFT,
                       N_max_right = N_MAX_RIGHT,
                       M_max_right = M_MAX_RIGHT,
                       r_min_right = R_MIN_RIGHT,
                       r_max_right = R_MAX_RIGHT
                ):
        spec.SpectralSystem.__init__(self, rho, N_max, M_max, r_min, r_max)
        self.N_max_right  = N_max_right
        self.N_max_total  = self.N_max + self.N_max_right + 2
        self.M_max_right  = M_max_right
        self.r_min_right  = r_min_right
        self.r_max_right  = r_max_right
        self.dxdr_right   = 2. / (self.r_max_right - self.r_min_right)
        self.x_right      = self.gauss_lobatto_grid(self.N_max_right)
        self.r_right      = self.r_of_x(self.x_right, self.r_min_right, self.r_max_right)
        self.source       = self.source_func2(self.r, self.r_right)
        self.coeffs       = self.solve(self.spectral_matrix2(), self.source)

        self.r_test_left  = np.linspace(self.r_min + 1.0e-4, self.r_max - 1e-4, self.N_max+1)
        self.r_test_right = np.linspace(self.r_min_right + 1.0e-4, self.r_max_right, self.N_max_right+1)
        self.x_test_left  = self.x_of_r(self.r_test_left, self.r_min, self.r_max)
        self.x_test_right = self.x_of_r(self.r_test_right, self.r_min_right, self.r_max_right)

        self.r            = np.concatenate((self.r_test_left, self.r_test_right))
        self.source_test  = self.source_func2(self.r_test_left, self.r_test_right)

        self.y_left       = self.y_func(self.x_test_left, self.N_max, self.coeffs[:self.N_max + 1])
        self.y_right      = self.y_func(self.x_test_right, self.N_max_right, self.coeffs[N_max + 1:])
        self.y            = np.concatenate((self.y_left, self.y_right))

        self.y_r_left     = self.y_r_func(self.x_test_left, self.N_max, self.dxdr, self.coeffs[:self.N_max + 1])
        self.y_r_right    = self.y_r_func(self.x_test_right, self.N_max_right, self.dxdr_right, self.coeffs[N_max + 1:])
        self.y_r          = np.concatenate((self.y_r_left, self.y_r_right))

        self.y_rr_left    = self.y_rr_func(self.x_test_left, self.N_max, self.dxdr, self.coeffs[:self.N_max + 1])
        self.y_rr_right   = self.y_rr_func(self.x_test_right, self.N_max_right, self.dxdr_right, self.coeffs[N_max + 1:])
        self.y_rr         = np.concatenate((self.y_rr_left, self.y_rr_right))


    ''' Rebuild the source vector '''
    def source_func2(self, rleft, rright):

        r_left_len  = self.N_max + 1
        r_right_len = self.N_max_right + 1
        r_len       = r_left_len + r_right_len
        source_vec  = np.zeros(r_len)

        for i in range(1, self.N_max):
            source_vec[i] = 4. * np.pi * self.rho(rleft[i])

        for i in range(1, self.N_max_right ):
            k = i + r_left_len
            source_vec[k] = 4. * np.pi * self.rho(rright[i])

        return source_vec


    ''' Rebuild spectral matrix '''
    def spectral_matrix2(self):

        r_left_len  = self.N_max + 1
        r_right_len = self.N_max_right + 1
        r_len       = r_left_len + r_right_len
        T           = np.zeros((r_len, r_len))


        # Top part of the matrix
        for j in range(r_left_len):
            T[0,j] = self.dxdr * ch.Cheby_poly_x(j, self.x[0])

        for i in range(1, r_left_len - 1):
            for j in range(r_left_len):
                inner  = 0.5 * self.dxdr * ch.Cheby_poly_xx(j, self.x[i]) + \
                         1./self.r[i] * ch.Cheby_poly_x(j, self.x[i])
                T[i,j] = 2. * self.dxdr * inner

        # Interface
        for j in range(r_left_len):
            T[self.N_max,j]   = ch.Cheby_poly(j, self.x[-1])
            T[self.N_max+1,j] = self.dxdr * ch.Cheby_poly_x(j, self.x[-1])

        for j in range(r_right_len):
            k = j + r_left_len
            T[self.N_max,k]   = - ch.Cheby_poly(j, self.x[0])
            T[self.N_max+1,k] = - self.dxdr_right * ch.Cheby_poly_x(j, self.x[0])


        # Bottom part of the matrix
        for j in range(r_right_len):
            k = j + r_left_len
            T[-1,k] = self.r_max_right * self.dxdr_right * \
                                ch.Cheby_poly_x(j, self.x_right[-1]) + \
                                ch.Cheby_poly(j, self.x_right[-1])

        for i in range(1, r_right_len - 1):
            for j in range(r_right_len):

                row = i + r_left_len
                col = j + r_left_len

                inner_right = 0.5 * self.dxdr_right * ch.Cheby_poly_xx(j, self.x_right[i]) + \
                              1./self.r_right[i] * ch.Cheby_poly_x(j, self.x_right[i])
                T[row,col]  = 2. * self.dxdr_right * inner_right

        return T