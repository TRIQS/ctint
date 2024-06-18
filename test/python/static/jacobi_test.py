################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011-2014 by M. Ferrero, O. Parcollet
# Copyright (C) 2018 by Simons Foundation
#   author : N. Wentzell
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
from triqs_ctint import Solver
from triqs.gf import *
from h5 import *
from triqs.operators import *
from numpy import matrix, array
import numpy as np
from triqs.utility.comparison_tests import *
import unittest

def jacobi_numerical(f,x,dx):
    n = len(x)
    jac = np.zeros((n,n))
    fi_x = f(x)
    for j in range(n):
        x_dx = x.copy()
        x_dx[j] = x[j] + dx 
        jac[:,j] = (f(x_dx) - fi_x)/dx
    return jac

class test_Gf_Base_Op(unittest.TestCase):

    def test_jacobi(self):
        dx = 1e-8
        beta = 10
        block_names = ['dn','up']
        orb_names = [0,1,2,3]
        mu = 0.5
        eps = matrix([[0.2,0.1,0.1,0.1],
                      [0.1,0.2,0.1,0.1],
                      [0.1,0.1,0.2,0.1],
                      [0.1,0.1,0.1,0.2]])
        gf_struct = [(bl, orb_names) for bl in block_names]
        S = Solver(beta = beta,
                   gf_struct = gf_struct,
                   n_iw = 1050,
                   n_tau = 100001)
        # --------- Initialize the non-interacting Green's function ----------
        for bl, g_bl in S.G0_iw: g_bl << inverse(iOmega_n + mu - inverse(iOmega_n - eps));

        Hint_1 = Operator()
        Hint_2 = Operator()
        for i in orb_names:
            Hint_1 += n('up',i)*n('dn',i)
        Hint_2 = Hint_1
        for i in orb_names:
            for j in orb_names:
                if i != j:
                    Hint_2 += n('up',i)*n('up',j)
                    Hint_2 += n('dn',i)*n('dn',j)

        for Hint in [Hint_1,Hint_2]:
            n_terms = len(list(Hint))
            solve_params = {
                'n_cycles': 1000,
                'h_int' : Hint,
                'n_s': 1
            }
            for i in range(10):
                np.random.seed(i)
                x = np.random.rand(3 * n_terms)
                assert_arrays_are_close(
                    S.jacobi(x,solve_params),
                    jacobi_numerical(lambda x: S.f(x,solve_params),x,dx),
                    1e-5,
                )


if __name__ == '__main__':
    unittest.main()
