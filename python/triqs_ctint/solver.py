from .solver_core import SolverCore
from triqs.gf import *
from triqs.utility import mpi

import numpy as np
from scipy.optimize import root


# === Some utility functions

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        po = np.get_printoptions()
        np.set_printoptions(precision=4)
        print(arg)
        np.set_printoptions(**po)

# === The SolverCore Wrapper

class Solver(SolverCore):

    def __init__(self, **constr_params):
        self.eval_count = 0
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : list of pairs [ (str,int), ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and the size of that block.
                    For example: ``[ ('up', 3), ('down', 3) ]``.
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
        n_tau : integer, optional
               Number of imaginary time points used for the Green's functions.
        use_D : bool, optional
               Use dynamic density-density interaction given via S.D0_iw[bl1, bl2][i,j]
        use_Jperp : bool, optional
               Use dynamic spin-spin interaction given via S.Jperp[i,j]
        n_tau_dynamical_interactions : int, optional
               Number of tau pts for D0_tau and jperp_tau (Default 10001)
        n_iw_dynamical_interactions : int, optional
               Number of matsubara freqs for D0_iw and jperp_iw (Default 200)
        """
        constr_params['gf_struct'] = fix_gf_struct_type(constr_params['gf_struct'])

        # Initialise the core solver
        SolverCore.__init__(self, **constr_params)



    def f(self, x,solve_params):
        #setting up needed parameters
        self.eval_count = self.eval_count + 1
        assert 'h_int' in solve_params
        if 'n_s' not in solve_params: solve_params['n_s'] = 1
        if 'delta' not in solve_params: solve_params['delta'] = 0.1
        solve_params.pop('delta')
        n_s = solve_params['n_s']
        h_int = solve_params['h_int']
        gf_struct = self.constr_params['gf_struct']
        n_terms = len(list(h_int))
        km = {}
        for bl, g_bl in self.G0_iw:
            self.G0_iw_inv[bl] << inverse(g_bl)
            km[bl] = make_zero_tail(g_bl, 2)
            km[bl][1] = np.eye(g_bl.target_shape[0])
        _alpha = x.reshape(n_terms, 2, 2, n_s)
        _params = solve_params.copy()
        _params['alpha'] = _alpha
        _params.update(self.constr_params)
        self.prepare_G0_shift_iw(**_params)
        _G0_shift_dens = { bl: g_bl.density(km[bl]).real for bl, g_bl in self.G0_shift_iw }
        # Evaluate the violation of the self-consistency condition
        _res = []
        for n, (term, coeff) in enumerate(h_int):
            bl_0, u_cdag_0 = term[0][1]
            bl_1, u_cdag_1 = term[1][1]
            bl_1, u_c_1 = term[2][1]
            bl_0, u_c_0 = term[3][1]
            idx_cdag_0 = dict(gf_struct)[bl_0].index(u_cdag_0)
            idx_cdag_1 = dict(gf_struct)[bl_1].index(u_cdag_1)
            idx_c_1 = dict(gf_struct)[bl_1].index(u_c_1)
            idx_c_0 = dict(gf_struct)[bl_0].index(u_c_0)
            _res.append(_G0_shift_dens[bl_0][idx_c_0,idx_cdag_0] - _alpha[n,0,0,0])
            if bl_0 == bl_1:
                _res.append(_G0_shift_dens[bl_0][idx_c_1,idx_cdag_0] - _alpha[n,0,1,0])
                _res.append(_G0_shift_dens[bl_0][idx_c_0,idx_cdag_1] - _alpha[n,1,0,0])
            else:
                _res.append(-_alpha[n,0,1,0])
                _res.append(-_alpha[n,1,0,0])
            _res.append(_G0_shift_dens[bl_1][idx_c_1,idx_cdag_1] - _alpha[n,1,1,0])
        return np.array(_res)

    def convert_alpha_to_ij(self,h_int,gf_struct,n,t):
        U_tilde = 0.0
        bl = None
        for m, (term, coeff) in enumerate(h_int):
            if n == m:
                bl_0, u_cdag_0 = term[0][1]
                bl_1, u_cdag_1 = term[1][1]
                bl_1, u_c_1 = term[2][1]
                bl_0, u_c_0 = term[3][1]
                idx_cdag_0 = dict(gf_struct)[bl_0].index(u_cdag_0)
                idx_cdag_1 = dict(gf_struct)[bl_1].index(u_cdag_1)
                idx_c_1 = dict(gf_struct)[bl_1].index(u_c_1)
                idx_c_0 = dict(gf_struct)[bl_0].index(u_c_0)
                if t == 0:
                    bl = bl_1
                    i  = idx_cdag_1
                    j  = idx_c_1
                    U_tilde = -coeff
                elif t == 1:
                    if bl_0 == bl_1:
                        bl = bl_1
                        i  = idx_cdag_1
                        j  = idx_c_0
                        U_tilde = coeff
                    else:
                        return None
                elif t == 2:
                    if bl_0 == bl_1:
                        bl = bl_1
                        i  = idx_cdag_0
                        j  = idx_c_1
                        U_tilde = coeff
                    else:
                        return None
                elif t == 3:
                    bl = bl_0
                    i  = idx_cdag_0
                    j  = idx_c_0
                    U_tilde = -coeff
        return [bl,i,j,U_tilde]
                        
            
    def jacobi(self,x,solve_params):
        #setting up needed parameters
        assert 'h_int' in solve_params
        if 'n_s' not in solve_params: solve_params['n_s'] = 1
        if 'delta' not in solve_params: solve_params['delta'] = 0.1
        solve_params.pop('delta')
        n_s = solve_params['n_s']
        h_int = solve_params['h_int']
        gf_struct = self.constr_params['gf_struct']
        n_terms = len(list(h_int))
        _alpha = x.reshape(n_terms, 2, 2, n_s)            
        _params = solve_params.copy()
        _params['alpha'] = _alpha
        _params.update(self.constr_params)
        km = {}
        for bl, g_bl in self.G0_iw:
            self.G0_iw_inv[bl] << inverse(g_bl)
            km[bl] = make_zero_tail(g_bl, 2)
            km[bl][1] = np.eye(g_bl.target_shape[0])
        #Set trivial, diagonal derivation:
        jac = np.zeros((len(x),len(x)))
        for i in range(len(x)):
            jac[i,i] = -1.0
        # Update the shifted Green function
        self.prepare_G0_shift_iw(**_params)
        bl= gf_struct[0][0]
        G =  self.G0_shift_iw[bl][0,0].copy()
        #create List with needed indices
        L  = []
        for n, (term, coeff) in enumerate(h_int):
            bl_0, u_cdag_0 = term[0][1]
            bl_1, u_cdag_1 = term[1][1]
            bl_1, u_c_1 = term[2][1]
            bl_0, u_c_0 = term[3][1]
            idx_cdag_0 = dict(gf_struct)[bl_0].index(u_cdag_0)
            idx_cdag_1 = dict(gf_struct)[bl_1].index(u_cdag_1)
            idx_c_1 = dict(gf_struct)[bl_1].index(u_c_1)
            idx_c_0 = dict(gf_struct)[bl_0].index(u_c_0)                              
            L.append([bl_0,idx_c_0,idx_cdag_0,n,0])
            if bl_0 == bl_1:
                L.append([bl_0,idx_c_0,idx_cdag_1,n,1])
                L.append([bl_1,idx_c_0,idx_cdag_1,n,2])
            L.append([bl_1,idx_c_1,idx_cdag_1,n,3])
        for K in L:
            i = K[3]*4+K[4]
            sigma = K[0]
            h = K[1]
            l = K[2]
            for n in range(n_terms):
                for t in range(4):
                    c = self.convert_alpha_to_ij(h_int,gf_struct,n,t)
                    if c == None:
                        continue;
                    else:
                        [sigma_2,z1,z2,U_tilde] = c
                        j = 4*n+t
                        G.zero()
                        if sigma == sigma_2:
                            G = -self.G0_shift_iw[sigma_2][h,z1]*self.G0_shift_iw[sigma_2][z2,l]*U_tilde
                    jac[i,j] += np.real(G.density())
        return jac


    
    def solve(self, **solve_params):
        """
        Solve the impurity problem.

        Parameters
        ----------
        solve_params : dict {'param':value} that is passed to the core solver.
                     The only two required parameters are
                        * `h_int`: The local interaction Hamiltonian
                        * `n_cycles`: The number of Monte-Carlo cycles
                     For the other optional parameters see documentation.
                     Note that in this Python Wrapper the alpha-tensor is optional.
                     If not given, it will be constructed automatically.
        """

        assert 'h_int' in solve_params
        assert 'n_cycles' in solve_params
        self.eval_count = 0
        delta = solve_params.pop('delta', 0.1)
        if 'n_s' not in solve_params: solve_params['n_s'] = 1

        h_int = solve_params['h_int']
        gf_struct = self.constr_params['gf_struct']
        n_s = solve_params['n_s']
        
        if 'alpha' not in solve_params:
            # --------- Determine the alpha tensor from SC Hartree Fock ----------
            mpi_print("Determine alpha-tensor")

            if n_s != 1:
                raise Exception("Automatic alpha determination not implemented for n_s != 1")

            def sign(a):
                if a >= 0: return 1
                if a < 0: return -1

            # Prepare the inverse of G0_iw and the known high-frequency moments
            km = {}
            for bl, g_bl in self.G0_iw:
                self.G0_iw_inv[bl] << inverse(g_bl)
                km[bl] = make_zero_tail(g_bl, 2)
                km[bl][1] = np.eye(g_bl.target_shape[0])
                
            if 'use_jacobi' not in solve_params: solve_params['use_jacobi'] = True
            use_jacobi = solve_params['use_jacobi']
            solve_params.pop('use_jacobi')
            # The numer of terms in h_int determines the leading dimension of alpha
            n_terms = len(list(h_int))

            # Calculate initial alpha guess from G0_iw.density(km)
            alpha_init = np.zeros((n_terms, 2, 2, n_s))
            G0_dens = { bl: g_bl.density(km[bl]).real for bl, g_bl in self.G0_iw }
            for n, (term, coeff) in enumerate(h_int):
                bl_0, u_cdag_0 = term[0][1]
                bl_1, u_cdag_1 = term[1][1]
                bl_1, u_c_1 = term[2][1]
                bl_0, u_c_0 = term[3][1]
                idx_cdag_0 = dict(gf_struct)[bl_0].index(u_cdag_0)
                idx_cdag_1 = dict(gf_struct)[bl_1].index(u_cdag_1)
                idx_c_1 = dict(gf_struct)[bl_1].index(u_c_1)
                idx_c_0 = dict(gf_struct)[bl_0].index(u_c_0)     
                alpha_init[n,0,0,0] = G0_dens[bl_0][idx_c_0,idx_cdag_0]
                alpha_init[n,1,1,0] = G0_dens[bl_1][idx_c_1,idx_cdag_1]
                if bl_0 == bl_1:
                    alpha_init[n,0,1,0] = G0_dens[bl_0][idx_c_1,idx_cdag_0]
                    alpha_init[n,1,0,0] = G0_dens[bl_0][idx_c_0,idx_cdag_1]
                else:
                    alpha_init[n,0,1,0] = 0
                    alpha_init[n,1,0,0] = 0

            
            # mpi_print("Init Alpha: " + str(alpha_init[...,0]))
            alpha_vec_init = alpha_init.reshape(4 * n_terms * n_s)
            alpha = alpha_init

            
            if mpi.is_master_node():
                # Find alpha on master_node

             
                self.eval_count = 0
                if use_jacobi == True:
                    mpi_print("Using Jacobi-Matrix for root search")
                    root_finder = root(self.f, alpha_vec_init,args=(solve_params),jac=self.jacobi,method="hybr")
                else:
                    root_finder = root(self.f, alpha_vec_init,args=(solve_params),method="hybr")

                    
                # Reshape result, Implement Fallback solution if unsuccessful
                if root_finder['success']:
                    alpha = root_finder['x'].reshape(n_terms, 2, 2, n_s)
                    solve_params['alpha'] = alpha
                else:
                    mpi_print("Could not determine alpha, falling back to G0_iw.density()")

            alpha = mpi.bcast(alpha, root=0)

            #_ Introduce alpha assymetry
            for n, (term, coeff) in enumerate(h_int):
                alpha[n,0,0,0] = alpha[n,0,0,0] - sign(coeff) * delta
                alpha[n,1,1,0] = alpha[n,1,1,0] + delta
                alpha[n,0,1,0] = alpha[n,0,1,0] + delta * (abs(alpha[n,0,1,0]) > 1e-6)
                alpha[n,1,0,0] = alpha[n,1,0,0] + sign(coeff) * delta * (abs(alpha[n,1,0,0]) > 1e-6)

            mpi_print(" --- Alpha Tensor : ")
            mpi_print(str(alpha[...,0]))
            solve_params['alpha'] = alpha

        solve_status = SolverCore.solve(self, **solve_params)
        return solve_status
