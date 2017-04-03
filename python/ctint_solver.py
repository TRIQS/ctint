from ctint import SolverCore
import pytriqs.utility.mpi as mpi
import numpy as np

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=100001):
        """
        :param beta: Inverse temperature.
        :param gf_struct: Structure of the Green's functions. It must be a
                          dictionary which maps the name of each block of the
                          Green's function as a string to a list of integer
                          indices.
                          For example: { 'up': [1,2,3], 'down', [1,2,3] }.
        :param n_iw: (optional, default = 1025) Number of Matsubara frequencies
                     used for the Green's functions.
        :param n_tau: (optional, default = 10001) Number of imaginary time points
                     used for the Green's functions.
        """

        # Initialise the core solver
        SolverCore.__init__(self, beta, gf_struct, n_tau, n_tau, n_iw)

    def solve(self, h_int, alpha, n_cycles, length_cycle = 50, n_warmup_cycles = 5000, random_seed = 34788, random_name = "", \
              verbosity = 3, max_time = -1, only_sign = False, measure_gw = True, measure_ft = False, measure_hist = False, correl_to_measure = []):
        """ Solve the impurity problem """

        # Call the core solver's core routine
        SolverCore.solve(self, h_int, alpha, n_cycles, length_cycle, n_warmup_cycles, random_seed, random_name, \
                         verbosity, max_time, only_sign, measure_gw, measure_ft, measure_hist, correl_to_measure)
