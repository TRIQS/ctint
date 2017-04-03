# Generated automatically using the command :
# c++2py.py ../c++/solver_core.hpp -p --members_read_only -m ctint -o ctint
from wrap_generator import *

# The module
module = module_(full_name = "ctint", doc = "", app_name = "ctint")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')
module.use_module('operators', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/arrays.hpp>
#include <triqs/python_tools/converters/optional.hpp>
#include <triqs/python_tools/converters/variant.hpp>
using namespace triqs::gfs;
using triqs::operators::many_body_operator;
using namespace triqs_ctint;
#include "./ctint_converters.hxx"
""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
        doc = r"The Solver class",   # doc of the C++ class
)

c.add_member(c_name = "average_sign",
             c_type = "double",
             read_only= True,
             doc = """Average sign of the CTINT """)

c.add_member(c_name = "M_tau",
             c_type = "std::optional<block_gf<imtime, matrix_valued> >",
             read_only= True,
             doc = """Building block for the Green function in imaginary time (Eq. (23) in Notes) """)

c.add_member(c_name = "M_iw_nfft",
             c_type = "std::optional<block_gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Same as M_tau, but measured directly in frequencies """)

c.add_member(c_name = "F_tau",
             c_type = "std::optional<block_gf<imtime, matrix_valued> >",
             read_only= True,
             doc = """The improved estimator F_tau """)

c.add_member(c_name = "M_iw",
             c_type = "std::optional<block_gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """The Fourier-transform of M_tau. Dependent on M_tau """)

c.add_member(c_name = "Giw",
             c_type = "std::optional<block_gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_tau """)

c.add_member(c_name = "Sigma_iw",
             c_type = "std::optional<block_gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Self-energy in Matsubara frequencies. Dependent on M_tau """)

c.add_member(c_name = "G0_iw",
             c_type = "block_gf<triqs::gfs::imfreq, triqs::gfs::matrix_valued>",
             read_only= True,
             doc = """Noninteracting Green Function in Matsubara frequencies """)

c.add_member(c_name = "D0_iw",
             c_type = "std::optional<block_gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Dynamic density-density interaction in Matsubara frequencies """)

c.add_member(c_name = "Jperp_iw",
             c_type = "std::optional<gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Dynamic spin-spin interaction in Matsubara frequencies """)

c.add_constructor("""(**triqs_ctint::constr_params_t)""",
                  doc = """+------------------------------+-------------+---------+----------------------------------------------------------------+
| Parameter Name               | Type        | Default | Documentation                                                  |
+==============================+=============+=========+================================================================+
| n_tau                        | int         | 10000   | Number of tau points for gf<imtime, matrix_valued>             |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| n_iw                         | int         | 500     | Number of Matsubara frequencies for gf<imfreq, matrix_valued>  |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| beta                         | double      |         | Inverse temperature                                            |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| gf_struct                    | gf_struct_t |         | block structure of the gf                                      |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| use_D                        | bool        | false   | Switch for dynamic density-density interaction                 |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| use_Jperp                    | bool        | false   | Switch for dynamic spin-spin interaction                       |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| n_tau_dynamical_interactions | int         | 10001   | Number of tau pts for D0_tau and jperp_tau                     |
+------------------------------+-------------+---------+----------------------------------------------------------------+
| n_iw_dynamical_interactions  | int         | 200     | Number of matsubara freqs for D0_iw and jperp_iw               |
+------------------------------+-------------+---------+----------------------------------------------------------------+ """)

c.add_method("""void solve (**triqs_ctint::solve_params_t)""",
             doc = """+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| Parameter Name       | Type                | Default                                        | Documentation                                                                |
+======================+=====================+================================================+==============================================================================+
| hartree_shift        | std::vector<double> | std::vector<double>{}                          | Shift of the chemical potential mu_sigma --> mu_sigma + hartree_shift_sigma  |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| h_int                | many_body_operator  |                                                | Interaction Hamiltonian                                                      |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| use_alpha            | bool                | false                                          | Switch for the use of the alpha function. Compare Sec. 1.3 in Notes.         |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| n_s                  | int                 | 2                                              | Number of auxiliary spins                                                    |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| alpha                | alpha_t             |                                                | Alpha parameter                                                              |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| n_cycles             | int                 |                                                | Number of MC cycles                                                          |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| length_cycle         | int                 | 50                                             | Length of a MC cycles                                                        |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| n_warmup_cycles      | int                 | 5000                                           | Number of warmup cycles                                                      |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| random_seed          | int                 | 34788+928374*triqs::mpi::communicator().rank() | Random seed of the random generator                                          |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| random_name          | std::string         | ""                                             | Name of the random generator                                                 |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| use_double_insertion | bool                | false                                          | Use double insertion                                                         |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| max_time             | int                 | -1                                             | Maximum running time in seconds (-1 : no limit)                              |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| verbosity            | int                 | triqs::mpi::communicator().rank()==0?3:0       | Verbosity                                                                    |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| measure_average_sign | bool                | true                                           | Measure the MC sign                                                          |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| measure_M_tau        | bool                | false                                          | Measure M(tau)                                                               |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| measure_M_iw         | bool                | false                                          | Measure M(iomega) using nfft                                                 |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| measure_F_tau        | bool                | false                                          | Measure F(tau)                                                               |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| nfft_buf_size        | int                 | 10000                                          | Size of the Nfft buffer                                                      |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+
| post_process         | bool                | true                                           | Perform post processing                                                      |
+----------------------+---------------------+------------------------------------------------+------------------------------------------------------------------------------+ """)

module.add_class(c)

module.generate_code()