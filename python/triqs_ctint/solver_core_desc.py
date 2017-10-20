# Generated automatically using the command :
# c++2py ../../c++/triqs_ctint/solver_core.hpp -p --members_read_only -m solver_core -o solver_core -C pytriqs --cxxflags="-std=c++17 -DHAS_OPTIONAL_HEADER"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = "", app_name = "solver_core")

# Imports
import pytriqs.gf
import pytriqs.operators

# Add here all includes
module.add_include("../triqs_ctint/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/variant.hpp>
#include "solver_core_converters.hxx"
""")
# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_ctint::solver_core",   # name of the C++ class
        doc = """The Solver class""",   # doc of the C++ class
)

c.add_member(c_name = "average_sign",
             c_type = "double",
             read_only= True,
             doc = """Average sign of the CTINT""")

c.add_member(c_name = "average_k",
             c_type = "double",
             read_only= True,
             doc = """Average perturbation order""")

c.add_member(c_name = "M_tau",
             c_type = "std::optional<g_tau_t>",
             read_only= True,
             doc = """Building block for the Green function in imaginary time (Eq. (23) in Notes)""")

c.add_member(c_name = "M_iw_nfft",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """Same as M_tau, but measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "F_tau",
             c_type = "std::optional<g_tau_t>",
             read_only= True,
             doc = """The improved estimator F_tau""")

c.add_member(c_name = "M4_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = """Same as M4_tau, but measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3pp_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (pp channel) in Matsubara frequencies""")

c.add_member(c_name = "M3ph_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (ph channel) in Matsubara frequencies""")

c.add_member(c_name = "M3pp_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (pp channel) in imaginary time""")

c.add_member(c_name = "M3ph_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (ph channel) in imaginary time""")

c.add_member(c_name = "M2pp_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """Building block for the susceptibility (pp channel) in imaginary time""")

c.add_member(c_name = "M2ph_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """Building block for the susceptibility (ph channel) in imaginary time""")

c.add_member(c_name = "M_iw",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """The Fourier-transform of M_tau. Dependent on M_tau""")

c.add_member(c_name = "G_iw",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_iw""")

c.add_member(c_name = "Sigma_iw",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """Self-energy in Matsubara frequencies. Dependent on M_iw""")

c.add_member(c_name = "M3pp_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (pp channel) in Matsubara frequencies""")

c.add_member(c_name = "M3ph_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (ph channel) in Matsubara frequencies""")

c.add_member(c_name = "M2pp_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """Building block for the susceptibility (pp channel) in Matsubara frequencies""")

c.add_member(c_name = "M2ph_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """Building block for the susceptibility (ph channel) in Matsubara frequencies""")

c.add_member(c_name = "F_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = """The two-particle vertex function in purely fermionic notation (iw1, iw2, iw3)""")

c.add_member(c_name = "G2c_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = """The connected part of the two-particle Green function""")

c.add_member(c_name = "G2_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = """The two-particle Green function""")

c.add_member(c_name = "chi2pp_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-particle channel in imaginary times""")

c.add_member(c_name = "chi2ph_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-hole channel in imaginary times""")

c.add_member(c_name = "chi2pp_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-particle channel in Matsubara frequencies""")

c.add_member(c_name = "chi2ph_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-hole channel in Matsubara frequencies""")

c.add_member(c_name = "chi3pp_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_3` in the particle-particle channel in Matsubara frequencies""")

c.add_member(c_name = "chi3ph_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_3` in the particle-hole channel in Matsubara frequencies""")

c.add_member(c_name = "chi3pp_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_3` in the particle-particle channel in Matsubara frequencies as obtained by the NFFT :math:`M_3` measurement""")

c.add_member(c_name = "chi3ph_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_3` in the particle-hole channel in Matsubara frequencies as obtained by the NFFT :math:`M_3` measurement""")

c.add_member(c_name = "G0_iw",
             c_type = "triqs_ctint::g_iw_t",
             read_only= True,
             doc = """Noninteracting Green Function in Matsubara frequencies""")

c.add_member(c_name = "D0_iw",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """Dynamic density-density interaction in Matsubara frequencies""")

c.add_member(c_name = "Jperp_iw",
             c_type = "std::optional<gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Dynamic spin-spin interaction in Matsubara frequencies""")

c.add_constructor("""(**triqs_ctint::constr_params_t)""", doc = """Construct a CTINT solver\n\n :param construct_parameters: Set of parameters specific to the CTINT solver
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| Parameter Name               | Type                              | Default | Documentation                                                  |
+==============================+===================================+=========+================================================================+
| n_tau                        | int                               | 10000   | Number of tau points for gf<imtime, matrix_valued>             |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| n_iw                         | int                               | 500     | Number of Matsubara frequencies for gf<imfreq, matrix_valued>  |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| beta                         | double                            |         | Inverse temperature                                            |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| gf_struct                    | triqs::hilbert_space::gf_struct_t |         | block structure of the gf                                      |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| use_D                        | bool                              | false   | Switch for dynamic density-density interaction                 |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| use_Jperp                    | bool                              | false   | Switch for dynamic spin-spin interaction                       |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| n_tau_dynamical_interactions | int                               | 10001   | Number of tau pts for D0_tau and jperp_tau                     |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| n_iw_dynamical_interactions  | int                               | 200     | Number of matsubara freqs for D0_iw and jperp_iw               |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+""")

c.add_method("""void solve (**triqs_ctint::solve_params_t)""",
             doc = """Solve method that performs CTINT calculation\n\n :param solve_params_t: Set of parameters specific to the CTINT run
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| Parameter Name       | Type                                 | Default                                        | Documentation                                                         |
+======================+======================================+================================================+=======================================================================+
| h_int                | triqs::operators::many_body_operator |                                                | Interaction Hamiltonian                                               |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| use_alpha            | bool                                 | false                                          | Switch for the use of the alpha function. Compare Sec. 1.3 in Notes.  |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_s                  | int                                  | 2                                              | Number of auxiliary spins                                             |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| alpha                | triqs_ctint::alpha_t                 |                                                | Alpha parameter                                                       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_cycles             | int                                  |                                                | Number of MC cycles                                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| length_cycle         | int                                  | 50                                             | Length of a MC cycles                                                 |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_warmup_cycles      | int                                  | 5000                                           | Number of warmup cycles                                               |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| random_seed          | int                                  | 34788+928374*triqs::mpi::communicator().rank() | Random seed of the random generator                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| random_name          | std::string                          | ""                                             | Name of the random generator                                          |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| use_double_insertion | bool                                 | false                                          | Use double insertion                                                  |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| max_time             | int                                  | -1                                             | Maximum running time in seconds (-1 : no limit)                       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| verbosity            | int                                  | triqs::mpi::communicator().rank()==0?3:0       | Verbosity                                                             |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_average_sign | bool                                 | true                                           | Measure the MC sign                                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_average_k    | bool                                 | true                                           | Measure the average perturbation order                                |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M_tau        | bool                                 | false                                          | Measure M(tau)                                                        |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M_iw         | bool                                 | false                                          | Measure M(iomega) using nfft                                          |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_F_tau        | bool                                 | false                                          | Measure F(tau)                                                        |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M4_iw        | bool                                 | false                                          | Measure M4(iw) NFFT                                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_iw_M4              | int                                  | 32                                             | Number of positive Matsubara frequencies in M4                        |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M3pp_iw      | bool                                 | false                                          | Measure M3pp(iw)                                                      |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M3ph_iw      | bool                                 | false                                          | Measure M3ph(iw)                                                      |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_iw_M3              | int                                  | 64                                             | Number of positive Matsubara frequencies in M3                        |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M3pp_tau     | bool                                 | false                                          | Measure M3pp(iw)                                                      |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M3ph_tau     | bool                                 | false                                          | Measure M3ph(iw)                                                      |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_tau_M3             | int                                  | 1000                                           | Number of imaginary time points in M3                                 |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M2pp_tau     | bool                                 | false                                          | Measure M2pp(tau)                                                     |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| measure_M2ph_tau     | bool                                 | false                                          | Measure M2ph(tau)                                                     |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_tau_M2             | int                                  | 10000                                          | Number of imaginary time points in M2                                 |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| n_iw_M2              | int                                  | 128                                            | Number of positive Matsubara frequencies in M2                        |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| nfft_buf_size        | int                                  | 500                                            | Size of the Nfft buffer                                               |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+
| post_process         | bool                                 | true                                           | Perform post processing                                               |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------------------+""")

c.add_property(name = "solve",
               getter = cfunction("void solve ()"),
               doc = """""")

module.add_class(c)

module.generate_code()