# Generated automatically using the command :
# c++2py ../../c++/triqs_ctint/solver_core.hpp --members_read_only -N triqs_ctint -a triqs_ctint -m solver_core -o solver_core -C pytriqs --moduledoc="The TRIQS ctint solver" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = "The TRIQS ctint solver", app_name = "triqs_ctint")

# Imports
module.add_imports(*['pytriqs.gf', 'pytriqs.operators'])

# Add here all includes
module.add_include("triqs_ctint/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>
#include <triqs/cpp2py_converters/h5.hpp>

using namespace triqs_ctint;
""")


# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_ctint::solver_core",   # name of the C++ class
        doc = """The Solver class""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "average_sign",
             c_type = "triqs_ctint::mc_weight_t",
             read_only= True,
             doc = """Average sign of the CTINT""")

c.add_member(c_name = "average_k",
             c_type = "double",
             read_only= True,
             doc = """Average perturbation order""")

c.add_member(c_name = "histogram",
             c_type = "std::optional<std::vector<double> >",
             read_only= True,
             doc = """Average perturbation order distribution""")

c.add_member(c_name = "density",
             c_type = "std::optional<std::vector<matrix<dcomplex> > >",
             read_only= True,
             doc = """The density matrix (measured by operator insertion)""")

c.add_member(c_name = "M_tau",
             c_type = "std::optional<block_gf<imtime, M_tau_target_t> >",
             read_only= True,
             doc = """Building block for the Green function in imaginary time (Eq. (23) in Notes)""")

c.add_member(c_name = "M_hartree",
             c_type = "std::optional<std::vector<matrix<M_tau_scalar_t> > >",
             read_only= True,
             doc = """Hartree-term of M_tau""")

c.add_member(c_name = "M_iw_nfft",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """Same as M_tau, but measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M4_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = """Building block for the full vertex function measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3pp_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (pp channel) in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3ph_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (ph channel) in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3pp_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (pp channel) in imaginary time""")

c.add_member(c_name = "M3ph_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = """Building block for the fermion boson vertex (ph channel) in imaginary time""")

c.add_member(c_name = "M3pp_delta",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """Equal-time peak in M3pp_tau""")

c.add_member(c_name = "M3ph_delta",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """Equal-time peak in M3ph_tau""")

c.add_member(c_name = "chi2pp_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-particle channel in imaginary times as obtained by operator insertion""")

c.add_member(c_name = "chi2ph_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-hole channel in imaginary times as obtained by operator insertion""")

c.add_member(c_name = "chiAB_tau",
             c_type = "std::optional<gf<imtime> >",
             read_only= True,
             doc = """The correlation function :math:`\\chi_AB` in imaginary times""")

c.add_member(c_name = "M_iw",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = """The Fourier-transform of M_tau. Dependent on M_tau""")

c.add_member(c_name = "G_iw",
             c_type = "triqs_ctint::g_iw_t",
             read_only= True,
             doc = """Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_iw""")

c.add_member(c_name = "Sigma_iw",
             c_type = "triqs_ctint::g_iw_t",
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

c.add_member(c_name = "chi2pp_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-particle channel in Matsubara frequencies""")

c.add_member(c_name = "chi2ph_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-hole channel in Matsubara frequencies""")

c.add_member(c_name = "M2pp_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """M2 in the particle-particle channel in imaginary time as obtained from M3""")

c.add_member(c_name = "M2ph_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """M2 in the particle-hole channel in imaginary time as obtained from M3""")

c.add_member(c_name = "chi2pp_new_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-particle channel in imaginary times as obtained from M3pp_tau""")

c.add_member(c_name = "chi2ph_new_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-hole channel in imaginary times as obtained from M3ph_tau""")

c.add_member(c_name = "chi2pp_new_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-particle channel in imaginary frequencies as obtained from M3pp_tau""")

c.add_member(c_name = "chi2ph_new_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = """The equal time correlator :math:`\\chi_2` in the particle-hole channel in imaginary frequencies as obtained from M3ph_tau""")

c.add_member(c_name = "chiAB_iw",
             c_type = "std::optional<gf<imfreq> >",
             read_only= True,
             doc = """The correlation function :math:`\\chi_AB` in imaginary frequencies""")

c.add_member(c_name = "M3pp_tau_conn",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = """The connected part of M3pp_tau""")

c.add_member(c_name = "M3ph_tau_conn",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = """The connected part of M3ph_tau""")

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
             c_type = "std::optional<block2_gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Dynamic density-density interaction in Matsubara frequencies""")

c.add_member(c_name = "Jperp_iw",
             c_type = "std::optional<gf<imfreq, matrix_valued> >",
             read_only= True,
             doc = """Dynamic spin-spin interaction in Matsubara frequencies""")

c.add_member(c_name = "G0_shift_iw",
             c_type = "triqs_ctint::g_iw_t",
             read_only= True,
             doc = """The shifted noninteracting Green Function in Matsubara frequencies""")

c.add_member(c_name = "G0_shift_tau",
             c_type = "triqs_ctint::g_tau_t",
             read_only= True,
             doc = """The shifted noninteracting Green Function in imaginary time""")

c.add_member(c_name = "constr_params",
             c_type = "triqs_ctint::constr_params_t",
             read_only= True,
             doc = """""")

c.add_member(c_name = "last_solve_params",
             c_type = "std::optional<solve_params_t>",
             read_only= True,
             doc = """""")

c.add_constructor("""(**triqs_ctint::constr_params_t)""", doc = """Construct a CTINT solver\n\n :param construct_parameters: Set of parameters specific to the CTINT solver
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| Parameter Name               | Type                              | Default | Documentation                                                  |
+==============================+===================================+=========+================================================================+
| n_tau                        | int                               | 5001    | Number of tau points for gf<imtime, matrix_valued>             |
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
| n_tau_dynamical_interactions | int                               | n_tau   | Number of tau pts for D0_tau and jperp_tau                     |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
| n_iw_dynamical_interactions  | int                               | n_iw    | Number of matsubara freqs for D0_iw and jperp_iw               |
+------------------------------+-----------------------------------+---------+----------------------------------------------------------------+
""")

c.add_method("""void solve (**triqs_ctint::solve_params_t)""",
             doc = """Solve method that performs CTINT calculation\n\n :param solve_params_t: Set of parameters specific to the CTINT run
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| Parameter Name       | Type                                 | Default                                        | Documentation                                             |
+======================+======================================+================================================+===========================================================+
| h_int                | triqs::operators::many_body_operator |                                                | Interaction Hamiltonian                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_s                  | int                                  | 1                                              | Number of auxiliary spins                                 |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| alpha                | triqs_ctint::alpha_t                 |                                                | Alpha parameter                                           |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_cycles             | int                                  |                                                | Number of MC cycles                                       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| length_cycle         | int                                  | 50                                             | Length of a MC cycles                                     |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_warmup_cycles      | int                                  | 5000                                           | Number of warmup cycles                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| random_seed          | int                                  | 34788+928374*triqs::mpi::communicator().rank() | Random seed of the random generator                       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| random_name          | std::string                          | ""                                             | Name of the random generator                              |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| use_double_insertion | bool                                 | false                                          | Use double insertion                                      |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| max_time             | int                                  | -1                                             | Maximum running time in seconds (-1 : no limit)           |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| verbosity            | int                                  | triqs::mpi::communicator().rank()==0?3:0       | Verbosity                                                 |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_average_sign | bool                                 | true                                           | Measure the MC sign                                       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_average_k    | bool                                 | true                                           | Measure the average perturbation order                    |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_histogram    | bool                                 | false                                          | Measure the average perturbation order distribution       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_density      | bool                                 | false                                          | Measure the density matrix by operator insertion          |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M_tau        | bool                                 | true                                           | Measure M(tau)                                            |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M_iw         | bool                                 | false                                          | Measure M(iomega) using nfft                              |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M4_iw        | bool                                 | false                                          | Measure M4(iw) NFFT                                       |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_iw_M4              | int                                  | 32                                             | Number of positive Matsubara frequencies in M4            |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M3pp_iw      | bool                                 | false                                          | Measure M3pp(iw)                                          |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M3ph_iw      | bool                                 | false                                          | Measure M3ph(iw)                                          |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_iw_M3              | int                                  | 64                                             | Number of positive fermionic Matsubara frequencies in M3  |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_iW_M3              | int                                  | 32                                             | Number of positive bosonic Matsubara frequencies in M3    |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M3pp_tau     | bool                                 | false                                          | Measure M3pp(tau)                                         |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_M3ph_tau     | bool                                 | false                                          | Measure M3ph(tau)                                         |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_tau_M3             | int                                  | 1000                                           | Number of imaginary time points in M3                     |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_chi2pp_tau   | bool                                 | false                                          | Measure of chi2pp by insertion                            |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_chi2ph_tau   | bool                                 | false                                          | Measure of chi2ph by insertion                            |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_tau_chi2           | int                                  | 10000                                          | Number of imaginary time points in chi2                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| n_iw_chi2            | int                                  | 128                                            | Number of positive Matsubara frequencies in chi2          |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| measure_chiAB_tau    | bool                                 | false                                          | Measure of chiAB by insertion                             |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| chi_A_vec            | std::vector<many_body_operator>      | {}                                             | The list of all operators A                               |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| chi_B_vec            | std::vector<many_body_operator>      | {}                                             | The list of all operators B                               |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| nfft_buf_size        | int                                  | 500                                            | Size of the Nfft buffer                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
| post_process         | bool                                 | true                                           | Perform post processing                                   |
+----------------------+--------------------------------------+------------------------------------------------+-----------------------------------------------------------+
""")

c.add_method("""void post_process ()""",
             doc = """""")

c.add_method("""std::string hdf5_scheme ()""",
             is_static = True,
             doc = """""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "triqs_ctint::solve_params_t",
        doc = """The parameters for the solve function""",
)
c.add_member(c_name = "h_int",
             c_type = "triqs::operators::many_body_operator",
             initializer = """  """,
             doc = """Interaction Hamiltonian""")

c.add_member(c_name = "n_s",
             c_type = "int",
             initializer = """ 1 """,
             doc = """Number of auxiliary spins""")

c.add_member(c_name = "alpha",
             c_type = "triqs_ctint::alpha_t",
             initializer = """  """,
             doc = """Alpha parameter""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = """Number of MC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 50 """,
             doc = """Length of a MC cycles""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = """Number of warmup cycles""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*triqs::mpi::communicator().rank() """,
             doc = """Random seed of the random generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = """Name of the random generator""")

c.add_member(c_name = "use_double_insertion",
             c_type = "bool",
             initializer = """ false """,
             doc = """Use double insertion""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = """Maximum running time in seconds (-1 : no limit)""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ triqs::mpi::communicator().rank()==0?3:0 """,
             doc = """Verbosity""")

c.add_member(c_name = "measure_average_sign",
             c_type = "bool",
             initializer = """ true """,
             doc = """Measure the MC sign""")

c.add_member(c_name = "measure_average_k",
             c_type = "bool",
             initializer = """ true """,
             doc = """Measure the average perturbation order""")

c.add_member(c_name = "measure_histogram",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure the average perturbation order distribution""")

c.add_member(c_name = "measure_density",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure the density matrix by operator insertion""")

c.add_member(c_name = "measure_M_tau",
             c_type = "bool",
             initializer = """ true """,
             doc = """Measure M(tau)""")

c.add_member(c_name = "measure_M_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure M(iomega) using nfft""")

c.add_member(c_name = "measure_M4_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure M4(iw) NFFT""")

c.add_member(c_name = "n_iw_M4",
             c_type = "int",
             initializer = """ 32 """,
             doc = """Number of positive Matsubara frequencies in M4""")

c.add_member(c_name = "measure_M3pp_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure M3pp(iw)""")

c.add_member(c_name = "measure_M3ph_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure M3ph(iw)""")

c.add_member(c_name = "n_iw_M3",
             c_type = "int",
             initializer = """ 64 """,
             doc = """Number of positive fermionic Matsubara frequencies in M3""")

c.add_member(c_name = "n_iW_M3",
             c_type = "int",
             initializer = """ 32 """,
             doc = """Number of positive bosonic Matsubara frequencies in M3""")

c.add_member(c_name = "measure_M3pp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure M3pp(tau)""")

c.add_member(c_name = "measure_M3ph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure M3ph(tau)""")

c.add_member(c_name = "n_tau_M3",
             c_type = "int",
             initializer = """ 1000 """,
             doc = """Number of imaginary time points in M3""")

c.add_member(c_name = "measure_chi2pp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure of chi2pp by insertion""")

c.add_member(c_name = "measure_chi2ph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure of chi2ph by insertion""")

c.add_member(c_name = "n_tau_chi2",
             c_type = "int",
             initializer = """ 10000 """,
             doc = """Number of imaginary time points in chi2""")

c.add_member(c_name = "n_iw_chi2",
             c_type = "int",
             initializer = """ 128 """,
             doc = """Number of positive Matsubara frequencies in chi2""")

c.add_member(c_name = "measure_chiAB_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure of chiAB by insertion""")

c.add_member(c_name = "chi_A_vec",
             c_type = "std::vector<many_body_operator>",
             initializer = """ {} """,
             doc = """The list of all operators A""")

c.add_member(c_name = "chi_B_vec",
             c_type = "std::vector<many_body_operator>",
             initializer = """ {} """,
             doc = """The list of all operators B""")

c.add_member(c_name = "nfft_buf_size",
             c_type = "int",
             initializer = """ 500 """,
             doc = """Size of the Nfft buffer""")

c.add_member(c_name = "post_process",
             c_type = "bool",
             initializer = """ true """,
             doc = """Perform post processing""")

module.add_converter(c)

# Converter for constr_params_t
c = converter_(
        c_type = "triqs_ctint::constr_params_t",
        doc = """The parameters for the solver construction""",
)
c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 5001 """,
             doc = """Number of tau points for gf<imtime, matrix_valued>""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 500 """,
             doc = """Number of Matsubara frequencies for gf<imfreq, matrix_valued>""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = """Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """  """,
             doc = """block structure of the gf""")

c.add_member(c_name = "use_D",
             c_type = "bool",
             initializer = """ false """,
             doc = """Switch for dynamic density-density interaction""")

c.add_member(c_name = "use_Jperp",
             c_type = "bool",
             initializer = """ false """,
             doc = """Switch for dynamic spin-spin interaction""")

c.add_member(c_name = "n_tau_dynamical_interactions",
             c_type = "int",
             initializer = """ res.n_tau """,
             doc = """Number of tau pts for D0_tau and jperp_tau""")

c.add_member(c_name = "n_iw_dynamical_interactions",
             c_type = "int",
             initializer = """ res.n_iw """,
             doc = """Number of matsubara freqs for D0_iw and jperp_iw""")

module.add_converter(c)


module.generate_code()