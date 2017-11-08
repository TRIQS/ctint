# Generated automatically using the command :
# c++2py ../../c++/triqs_ctint/solver_core.hpp --members_read_only -N triqs_ctint -a triqs_ctint -m solver_core -o solver_core -C triqs -C nda_py --moduledoc="The TRIQS ctint solver" --includes=../../c++ --cxxflags="-std=c++20 $(triqs++ -cxxflags)" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The TRIQS ctint solver", app_name = "triqs_ctint")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes', 'triqs.operators', 'h5._h5py'])

# Add here all includes
module.add_include("triqs_ctint/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/std_array.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace triqs_ctint;
""")


# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_ctint::solver_core",   # name of the C++ class
        doc = r"""The Solver class""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "average_sign",
             c_type = "mc_weight_t",
             read_only= True,
             doc = r"""Average sign of the CTINT""")

c.add_member(c_name = "nmeasures",
             c_type = "uint64_t",
             read_only= True,
             doc = r"""Total number of measures""")

c.add_member(c_name = "average_k",
             c_type = "double",
             read_only= True,
             doc = r"""Average perturbation order""")

c.add_member(c_name = "auto_corr_time",
             c_type = "double",
             read_only= True,
             doc = r"""Auto-correlation time""")

c.add_member(c_name = "warmup_time",
             c_type = "double",
             read_only= True,
             doc = r"""Warmup time in seconds""")

c.add_member(c_name = "accumulation_time",
             c_type = "double",
             read_only= True,
             doc = r"""Accumulation time in seconds""")

c.add_member(c_name = "histogram",
             c_type = "std::optional<std::vector<double>>",
             read_only= True,
             doc = r"""Average perturbation order distribution""")

c.add_member(c_name = "density",
             c_type = "std::optional<block_matrix_t>",
             read_only= True,
             doc = r"""The density matrix (measured by operator insertion)""")

c.add_member(c_name = "M_tau",
             c_type = "std::optional<block_gf<imtime, M_tau_target_t>>",
             read_only= True,
             doc = r"""Building block for the Green function in imaginary time (Eq. (23) in Notes)""")

c.add_member(c_name = "M_hartree",
             c_type = "std::optional<block_matrix_t>",
             read_only= True,
             doc = r"""Hartree-term of M_tau""")

c.add_member(c_name = "M_iw_nfft",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = r"""Same as M_tau, but measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M4_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""Building block for the full vertex function measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M4pp_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""Building block for the full vertex function (pp channel) measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M4ph_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""Building block for the full vertex function (ph channel) measured directly in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3pp_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (pp channel) in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3ph_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (ph channel) in Matsubara frequencies using NFFT""")

c.add_member(c_name = "M3pp_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (pp channel) in imaginary time""")

c.add_member(c_name = "M3ph_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (ph channel) in imaginary time""")

c.add_member(c_name = "M3xph_tau",
             c_type = "std::optional<chi3_tau_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (xph channel) in imaginary time""")

c.add_member(c_name = "M3pp_delta",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""Equal-time peak in M3pp_tau""")

c.add_member(c_name = "M3ph_delta",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""Equal-time peak in M3ph_tau""")

c.add_member(c_name = "M3xph_delta",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""Equal-time peak in M3ph_tau""")

c.add_member(c_name = "chi2pp_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-particle channel in imaginary times as obtained by operator insertion""")

c.add_member(c_name = "chi2ph_tau",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-hole channel in imaginary times as obtained by operator insertion""")

c.add_member(c_name = "chiAB_tau",
             c_type = "std::optional<gf<imtime>>",
             read_only= True,
             doc = r"""The correlation function $\chi_AB$ in imaginary times""")

c.add_member(c_name = "M_iw",
             c_type = "std::optional<g_iw_t>",
             read_only= True,
             doc = r"""The Fourier-transform of M_tau. Dependent on M_tau""")

c.add_member(c_name = "G_iw",
             c_type = "g_iw_t",
             read_only= True,
             doc = r"""Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_iw""")

c.add_member(c_name = "Sigma_iw",
             c_type = "g_iw_t",
             read_only= True,
             doc = r"""Self-energy in Matsubara frequencies. Dependent on M_iw""")

c.add_member(c_name = "M3pp_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (pp channel) in Matsubara frequencies""")

c.add_member(c_name = "M3ph_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (ph channel) in Matsubara frequencies""")

c.add_member(c_name = "M3xph_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""Building block for the fermion boson vertex (xph channel) in Matsubara frequencies""")

c.add_member(c_name = "F_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The two-particle vertex function in purely fermionic notation (iw1, iw2, iw3)""")

c.add_member(c_name = "Fpp_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The two-particle vertex function (pp channel)""")

c.add_member(c_name = "Fph_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The two-particle vertex function (ph channel)""")

c.add_member(c_name = "G2_conn_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The connected part of the two-particle Green function""")

c.add_member(c_name = "G2pp_conn_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The connected part of the two-particle Green function (pp channel)""")

c.add_member(c_name = "G2ph_conn_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The connected part of the two-particle Green function (ph channel)""")

c.add_member(c_name = "G2_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The two-particle Green function""")

c.add_member(c_name = "G2pp_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The two-particle Green function (pp channel)""")

c.add_member(c_name = "G2ph_iw",
             c_type = "std::optional<chi4_iw_t>",
             read_only= True,
             doc = r"""The two-particle Green function (ph channel)""")

c.add_member(c_name = "chi2pp_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-particle channel in Matsubara frequencies""")

c.add_member(c_name = "chi2ph_iw",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-hole channel in Matsubara frequencies""")

c.add_member(c_name = "chi2pp_conn_tau_from_M3",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""M2 in the particle-particle channel in imaginary time as obtained from M3""")

c.add_member(c_name = "chi2ph_conn_tau_from_M3",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""M2 in the particle-hole channel in imaginary time as obtained from M3""")

c.add_member(c_name = "chi2xph_conn_tau_from_M3",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""M2 in the particle-hole-cross channel in imaginary time as obtained from M3""")

c.add_member(c_name = "chi2pp_tau_from_M3",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-particle channel in imaginary times as obtained from M3pp_tau""")

c.add_member(c_name = "chi2ph_tau_from_M3",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-hole channel in imaginary times as obtained from M3ph_tau""")

c.add_member(c_name = "chi2xph_tau_from_M3",
             c_type = "std::optional<chi2_tau_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-hole-cross channel in imaginary times as obtained from M3ph_tau""")

c.add_member(c_name = "chi2pp_iw_from_M3",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-particle channel in imaginary frequencies as obtained from M3pp_tau""")

c.add_member(c_name = "chi2ph_iw_from_M3",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-hole channel in imaginary frequencies as obtained from M3ph_tau""")

c.add_member(c_name = "chi2xph_iw_from_M3",
             c_type = "std::optional<chi2_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_2$ in the particle-hole-cross channel in imaginary frequencies as obtained from M3ph_tau""")

c.add_member(c_name = "chiAB_iw",
             c_type = "std::optional<gf<imfreq>>",
             read_only= True,
             doc = r"""The correlation function $\chi_AB$ in imaginary frequencies""")

c.add_member(c_name = "chi3pp_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_3$ in the particle-particle channel in Matsubara frequencies""")

c.add_member(c_name = "chi3ph_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_3$ in the particle-hole channel in Matsubara frequencies""")

c.add_member(c_name = "chi3xph_iw",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_3$ in the particle-hole-cross channel in Matsubara frequencies""")

c.add_member(c_name = "chi3pp_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_3$ in the particle-particle channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement""")

c.add_member(c_name = "chi3ph_iw_nfft",
             c_type = "std::optional<chi3_iw_t>",
             read_only= True,
             doc = r"""The equal time correlator $\chi_3$ in the particle-hole channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement""")

c.add_member(c_name = "G0_iw",
             c_type = "g_iw_t",
             read_only= True,
             doc = r"""Noninteracting Green Function in Matsubara frequencies""")

c.add_member(c_name = "G0_iw_inv",
             c_type = "g_iw_t",
             read_only= True,
             doc = r"""The inverse of the noninteracting Green Function""")

c.add_member(c_name = "D0_iw",
             c_type = "std::optional<block2_gf<imfreq, matrix_valued>>",
             read_only= True,
             doc = r"""Dynamic density-density interaction in Matsubara frequencies""")

c.add_member(c_name = "Jperp_iw",
             c_type = "std::optional<gf<imfreq, matrix_valued>>",
             read_only= True,
             doc = r"""Dynamic spin-spin interaction in Matsubara frequencies""")

c.add_member(c_name = "G0_shift_iw",
             c_type = "g_iw_t",
             read_only= True,
             doc = r"""The shifted noninteracting Green Function in Matsubara frequencies""")

c.add_member(c_name = "G0_shift_tau",
             c_type = "g_tau_t",
             read_only= True,
             doc = r"""The shifted noninteracting Green Function in imaginary time""")

c.add_member(c_name = "constr_params",
             c_type = "constr_params_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "last_solve_params",
             c_type = "std::optional<solve_params_t>",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(**constr_params_t)""", doc = r"""Construct a CTINT solver



+------------------------------+-------------+-----------+----------------------------------------------------------------+
| Parameter Name               | Type        | Default   | Documentation                                                  |
+==============================+=============+===========+================================================================+
| n_tau                        | int         | 5001      | Number of tau points for gf<imtime, matrix_valued>             |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| n_iw                         | int         | 500       | Number of Matsubara frequencies for gf<imfreq, matrix_valued>  |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| beta                         | double      | --        | Inverse temperature                                            |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| gf_struct                    | gf_struct_t | --        | block structure of the gf                                      |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| use_D                        | bool        | false     | Switch for dynamic density-density interaction                 |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| use_Jperp                    | bool        | false     | Switch for dynamic spin-spin interaction                       |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| n_tau_dynamical_interactions | int         | n_tau     | Number of tau pts for D0_tau and jperp_tau                     |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
| n_iw_dynamical_interactions  | int         | n_iw      | Number of matsubara freqs for D0_iw and jperp_iw               |
+------------------------------+-------------+-----------+----------------------------------------------------------------+
""")

c.add_method("""void solve (**solve_params_t)""",
             doc = r"""Solve method that performs CTINT calculation



+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                            | Default                                 | Documentation                                                                                                                         |
+===============================+=================================+=========================================+=======================================================================================================================================+
| h_int                         | many_body_operator              | --                                      | Interaction Hamiltonian                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_s                           | int                             | 2                                       | Number of auxiliary spins                                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| alpha                         | alpha_t                         | --                                      | Alpha tensor                                                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                             | --                                      | Number of MC cycles                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                             | 50                                      | Length of a MC cycles                                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                             | 5000                                    | Number of warmup cycles                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| random_seed                   | int                             | 34788+928374*mpi::communicator().rank() | Random seed of the random generator                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| random_name                   | std::string                     | ""                                      | Name of the random generator                                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| use_double_insertion          | bool                            | false                                   | Use double insertion                                                                                                                  |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| max_time                      | int                             | -1                                      | Maximum running time in seconds (-1 : no limit)                                                                                       |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| max_order                     | int                             | -1                                      | Maximum pertubation order which is accepted in move::insert/remove (<0 : unlimited)                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| verbosity                     | int                             | mpi::communicator().rank()==0?3:0       | Verbosity                                                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_sign_only             | bool                            | false                                   | Measure Sign only mode                                                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_average_sign          | bool                            | true                                    | Measure the MC sign                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_average_k             | bool                            | true                                    | Measure the average perturbation order                                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_auto_corr_time        | bool                            | true                                    | Measure the auto-correlation time                                                                                                     |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_histogram             | bool                            | false                                   | Measure the average perturbation order distribution                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_density               | bool                            | false                                   | Measure the density matrix by operator insertion                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M_tau                 | bool                            | true                                    | Measure M(tau)                                                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M_iw                  | bool                            | false                                   | Measure M(iomega) using nfft                                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M4_iw                 | bool                            | false                                   | Measure M4(iw) NFFT                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M4pp_iw               | bool                            | false                                   | Measure M4pp(iw) NFFT                                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M4ph_iw               | bool                            | false                                   | Measure M4ph(iw) NFFT                                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iW_M4                       | int                             | 32                                      | Number of positive bosonic Matsubara frequencies in M4                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_M4                       | int                             | 32                                      | Number of positive fermionic Matsubara frequencies in M4                                                                              |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3pp_iw               | bool                            | false                                   | Measure M3pp(iw)                                                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3ph_iw               | bool                            | false                                   | Measure M3ph(iw)                                                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_M3                       | int                             | 64                                      | Number of positive fermionic Matsubara frequencies in M3                                                                              |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iW_M3                       | int                             | 32                                      | Number of positive bosonic Matsubara frequencies in M3                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3pp_tau              | bool                            | false                                   | Measure M3pp(tau)                                                                                                                     |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3ph_tau              | bool                            | false                                   | Measure M3ph(tau)                                                                                                                     |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3xph_tau             | bool                            | false                                   | Measure M3xph(tau)                                                                                                                    |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_tau_M3                      | int                             | 201                                     | Number of imaginary time points in M3                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_chi2pp_tau            | bool                            | false                                   | Measure of chi2pp by insertion                                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_chi2ph_tau            | bool                            | false                                   | Measure of chi2ph by insertion                                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_tau_chi2                    | int                             | 201                                     | Number of imaginary time points in chi2                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_chi2                     | int                             | 32                                      | Number of positive Matsubara frequencies in chi2                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_chiAB_tau             | bool                            | false                                   | Measure of chiAB by insertion                                                                                                         |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| chi_A_vec                     | std::vector<many_body_operator> | {}                                      | The list of all operators A                                                                                                           |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| chi_B_vec                     | std::vector<many_body_operator> | {}                                      | The list of all operators B                                                                                                           |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| nfft_buf_size                 | int                             | 500                                     | Size of the Nfft buffer                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| post_process                  | bool                            | true                                    | Perform post processing                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_init_size                 | int                             | 1000                                    | The maximum size of the determinant matrix before a resize                                                                            |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_n_operations_before_check | int                             | 100                                     | Max number of ops before the test of deviation of the det, M^-1 is performed.                                                         |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_precision_warning         | double                          | 1.e-8                                   | Threshold for determinant precision warnings                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_precision_error           | double                          | 1.e-5                                   | Threshold for determinant precision error                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_singular_threshold        | double                          | -1                                      | Bound for the determinant matrix being singular: abs(det) < singular_threshold. For negative threshold check if !isnormal(abs(det)).  |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
""")

c.add_method("""void prepare_G0_shift_iw (**params_t)""",
             doc = r"""



+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                            | Default                                 | Documentation                                                                                                                         |
+===============================+=================================+=========================================+=======================================================================================================================================+
| n_tau                         | int                             | 5001                                    | Number of tau points for gf<imtime, matrix_valued>                                                                                    |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw                          | int                             | 500                                     | Number of Matsubara frequencies for gf<imfreq, matrix_valued>                                                                         |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| beta                          | double                          | --                                      | Inverse temperature                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| gf_struct                     | gf_struct_t                     | --                                      | block structure of the gf                                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| use_D                         | bool                            | false                                   | Switch for dynamic density-density interaction                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| use_Jperp                     | bool                            | false                                   | Switch for dynamic spin-spin interaction                                                                                              |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_tau_dynamical_interactions  | int                             | n_tau                                   | Number of tau pts for D0_tau and jperp_tau                                                                                            |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_dynamical_interactions   | int                             | n_iw                                    | Number of matsubara freqs for D0_iw and jperp_iw                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| h_int                         | many_body_operator              | --                                      | Interaction Hamiltonian                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_s                           | int                             | 1                                       | Number of auxiliary spins                                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| alpha                         | alpha_t                         | --                                      | Alpha tensor                                                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                             | --                                      | Number of MC cycles                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                             | 100                                     | Length of a MC cycles                                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                             | 5000                                    | Number of warmup cycles                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| random_seed                   | int                             | 34788+928374*mpi::communicator().rank() | Random seed of the random generator                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| random_name                   | std::string                     | ""                                      | Name of the random generator                                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| use_double_insertion          | bool                            | true                                    | Use double insertion                                                                                                                  |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| insertion_types               | std::vector<int>                | {}                                      | Types of insertions to use                                                                                                            |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| max_time                      | int                             | -1                                      | Maximum running time in seconds (-1 : no limit)                                                                                       |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| max_order                     | int                             | -1                                      | Maximum pertubation order which is accepted in move::insert/remove (<0 : unlimited)                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| verbosity                     | int                             | mpi::communicator().rank()==0?3:0       | Verbosity                                                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_sign_only             | bool                            | false                                   | Measure Sign only mode                                                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_average_sign          | bool                            | true                                    | Measure the MC sign                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_average_k             | bool                            | true                                    | Measure the average perturbation order                                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_auto_corr_time        | bool                            | true                                    | Measure the auto-correlation time                                                                                                     |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_histogram             | bool                            | false                                   | Measure the average perturbation order distribution                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_density               | bool                            | false                                   | Measure the density matrix by operator insertion                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M_tau                 | bool                            | true                                    | Measure M(tau)                                                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M_iw                  | bool                            | false                                   | Measure M(iomega) using nfft                                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M4_iw                 | bool                            | false                                   | Measure M4(iw) NFFT                                                                                                                   |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M4pp_iw               | bool                            | false                                   | Measure M4pp(iw) NFFT                                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M4ph_iw               | bool                            | false                                   | Measure M4ph(iw) NFFT                                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iW_M4                       | int                             | 32                                      | Number of positive bosonic Matsubara frequencies in M4                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_M4                       | int                             | 32                                      | Number of positive fermionic Matsubara frequencies in M4                                                                              |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3pp_iw               | bool                            | false                                   | Measure M3pp(iw)                                                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3ph_iw               | bool                            | false                                   | Measure M3ph(iw)                                                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_M3                       | int                             | 64                                      | Number of positive fermionic Matsubara frequencies in M3                                                                              |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iW_M3                       | int                             | 32                                      | Number of positive bosonic Matsubara frequencies in M3                                                                                |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3pp_tau              | bool                            | false                                   | Measure M3pp(tau)                                                                                                                     |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3ph_tau              | bool                            | false                                   | Measure M3ph(tau)                                                                                                                     |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_M3xph_tau             | bool                            | false                                   | Measure M3xph(tau)                                                                                                                    |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_tau_M3                      | int                             | 201                                     | Number of imaginary time points in M3                                                                                                 |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_chi2pp_tau            | bool                            | false                                   | Measure of chi2pp by insertion                                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_chi2ph_tau            | bool                            | false                                   | Measure of chi2ph by insertion                                                                                                        |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_tau_chi2                    | int                             | 201                                     | Number of imaginary time points in chi2                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| n_iw_chi2                     | int                             | 32                                      | Number of positive Matsubara frequencies in chi2                                                                                      |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| measure_chiAB_tau             | bool                            | false                                   | Measure of chiAB by insertion                                                                                                         |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| chi_A_vec                     | std::vector<many_body_operator> | {}                                      | The list of all operators A                                                                                                           |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| chi_B_vec                     | std::vector<many_body_operator> | {}                                      | The list of all operators B                                                                                                           |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| nfft_buf_size                 | int                             | 100000                                  | Size of the Nfft buffer                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| post_process                  | bool                            | true                                    | Perform post processing                                                                                                               |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_init_size                 | int                             | 1000                                    | The maximum size of the determinant matrix before a resize                                                                            |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_n_operations_before_check | int                             | 100                                     | Max number of ops before the test of deviation of the det, M^-1 is performed.                                                         |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_precision_warning         | double                          | 1.e-8                                   | Threshold for determinant precision warnings                                                                                          |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_precision_error           | double                          | 1.e-5                                   | Threshold for determinant precision error                                                                                             |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| det_singular_threshold        | double                          | -1                                      | Bound for the determinant matrix being singular: abs(det) < singular_threshold. For negative threshold check if !isnormal(abs(det)).  |
+-------------------------------+---------------------------------+-----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
""")

c.add_method("""void post_process ()""",
             doc = r"""""")

c.add_method("""std::string hdf5_format ()""",
             is_static = True,
             doc = r"""""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "triqs_ctint::solve_params_t",
        doc = r"""The parameters for the solve function""",
)
c.add_member(c_name = "h_int",
             c_type = "many_body_operator",
             initializer = """  """,
             doc = r"""Interaction Hamiltonian""")

c.add_member(c_name = "n_s",
             c_type = "int",
             initializer = """ 2 """,
             doc = r"""Number of auxiliary spins""")

c.add_member(c_name = "alpha",
             c_type = "alpha_t",
             initializer = """  """,
             doc = r"""Alpha tensor""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = r"""Number of MC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 50 """,
             doc = r"""Length of a MC cycles""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = r"""Number of warmup cycles""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*mpi::communicator().rank() """,
             doc = r"""Random seed of the random generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Name of the random generator""")

c.add_member(c_name = "use_double_insertion",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Use double insertion""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum running time in seconds (-1 : no limit)""")

c.add_member(c_name = "max_order",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum pertubation order which is accepted in move::insert/remove
     (<0 : unlimited)""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ mpi::communicator().rank()==0?3:0 """,
             doc = r"""Verbosity""")

c.add_member(c_name = "measure_sign_only",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure Sign only mode""")

c.add_member(c_name = "measure_average_sign",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the MC sign""")

c.add_member(c_name = "measure_average_k",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the average perturbation order""")

c.add_member(c_name = "measure_auto_corr_time",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the auto-correlation time""")

c.add_member(c_name = "measure_histogram",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the average perturbation order distribution""")

c.add_member(c_name = "measure_density",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the density matrix by operator insertion""")

c.add_member(c_name = "measure_M_tau",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure M(tau)""")

c.add_member(c_name = "measure_M_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M(iomega) using nfft""")

c.add_member(c_name = "measure_M4_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M4(iw) NFFT""")

c.add_member(c_name = "measure_M4pp_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M4pp(iw) NFFT""")

c.add_member(c_name = "measure_M4ph_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M4ph(iw) NFFT""")

c.add_member(c_name = "n_iW_M4",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive bosonic Matsubara frequencies in M4""")

c.add_member(c_name = "n_iw_M4",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive fermionic Matsubara frequencies in M4""")

c.add_member(c_name = "measure_M3pp_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3pp(iw)""")

c.add_member(c_name = "measure_M3ph_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3ph(iw)""")

c.add_member(c_name = "n_iw_M3",
             c_type = "int",
             initializer = """ 64 """,
             doc = r"""Number of positive fermionic Matsubara frequencies in M3""")

c.add_member(c_name = "n_iW_M3",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive bosonic Matsubara frequencies in M3""")

c.add_member(c_name = "measure_M3pp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3pp(tau)""")

c.add_member(c_name = "measure_M3ph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3ph(tau)""")

c.add_member(c_name = "measure_M3xph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3xph(tau)""")

c.add_member(c_name = "n_tau_M3",
             c_type = "int",
             initializer = """ 201 """,
             doc = r"""Number of imaginary time points in M3""")

c.add_member(c_name = "measure_chi2pp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure of chi2pp by insertion""")

c.add_member(c_name = "measure_chi2ph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure of chi2ph by insertion""")

c.add_member(c_name = "n_tau_chi2",
             c_type = "int",
             initializer = """ 201 """,
             doc = r"""Number of imaginary time points in chi2""")

c.add_member(c_name = "n_iw_chi2",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive Matsubara frequencies in chi2""")

c.add_member(c_name = "measure_chiAB_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure of chiAB by insertion""")

c.add_member(c_name = "chi_A_vec",
             c_type = "std::vector<many_body_operator>",
             initializer = """ {} """,
             doc = r"""The list of all operators A""")

c.add_member(c_name = "chi_B_vec",
             c_type = "std::vector<many_body_operator>",
             initializer = """ {} """,
             doc = r"""The list of all operators B""")

c.add_member(c_name = "nfft_buf_size",
             c_type = "int",
             initializer = """ 500 """,
             doc = r"""Size of the Nfft buffer""")

c.add_member(c_name = "post_process",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Perform post processing""")

c.add_member(c_name = "det_init_size",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""The maximum size of the determinant matrix before a resize""")

c.add_member(c_name = "det_n_operations_before_check",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Max number of ops before the test of deviation of the det, M^-1 is performed.""")

c.add_member(c_name = "det_precision_warning",
             c_type = "double",
             initializer = """ 1.e-8 """,
             doc = r"""Threshold for determinant precision warnings""")

c.add_member(c_name = "det_precision_error",
             c_type = "double",
             initializer = """ 1.e-5 """,
             doc = r"""Threshold for determinant precision error""")

c.add_member(c_name = "det_singular_threshold",
             c_type = "double",
             initializer = """ -1 """,
             doc = r"""Bound for the determinant matrix being singular: abs(det) < singular_threshold.
     For negative threshold check if !isnormal(abs(det)).""")

module.add_converter(c)

# Converter for params_t
c = converter_(
        c_type = "triqs_ctint::params_t",
        doc = r"""A struct combining both constr_params_t and solve_params_t""",
)
c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 5001 """,
             doc = r"""Number of tau points for gf<imtime, matrix_valued>""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 500 """,
             doc = r"""Number of Matsubara frequencies for gf<imfreq, matrix_valued>""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "gf_struct_t",
             initializer = """  """,
             doc = r"""block structure of the gf""")

c.add_member(c_name = "use_D",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Switch for dynamic density-density interaction""")

c.add_member(c_name = "use_Jperp",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Switch for dynamic spin-spin interaction""")

c.add_member(c_name = "n_tau_dynamical_interactions",
             c_type = "int",
             initializer = """ res.n_tau """,
             doc = r"""Number of tau pts for D0_tau and jperp_tau""")

c.add_member(c_name = "n_iw_dynamical_interactions",
             c_type = "int",
             initializer = """ res.n_iw """,
             doc = r"""Number of matsubara freqs for D0_iw and jperp_iw""")

c.add_member(c_name = "h_int",
             c_type = "many_body_operator",
             initializer = """  """,
             doc = r"""Interaction Hamiltonian""")

c.add_member(c_name = "n_s",
             c_type = "int",
             initializer = """ 1 """,
             doc = r"""Number of auxiliary spins""")

c.add_member(c_name = "alpha",
             c_type = "alpha_t",
             initializer = """  """,
             doc = r"""Alpha tensor""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = r"""Number of MC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Length of a MC cycles""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = r"""Number of warmup cycles""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*mpi::communicator().rank() """,
             doc = r"""Random seed of the random generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Name of the random generator""")

c.add_member(c_name = "use_double_insertion",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Use double insertion""")

c.add_member(c_name = "insertion_types",
             c_type = "std::vector<int>",
             initializer = """ {} """,
             doc = r"""Types of insertions to use""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum running time in seconds (-1 : no limit)""")

c.add_member(c_name = "max_order",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum pertubation order which is accepted in move::insert/remove
     (<0 : unlimited)""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ mpi::communicator().rank()==0?3:0 """,
             doc = r"""Verbosity""")

c.add_member(c_name = "measure_sign_only",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure Sign only mode""")

c.add_member(c_name = "measure_average_sign",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the MC sign""")

c.add_member(c_name = "measure_average_k",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the average perturbation order""")

c.add_member(c_name = "measure_auto_corr_time",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the auto-correlation time""")

c.add_member(c_name = "measure_histogram",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the average perturbation order distribution""")

c.add_member(c_name = "measure_density",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the density matrix by operator insertion""")

c.add_member(c_name = "measure_M_tau",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure M(tau)""")

c.add_member(c_name = "measure_M_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M(iomega) using nfft""")

c.add_member(c_name = "measure_M4_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M4(iw) NFFT""")

c.add_member(c_name = "measure_M4pp_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M4pp(iw) NFFT""")

c.add_member(c_name = "measure_M4ph_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M4ph(iw) NFFT""")

c.add_member(c_name = "n_iW_M4",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive bosonic Matsubara frequencies in M4""")

c.add_member(c_name = "n_iw_M4",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive fermionic Matsubara frequencies in M4""")

c.add_member(c_name = "measure_M3pp_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3pp(iw)""")

c.add_member(c_name = "measure_M3ph_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3ph(iw)""")

c.add_member(c_name = "n_iw_M3",
             c_type = "int",
             initializer = """ 64 """,
             doc = r"""Number of positive fermionic Matsubara frequencies in M3""")

c.add_member(c_name = "n_iW_M3",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive bosonic Matsubara frequencies in M3""")

c.add_member(c_name = "measure_M3pp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3pp(tau)""")

c.add_member(c_name = "measure_M3ph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3ph(tau)""")

c.add_member(c_name = "measure_M3xph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure M3xph(tau)""")

c.add_member(c_name = "n_tau_M3",
             c_type = "int",
             initializer = """ 201 """,
             doc = r"""Number of imaginary time points in M3""")

c.add_member(c_name = "measure_chi2pp_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure of chi2pp by insertion""")

c.add_member(c_name = "measure_chi2ph_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure of chi2ph by insertion""")

c.add_member(c_name = "n_tau_chi2",
             c_type = "int",
             initializer = """ 201 """,
             doc = r"""Number of imaginary time points in chi2""")

c.add_member(c_name = "n_iw_chi2",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of positive Matsubara frequencies in chi2""")

c.add_member(c_name = "measure_chiAB_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure of chiAB by insertion""")

c.add_member(c_name = "chi_A_vec",
             c_type = "std::vector<many_body_operator>",
             initializer = """ {} """,
             doc = r"""The list of all operators A""")

c.add_member(c_name = "chi_B_vec",
             c_type = "std::vector<many_body_operator>",
             initializer = """ {} """,
             doc = r"""The list of all operators B""")

c.add_member(c_name = "nfft_buf_size",
             c_type = "int",
             initializer = """ 100000 """,
             doc = r"""Size of the Nfft buffer""")

c.add_member(c_name = "post_process",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Perform post processing""")

c.add_member(c_name = "det_init_size",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""The maximum size of the determinant matrix before a resize""")

c.add_member(c_name = "det_n_operations_before_check",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Max number of ops before the test of deviation of the det, M^-1 is performed.""")

c.add_member(c_name = "det_precision_warning",
             c_type = "double",
             initializer = """ 1.e-8 """,
             doc = r"""Threshold for determinant precision warnings""")

c.add_member(c_name = "det_precision_error",
             c_type = "double",
             initializer = """ 1.e-5 """,
             doc = r"""Threshold for determinant precision error""")

c.add_member(c_name = "det_singular_threshold",
             c_type = "double",
             initializer = """ -1 """,
             doc = r"""Bound for the determinant matrix being singular: abs(det) < singular_threshold.
     For negative threshold check if !isnormal(abs(det)).""")

module.add_converter(c)

# Converter for constr_params_t
c = converter_(
        c_type = "triqs_ctint::constr_params_t",
        doc = r"""The parameters for the solver construction""",
)
c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 5001 """,
             doc = r"""Number of tau points for gf<imtime, matrix_valued>""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 500 """,
             doc = r"""Number of Matsubara frequencies for gf<imfreq, matrix_valued>""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "gf_struct_t",
             initializer = """  """,
             doc = r"""block structure of the gf""")

c.add_member(c_name = "use_D",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Switch for dynamic density-density interaction""")

c.add_member(c_name = "use_Jperp",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Switch for dynamic spin-spin interaction""")

c.add_member(c_name = "n_tau_dynamical_interactions",
             c_type = "int",
             initializer = """ res.n_tau """,
             doc = r"""Number of tau pts for D0_tau and jperp_tau""")

c.add_member(c_name = "n_iw_dynamical_interactions",
             c_type = "int",
             initializer = """ res.n_iw """,
             doc = r"""Number of matsubara freqs for D0_iw and jperp_iw""")

module.add_converter(c)


module.generate_code()
