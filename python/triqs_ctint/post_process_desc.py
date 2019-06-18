# Generated automatically using the command :
# c++2py ../../c++/triqs_ctint/post_process.hpp --members_read_only -N triqs_ctint -a triqs_ctint -m post_process -o post_process -C pytriqs --moduledoc="The TRIQS ctint postprocess functionality" --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "post_process", doc = r"The TRIQS ctint postprocess functionality", app_name = "triqs_ctint")

# Imports
module.add_imports(*['pytriqs.gf', 'pytriqs.operators'])

# Add here all includes
module.add_include("triqs_ctint/post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace triqs_ctint;
""")


module.add_function ("triqs_ctint::chi4_iw_t triqs_ctint::G2c_from_M4 (chi4_iw_t::view_type M4_iw, triqs_ctint::g_iw_cv_t M_iw, triqs_ctint::g_iw_cv_t G0_iw)", doc = r"""Calculate the connected part of the two-particle Green function from M4_iw and M_iw""")

module.add_function ("triqs_ctint::chi4_iw_t triqs_ctint::F_from_G2c (chi4_iw_t::view_type G2c_iw, triqs_ctint::g_iw_cv_t G_iw)", doc = r"""Calculate the vertex function :math:`F` from G2c_iw and G_iw""")

module.add_function ("triqs_ctint::chi4_iw_t triqs_ctint::G2_from_G2c (chi4_iw_t::view_type G2c_iw, triqs_ctint::g_iw_cv_t G_iw)", doc = r"""Calculate the two-particle Green function from G2c_iw and G_iw""")

module.add_function ("triqs_ctint::chi4_iw_t triqs_ctint::chi_tilde_ph_from_G2c (chi4_iw_t::view_type G2c_iw, triqs_ctint::g_iw_cv_t G_iw, triqs::hilbert_space::gf_struct_t gf_struct)", doc = r"""Calculate the generalized ph susceptibility from G2c_iw and G_iw""")

module.add_function ("triqs_ctint::chi3_iw_t triqs_ctint::chi3_from_M3_PP (triqs_ctint::chi3_iw_cv_t M3_iw, triqs_ctint::g_iw_cv_t M_iw, triqs_ctint::g_iw_cv_t G0_iw, triqs_ctint::block_matrix_t dens_G)", doc = r"""""")

module.add_function ("triqs_ctint::chi3_iw_t triqs_ctint::chi3_from_M3_PH (triqs_ctint::chi3_iw_cv_t M3_iw, triqs_ctint::g_iw_cv_t M_iw, triqs_ctint::g_iw_cv_t G0_iw, triqs_ctint::block_matrix_t dens_G)", doc = r"""""")

module.add_function ("triqs_ctint::chi2_tau_t triqs_ctint::chi2_from_chi2_conn_PP (triqs_ctint::chi2_tau_cv_t chi2_conn_tau, triqs_ctint::g_iw_cv_t G_iw, triqs_ctint::block_matrix_t dens_G)", doc = r"""""")

module.add_function ("triqs_ctint::chi2_tau_t triqs_ctint::chi2_from_chi2_conn_PH (triqs_ctint::chi2_tau_cv_t chi2_conn_tau, triqs_ctint::g_iw_cv_t G_iw, triqs_ctint::block_matrix_t dens_G)", doc = r"""""")

module.add_function ("gf<triqs::gfs::imtime,triqs::gfs::matrix_valued> triqs_ctint::chiAB_from_chi2_PP (triqs_ctint::chi2_tau_cv_t chi2pp_tau, triqs::hilbert_space::gf_struct_t gf_struct, std::vector<many_body_operator> A_op_vec, std::vector<many_body_operator> B_op_vec)", doc = r"""""")

module.add_function ("gf<triqs::gfs::imtime,triqs::gfs::matrix_valued> triqs_ctint::chiAB_from_chi2_PH (triqs_ctint::chi2_tau_cv_t chi2ph_tau, triqs::hilbert_space::gf_struct_t gf_struct, std::vector<many_body_operator> A_op_vec, std::vector<many_body_operator> B_op_vec)", doc = r"""""")

module.add_function ("triqs_ctint::chi2_tau_t triqs_ctint::chi2_conn_from_M3_PP (triqs_ctint::chi3_tau_t M3pp_tau, triqs_ctint::chi2_tau_t M3pp_delta, triqs_ctint::g_iw_cv_t M_iw, triqs_ctint::g_iw_cv_t G0_iw, triqs_ctint::g_tau_cv_t M_tau, triqs_ctint::block_matrix_t M_hartree, triqs_ctint::g_tau_cv_t G0_tau)", doc = r"""""")

module.add_function ("triqs_ctint::chi2_tau_t triqs_ctint::chi2_conn_from_M3_PH (triqs_ctint::chi3_tau_t M3ph_tau, triqs_ctint::chi2_tau_t M3ph_delta, triqs_ctint::g_iw_cv_t M_iw, triqs_ctint::g_iw_cv_t G0_iw, triqs_ctint::g_tau_cv_t M_tau, triqs_ctint::block_matrix_t M_hartree, triqs_ctint::g_tau_cv_t G0_tau)", doc = r"""""")



module.generate_code()