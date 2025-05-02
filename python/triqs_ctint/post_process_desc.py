# Generated automatically using the command :
# c++2py ../../c++/triqs_ctint/post_process.hpp --members_read_only -N triqs_ctint -a triqs_ctint -m post_process -o post_process -C triqs --moduledoc="The TRIQS ctint postprocess functionality" --includes=../../c++ --cxxflags="-std=c++20 $(triqs++ -cxxflags)" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "post_process", doc = r"The TRIQS ctint postprocess functionality", app_name = "triqs_ctint")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes', 'triqs.operators'])

# Add here all includes
module.add_include("triqs_ctint/post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/std_array.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace triqs_ctint;
""")


module.add_function ("chi4_iw_t triqs_ctint::G2c_from_M4 (chi4_iw_t::view_type M4_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw)", doc = r"""Calculate the connected part of the two-particle Green function from M4_iw and M_iw""")

module.add_function ("chi4_iw_t triqs_ctint::G2ppc_from_M4pp (chi4_iw_t::view_type M4pp_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw)", doc = r"""Calculate the connected part of the two-particle Green function from M4pp_iw and M_iw""")

module.add_function ("chi4_iw_t triqs_ctint::G2phc_from_M4ph (chi4_iw_t::view_type M4ph_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw)", doc = r"""Calculate the connected part of the two-particle Green function from M4pp_iw and M_iw""")

module.add_function ("chi4_iw_t triqs_ctint::F_from_G2c (chi4_iw_t::view_type G2c_iw, g_iw_cv_t G_iw)", doc = r"""Calculate the vertex function :math:`F` from G2c_iw and G_iw""")

module.add_function ("chi4_iw_t triqs_ctint::Fpp_from_G2ppc (chi4_iw_t::view_type G2ppc_iw, g_iw_cv_t G_iw)", doc = r"""Calculate the vertex function :math:`Fpp` from G2ppc_iw and G_iw""")

module.add_function ("chi4_iw_t triqs_ctint::Fph_from_G2phc (chi4_iw_t::view_type G2phc_iw, g_iw_cv_t G_iw)", doc = r"""Calculate the vertex function :math:`Fph` from G2phc_iw and G_iw""")

module.add_function ("chi4_iw_t triqs_ctint::G2_from_G2c (chi4_iw_t::view_type G2c_iw, g_iw_cv_t G_iw)", doc = r"""Calculate the two-particle Green function from G2c_iw and G_iw""")

module.add_function ("chi4_iw_t triqs_ctint::G2pp_from_G2ppc (chi4_iw_t::view_type G2ppc_iw, g_iw_cv_t G_iw)", doc = r"""Calculate the two-particle Green function from G2ppc_iw and G_iw""")

module.add_function ("chi4_iw_t triqs_ctint::G2ph_from_G2phc (chi4_iw_t::view_type G2phc_iw, g_iw_cv_t G_iw)", doc = r"""Calculate the two-particle Green function from G2phc_iw and G_iw""")

module.add_function ("chi4_iw_t triqs_ctint::chi_tilde_ph_from_G2phc (chi4_iw_t::view_type G2phc_iw, g_iw_cv_t G_iw, gf_struct_t gf_struct)", doc = r"""Calculate the generalized ph susceptibility from G2phc_iw and G_iw""")

module.add_function ("chi3_iw_t triqs_ctint::chi3_from_M3_PP (chi3_iw_cv_t M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, block_matrix_t dens_G, block_matrix_t M_hartree)", doc = r"""""")

module.add_function ("chi3_iw_t triqs_ctint::chi3_from_M3_PH (chi3_iw_cv_t M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, block_matrix_t dens_G, block_matrix_t M_hartree)", doc = r"""""")

module.add_function ("chi2_tau_t triqs_ctint::chi2_from_chi2_conn_PP (chi2_tau_cv_t chi2_conn_tau, g_iw_cv_t G_iw, block_matrix_t dens_G)", doc = r"""""")

module.add_function ("chi2_tau_t triqs_ctint::chi2_from_chi2_conn_PH (chi2_tau_cv_t chi2_conn_tau, g_iw_cv_t G_iw, block_matrix_t dens_G)", doc = r"""""")

module.add_function ("gf<imtime, matrix_valued> triqs_ctint::chiAB_from_chi2_PP (chi2_tau_cv_t chi2pp_tau, gf_struct_t gf_struct, std::vector<many_body_operator> A_op_vec, std::vector<many_body_operator> B_op_vec)", doc = r"""""")

module.add_function ("gf<imtime, matrix_valued> triqs_ctint::chiAB_from_chi2_PH (chi2_tau_cv_t chi2ph_tau, gf_struct_t gf_struct, std::vector<many_body_operator> A_op_vec, std::vector<many_body_operator> B_op_vec)", doc = r"""""")

module.add_function ("chi2_tau_t triqs_ctint::chi2_conn_from_M3_PP (chi3_tau_t M3pp_tau, chi2_tau_t M3pp_delta, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau, block_matrix_t M_hartree, g_tau_cv_t G0_tau)", doc = r"""""")

module.add_function ("chi2_tau_t triqs_ctint::chi2_conn_from_M3_PH (chi3_tau_t M3ph_tau, chi2_tau_t M3ph_delta, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau, block_matrix_t M_hartree, g_tau_cv_t G0_tau)", doc = r"""""")



module.generate_code()
