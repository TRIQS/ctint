# Generated automatically using the command :
# c++2py.py ../c++/post_process.hpp -p -m postproc -o postproc
from wrap_generator import *

# The module
module = module_(full_name = "postproc", doc = "", app_name = "postproc")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/vector.hpp>
using namespace triqs::gfs;
""")
module.add_function ("void G_iw_from_M_iw (block_gf<imfreq> M_iw, block_gf<imfreq> G0_shift_iw, block_gf<imfreq> G_iw)", doc = """""")

module.add_function ("void Sigma_iw_from_M_iw (block_gf<imfreq> M_iw, block_gf<imfreq> G0_shift_iw, block_gf<imfreq> G_iw, block_gf<imfreq> Sigmaw, std::vector<std::vector<double>> fact)", doc = """""")

module.generate_code()
