
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// ==================== module functions ====================

// F_from_G2c
static auto const _c2py_fun_0 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::F_from_G2c(G2_conn_iw, G_iw); },
                                      "G2_conn_iw", "G_iw")};

// Fph_from_G2ph_conn
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2ph_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::Fph_from_G2ph_conn(G2ph_conn_iw, G_iw); },
                                      "G2ph_conn_iw", "G_iw")};

// Fpp_from_G2pp_conn
static auto const _c2py_fun_2 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2pp_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::Fpp_from_G2pp_conn(G2pp_conn_iw, G_iw); },
                                      "G2pp_conn_iw", "G_iw")};

// G2_conn_from_M4
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>, triqs::gfs::tensor_valued<4>,
                                      nda::C_layout, 2>::const_view_type M4_iw,
                 triqs_ctint::g_reg_iw_cv_t M_iw, triqs_ctint::g_reg_iw_cv_t G0_iw) { return triqs_ctint::G2_conn_from_M4(M4_iw, M_iw, G0_iw); },
              "M4_iw", "M_iw", "G0_iw")};

// G2_from_G2c
static auto const _c2py_fun_4 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::G2_from_G2c(G2_conn_iw, G_iw); },
                                      "G2_conn_iw", "G_iw")};

// G2ph_conn_from_M4ph
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>, triqs::gfs::tensor_valued<4>,
                           nda::C_layout, 2>::const_view_type M4ph_iw,
      triqs_ctint::g_reg_iw_cv_t M_iw, triqs_ctint::g_reg_iw_cv_t G0_iw) { return triqs_ctint::G2ph_conn_from_M4ph(M4ph_iw, M_iw, G0_iw); },
   "M4ph_iw", "M_iw", "G0_iw")};

// G2ph_from_G2ph_conn
static auto const _c2py_fun_6 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2ph_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::G2ph_from_G2ph_conn(G2ph_conn_iw, G_iw); },
                                      "G2ph_conn_iw", "G_iw")};

// G2pp_conn_from_M4pp
static auto const _c2py_fun_7 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>, triqs::gfs::tensor_valued<4>,
                           nda::C_layout, 2>::const_view_type M4pp_iw,
      triqs_ctint::g_reg_iw_cv_t M_iw, triqs_ctint::g_reg_iw_cv_t G0_iw) { return triqs_ctint::G2pp_conn_from_M4pp(M4pp_iw, M_iw, G0_iw); },
   "M4pp_iw", "M_iw", "G0_iw")};

// G2pp_from_G2pp_conn
static auto const _c2py_fun_8 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2pp_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::G2pp_from_G2pp_conn(G2pp_conn_iw, G_iw); },
                                      "G2pp_conn_iw", "G_iw")};

// chi3_from_M3_PH
static auto const _c2py_fun_9 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs_ctint::chi3_iw_cv_t M3_iw, triqs_ctint::g_reg_iw_cv_t M_iw, triqs_ctint::g_reg_iw_cv_t G0_iw, const triqs_ctint::block_matrix_t &dens_G,
      const triqs_ctint::block_matrix_t &M_hartree) { return triqs_ctint::chi3_from_M3_PH(M3_iw, M_iw, G0_iw, dens_G, M_hartree); },
   "M3_iw", "M_iw", "G0_iw", "dens_G", "M_hartree")};

// chi3_from_M3_PP
static auto const _c2py_fun_10 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs_ctint::chi3_iw_cv_t M3_iw, triqs_ctint::g_reg_iw_cv_t M_iw, triqs_ctint::g_reg_iw_cv_t G0_iw, const triqs_ctint::block_matrix_t &dens_G,
      const triqs_ctint::block_matrix_t &M_hartree) { return triqs_ctint::chi3_from_M3_PP(M3_iw, M_iw, G0_iw, dens_G, M_hartree); },
   "M3_iw", "M_iw", "G0_iw", "dens_G", "M_hartree")};

// chiAB_from_chi2_PH
static auto const _c2py_fun_11 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs_ctint::chi2_tau_cv_t chi2ph_tau, const triqs::gfs::gf_struct_t &gf_struct,
      const std::vector<triqs::operators::many_body_operator> &A_op_vec, const std::vector<triqs::operators::many_body_operator> &B_op_vec) {
     return triqs_ctint::chiAB_from_chi2_PH(chi2ph_tau, gf_struct, A_op_vec, B_op_vec);
   },
   "chi2ph_tau", "gf_struct", "A_op_vec", "B_op_vec")};

// chiAB_from_chi2_PP
static auto const _c2py_fun_12 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs_ctint::chi2_tau_cv_t chi2pp_tau, const triqs::gfs::gf_struct_t &gf_struct,
      const std::vector<triqs::operators::many_body_operator> &A_op_vec, const std::vector<triqs::operators::many_body_operator> &B_op_vec) {
     return triqs_ctint::chiAB_from_chi2_PP(chi2pp_tau, gf_struct, A_op_vec, B_op_vec);
   },
   "chi2pp_tau", "gf_struct", "A_op_vec", "B_op_vec")};

// chi_tilde_ph_from_G2ph_conn
static auto const _c2py_fun_13 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                                                              triqs::gfs::tensor_valued<4>, nda::C_layout, 2>::const_view_type G2ph_conn_iw,
                                         triqs_ctint::g_reg_iw_cv_t G_iw) { return triqs_ctint::chi_tilde_ph_from_G2ph_conn(G2ph_conn_iw, G_iw); },
                                      "G2ph_conn_iw", "G_iw")};

static const auto _c2py_doc_0  = _c2py_fun_0.doc(R"DOC(
Calculate the vertex function $F$ from G2_conn_iw and G_iw
)DOC");
static const auto _c2py_doc_1  = _c2py_fun_1.doc(R"DOC(
Calculate the vertex function $Fph$ from G2ph_conn_iw and G_iw
)DOC");
static const auto _c2py_doc_2  = _c2py_fun_2.doc(R"DOC(
Calculate the vertex function $Fpp$ from G2pp_conn_iw and G_iw
)DOC");
static const auto _c2py_doc_3  = _c2py_fun_3.doc(R"DOC(
Calculate the connected part of the two-particle Green function from M4_iw and M_iw
)DOC");
static const auto _c2py_doc_4  = _c2py_fun_4.doc(R"DOC(
Calculate the two-particle Green function from G2_conn_iw and G_iw
)DOC");
static const auto _c2py_doc_5  = _c2py_fun_5.doc(R"DOC(
Calculate the connected part of the two-particle Green function from M4pp_iw and M_iw
)DOC");
static const auto _c2py_doc_6  = _c2py_fun_6.doc(R"DOC(
Calculate the two-particle Green function from G2ph_conn_iw and G_iw
)DOC");
static const auto _c2py_doc_7  = _c2py_fun_7.doc(R"DOC(
Calculate the connected part of the two-particle Green function from M4pp_iw and M_iw
)DOC");
static const auto _c2py_doc_8  = _c2py_fun_8.doc(R"DOC(
Calculate the two-particle Green function from G2pp_conn_iw and G_iw
)DOC");
static const auto _c2py_doc_9  = _c2py_fun_9.doc(R"DOC()DOC");
static const auto _c2py_doc_10 = _c2py_fun_10.doc(R"DOC()DOC");
static const auto _c2py_doc_11 = _c2py_fun_11.doc(R"DOC()DOC");
static const auto _c2py_doc_12 = _c2py_fun_12.doc(R"DOC()DOC");
static const auto _c2py_doc_13 = _c2py_fun_13.doc(R"DOC(
Calculate the generalized ph susceptibility from G2ph_conn_iw and G_iw
)DOC");
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"F_from_G2c", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"Fph_from_G2ph_conn", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"Fpp_from_G2pp_conn", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"G2_conn_from_M4", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"G2_from_G2c", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"G2ph_conn_from_M4ph", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"G2ph_from_G2ph_conn", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {"G2pp_conn_from_M4pp", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"G2pp_from_G2pp_conn", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {"chi3_from_M3_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_9>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_9.c_str()},
   {"chi3_from_M3_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_10>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_10.c_str()},
   {"chiAB_from_chi2_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_11>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_11.c_str()},
   {"chiAB_from_chi2_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_12>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_12.c_str()},
   {"chi_tilde_ph_from_G2ph_conn", (PyCFunction)c2py::pyfkw<_c2py_fun_13>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_13.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "post_process",                          /* name of module */
                                        R"RAWDOC(Postprocess utilities.)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_post_process() {

  if (not c2py::check_python_version("post_process")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)

#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
