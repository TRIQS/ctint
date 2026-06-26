
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
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = triqs_ctint::constr_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "triqs_ctint.solver_core.ConstrParamsT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_ctint::constr_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing triqs_ctint::constr_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("n_tau", self_c.n_tau, true);
  de("dlr_wmax", self_c.dlr_wmax, false);
  de("dlr_eps", self_c.dlr_eps, true);
  de("beta", self_c.beta, false);
  de("gf_struct", self_c.gf_struct, false);
  de("use_D", self_c.use_D, true);
  de("use_Jperp", self_c.use_Jperp, true);
  de("n_tau_dynamical_interactions", self_c.n_tau_dynamical_interactions, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
dlr_wmax : {par_0}

beta : {par_1}

gf_struct : {par_2}

n_tau : {par_3}, default=5001

dlr_eps : {par_4}, default=1e-10

use_D : {par_5}, default=false

use_Jperp : {par_6}, default=false

n_tau_dynamical_interactions : {par_7}, default=this->n_tau

)DOC",
   "par",
   {c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<triqs::gfs::gf_struct_t>(), c2py::python_typename<int>(),
    c2py::python_typename<double>(), c2py::python_typename<bool>(), c2py::python_typename<bool>(), c2py::python_typename<int>()});
// block_names
static auto const _c2py_fun_0 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_0 const &self) -> decltype(auto) { return self.block_names(); }, "self")};

// n_blocks
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_0 const &self) -> decltype(auto) { return self.n_blocks(); }, "self")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
Names of the blocks of the Green's function.
)DOC");
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC(
Number of blocks of the Green's function.
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {
   {"block_names", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"n_blocks", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0 = R"DOC(Number of imaginary-time points for the single-particle quantities.)DOC";
constexpr auto _c2py_doc_member_1 = R"DOC(DLR bandwidth cutoff :math:`w_{max} = \Lambda / \beta` for the single-particle quantities.)DOC";
constexpr auto _c2py_doc_member_2 = R"DOC(DLR error tolerance :math:`\epsilon` for the single-particle quantities.)DOC";
constexpr auto _c2py_doc_member_3 = R"DOC(Inverse temperature :math:`\beta`.)DOC";
constexpr auto _c2py_doc_member_4 = R"DOC(Structure of the Green's function (names and sizes of blocks).)DOC";
constexpr auto _c2py_doc_member_5 = R"DOC(Use a dynamic density-density interaction?)DOC";
constexpr auto _c2py_doc_member_6 = R"DOC(Use a dynamic spin-spin interaction?)DOC";
constexpr auto _c2py_doc_member_7 = R"DOC(Number of imaginary-time points for :math:`D_0(\tau)` and :math:`J_\perp(\tau)`.)DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["n_tau"]                        = self_c.n_tau;
  dic["dlr_wmax"]                     = self_c.dlr_wmax;
  dic["dlr_eps"]                      = self_c.dlr_eps;
  dic["beta"]                         = self_c.beta;
  dic["gf_struct"]                    = self_c.gf_struct;
  dic["use_D"]                        = self_c.use_D;
  dic["use_Jperp"]                    = self_c.use_Jperp;
  dic["n_tau_dynamical_interactions"] = self_c.n_tau_dynamical_interactions;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::n_tau, _c2py_cls_0>("n_tau", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::dlr_wmax, _c2py_cls_0>("dlr_wmax", _c2py_doc_member_1),
   c2py::getsetdef_from_member<&_c2py_cls_0::dlr_eps, _c2py_cls_0>("dlr_eps", _c2py_doc_member_2),
   c2py::getsetdef_from_member<&_c2py_cls_0::beta, _c2py_cls_0>("beta", _c2py_doc_member_3),
   c2py::getsetdef_from_member<&_c2py_cls_0::gf_struct, _c2py_cls_0>("gf_struct", _c2py_doc_member_4),
   c2py::getsetdef_from_member<&_c2py_cls_0::use_D, _c2py_cls_0>("use_D", _c2py_doc_member_5),
   c2py::getsetdef_from_member<&_c2py_cls_0::use_Jperp, _c2py_cls_0>("use_Jperp", _c2py_doc_member_6),
   c2py::getsetdef_from_member<&_c2py_cls_0::n_tau_dynamical_interactions, _c2py_cls_0>("n_tau_dynamical_interactions", _c2py_doc_member_7),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> =
   R"DOC(Parameters used for constructing the solver class.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = triqs_ctint::solve_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "triqs_ctint.solver_core.SolveParamsT";

static int synth_constructor_1(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_ctint::solve_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_1> *)self)->_c = new _c2py_cls_1{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing triqs_ctint::solve_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  de("h_int", self_c.h_int, false);
  de("n_s", self_c.n_s, true);
  de("alpha", self_c.alpha, false);
  de("n_cycles", self_c.n_cycles, false);
  de("length_cycle", self_c.length_cycle, true);
  de("max_length_cycle", self_c.max_length_cycle, true);
  de("target_auto_corr_time", self_c.target_auto_corr_time, true);
  de("n_warmup_cycles", self_c.n_warmup_cycles, true);
  de("max_warmup_cycles", self_c.max_warmup_cycles, true);
  de("random_seed", self_c.random_seed, true);
  de("random_name", self_c.random_name, true);
  de("use_double_insertion", self_c.use_double_insertion, true);
  de("insertion_types", self_c.insertion_types, true);
  de("use_auxiliary_spin_flip", self_c.use_auxiliary_spin_flip, true);
  de("max_time", self_c.max_time, true);
  de("max_order", self_c.max_order, true);
  de("verbosity", self_c.verbosity, true);
  de("rethrow_exception", self_c.rethrow_exception, true);
  de("measure_sign_only", self_c.measure_sign_only, true);
  de("measure_average_sign", self_c.measure_average_sign, true);
  de("measure_average_k", self_c.measure_average_k, true);
  de("measure_histogram", self_c.measure_histogram, true);
  de("measure_densities", self_c.measure_densities, true);
  de("measure_density_matrix", self_c.measure_density_matrix, true);
  de("measure_M_tau", self_c.measure_M_tau, true);
  de("measure_M_iw", self_c.measure_M_iw, true);
  de("measure_M4_iw", self_c.measure_M4_iw, true);
  de("measure_M4pp_iw", self_c.measure_M4pp_iw, true);
  de("measure_M4ph_iw", self_c.measure_M4ph_iw, true);
  de("n_iW_M4", self_c.n_iW_M4, true);
  de("n_iw_M4", self_c.n_iw_M4, true);
  de("measure_M3pp_iw", self_c.measure_M3pp_iw, true);
  de("measure_M3ph_iw", self_c.measure_M3ph_iw, true);
  de("measure_M3pp_iw_full", self_c.measure_M3pp_iw_full, true);
  de("measure_M3ph_iw_full", self_c.measure_M3ph_iw_full, true);
  de("n_iw_M3", self_c.n_iw_M3, true);
  de("n_iW_M3", self_c.n_iW_M3, true);
  de("dlr2d_compress_grid", self_c.dlr2d_compress_grid, true);
  de("measure_M3pp_tau", self_c.measure_M3pp_tau, true);
  de("measure_M3ph_tau", self_c.measure_M3ph_tau, true);
  de("measure_M3xph_tau", self_c.measure_M3xph_tau, true);
  de("n_tau_M3", self_c.n_tau_M3, true);
  de("measure_chi2pp_tau", self_c.measure_chi2pp_tau, true);
  de("measure_chi2ph_tau", self_c.measure_chi2ph_tau, true);
  de("measure_chiAB_tau", self_c.measure_chiAB_tau, true);
  de("chi_ops", self_c.chi_ops, true);
  de("measure_static_obs", self_c.measure_static_obs, true);
  de("static_obs", self_c.static_obs, true);
  de("n_tau_static_obs", self_c.n_tau_static_obs, true);
  de("nfft_buf_size", self_c.nfft_buf_size, true);
  de("nfft_tol", self_c.nfft_tol, true);
  de("post_process", self_c.post_process, true);
  de("det_init_size", self_c.det_init_size, true);
  de("det_n_operations_before_check", self_c.det_n_operations_before_check, true);
  de("det_precision_warning", self_c.det_precision_warning, true);
  de("det_precision_error", self_c.det_precision_error, true);
  de("det_singular_threshold", self_c.det_singular_threshold, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = synth_constructor_1;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
h_int : {par_0}

alpha : {par_1}

n_cycles : {par_2}

n_s : {par_3}, default=2

length_cycle : {par_4}, default=-1

max_length_cycle : {par_5}, default=5000

target_auto_corr_time : {par_6}, default=2.0

n_warmup_cycles : {par_7}, default=-1

max_warmup_cycles : {par_8}, default=100000

random_seed : {par_9}, default=34788

random_name : {par_10}, default=""

use_double_insertion : {par_11}, default=true

insertion_types : {par_12}, default={}

use_auxiliary_spin_flip : {par_13}, default=false

max_time : {par_14}, default=-1

max_order : {par_15}, default=-1

verbosity : {par_16}, default== 0 ? 3 : 0

rethrow_exception : {par_17}, default=true

measure_sign_only : {par_18}, default=false

measure_average_sign : {par_19}, default=true

measure_average_k : {par_20}, default=true

measure_histogram : {par_21}, default=false

measure_densities : {par_22}, default=true

measure_density_matrix : {par_23}, default=false

measure_M_tau : {par_24}, default=true

measure_M_iw : {par_25}, default=false

measure_M4_iw : {par_26}, default=false

measure_M4pp_iw : {par_27}, default=false

measure_M4ph_iw : {par_28}, default=false

n_iW_M4 : {par_29}, default=32

n_iw_M4 : {par_30}, default=32

measure_M3pp_iw : {par_31}, default=false

measure_M3ph_iw : {par_32}, default=false

measure_M3pp_iw_full : {par_33}, default=false

measure_M3ph_iw_full : {par_34}, default=false

n_iw_M3 : {par_35}, default=64

n_iW_M3 : {par_36}, default=32

dlr2d_compress_grid : {par_37}, default=false

measure_M3pp_tau : {par_38}, default=false

measure_M3ph_tau : {par_39}, default=false

measure_M3xph_tau : {par_40}, default=false

n_tau_M3 : {par_41}, default=201

measure_chi2pp_tau : {par_42}, default=false

measure_chi2ph_tau : {par_43}, default=false

measure_chiAB_tau : {par_44}, default=false

chi_ops : {par_45}, default={}

measure_static_obs : {par_46}, default=false

static_obs : {par_47}, default={}

n_tau_static_obs : {par_48}, default=10

nfft_buf_size : {par_49}, default=100000

nfft_tol : {par_50}, default=1e-8

post_process : {par_51}, default=true

det_init_size : {par_52}, default=1000

det_n_operations_before_check : {par_53}, default=100

det_precision_warning : {par_54}, default=1.e-8

det_precision_error : {par_55}, default=1.e-5

det_singular_threshold : {par_56}, default=-1

)DOC",
                      "par",
                      {c2py::python_typename<triqs::operators::many_body_operator>(),
                       c2py::python_typename<triqs_ctint::alpha_t>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<std::string>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<std::vector<int>>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<std::vector<std::pair<triqs::operators::many_body_operator, triqs::operators::many_body_operator>>>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<std::vector<triqs::operators::many_body_operator>>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_8  = R"DOC(Interacting part of the local Hamiltonian.)DOC";
constexpr auto _c2py_doc_member_9  = R"DOC(Number of auxiliary spins.)DOC";
constexpr auto _c2py_doc_member_10 = R"DOC(The :math:`\alpha` tensor used in the determinantal expansion.)DOC";
constexpr auto _c2py_doc_member_11 = R"DOC(Number of QMC cycles.)DOC";
constexpr auto _c2py_doc_member_12 = R"DOC(Length of a single QMC cycle (-1: automatically determined from the autocorrelation time).)DOC";
constexpr auto _c2py_doc_member_13 = R"DOC(Maximum allowed length_cycle when auto-determined (safety cap).)DOC";
constexpr auto _c2py_doc_member_14 = R"DOC(Target autocorrelation time in units of length_cycle (used when length_cycle=-1).)DOC";
constexpr auto _c2py_doc_member_15 = R"DOC(Number of cycles for thermalization (-1: automatic convergence detection).)DOC";
constexpr auto _c2py_doc_member_16 = R"DOC(Maximum number of warmup cycles when using automatic warmup (safety cap).)DOC";
constexpr auto _c2py_doc_member_17 = R"DOC(Seed for the random number generator (shared by all MPI ranks; the rank is used as the
spawn key to derive an independent stream per rank, see triqs::mc_tools::random_generator).)DOC";
constexpr auto _c2py_doc_member_18 = R"DOC(Name of the random number generator.)DOC";
constexpr auto _c2py_doc_member_19 = R"DOC(Use double insertion?)DOC";
constexpr auto _c2py_doc_member_20 = R"DOC(Types of insertions to use.)DOC";
constexpr auto _c2py_doc_member_21 = R"DOC(Use auxiliary spin-flip insertion (requires :math:`n_s = 2`)?)DOC";
constexpr auto _c2py_doc_member_22 = R"DOC(Maximum runtime in seconds, use -1 to set infinite.)DOC";
constexpr auto _c2py_doc_member_23 = R"DOC(Maximum perturbation order accepted during insertion and removal moves (use -1 for unlimited).)DOC";
constexpr auto _c2py_doc_member_24 = R"DOC(Verbosity level.)DOC";
constexpr auto _c2py_doc_member_25 = R"DOC(Catch exceptions on the nodes and rethrow them on rank 0?)DOC";
constexpr auto _c2py_doc_member_26 = R"DOC(Measure the sign only?)DOC";
constexpr auto _c2py_doc_member_27 = R"DOC(Measure the Monte-Carlo sign?)DOC";
constexpr auto _c2py_doc_member_28 = R"DOC(Measure the average perturbation order?)DOC";
constexpr auto _c2py_doc_member_29 = R"DOC(Measure the perturbation-order distribution?)DOC";
constexpr auto _c2py_doc_member_30 = R"DOC(Measure the diagonal densities by operator insertion?)DOC";
constexpr auto _c2py_doc_member_31 = R"DOC(Measure the full density matrix by operator insertion (needed for :math:`\chi^{(3)}`)?)DOC";
constexpr auto _c2py_doc_member_32 = R"DOC(Measure :math:`M(\tau)`?)DOC";
constexpr auto _c2py_doc_member_33 = R"DOC(Measure :math:`M(i\omega)` using NFFT?)DOC";
constexpr auto _c2py_doc_member_34 = R"DOC(Measure :math:`M^{(4)}(i\omega)` using NFFT?)DOC";
constexpr auto _c2py_doc_member_35 = R"DOC(Measure :math:`M^{(4)}_{pp}(i\omega)` using NFFT?)DOC";
constexpr auto _c2py_doc_member_36 = R"DOC(Measure :math:`M^{(4)}_{ph}(i\omega)` using NFFT?)DOC";
constexpr auto _c2py_doc_member_37 = R"DOC(Number of positive bosonic Matsubara frequencies in :math:`M^{(4)}`.)DOC";
constexpr auto _c2py_doc_member_38 = R"DOC(Number of positive fermionic Matsubara frequencies in :math:`M^{(4)}`.)DOC";
constexpr auto _c2py_doc_member_39 = R"DOC(Measure :math:`M^{(3)}_{pp}(i\omega)`?)DOC";
constexpr auto _c2py_doc_member_40 = R"DOC(Measure :math:`M^{(3)}_{ph}(i\omega)`?)DOC";
constexpr auto _c2py_doc_member_41 = R"DOC(Measure :math:`M^{(3)}_{pp}(i\omega)` on the full frequency grid?)DOC";
constexpr auto _c2py_doc_member_42 = R"DOC(Measure :math:`M^{(3)}_{ph}(i\omega)` on the full frequency grid?)DOC";
constexpr auto _c2py_doc_member_43 = R"DOC(Number of positive fermionic Matsubara frequencies in :math:`M^{(3)}`.)DOC";
constexpr auto _c2py_doc_member_44 = R"DOC(Number of positive bosonic Matsubara frequencies in :math:`M^{(3)}`.)DOC";
constexpr auto _c2py_doc_member_45 = R"DOC(Compress the DLR2D imaginary-frequency grid?)DOC";
constexpr auto _c2py_doc_member_46 = R"DOC(Measure :math:`M^{(3)}_{pp}(\tau)`?)DOC";
constexpr auto _c2py_doc_member_47 = R"DOC(Measure :math:`M^{(3)}_{ph}(\tau)`?)DOC";
constexpr auto _c2py_doc_member_48 = R"DOC(Measure :math:`M^{(3)}_{xph}(\tau)`?)DOC";
constexpr auto _c2py_doc_member_49 = R"DOC(Number of imaginary-time points in :math:`M^{(3)}`.)DOC";
constexpr auto _c2py_doc_member_50 = R"DOC(Measure :math:`\chi^{(2)}_{pp}(\tau)` by insertion?)DOC";
constexpr auto _c2py_doc_member_51 = R"DOC(Measure :math:`\chi^{(2)}_{ph}(\tau)` by insertion?)DOC";
constexpr auto _c2py_doc_member_52 = R"DOC(Measure :math:`\chi_{AB}(\tau)` by insertion?)DOC";
constexpr auto _c2py_doc_member_53 = R"DOC(List of operator pairs :math:`(A, B)` for the :math:`\chi_{AB}` measurement.)DOC";
constexpr auto _c2py_doc_member_54 = R"DOC(Measure static expectation values of arbitrary operators by insertion with tau-averaging?)DOC";
constexpr auto _c2py_doc_member_55 =
   R"DOC(List of operators :math:`C_i` for the static-observable measurement (measures :math:`\langle C_i \rangle`).)DOC";
constexpr auto _c2py_doc_member_56 = R"DOC(Number of tau points for tau-averaging in the static_obs measurement.)DOC";
constexpr auto _c2py_doc_member_57 = R"DOC(Size of the NFFT buffer.)DOC";
constexpr auto _c2py_doc_member_58 = R"DOC(Tolerance for the NFFT transform.)DOC";
constexpr auto _c2py_doc_member_59 = R"DOC(Perform post-processing?)DOC";
constexpr auto _c2py_doc_member_60 = R"DOC(The maximum size of the determinant matrix before a resize.)DOC";
constexpr auto _c2py_doc_member_61 = R"DOC(Maximum number of operations before testing the accuracy of :math:`\det(M)` and :math:`M^{-1}`.)DOC";
constexpr auto _c2py_doc_member_62 = R"DOC(Threshold for determinant precision warnings.)DOC";
constexpr auto _c2py_doc_member_63 = R"DOC(Threshold for determinant precision errors.)DOC";
constexpr auto _c2py_doc_member_64 =
   R"DOC(Bound for the determinant matrix being singular (if :math:`< 0`, checks for subnormal numbers instead).)DOC";
static PyObject *prop_get_dict_1(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  c2py::pydict dic;
  dic["h_int"]                         = self_c.h_int;
  dic["n_s"]                           = self_c.n_s;
  dic["alpha"]                         = self_c.alpha;
  dic["n_cycles"]                      = self_c.n_cycles;
  dic["length_cycle"]                  = self_c.length_cycle;
  dic["max_length_cycle"]              = self_c.max_length_cycle;
  dic["target_auto_corr_time"]         = self_c.target_auto_corr_time;
  dic["n_warmup_cycles"]               = self_c.n_warmup_cycles;
  dic["max_warmup_cycles"]             = self_c.max_warmup_cycles;
  dic["random_seed"]                   = self_c.random_seed;
  dic["random_name"]                   = self_c.random_name;
  dic["use_double_insertion"]          = self_c.use_double_insertion;
  dic["insertion_types"]               = self_c.insertion_types;
  dic["use_auxiliary_spin_flip"]       = self_c.use_auxiliary_spin_flip;
  dic["max_time"]                      = self_c.max_time;
  dic["max_order"]                     = self_c.max_order;
  dic["verbosity"]                     = self_c.verbosity;
  dic["rethrow_exception"]             = self_c.rethrow_exception;
  dic["measure_sign_only"]             = self_c.measure_sign_only;
  dic["measure_average_sign"]          = self_c.measure_average_sign;
  dic["measure_average_k"]             = self_c.measure_average_k;
  dic["measure_histogram"]             = self_c.measure_histogram;
  dic["measure_densities"]             = self_c.measure_densities;
  dic["measure_density_matrix"]        = self_c.measure_density_matrix;
  dic["measure_M_tau"]                 = self_c.measure_M_tau;
  dic["measure_M_iw"]                  = self_c.measure_M_iw;
  dic["measure_M4_iw"]                 = self_c.measure_M4_iw;
  dic["measure_M4pp_iw"]               = self_c.measure_M4pp_iw;
  dic["measure_M4ph_iw"]               = self_c.measure_M4ph_iw;
  dic["n_iW_M4"]                       = self_c.n_iW_M4;
  dic["n_iw_M4"]                       = self_c.n_iw_M4;
  dic["measure_M3pp_iw"]               = self_c.measure_M3pp_iw;
  dic["measure_M3ph_iw"]               = self_c.measure_M3ph_iw;
  dic["measure_M3pp_iw_full"]          = self_c.measure_M3pp_iw_full;
  dic["measure_M3ph_iw_full"]          = self_c.measure_M3ph_iw_full;
  dic["n_iw_M3"]                       = self_c.n_iw_M3;
  dic["n_iW_M3"]                       = self_c.n_iW_M3;
  dic["dlr2d_compress_grid"]           = self_c.dlr2d_compress_grid;
  dic["measure_M3pp_tau"]              = self_c.measure_M3pp_tau;
  dic["measure_M3ph_tau"]              = self_c.measure_M3ph_tau;
  dic["measure_M3xph_tau"]             = self_c.measure_M3xph_tau;
  dic["n_tau_M3"]                      = self_c.n_tau_M3;
  dic["measure_chi2pp_tau"]            = self_c.measure_chi2pp_tau;
  dic["measure_chi2ph_tau"]            = self_c.measure_chi2ph_tau;
  dic["measure_chiAB_tau"]             = self_c.measure_chiAB_tau;
  dic["chi_ops"]                       = self_c.chi_ops;
  dic["measure_static_obs"]            = self_c.measure_static_obs;
  dic["static_obs"]                    = self_c.static_obs;
  dic["n_tau_static_obs"]              = self_c.n_tau_static_obs;
  dic["nfft_buf_size"]                 = self_c.nfft_buf_size;
  dic["nfft_tol"]                      = self_c.nfft_tol;
  dic["post_process"]                  = self_c.post_process;
  dic["det_init_size"]                 = self_c.det_init_size;
  dic["det_n_operations_before_check"] = self_c.det_n_operations_before_check;
  dic["det_precision_warning"]         = self_c.det_precision_warning;
  dic["det_precision_error"]           = self_c.det_precision_error;
  dic["det_singular_threshold"]        = self_c.det_singular_threshold;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_1::h_int, _c2py_cls_1>("h_int", _c2py_doc_member_8),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_s, _c2py_cls_1>("n_s", _c2py_doc_member_9),
   c2py::getsetdef_from_member<&_c2py_cls_1::alpha, _c2py_cls_1>("alpha", _c2py_doc_member_10),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_cycles, _c2py_cls_1>("n_cycles", _c2py_doc_member_11),
   c2py::getsetdef_from_member<&_c2py_cls_1::length_cycle, _c2py_cls_1>("length_cycle", _c2py_doc_member_12),
   c2py::getsetdef_from_member<&_c2py_cls_1::max_length_cycle, _c2py_cls_1>("max_length_cycle", _c2py_doc_member_13),
   c2py::getsetdef_from_member<&_c2py_cls_1::target_auto_corr_time, _c2py_cls_1>("target_auto_corr_time", _c2py_doc_member_14),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_warmup_cycles, _c2py_cls_1>("n_warmup_cycles", _c2py_doc_member_15),
   c2py::getsetdef_from_member<&_c2py_cls_1::max_warmup_cycles, _c2py_cls_1>("max_warmup_cycles", _c2py_doc_member_16),
   c2py::getsetdef_from_member<&_c2py_cls_1::random_seed, _c2py_cls_1>("random_seed", _c2py_doc_member_17),
   c2py::getsetdef_from_member<&_c2py_cls_1::random_name, _c2py_cls_1>("random_name", _c2py_doc_member_18),
   c2py::getsetdef_from_member<&_c2py_cls_1::use_double_insertion, _c2py_cls_1>("use_double_insertion", _c2py_doc_member_19),
   c2py::getsetdef_from_member<&_c2py_cls_1::insertion_types, _c2py_cls_1>("insertion_types", _c2py_doc_member_20),
   c2py::getsetdef_from_member<&_c2py_cls_1::use_auxiliary_spin_flip, _c2py_cls_1>("use_auxiliary_spin_flip", _c2py_doc_member_21),
   c2py::getsetdef_from_member<&_c2py_cls_1::max_time, _c2py_cls_1>("max_time", _c2py_doc_member_22),
   c2py::getsetdef_from_member<&_c2py_cls_1::max_order, _c2py_cls_1>("max_order", _c2py_doc_member_23),
   c2py::getsetdef_from_member<&_c2py_cls_1::verbosity, _c2py_cls_1>("verbosity", _c2py_doc_member_24),
   c2py::getsetdef_from_member<&_c2py_cls_1::rethrow_exception, _c2py_cls_1>("rethrow_exception", _c2py_doc_member_25),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_sign_only, _c2py_cls_1>("measure_sign_only", _c2py_doc_member_26),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_average_sign, _c2py_cls_1>("measure_average_sign", _c2py_doc_member_27),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_average_k, _c2py_cls_1>("measure_average_k", _c2py_doc_member_28),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_histogram, _c2py_cls_1>("measure_histogram", _c2py_doc_member_29),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_densities, _c2py_cls_1>("measure_densities", _c2py_doc_member_30),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_density_matrix, _c2py_cls_1>("measure_density_matrix", _c2py_doc_member_31),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M_tau, _c2py_cls_1>("measure_M_tau", _c2py_doc_member_32),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M_iw, _c2py_cls_1>("measure_M_iw", _c2py_doc_member_33),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M4_iw, _c2py_cls_1>("measure_M4_iw", _c2py_doc_member_34),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M4pp_iw, _c2py_cls_1>("measure_M4pp_iw", _c2py_doc_member_35),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M4ph_iw, _c2py_cls_1>("measure_M4ph_iw", _c2py_doc_member_36),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_iW_M4, _c2py_cls_1>("n_iW_M4", _c2py_doc_member_37),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_iw_M4, _c2py_cls_1>("n_iw_M4", _c2py_doc_member_38),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3pp_iw, _c2py_cls_1>("measure_M3pp_iw", _c2py_doc_member_39),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3ph_iw, _c2py_cls_1>("measure_M3ph_iw", _c2py_doc_member_40),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3pp_iw_full, _c2py_cls_1>("measure_M3pp_iw_full", _c2py_doc_member_41),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3ph_iw_full, _c2py_cls_1>("measure_M3ph_iw_full", _c2py_doc_member_42),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_iw_M3, _c2py_cls_1>("n_iw_M3", _c2py_doc_member_43),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_iW_M3, _c2py_cls_1>("n_iW_M3", _c2py_doc_member_44),
   c2py::getsetdef_from_member<&_c2py_cls_1::dlr2d_compress_grid, _c2py_cls_1>("dlr2d_compress_grid", _c2py_doc_member_45),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3pp_tau, _c2py_cls_1>("measure_M3pp_tau", _c2py_doc_member_46),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3ph_tau, _c2py_cls_1>("measure_M3ph_tau", _c2py_doc_member_47),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_M3xph_tau, _c2py_cls_1>("measure_M3xph_tau", _c2py_doc_member_48),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_tau_M3, _c2py_cls_1>("n_tau_M3", _c2py_doc_member_49),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_chi2pp_tau, _c2py_cls_1>("measure_chi2pp_tau", _c2py_doc_member_50),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_chi2ph_tau, _c2py_cls_1>("measure_chi2ph_tau", _c2py_doc_member_51),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_chiAB_tau, _c2py_cls_1>("measure_chiAB_tau", _c2py_doc_member_52),
   c2py::getsetdef_from_member<&_c2py_cls_1::chi_ops, _c2py_cls_1>("chi_ops", _c2py_doc_member_53),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_static_obs, _c2py_cls_1>("measure_static_obs", _c2py_doc_member_54),
   c2py::getsetdef_from_member<&_c2py_cls_1::static_obs, _c2py_cls_1>("static_obs", _c2py_doc_member_55),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_tau_static_obs, _c2py_cls_1>("n_tau_static_obs", _c2py_doc_member_56),
   c2py::getsetdef_from_member<&_c2py_cls_1::nfft_buf_size, _c2py_cls_1>("nfft_buf_size", _c2py_doc_member_57),
   c2py::getsetdef_from_member<&_c2py_cls_1::nfft_tol, _c2py_cls_1>("nfft_tol", _c2py_doc_member_58),
   c2py::getsetdef_from_member<&_c2py_cls_1::post_process, _c2py_cls_1>("post_process", _c2py_doc_member_59),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_init_size, _c2py_cls_1>("det_init_size", _c2py_doc_member_60),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_n_operations_before_check, _c2py_cls_1>("det_n_operations_before_check", _c2py_doc_member_61),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_precision_warning, _c2py_cls_1>("det_precision_warning", _c2py_doc_member_62),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_precision_error, _c2py_cls_1>("det_precision_error", _c2py_doc_member_63),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_singular_threshold, _c2py_cls_1>("det_singular_threshold", _c2py_doc_member_64),
   {"__dict__", (getter)prop_get_dict_1, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> =
   R"DOC(Parameters passed to the solve method of the solver class.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_1>;
// --------- class _c2py_cls_2 -----------
using _c2py_cls_2                                            = triqs_ctint::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_2>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_2> = "triqs_ctint.solver_core.SolverCore";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_2, const triqs_ctint::constr_params_t &>("constr_params_")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_2> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_2> = _c2py_init_0.doc(R"DOC(
Construct a CT-INT solver.

Parameters
----------
constr_params_ : {par_0}
   Set of parameters used to construct the solver.
)DOC",
                                                                    {{c2py::python_typename<const triqs_ctint::constr_params_t &>()}});
// post_process
static auto const _c2py_fun_2 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_2 &self) -> decltype(auto) { return self.post_process(); }, "self")};

// solve
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_2 &self, const triqs_ctint::solve_params_t &solve_params) -> decltype(auto) { return self.solve(solve_params); },
                 "self", "solve_params")};

static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC(
Retrigger post-processing with the last set of parameters.
)DOC");
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
Solve the impurity problem with a CT-INT calculation.

Parameters
----------
solve_params : {par_0}
   Set of parameters used for the solve.
)DOC",
                                                {{c2py::python_typename<const triqs_ctint::solve_params_t &>()}});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_2>[] = {
   {"post_process", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"solve", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_2>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_2>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_2>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_65 = R"DOC(Non-interacting Green's function :math:`G_0(i\omega)` in Matsubara frequencies.)DOC";
constexpr auto _c2py_doc_member_66 = R"DOC(Inverse of the non-interacting Green's function :math:`G_0^{-1}(i\omega)`.)DOC";
constexpr auto _c2py_doc_member_67 = R"DOC(Dynamic density-density interaction :math:`D_0(i\omega)` in Matsubara frequencies (DLR mesh).)DOC";
constexpr auto _c2py_doc_member_68 = R"DOC(Dynamic spin-spin interaction :math:`J_\perp(i\omega)` in Matsubara frequencies (DLR mesh).)DOC";
constexpr auto _c2py_doc_member_69 = R"DOC(The shifted non-interacting Green's function in Matsubara frequencies.)DOC";
constexpr auto _c2py_doc_member_70 = R"DOC(The shifted non-interacting Green's function in imaginary time.)DOC";
constexpr auto _c2py_doc_member_71 = R"DOC(Parameters used to construct the solver.)DOC";
constexpr auto _c2py_doc_member_72 = R"DOC(Parameters used in the last solve (empty until the solver has been run).)DOC";
constexpr auto _c2py_doc_member_73 = R"DOC(Average sign of the CTINT)DOC";
constexpr auto _c2py_doc_member_74 = R"DOC(Total number of measures)DOC";
constexpr auto _c2py_doc_member_75 = R"DOC(Average perturbation order)DOC";
constexpr auto _c2py_doc_member_76 = R"DOC(Error bar for average sign)DOC";
constexpr auto _c2py_doc_member_77 = R"DOC(Error bar for average perturbation order)DOC";
constexpr auto _c2py_doc_member_78 = R"DOC(Auto-correlation time)DOC";
constexpr auto _c2py_doc_member_79 = R"DOC(Number of warmup cycles actually performed)DOC";
constexpr auto _c2py_doc_member_80 = R"DOC(The length_cycle value used during accumulation (after auto-determination))DOC";
constexpr auto _c2py_doc_member_81 = R"DOC(Warmup time in seconds)DOC";
constexpr auto _c2py_doc_member_82 = R"DOC(Accumulation time in seconds)DOC";
constexpr auto _c2py_doc_member_83 = R"DOC(Average perturbation order distribution)DOC";
constexpr auto _c2py_doc_member_84 = R"DOC(The diagonal densities (measured by operator insertion))DOC";
constexpr auto _c2py_doc_member_85 = R"DOC(Error bars for densities from linear binning)DOC";
constexpr auto _c2py_doc_member_86 = R"DOC(The full density matrix (measured by operator insertion, needed for chi3))DOC";
constexpr auto _c2py_doc_member_87 = R"DOC(Error bars for density_matrix from linear binning)DOC";
constexpr auto _c2py_doc_member_88 = R"DOC(Building block for the Green function in imaginary time (Eq. (23) in Notes))DOC";
constexpr auto _c2py_doc_member_89 = R"DOC(Hartree-term of M_tau)DOC";
constexpr auto _c2py_doc_member_90 = R"DOC(Same as M_tau, but measured directly in Matsubara frequencies using NFFT on DLR grid)DOC";
constexpr auto _c2py_doc_member_91 = R"DOC(Building block for the full vertex function measured directly in Matsubara frequencies using NFFT)DOC";
constexpr auto _c2py_doc_member_92 =
   R"DOC(Building block for the full vertex function (pp channel) measured directly in Matsubara frequencies using NFFT)DOC";
constexpr auto _c2py_doc_member_93 =
   R"DOC(Building block for the full vertex function (ph channel) measured directly in Matsubara frequencies using NFFT)DOC";
constexpr auto _c2py_doc_member_94 =
   R"DOC(Building block for the fermion boson vertex (pp channel) in Matsubara frequencies using NFFT on DLR2D grid)DOC";
constexpr auto _c2py_doc_member_95 =
   R"DOC(Building block for the fermion boson vertex (ph channel) in Matsubara frequencies using NFFT on DLR2D grid)DOC";
constexpr auto _c2py_doc_member_96 =
   R"DOC(Building block for the fermion boson vertex (pp channel) in Matsubara frequencies using NFFT on full grid)DOC";
constexpr auto _c2py_doc_member_97 =
   R"DOC(Building block for the fermion boson vertex (ph channel) in Matsubara frequencies using NFFT on full grid)DOC";
constexpr auto _c2py_doc_member_98  = R"DOC(Building block for the fermion boson vertex (pp channel) in imaginary time)DOC";
constexpr auto _c2py_doc_member_99  = R"DOC(Building block for the fermion boson vertex (ph channel) in imaginary time)DOC";
constexpr auto _c2py_doc_member_100 = R"DOC(Building block for the fermion boson vertex (xph channel) in imaginary time)DOC";
constexpr auto _c2py_doc_member_101 = R"DOC(Equal-time peak in M3pp_tau)DOC";
constexpr auto _c2py_doc_member_102 = R"DOC(Equal-time peak in M3ph_tau)DOC";
constexpr auto _c2py_doc_member_103 = R"DOC(Equal-time peak in M3xph_tau)DOC";
constexpr auto _c2py_doc_member_104 =
   R"DOC(The equal time correlator $$ in the particle-particle channel in imaginary times as obtained by operator insertion)DOC";
constexpr auto _c2py_doc_member_105 =
   R"DOC(The equal time correlator $$ in the particle-hole channel in imaginary times as obtained by operator insertion)DOC";
constexpr auto _c2py_doc_member_106 = R"DOC(The correlation function $$ in imaginary times)DOC";
constexpr auto _c2py_doc_member_107 = R"DOC(Static expectation values $ C_i $ measured by operator insertion with tau-averaging)DOC";
constexpr auto _c2py_doc_member_108 = R"DOC(Error bars for static_obs from linear binning)DOC";
constexpr auto _c2py_doc_member_109 = R"DOC(The Fourier-transform of M_tau. Dependent on M_tau)DOC";
constexpr auto _c2py_doc_member_110 = R"DOC(Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_iw)DOC";
constexpr auto _c2py_doc_member_111 = R"DOC(Dynamic self-energy in Matsubara frequencies (DLR, decays to zero). Dependent on M_iw)DOC";
constexpr auto _c2py_doc_member_112 = R"DOC(Static (Hartree) part of the self-energy. Sigma = Sigma_dyn + Sigma_hartree)DOC";
constexpr auto _c2py_doc_member_113 = R"DOC(Building block for the fermion boson vertex (pp channel) in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_114 = R"DOC(Building block for the fermion boson vertex (ph channel) in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_115 = R"DOC(Building block for the fermion boson vertex (xph channel) in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_116 = R"DOC(The two-particle vertex function in purely fermionic notation (iw1, iw2, iw3))DOC";
constexpr auto _c2py_doc_member_117 = R"DOC(The two-particle vertex function (pp channel))DOC";
constexpr auto _c2py_doc_member_118 = R"DOC(The two-particle vertex function (ph channel))DOC";
constexpr auto _c2py_doc_member_119 = R"DOC(The connected part of the two-particle Green function)DOC";
constexpr auto _c2py_doc_member_120 = R"DOC(The connected part of the two-particle Green function (pp channel))DOC";
constexpr auto _c2py_doc_member_121 = R"DOC(The connected part of the two-particle Green function (ph channel))DOC";
constexpr auto _c2py_doc_member_122 = R"DOC(The two-particle Green function)DOC";
constexpr auto _c2py_doc_member_123 = R"DOC(The two-particle Green function (pp channel))DOC";
constexpr auto _c2py_doc_member_124 = R"DOC(The two-particle Green function (ph channel))DOC";
constexpr auto _c2py_doc_member_125 = R"DOC(The equal time correlator $$ in the particle-particle channel in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_126 = R"DOC(The equal time correlator $$ in the particle-hole channel in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_127 = R"DOC(The correlation function $$ in imaginary frequencies)DOC";
constexpr auto _c2py_doc_member_128 = R"DOC(The equal time correlator $$ in the particle-particle channel in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_129 = R"DOC(The equal time correlator $$ in the particle-hole channel in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_130 = R"DOC(The equal time correlator $$ in the particle-hole-cross channel in Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_131 =
   R"DOC(The equal time correlator $$ in the particle-particle channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement on DLR2D grid)DOC";
constexpr auto _c2py_doc_member_132 =
   R"DOC(The equal time correlator $$ in the particle-hole channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement on DLR2D grid)DOC";
constexpr auto _c2py_doc_member_133 = R"DOC(chi3 pp channel from full-grid NFFT M3 measurement)DOC";
constexpr auto _c2py_doc_member_134 = R"DOC(chi3 ph channel from full-grid NFFT M3 measurement)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_2>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_2::G0_iw, _c2py_cls_2>("G0_iw", _c2py_doc_member_65),
   c2py::getsetdef_from_member<&_c2py_cls_2::G0_iw_inv, _c2py_cls_2>("G0_iw_inv", _c2py_doc_member_66),
   c2py::getsetdef_from_member<&_c2py_cls_2::D0_iw, _c2py_cls_2>("D0_iw", _c2py_doc_member_67),
   c2py::getsetdef_from_member<&_c2py_cls_2::Jperp_iw, _c2py_cls_2>("Jperp_iw", _c2py_doc_member_68),
   c2py::getsetdef_from_member<&_c2py_cls_2::G0_shift_iw, _c2py_cls_2>("G0_shift_iw", _c2py_doc_member_69),
   c2py::getsetdef_from_member<&_c2py_cls_2::G0_shift_tau, _c2py_cls_2>("G0_shift_tau", _c2py_doc_member_70),
   c2py::getsetdef_from_member<&_c2py_cls_2::constr_params, _c2py_cls_2>("constr_params", _c2py_doc_member_71),
   c2py::getsetdef_from_member<&_c2py_cls_2::last_solve_params, _c2py_cls_2>("last_solve_params", _c2py_doc_member_72),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_sign, _c2py_cls_2>("average_sign", _c2py_doc_member_73),
   c2py::getsetdef_from_member<&_c2py_cls_2::nmeasures, _c2py_cls_2>("nmeasures", _c2py_doc_member_74),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_k, _c2py_cls_2>("average_k", _c2py_doc_member_75),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_sign_error, _c2py_cls_2>("average_sign_error", _c2py_doc_member_76),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_k_error, _c2py_cls_2>("average_k_error", _c2py_doc_member_77),
   c2py::getsetdef_from_member<&_c2py_cls_2::auto_corr_time, _c2py_cls_2>("auto_corr_time", _c2py_doc_member_78),
   c2py::getsetdef_from_member<&_c2py_cls_2::warmup_cycles_done, _c2py_cls_2>("warmup_cycles_done", _c2py_doc_member_79),
   c2py::getsetdef_from_member<&_c2py_cls_2::length_cycle_used, _c2py_cls_2>("length_cycle_used", _c2py_doc_member_80),
   c2py::getsetdef_from_member<&_c2py_cls_2::warmup_time, _c2py_cls_2>("warmup_time", _c2py_doc_member_81),
   c2py::getsetdef_from_member<&_c2py_cls_2::accumulation_time, _c2py_cls_2>("accumulation_time", _c2py_doc_member_82),
   c2py::getsetdef_from_member<&_c2py_cls_2::histogram, _c2py_cls_2>("histogram", _c2py_doc_member_83),
   c2py::getsetdef_from_member<&_c2py_cls_2::densities, _c2py_cls_2>("densities", _c2py_doc_member_84),
   c2py::getsetdef_from_member<&_c2py_cls_2::densities_errors, _c2py_cls_2>("densities_errors", _c2py_doc_member_85),
   c2py::getsetdef_from_member<&_c2py_cls_2::density_matrix, _c2py_cls_2>("density_matrix", _c2py_doc_member_86),
   c2py::getsetdef_from_member<&_c2py_cls_2::density_matrix_errors, _c2py_cls_2>("density_matrix_errors", _c2py_doc_member_87),
   c2py::getsetdef_from_member<&_c2py_cls_2::M_tau, _c2py_cls_2>("M_tau", _c2py_doc_member_88),
   c2py::getsetdef_from_member<&_c2py_cls_2::M_hartree, _c2py_cls_2>("M_hartree", _c2py_doc_member_89),
   c2py::getsetdef_from_member<&_c2py_cls_2::M_iw_nfft, _c2py_cls_2>("M_iw_nfft", _c2py_doc_member_90),
   c2py::getsetdef_from_member<&_c2py_cls_2::M4_iw, _c2py_cls_2>("M4_iw", _c2py_doc_member_91),
   c2py::getsetdef_from_member<&_c2py_cls_2::M4pp_iw, _c2py_cls_2>("M4pp_iw", _c2py_doc_member_92),
   c2py::getsetdef_from_member<&_c2py_cls_2::M4ph_iw, _c2py_cls_2>("M4ph_iw", _c2py_doc_member_93),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3pp_iw_nfft, _c2py_cls_2>("M3pp_iw_nfft", _c2py_doc_member_94),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3ph_iw_nfft, _c2py_cls_2>("M3ph_iw_nfft", _c2py_doc_member_95),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3pp_iw_nfft_full, _c2py_cls_2>("M3pp_iw_nfft_full", _c2py_doc_member_96),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3ph_iw_nfft_full, _c2py_cls_2>("M3ph_iw_nfft_full", _c2py_doc_member_97),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3pp_tau, _c2py_cls_2>("M3pp_tau", _c2py_doc_member_98),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3ph_tau, _c2py_cls_2>("M3ph_tau", _c2py_doc_member_99),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3xph_tau, _c2py_cls_2>("M3xph_tau", _c2py_doc_member_100),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3pp_delta, _c2py_cls_2>("M3pp_delta", _c2py_doc_member_101),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3ph_delta, _c2py_cls_2>("M3ph_delta", _c2py_doc_member_102),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3xph_delta, _c2py_cls_2>("M3xph_delta", _c2py_doc_member_103),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi2pp_tau, _c2py_cls_2>("chi2pp_tau", _c2py_doc_member_104),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi2ph_tau, _c2py_cls_2>("chi2ph_tau", _c2py_doc_member_105),
   c2py::getsetdef_from_member<&_c2py_cls_2::chiAB_tau, _c2py_cls_2>("chiAB_tau", _c2py_doc_member_106),
   c2py::getsetdef_from_member<&_c2py_cls_2::static_obs, _c2py_cls_2>("static_obs", _c2py_doc_member_107),
   c2py::getsetdef_from_member<&_c2py_cls_2::static_obs_errors, _c2py_cls_2>("static_obs_errors", _c2py_doc_member_108),
   c2py::getsetdef_from_member<&_c2py_cls_2::M_iw, _c2py_cls_2>("M_iw", _c2py_doc_member_109),
   c2py::getsetdef_from_member<&_c2py_cls_2::G_iw, _c2py_cls_2>("G_iw", _c2py_doc_member_110),
   c2py::getsetdef_from_member<&_c2py_cls_2::Sigma_dyn_iw, _c2py_cls_2>("Sigma_dyn_iw", _c2py_doc_member_111),
   c2py::getsetdef_from_member<&_c2py_cls_2::Sigma_hartree, _c2py_cls_2>("Sigma_hartree", _c2py_doc_member_112),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3pp_iw, _c2py_cls_2>("M3pp_iw", _c2py_doc_member_113),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3ph_iw, _c2py_cls_2>("M3ph_iw", _c2py_doc_member_114),
   c2py::getsetdef_from_member<&_c2py_cls_2::M3xph_iw, _c2py_cls_2>("M3xph_iw", _c2py_doc_member_115),
   c2py::getsetdef_from_member<&_c2py_cls_2::F_iw, _c2py_cls_2>("F_iw", _c2py_doc_member_116),
   c2py::getsetdef_from_member<&_c2py_cls_2::Fpp_iw, _c2py_cls_2>("Fpp_iw", _c2py_doc_member_117),
   c2py::getsetdef_from_member<&_c2py_cls_2::Fph_iw, _c2py_cls_2>("Fph_iw", _c2py_doc_member_118),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_conn_iw, _c2py_cls_2>("G2_conn_iw", _c2py_doc_member_119),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2pp_conn_iw, _c2py_cls_2>("G2pp_conn_iw", _c2py_doc_member_120),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2ph_conn_iw, _c2py_cls_2>("G2ph_conn_iw", _c2py_doc_member_121),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw, _c2py_cls_2>("G2_iw", _c2py_doc_member_122),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2pp_iw, _c2py_cls_2>("G2pp_iw", _c2py_doc_member_123),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2ph_iw, _c2py_cls_2>("G2ph_iw", _c2py_doc_member_124),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi2pp_iw, _c2py_cls_2>("chi2pp_iw", _c2py_doc_member_125),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi2ph_iw, _c2py_cls_2>("chi2ph_iw", _c2py_doc_member_126),
   c2py::getsetdef_from_member<&_c2py_cls_2::chiAB_iw, _c2py_cls_2>("chiAB_iw", _c2py_doc_member_127),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3pp_iw, _c2py_cls_2>("chi3pp_iw", _c2py_doc_member_128),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3ph_iw, _c2py_cls_2>("chi3ph_iw", _c2py_doc_member_129),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3xph_iw, _c2py_cls_2>("chi3xph_iw", _c2py_doc_member_130),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3pp_iw_nfft, _c2py_cls_2>("chi3pp_iw_nfft", _c2py_doc_member_131),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3ph_iw_nfft, _c2py_cls_2>("chi3ph_iw_nfft", _c2py_doc_member_132),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3pp_iw_nfft_full, _c2py_cls_2>("chi3pp_iw_nfft_full", _c2py_doc_member_133),
   c2py::getsetdef_from_member<&_c2py_cls_2::chi3ph_iw_nfft_full, _c2py_cls_2>("chi3ph_iw_nfft_full", _c2py_doc_member_134),

   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_2> = R"DOC(The CT-INT solver.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_2>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "solver_core",                           /* name of module */
                                        R"RAWDOC(The TRIQS ctint solver)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_solver_core() {

  if (not c2py::check_python_version("solver_core")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_2>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "ConstrParamsT");
  _add_type(_c2py_cls_1, "SolveParamsT");
  _add_type(_c2py_cls_2, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_2>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
