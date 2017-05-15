// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/solver_core.hpp -p --members_read_only -m ctint -o ctint


// --- C++ Python converter for constr_params_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<constr_params_t> {
 static PyObject *c2py(constr_params_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "n_tau"                       , convert_to_python(x.n_tau));
  PyDict_SetItemString( d, "n_iw"                        , convert_to_python(x.n_iw));
  PyDict_SetItemString( d, "beta"                        , convert_to_python(x.beta));
  PyDict_SetItemString( d, "gf_struct"                   , convert_to_python(x.gf_struct));
  PyDict_SetItemString( d, "use_D"                       , convert_to_python(x.use_D));
  PyDict_SetItemString( d, "use_Jperp"                   , convert_to_python(x.use_Jperp));
  PyDict_SetItemString( d, "n_tau_dynamical_interactions", convert_to_python(x.n_tau_dynamical_interactions));
  PyDict_SetItemString( d, "n_iw_dynamical_interactions" , convert_to_python(x.n_iw_dynamical_interactions));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static constr_params_t py2c(PyObject *dic) {
  constr_params_t res;
  _get_optional(dic, "n_tau"                       , res.n_tau                          ,10000);
  _get_optional(dic, "n_iw"                        , res.n_iw                           ,500);
  res.beta = convert_from_python<double>(PyDict_GetItemString(dic, "beta"));
  res.gf_struct = convert_from_python<triqs::hilbert_space::gf_struct_t>(PyDict_GetItemString(dic, "gf_struct"));
  _get_optional(dic, "use_D"                       , res.use_D                          ,false);
  _get_optional(dic, "use_Jperp"                   , res.use_Jperp                      ,false);
  _get_optional(dic, "n_tau_dynamical_interactions", res.n_tau_dynamical_interactions   ,10001);
  _get_optional(dic, "n_iw_dynamical_interactions" , res.n_iw_dynamical_interactions    ,200);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"n_tau","n_iw","beta","gf_struct","use_D","use_Jperp","n_tau_dynamical_interactions","n_iw_dynamical_interactions"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_optional <int                              >(dic, fs, err, "n_tau"                       , "int");
  _check_optional <int                              >(dic, fs, err, "n_iw"                        , "int");
  _check_mandatory<double                           >(dic, fs, err, "beta"                        , "double");
  _check_mandatory<triqs::hilbert_space::gf_struct_t>(dic, fs, err, "gf_struct"                   , "triqs::hilbert_space::gf_struct_t");
  _check_optional <bool                             >(dic, fs, err, "use_D"                       , "bool");
  _check_optional <bool                             >(dic, fs, err, "use_Jperp"                   , "bool");
  _check_optional <int                              >(dic, fs, err, "n_tau_dynamical_interactions", "int");
  _check_optional <int                              >(dic, fs, err, "n_iw_dynamical_interactions" , "int");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class constr_params_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}


// --- C++ Python converter for solve_params_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<solve_params_t> {
 static PyObject *c2py(solve_params_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "hartree_shift"       , convert_to_python(x.hartree_shift));
  PyDict_SetItemString( d, "h_int"               , convert_to_python(x.h_int));
  PyDict_SetItemString( d, "use_alpha"           , convert_to_python(x.use_alpha));
  PyDict_SetItemString( d, "n_s"                 , convert_to_python(x.n_s));
  PyDict_SetItemString( d, "alpha"               , convert_to_python(x.alpha));
  PyDict_SetItemString( d, "n_cycles"            , convert_to_python(x.n_cycles));
  PyDict_SetItemString( d, "length_cycle"        , convert_to_python(x.length_cycle));
  PyDict_SetItemString( d, "n_warmup_cycles"     , convert_to_python(x.n_warmup_cycles));
  PyDict_SetItemString( d, "random_seed"         , convert_to_python(x.random_seed));
  PyDict_SetItemString( d, "random_name"         , convert_to_python(x.random_name));
  PyDict_SetItemString( d, "use_double_insertion", convert_to_python(x.use_double_insertion));
  PyDict_SetItemString( d, "max_time"            , convert_to_python(x.max_time));
  PyDict_SetItemString( d, "verbosity"           , convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "measure_average_sign", convert_to_python(x.measure_average_sign));
  PyDict_SetItemString( d, "measure_M_tau"       , convert_to_python(x.measure_M_tau));
  PyDict_SetItemString( d, "measure_M_iw"        , convert_to_python(x.measure_M_iw));
  PyDict_SetItemString( d, "measure_F_tau"       , convert_to_python(x.measure_F_tau));
  PyDict_SetItemString( d, "measure_M4_iw"       , convert_to_python(x.measure_M4_iw));
  PyDict_SetItemString( d, "n_iw_M4"             , convert_to_python(x.n_iw_M4));
  PyDict_SetItemString( d, "measure_M3pp_iw"     , convert_to_python(x.measure_M3pp_iw));
  PyDict_SetItemString( d, "measure_M3ph_iw"     , convert_to_python(x.measure_M3ph_iw));
  PyDict_SetItemString( d, "n_iw_M3"             , convert_to_python(x.n_iw_M3));
  PyDict_SetItemString( d, "measure_M3pp_tau"    , convert_to_python(x.measure_M3pp_tau));
  PyDict_SetItemString( d, "measure_M3ph_tau"    , convert_to_python(x.measure_M3ph_tau));
  PyDict_SetItemString( d, "n_tau_M3"            , convert_to_python(x.n_tau_M3));
  PyDict_SetItemString( d, "measure_M2pp_tau"    , convert_to_python(x.measure_M2pp_tau));
  PyDict_SetItemString( d, "measure_M2ph_tau"    , convert_to_python(x.measure_M2ph_tau));
  PyDict_SetItemString( d, "measure_M2xph_tau"   , convert_to_python(x.measure_M2xph_tau));
  PyDict_SetItemString( d, "n_tau_M2"            , convert_to_python(x.n_tau_M2));
  PyDict_SetItemString( d, "n_iw_M2"             , convert_to_python(x.n_iw_M2));
  PyDict_SetItemString( d, "nfft_buf_size"       , convert_to_python(x.nfft_buf_size));
  PyDict_SetItemString( d, "post_process"        , convert_to_python(x.post_process));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static solve_params_t py2c(PyObject *dic) {
  solve_params_t res;
  _get_optional(dic, "hartree_shift"       , res.hartree_shift          ,std::vector<double>{});
  res.h_int = convert_from_python<triqs::operators::many_body_operator>(PyDict_GetItemString(dic, "h_int"));
  _get_optional(dic, "use_alpha"           , res.use_alpha              ,false);
  _get_optional(dic, "n_s"                 , res.n_s                    ,2);
  res.alpha = convert_from_python<triqs_ctint::alpha_t>(PyDict_GetItemString(dic, "alpha"));
  res.n_cycles = convert_from_python<int>(PyDict_GetItemString(dic, "n_cycles"));
  _get_optional(dic, "length_cycle"        , res.length_cycle           ,50);
  _get_optional(dic, "n_warmup_cycles"     , res.n_warmup_cycles        ,5000);
  _get_optional(dic, "random_seed"         , res.random_seed            ,34788+928374*triqs::mpi::communicator().rank());
  _get_optional(dic, "random_name"         , res.random_name            ,"");
  _get_optional(dic, "use_double_insertion", res.use_double_insertion   ,false);
  _get_optional(dic, "max_time"            , res.max_time               ,-1);
  _get_optional(dic, "verbosity"           , res.verbosity              ,triqs::mpi::communicator().rank()==0?3:0);
  _get_optional(dic, "measure_average_sign", res.measure_average_sign   ,true);
  _get_optional(dic, "measure_M_tau"       , res.measure_M_tau          ,false);
  _get_optional(dic, "measure_M_iw"        , res.measure_M_iw           ,false);
  _get_optional(dic, "measure_F_tau"       , res.measure_F_tau          ,false);
  _get_optional(dic, "measure_M4_iw"       , res.measure_M4_iw          ,false);
  _get_optional(dic, "n_iw_M4"             , res.n_iw_M4                ,32);
  _get_optional(dic, "measure_M3pp_iw"     , res.measure_M3pp_iw        ,false);
  _get_optional(dic, "measure_M3ph_iw"     , res.measure_M3ph_iw        ,false);
  _get_optional(dic, "n_iw_M3"             , res.n_iw_M3                ,64);
  _get_optional(dic, "measure_M3pp_tau"    , res.measure_M3pp_tau       ,false);
  _get_optional(dic, "measure_M3ph_tau"    , res.measure_M3ph_tau       ,false);
  _get_optional(dic, "n_tau_M3"            , res.n_tau_M3               ,1000);
  _get_optional(dic, "measure_M2pp_tau"    , res.measure_M2pp_tau       ,false);
  _get_optional(dic, "measure_M2ph_tau"    , res.measure_M2ph_tau       ,false);
  _get_optional(dic, "measure_M2xph_tau"   , res.measure_M2xph_tau      ,false);
  _get_optional(dic, "n_tau_M2"            , res.n_tau_M2               ,10000);
  _get_optional(dic, "n_iw_M2"             , res.n_iw_M2                ,128);
  _get_optional(dic, "nfft_buf_size"       , res.nfft_buf_size          ,500);
  _get_optional(dic, "post_process"        , res.post_process           ,true);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"hartree_shift","h_int","use_alpha","n_s","alpha","n_cycles","length_cycle","n_warmup_cycles","random_seed","random_name","use_double_insertion","max_time","verbosity","measure_average_sign","measure_M_tau","measure_M_iw","measure_F_tau","measure_M4_iw","n_iw_M4","measure_M3pp_iw","measure_M3ph_iw","n_iw_M3","measure_M3pp_tau","measure_M3ph_tau","n_tau_M3","measure_M2pp_tau","measure_M2ph_tau","measure_M2xph_tau","n_tau_M2","n_iw_M2","nfft_buf_size","post_process"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_optional <std::vector<double>                 >(dic, fs, err, "hartree_shift"       , "std::vector<double>");
  _check_mandatory<triqs::operators::many_body_operator>(dic, fs, err, "h_int"               , "triqs::operators::many_body_operator");
  _check_optional <bool                                >(dic, fs, err, "use_alpha"           , "bool");
  _check_optional <int                                 >(dic, fs, err, "n_s"                 , "int");
  _check_mandatory<triqs_ctint::alpha_t                >(dic, fs, err, "alpha"               , "triqs_ctint::alpha_t");
  _check_mandatory<int                                 >(dic, fs, err, "n_cycles"            , "int");
  _check_optional <int                                 >(dic, fs, err, "length_cycle"        , "int");
  _check_optional <int                                 >(dic, fs, err, "n_warmup_cycles"     , "int");
  _check_optional <int                                 >(dic, fs, err, "random_seed"         , "int");
  _check_optional <std::string                         >(dic, fs, err, "random_name"         , "std::string");
  _check_optional <bool                                >(dic, fs, err, "use_double_insertion", "bool");
  _check_optional <int                                 >(dic, fs, err, "max_time"            , "int");
  _check_optional <int                                 >(dic, fs, err, "verbosity"           , "int");
  _check_optional <bool                                >(dic, fs, err, "measure_average_sign", "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M_tau"       , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M_iw"        , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_F_tau"       , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M4_iw"       , "bool");
  _check_optional <int                                 >(dic, fs, err, "n_iw_M4"             , "int");
  _check_optional <bool                                >(dic, fs, err, "measure_M3pp_iw"     , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M3ph_iw"     , "bool");
  _check_optional <int                                 >(dic, fs, err, "n_iw_M3"             , "int");
  _check_optional <bool                                >(dic, fs, err, "measure_M3pp_tau"    , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M3ph_tau"    , "bool");
  _check_optional <int                                 >(dic, fs, err, "n_tau_M3"            , "int");
  _check_optional <bool                                >(dic, fs, err, "measure_M2pp_tau"    , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M2ph_tau"    , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_M2xph_tau"   , "bool");
  _check_optional <int                                 >(dic, fs, err, "n_tau_M2"            , "int");
  _check_optional <int                                 >(dic, fs, err, "n_iw_M2"             , "int");
  _check_optional <int                                 >(dic, fs, err, "nfft_buf_size"       , "int");
  _check_optional <bool                                >(dic, fs, err, "post_process"        , "bool");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class solve_params_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}