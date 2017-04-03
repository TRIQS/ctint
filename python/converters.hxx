// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/post_processor.hpp -p -m postproc -o postproc

// --- C++ Python converter for constr_params_t

namespace triqs {
  namespace py_tools {

    template <> struct py_converter<constr_params_t> {
      static PyObject *c2py(constr_params_t const &x) {
        PyObject *d = PyDict_New();
        PyDict_SetItemString(d, "beta", convert_to_python(x.beta));
        PyDict_SetItemString(d, "gf_struct", convert_to_python(x.gf_struct));
        PyDict_SetItemString(d, "n_tau_g0", convert_to_python(x.n_tau_g0));
        PyDict_SetItemString(d, "n_tau_f", convert_to_python(x.n_tau_f));
        PyDict_SetItemString(d, "n_iw", convert_to_python(x.n_iw));
        PyDict_SetItemString(d, "n_tau_dynamical_interactions", convert_to_python(x.n_tau_dynamical_interactions));
        PyDict_SetItemString(d, "n_iw_dynamical_interactions", convert_to_python(x.n_iw_dynamical_interactions));
        PyDict_SetItemString(d, "n_tau_nnt", convert_to_python(x.n_tau_nnt));
        PyDict_SetItemString(d, "n_tau_g2t", convert_to_python(x.n_tau_g2t));
        PyDict_SetItemString(d, "n_w_f_g2w", convert_to_python(x.n_w_f_g2w));
        PyDict_SetItemString(d, "n_w_b_g2w", convert_to_python(x.n_w_b_g2w));
        PyDict_SetItemString(d, "n_tau_M4t", convert_to_python(x.n_tau_M4t));
        PyDict_SetItemString(d, "n_w_f_M4w", convert_to_python(x.n_w_f_M4w));
        PyDict_SetItemString(d, "n_w_b_M4w", convert_to_python(x.n_w_b_M4w));
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
        res.beta      = convert_from_python<double>(PyDict_GetItemString(dic, "beta"));
        res.gf_struct = convert_from_python<gf_struct_t>(PyDict_GetItemString(dic, "gf_struct"));
        _get_optional(dic, "n_tau_g0", res.n_tau_g0, 10001);
        _get_optional(dic, "n_tau_f", res.n_tau_f, 10001);
        _get_optional(dic, "n_iw", res.n_iw, 500);
        _get_optional(dic, "n_tau_dynamical_interactions", res.n_tau_dynamical_interactions, 10001);
        _get_optional(dic, "n_iw_dynamical_interactions", res.n_iw_dynamical_interactions, 200);
        _get_optional(dic, "n_tau_nnt", res.n_tau_nnt, 241);
        _get_optional(dic, "n_tau_g2t", res.n_tau_g2t, 121);
        _get_optional(dic, "n_w_f_g2w", res.n_w_f_g2w, 60);
        _get_optional(dic, "n_w_b_g2w", res.n_w_b_g2w, 60);
        _get_optional(dic, "n_tau_M4t", res.n_tau_M4t, 4);
        _get_optional(dic, "n_w_f_M4w", res.n_w_f_M4w, 4);
        _get_optional(dic, "n_w_b_M4w", res.n_w_b_M4w, 4);
        return res;
      }

      template <typename T> static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
        if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
          fs << "\n"
             << ++err << " The parameter " << name << " does not have the right type : expecting " << tname << " in C++, but got '"
             << PyDict_GetItemString(dic, name)->ob_t->tp_name << "' in Python.";
      }

      template <typename T> static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
        if (!PyDict_Contains(dic, pyref::string(name)))
          fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
        else
          _check<T>(dic, fs, err, name, tname);
      }

      template <typename T> static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
        if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
      }

      static bool is_convertible(PyObject *dic, bool raise_exception) {
        if (!PyDict_Check(dic)) {
          if (raise_exception) { PyErr_SetString(PyExc_TypeError, "Not a python dict"); }
          return false;
        }
        std::stringstream fs, fs2;
        int err = 0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
        std::vector<std::string> ks, all_keys = {"beta",
                                                 "gf_struct",
                                                 "n_tau_g0",
                                                 "n_tau_f",
                                                 "n_iw",
                                                 "n_tau_dynamical_interactions",
                                                 "n_iw_dynamical_interactions",
                                                 "n_tau_nnt",
                                                 "n_tau_g2t",
                                                 "n_w_f_g2w",
                                                 "n_w_b_g2w",
                                                 "n_tau_M4t",
                                                 "n_w_f_M4w",
                                                 "n_w_b_M4w"};
        pyref keys = PyDict_Keys(dic);
        if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
          fs << "\nThe dict keys are not strings";
          goto _error;
        }
        ks = convert_from_python<std::vector<std::string>>(keys);
        for (auto &k : ks)
          if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
            fs << "\n" << ++err << " The parameter '" << k << "' is not recognized.";
#endif

        _check_mandatory<double>(dic, fs, err, "beta", "double");
        _check_mandatory<gf_struct_t>(dic, fs, err, "gf_struct", "gf_struct_t");
        _check_optional<int>(dic, fs, err, "n_tau_g0", "int");
        _check_optional<int>(dic, fs, err, "n_tau_f", "int");
        _check_optional<int>(dic, fs, err, "n_iw", "int");
        _check_optional<int>(dic, fs, err, "n_tau_dynamical_interactions", "int");
        _check_optional<int>(dic, fs, err, "n_iw_dynamical_interactions", "int");
        _check_optional<int>(dic, fs, err, "n_tau_nnt", "int");
        _check_optional<int>(dic, fs, err, "n_tau_g2t", "int");
        _check_optional<int>(dic, fs, err, "n_w_f_g2w", "int");
        _check_optional<int>(dic, fs, err, "n_w_b_g2w", "int");
        _check_optional<int>(dic, fs, err, "n_tau_M4t", "int");
        _check_optional<int>(dic, fs, err, "n_w_f_M4w", "int");
        _check_optional<int>(dic, fs, err, "n_w_b_M4w", "int");
        if (err) goto _error;
        return true;

      _error:
        fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err << " error" << (err > 1 ? "s" : "")
            << " in Python -> C++ transcription for the class constr_params_t\n"
            << fs.str();
        if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
        return false;
      }
    };
  }
}
