#include "./params.hpp"

namespace triqs_ctint {

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "n_tau", cp.n_tau);
    h5_write(grp, "n_iw", cp.n_iw);
    h5_write(grp, "beta", cp.beta);
    h5_write(grp, "gf_struct", cp.gf_struct);
    h5_write(grp, "use_D", cp.use_D);
    h5_write(grp, "use_Jperp", cp.use_Jperp);
    h5_write(grp, "n_tau_dynamical_interactions", cp.n_tau_dynamical_interactions);
    h5_write(grp, "n_iw_dynamical_interactions", cp.n_iw_dynamical_interactions);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "n_tau", cp.n_tau);
    h5_read(grp, "n_iw", cp.n_iw);
    h5_read(grp, "beta", cp.beta);
    h5_read(grp, "gf_struct", cp.gf_struct);
    h5_read(grp, "use_D", cp.use_D);
    h5_read(grp, "use_Jperp", cp.use_Jperp);
    h5_read(grp, "n_tau_dynamical_interactions", cp.n_tau_dynamical_interactions);
    h5_read(grp, "n_iw_dynamical_interactions", cp.n_iw_dynamical_interactions);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "h_int", sp.h_int);
    h5_write(grp, "n_s", sp.n_s);
    h5_write(grp, "alpha", sp.alpha);
    h5_write(grp, "n_cycles", sp.n_cycles);
    h5_write(grp, "length_cycle", sp.length_cycle);
    h5_write(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5_write(grp, "random_seed", sp.random_seed);
    h5_write(grp, "random_name", sp.random_name);
    h5_write(grp, "use_double_insertion", sp.use_double_insertion);
    h5_write(grp, "max_time", sp.max_time);
    h5_write(grp, "verbosity", sp.verbosity);
    h5_write(grp, "measure_average_sign", sp.measure_average_sign);
    h5_write(grp, "measure_average_k", sp.measure_average_k);
    h5_write(grp, "measure_histogram", sp.measure_histogram);
    h5_write(grp, "measure_M_tau", sp.measure_M_tau);
    h5_write(grp, "measure_M_iw", sp.measure_M_iw);
    h5_write(grp, "measure_M4_iw", sp.measure_M4_iw);
    h5_write(grp, "n_iw_M4", sp.n_iw_M4);
    h5_write(grp, "measure_M3pp_iw", sp.measure_M3pp_iw);
    h5_write(grp, "measure_M3ph_iw", sp.measure_M3ph_iw);
    h5_write(grp, "n_iw_M3", sp.n_iw_M3);
    h5_write(grp, "n_iW_M3", sp.n_iW_M3);
    h5_write(grp, "measure_M3pp_tau", sp.measure_M3pp_tau);
    h5_write(grp, "measure_M3ph_tau", sp.measure_M3ph_tau);
    h5_write(grp, "n_tau_M3", sp.n_tau_M3);
    h5_write(grp, "measure_chi2pp_tau", sp.measure_chi2pp_tau);
    h5_write(grp, "measure_chi2ph_tau", sp.measure_chi2ph_tau);
    h5_write(grp, "n_tau_chi2", sp.n_tau_chi2);
    h5_write(grp, "n_iw_chi2", sp.n_iw_chi2);
    h5_write(grp, "measure_chiAB_tau", sp.measure_chiAB_tau);
    h5_write(grp, "chi_A_vec", sp.chi_A_vec);
    h5_write(grp, "chi_B_vec", sp.chi_B_vec);
    h5_write(grp, "nfft_buf_size", sp.nfft_buf_size);
    h5_write(grp, "post_process", sp.post_process);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    auto grp = h5group.open_group(subgroup_name);
    // Take care! Do not read random_seed and verbosity as they should be different based on mpi rank
    h5_read(grp, "h_int", sp.h_int);
    h5_read(grp, "n_s", sp.n_s);
    h5_read(grp, "alpha", sp.alpha);
    h5_read(grp, "n_cycles", sp.n_cycles);
    h5_read(grp, "length_cycle", sp.length_cycle);
    h5_read(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5_read(grp, "random_name", sp.random_name);
    h5_read(grp, "use_double_insertion", sp.use_double_insertion);
    h5_read(grp, "max_time", sp.max_time);
    h5_read(grp, "measure_average_sign", sp.measure_average_sign);
    h5_read(grp, "measure_average_k", sp.measure_average_k);
    h5_read(grp, "measure_histogram", sp.measure_histogram);
    h5_read(grp, "measure_M_tau", sp.measure_M_tau);
    h5_read(grp, "measure_M_iw", sp.measure_M_iw);
    h5_read(grp, "measure_M4_iw", sp.measure_M4_iw);
    h5_read(grp, "n_iw_M4", sp.n_iw_M4);
    h5_read(grp, "measure_M3pp_iw", sp.measure_M3pp_iw);
    h5_read(grp, "measure_M3ph_iw", sp.measure_M3ph_iw);
    h5_read(grp, "n_iw_M3", sp.n_iw_M3);
    h5_read(grp, "n_iW_M3", sp.n_iW_M3);
    h5_read(grp, "measure_M3pp_tau", sp.measure_M3pp_tau);
    h5_read(grp, "measure_M3ph_tau", sp.measure_M3ph_tau);
    h5_read(grp, "n_tau_M3", sp.n_tau_M3);
    h5_read(grp, "measure_chi2pp_tau", sp.measure_chi2pp_tau);
    h5_read(grp, "measure_chi2ph_tau", sp.measure_chi2ph_tau);
    h5_read(grp, "n_tau_chi2", sp.n_tau_chi2);
    h5_read(grp, "n_iw_chi2", sp.n_iw_chi2);
    h5_read(grp, "measure_chiAB_tau", sp.measure_chiAB_tau);
    h5_read(grp, "chi_A_vec", sp.chi_A_vec);
    h5_read(grp, "chi_B_vec", sp.chi_B_vec);
    h5_read(grp, "nfft_buf_size", sp.nfft_buf_size);
    h5_read(grp, "post_process", sp.post_process);
  }

} // namespace triqs_ctint
