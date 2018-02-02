#include "./params.hpp"

namespace triqs_ctint {

  std::pair<int, int> get_int_indices( canonical_ops_t const &op, gf_struct_t const &gf_struct) {

    // The Fundamental operator-set allows for easy check of index validity
    triqs::hilbert_space::fundamental_operator_set fs(gf_struct);
    if (!fs.has_indices(op.indices)) TRIQS_RUNTIME_ERROR << " Index of c/c^+ operator not compatible with Green Function structure ";

    // Get block-name with apply visitor, lambda(0) is called to determine return type ...
    std::string bl_name = visit([](auto idx) { return std::to_string(idx); }, op.indices[0]);

    // Capture positions in block and nonblock list
    int bl_int_idx    = std::distance(gf_struct.cbegin(), gf_struct.find(bl_name));
    auto idx_lst      = gf_struct.at(bl_name);
    int nonbl_int_idx = std::distance(idx_lst.cbegin(), std::find(idx_lst.cbegin(), idx_lst.cend(), op.indices[1]));

    return std::make_pair(bl_int_idx, nonbl_int_idx);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
    h5_write(grp, "n_tau", cp.n_tau);
    h5_write(grp, "n_iw", cp.n_iw);
    h5_write(grp, "beta", cp.beta);
    //h5_write(grp, "gf_struct", cp.gf_struct);
    h5_write(grp, "use_D", cp.use_D);
    h5_write(grp, "use_Jperp", cp.use_Jperp);
    h5_write(grp, "n_tau_dynamical_interactions", cp.n_tau_dynamical_interactions);
    h5_write(grp, "n_iw_dynamical_interactions", cp.n_iw_dynamical_interactions);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
    h5_read(grp, "n_tau", cp.n_tau);
    h5_read(grp, "n_iw", cp.n_iw);
    h5_read(grp, "beta", cp.beta);
    //h5_read(grp, "gf_struct", cp.gf_struct);
    h5_read(grp, "use_D", cp.use_D);
    h5_read(grp, "use_Jperp", cp.use_Jperp);
    h5_read(grp, "n_tau_dynamical_interactions", cp.n_tau_dynamical_interactions);
    h5_read(grp, "n_iw_dynamical_interactions", cp.n_iw_dynamical_interactions);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
    h5_write(grp, "h_int", sp.h_int);
    h5_write(grp, "use_alpha", sp.use_alpha);
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
    h5_write(grp, "measure_M_tau", sp.measure_M_tau);
    h5_write(grp, "measure_M_iw", sp.measure_M_iw);
    h5_write(grp, "measure_M4_iw", sp.measure_M4_iw);
    h5_write(grp, "n_iw_M4", sp.n_iw_M4);
    h5_write(grp, "measure_M3pp_iw", sp.measure_M3pp_iw);
    h5_write(grp, "measure_M3ph_iw ", sp.measure_M3ph_iw);
    h5_write(grp, "n_iw_M3", sp.n_iw_M3);
    h5_write(grp, "n_iW_M3", sp.n_iW_M3);
    h5_write(grp, "measure_M3pp_tau", sp.measure_M3pp_tau);
    h5_write(grp, "measure_M3ph_tau", sp.measure_M3ph_tau);
    h5_write(grp, "n_tau_M3", sp.n_tau_M3);
    h5_write(grp, "measure_chi2pp_tau", sp.measure_chi2pp_tau);
    h5_write(grp, "measure_chi2ph_tau", sp.measure_chi2ph_tau);
    h5_write(grp, "n_tau_chi2", sp.n_tau_chi2);
    h5_write(grp, "n_iw_chi2", sp.n_iw_chi2);
    h5_write(grp, "nfft_buf_size", sp.nfft_buf_size);
    h5_write(grp, "post_process", sp.post_process);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
    h5_read(grp, "h_int", sp.h_int);
    h5_read(grp, "use_alpha", sp.use_alpha);
    h5_read(grp, "n_s", sp.n_s);
    h5_read(grp, "alpha", sp.alpha);
    h5_read(grp, "n_cycles", sp.n_cycles);
    h5_read(grp, "length_cycle", sp.length_cycle);
    h5_read(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5_read(grp, "random_seed", sp.random_seed);
    h5_read(grp, "random_name", sp.random_name);
    h5_read(grp, "use_double_insertion", sp.use_double_insertion);
    h5_read(grp, "max_time", sp.max_time);
    h5_read(grp, "verbosity", sp.verbosity);
    h5_read(grp, "measure_average_sign", sp.measure_average_sign);
    h5_read(grp, "measure_M_tau", sp.measure_M_tau);
    h5_read(grp, "measure_M_iw", sp.measure_M_iw);
    h5_read(grp, "measure_M4_iw", sp.measure_M4_iw);
    h5_read(grp, "n_iw_M4", sp.n_iw_M4);
    h5_read(grp, "measure_M3pp_iw", sp.measure_M3pp_iw);
    h5_read(grp, "measure_M3ph_iw ", sp.measure_M3ph_iw);
    h5_read(grp, "n_iw_M3", sp.n_iw_M3);
    h5_read(grp, "n_iW_M3", sp.n_iW_M3);
    h5_read(grp, "measure_M3pp_tau", sp.measure_M3pp_tau);
    h5_read(grp, "measure_M3ph_tau", sp.measure_M3ph_tau);
    h5_read(grp, "n_tau_M3", sp.n_tau_M3);
    h5_read(grp, "measure_chi2pp_tau", sp.measure_chi2pp_tau);
    h5_read(grp, "measure_chi2ph_tau", sp.measure_chi2ph_tau);
    h5_read(grp, "n_tau_chi2", sp.n_tau_chi2);
    h5_read(grp, "n_iw_chi2", sp.n_iw_chi2);
    h5_read(grp, "nfft_buf_size", sp.nfft_buf_size);
    h5_read(grp, "post_process", sp.post_process);
  }

} // namespace triqs_ctint
