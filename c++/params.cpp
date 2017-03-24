#include "./params.hpp"

namespace triqs_ctint {

  array<array<dcomplex, 2>, 2> params_t::get_U() const { // TODO Clean and Rework!

    // Calculate U-Matrix U(block1,block2)(i,j)
    array<array<dcomplex, 2>, 2> U(n_blocks(), n_blocks());

    // Extract flat U-Matrix, where indeces combine the block (non-leading) and non-block (leading) index
    array<dcomplex, 2> Uall = dict_to_matrix(extract_U_dict2(h_int), gf_struct);

    int ibl = 0, jbl = 0, offset_col = 0, offset_row = 0;
    for (auto const &bl1 : gf_struct) {
      jbl        = 0;
      offset_col = 0;
      for (auto const &bl2 : gf_struct) {
        U(ibl, jbl).resize(bl1.second.size(), bl2.second.size());
        for (int i = 0; i < bl1.second.size(); i++)
          for (int j = 0; j < bl2.second.size(); j++) { U(ibl, jbl)(i, j) = Uall(offset_row + i, offset_col + j); } // j
        offset_col += bl2.second.size();
        jbl++;
      } // bl2
      offset_row += bl1.second.size();
      ibl++;
    } // bl1

    return U;
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
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
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
    h5_write(grp, "hartree_shift", sp.hartree_shift);
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
    h5_write(grp, "measure_F_tau", sp.measure_F_tau);
    h5_write(grp, "measure_M4_iw", sp.measure_M4_iw);
    h5_write(grp, "n_iw_M4", sp.n_iw_M4);
    h5_write(grp, "measure_M3pp_iw", sp.measure_M3pp_iw);
    h5_write(grp, "measure_M3ph_iw ", sp.measure_M3ph_iw);
    h5_write(grp, "n_iw_M3", sp.n_iw_M3);
    h5_write(grp, "measure_M2pp_tau", sp.measure_M2pp_tau);
    h5_write(grp, "measure_M2ph_tau", sp.measure_M2ph_tau);
    h5_write(grp, "measure_M2xph_tau", sp.measure_M2xph_tau);
    h5_write(grp, "n_tau_M2", sp.n_tau_M2);
    h5_write(grp, "n_iw_M2", sp.n_iw_M2);
    h5_write(grp, "nfft_buf_size", sp.nfft_buf_size);
    h5_write(grp, "post_process", sp.post_process);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
    h5_read(grp, "hartree_shift", sp.hartree_shift);
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
    h5_read(grp, "measure_F_tau", sp.measure_F_tau);
    h5_read(grp, "measure_M4_iw", sp.measure_M4_iw);
    h5_read(grp, "n_iw_M4", sp.n_iw_M4);
    h5_read(grp, "measure_M3pp_iw", sp.measure_M3pp_iw);
    h5_read(grp, "measure_M3ph_iw ", sp.measure_M3ph_iw);
    h5_read(grp, "n_iw_M3", sp.n_iw_M3);
    h5_read(grp, "measure_M2pp_tau", sp.measure_M2pp_tau);
    h5_read(grp, "measure_M2ph_tau", sp.measure_M2ph_tau);
    h5_read(grp, "measure_M2xph_tau", sp.measure_M2xph_tau);
    h5_read(grp, "n_tau_M2", sp.n_tau_M2);
    h5_read(grp, "n_iw_M2", sp.n_iw_M2);
    h5_read(grp, "nfft_buf_size", sp.nfft_buf_size);
    h5_read(grp, "post_process", sp.post_process);
  }

} // namespace triqs_ctint
