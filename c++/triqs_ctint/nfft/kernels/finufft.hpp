// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"

namespace triqs::utility::nfft {

  template <int Rank> struct kernel_finufft_t {

    kernel_finufft_t() = default;

    /// Initialize for type 1 (non-uniform tau -> uniform Matsubara grid)
    void init_type1(std::array<int64_t, Rank> const &niws, int /*buf_size*/, double tol) {
      finufft_default_opts(&opts);
      opts.nthreads         = 1;
      auto Ns               = std::vector(niws.rbegin(), niws.rend());
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(1, Rank, Ns.data(), 1, 1, tol, &raw_plan, &opts));
      plan.reset(raw_plan);
    }

    /// Initialize for type 1 + gather (non-uniform targets via bounding uniform grid)
    void init_type1_gather(std::vector<std::array<target_mf_t, Rank>> const &target_mf, int64_t n_targets, double tol) {
      // Find max |n| per dimension to determine bounding uniform grid
      std::array<int64_t, Rank> max_abs_n{};
      for (int64_t d = 0; d < n_targets; ++d)
        for (int r = 0; r < Rank; ++r) max_abs_n[r] = std::max(max_abs_n[r], std::abs(static_cast<int64_t>(target_mf[d][r].n)));

      // Grid size per dimension: covers modes [-N/2, N/2-1], need N/2 > max_abs_n
      for (int r = 0; r < Rank; ++r) gather_niws[r] = 2 * (max_abs_n[r] + 1);

      // Precompute gather indices and signs for each target point
      // fk_arr is row-major with mode k in dim r at array index (k + N_r/2)
      // Sign: (-1)^(sum of array indices) * common_factor, where common_factor = prod_r (-1)^(N_r/2)
      // This simplifies to (-1)^(sum_r(n_r + N_r/2)) * common_factor = (-1)^(sum_r n_r) * (-1)^(sum_r N_r/2) * common_factor
      // Since common_factor = (-1)^(sum_r N_r/2), sign = (-1)^(sum_r n_r) * common_factor^2 = (-1)^(sum_r n_r)
      // Wait: common_factor = prod_r [(-1)^(N_r/2)] which equals (-1)^(sum_r N_r/2).
      // sign = (-1)^(sum_r(n_r + N_r/2)) = (-1)^(sum_r n_r) * (-1)^(sum_r N_r/2) = (-1)^(sum_r n_r) * common_factor
      // But that's not right either. Let me just follow execute_type1 exactly:
      // factor = common_factor * (idx_sum % 2 ? -1 : 1) where idx_sum = sum of array indices
      // array index for target d, dim r = target_mf[d][r].n + gather_niws[r]/2
      gather_indices.resize(n_targets);
      gather_signs.resize(n_targets);
      int common_factor = 1;
      for (int r = 0; r < Rank; ++r) common_factor *= (gather_niws[r] / 2) % 2 ? -1 : 1;

      for (int64_t d = 0; d < n_targets; ++d) {
        int64_t flat    = 0;
        int64_t idx_sum = 0;
        for (int r = 0; r < Rank; ++r) {
          int64_t arr_idx = target_mf[d][r].n + gather_niws[r] / 2;
          if (r > 0) flat *= gather_niws[r];
          flat += arr_idx;
          idx_sum += arr_idx;
        }
        gather_indices[d] = flat;
        gather_signs[d]   = common_factor * (idx_sum % 2 ? -1 : 1);
      }

      // Allocate output array for type1 transform
      gather_fk_arr.resize(nda::stdutil::make_std_array<long>(gather_niws));

      // Init FINUFFT type1 plan with the bounding grid
      finufft_default_opts(&opts);
      opts.nthreads         = 1;
      auto Ns               = std::vector(gather_niws.rbegin(), gather_niws.rend());
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(1, Rank, Ns.data(), 1, 1, tol, &raw_plan, &opts));
      plan.reset(raw_plan);
    }

    /// Initialize for type 3 (non-uniform tau -> non-uniform Matsubara)
    void init_type3(std::vector<std::array<target_mf_t, Rank>> const &target_mf, int64_t n_targets, double tol) {
      s_arr.resize(Rank, n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < n_targets; ++d) s_arr(r, d) = std::imag(dcomplex(target_mf[d][r]));
      finufft_default_opts(&opts);
      opts.nthreads         = 1;
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(3, Rank, nullptr, 1, 1, tol, &raw_plan, &opts));
      plan.reset(raw_plan);
    }

    /// Initialize both type1_gather and type3 paths (for 3-way automatic dispatch)
    void init_type1_gather_and_type3(std::vector<std::array<target_mf_t, Rank>> const &target_mf, int64_t n_targets, double tol) {
      init_type1_gather(target_mf, n_targets, tol);

      // Set up type3 target frequencies
      s_arr.resize(Rank, n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < n_targets; ++d) s_arr(r, d) = std::imag(dcomplex(target_mf[d][r]));

      // Create separate type3 plan
      finufft_opts opts_t3{};
      finufft_default_opts(&opts_t3);
      opts_t3.nthreads      = 1;
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(3, Rank, nullptr, 1, 1, tol, &raw_plan, &opts_t3));
      plan_t3.reset(raw_plan);
    }

    void execute_type1(shared_state_t<Rank> &state, nda::array_view<dcomplex, Rank> fiw_arr, nda::array<dcomplex, Rank> &fk_arr,
                       int common_factor) {
      set_pts(state, nullptr, plan);
      check_finufft(finufft_execute(plan.get(), state.fx_arr.data(), fk_arr.data()));
      for (auto idx_tpl : fiw_arr.indices()) {
        auto idx_sum = std::apply([](auto... idx) { return (idx + ... + 0); }, idx_tpl);
        int factor   = common_factor * (idx_sum % 2 ? -1 : 1);
        std::apply(fiw_arr, idx_tpl) += std::apply(fk_arr, idx_tpl) * factor;
      }
    }

    void execute_type1_gather(shared_state_t<Rank> &state, nda::array_view<dcomplex, 1> fiw_vec) {
      set_pts(state, nullptr, plan);
      check_finufft(finufft_execute(plan.get(), state.fx_arr.data(), gather_fk_arr.data()));
      // Gather selected modes with precomputed signs
      auto const *fk_data = gather_fk_arr.data();
      for (int64_t d = 0; d < static_cast<int64_t>(gather_indices.size()); ++d)
        fiw_vec(d) += fk_data[gather_indices[d]] * static_cast<double>(gather_signs[d]);
    }

    void execute_type3(shared_state_t<Rank> &state, nda::array_view<dcomplex, 1> fiw_vec) {
      auto &p = plan_t3 ? plan_t3 : plan;
      set_pts(state, &s_arr, p);
      check_finufft(finufft_execute(p.get(), state.fx_arr.data(), state.fk_vec.data()));
      fiw_vec += state.fk_vec;
    }

    // Expose for calibration
    void set_pts_type1(shared_state_t<Rank> &state) { set_pts(state, nullptr, plan); }
    void set_pts_type3(shared_state_t<Rank> &state) { set_pts(state, &s_arr, plan_t3 ? plan_t3 : plan); }
    finufft_plan get_plan() const { return plan.get(); }
    finufft_plan get_plan_t3() const { return plan_t3 ? plan_t3.get() : plan.get(); }
    nda::array<dcomplex, Rank> &get_gather_fk_arr() { return gather_fk_arr; }
    std::array<int64_t, Rank> const &get_gather_niws() const { return gather_niws; }

    /// Release the type3 plan (after calibration picks type1_gather)
    void release_type3() {
      plan_t3.reset();
      s_arr = {};
    }

    /// Release the type1_gather plan and data (after calibration picks type3).
    /// Moves plan_t3 into plan so execute_type3 works with the primary plan.
    void release_type1_gather() {
      if (plan_t3) plan = std::move(plan_t3);
      gather_indices.clear();
      gather_signs.clear();
      gather_fk_arr = {};
      gather_niws   = {};
    }

    private:
    finufft_plan_ptr plan;
    finufft_plan_ptr plan_t3; // separate type3 plan (when both paths coexist)
    finufft_opts opts{};
    nda::array<double, 2> s_arr;              // type3 target frequencies
    std::array<int64_t, Rank> gather_niws{};  // type1_gather grid sizes
    std::vector<int64_t> gather_indices;      // type1_gather: flat index per target
    std::vector<int> gather_signs;            // type1_gather: sign factor per target
    nda::array<dcomplex, Rank> gather_fk_arr; // type1_gather: FINUFFT output buffer

    // FINUFFT expects coordinates in reverse rank order
    void set_pts(shared_state_t<Rank> &state, nda::array<double, 2> *tgt, finufft_plan_ptr const &p) {
      auto _ = nda::range::all;
      auto n_tgt = tgt ? state.n_targets : int64_t{0};
      auto t     = [&](int r) -> double * { return tgt ? (*tgt)(r, _).data() : nullptr; };
      if constexpr (Rank == 1)
        check_finufft(finufft_setpts(p.get(), state.buf_counter, state.x_arr(0, _).data(), nullptr, nullptr, n_tgt, t(0), nullptr, nullptr));
      else if constexpr (Rank == 2)
        check_finufft(
           finufft_setpts(p.get(), state.buf_counter, state.x_arr(1, _).data(), state.x_arr(0, _).data(), nullptr, n_tgt, t(1), t(0), nullptr));
      else
        check_finufft(finufft_setpts(p.get(), state.buf_counter, state.x_arr(2, _).data(), state.x_arr(1, _).data(), state.x_arr(0, _).data(), n_tgt,
                                     t(2), t(1), t(0)));
    }
  };

} // namespace triqs::utility::nfft
