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

    void execute_type1(shared_state_t<Rank> &state, nda::array_view<dcomplex, Rank> fiw_arr, nda::array<dcomplex, Rank> &fk_arr,
                       int common_factor) {
      set_pts(state, nullptr);
      check_finufft(finufft_execute(plan.get(), state.fx_arr.data(), fk_arr.data()));
      for (auto idx_tpl : fiw_arr.indices()) {
        auto idx_sum = std::apply([](auto... idx) { return (idx + ... + 0); }, idx_tpl);
        int factor   = common_factor * (idx_sum % 2 ? -1 : 1);
        std::apply(fiw_arr, idx_tpl) += std::apply(fk_arr, idx_tpl) * factor;
      }
    }

    void execute_type3(shared_state_t<Rank> &state, nda::array_view<dcomplex, 1> fiw_vec) {
      set_pts(state, &s_arr);
      check_finufft(finufft_execute(plan.get(), state.fx_arr.data(), state.fk_vec.data()));
      fiw_vec += state.fk_vec;
    }

    // Expose for calibration
    void set_pts_type3(shared_state_t<Rank> &state) { set_pts(state, &s_arr); }
    finufft_plan get_plan() const { return plan.get(); }

    private:
    finufft_plan_ptr plan;
    finufft_opts opts{};
    nda::array<double, 2> s_arr; // type3 target frequencies

    // FINUFFT expects coordinates in reverse rank order
    void set_pts(shared_state_t<Rank> &state, nda::array<double, 2> *tgt) {
      auto _ = nda::range::all;
      auto n_tgt = tgt ? state.n_targets : int64_t{0};
      auto t     = [&](int r) -> double * { return tgt ? (*tgt)(r, _).data() : nullptr; };
      if constexpr (Rank == 1)
        check_finufft(
           finufft_setpts(plan.get(), state.buf_counter, state.x_arr(0, _).data(), nullptr, nullptr, n_tgt, t(0), nullptr, nullptr));
      else if constexpr (Rank == 2)
        check_finufft(finufft_setpts(plan.get(), state.buf_counter, state.x_arr(1, _).data(), state.x_arr(0, _).data(), nullptr, n_tgt,
                                     t(1), t(0), nullptr));
      else
        check_finufft(finufft_setpts(plan.get(), state.buf_counter, state.x_arr(2, _).data(), state.x_arr(1, _).data(),
                                     state.x_arr(0, _).data(), n_tgt, t(2), t(1), t(0)));
    }
  };

} // namespace triqs::utility::nfft
