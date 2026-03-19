// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "common.hpp"
#include "kernels.hpp"
#include "calibration.hpp"
#include <optional>
#include <variant>

namespace triqs::utility::nfft {

  template <int Rank> struct buffer_t {

    static_assert(Rank >= 1 and Rank <= 3, "buffer_t only supports Rank 1, 2, and 3");

    template <int TolDigits>
    using finufft_kernel_t = kernel_finufft_t<Rank, TolDigits>;

    using finufft_variant_t = std::variant<std::monostate, finufft_kernel_t<6>, finufft_kernel_t<8>, finufft_kernel_t<10>, finufft_kernel_t<12>>;

    using do_nfft_fn_t = void (*)(buffer_t &);

    buffer_t() = default;

    /// Type 1: non-uniform tau -> uniform Matsubara grid (with automatic NAF/FINUFFT dispatch)
    buffer_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, double tol_ = 1e-8)
       : fiw_arr(std::move(fiw_arr_)),
         niws(nda::stdutil::make_std_array<int64_t>(fiw_arr.shape())),
         fk_arr(fiw_arr.shape()),
         tol(tol_) {

      state_.buf_size = buf_size_;
      state_.beta     = beta_;
      state_.x_arr.resize(Rank, buf_size_);
      state_.fx_arr.resize(buf_size_);

      for (int n : niws) {
        if (n % 2 != 0) NDA_RUNTIME_ERROR << " dimension with uneven frequency count not allowed in NFFT Buffer \n";
        common_factor *= (n / 2) % 2 ? -1 : 1;
      }

      dispatch_tol_digits([&]<int TolDigits>() { init_type1_impl<TolDigits>(beta_); });
    }

    /// Non-uniform target constructor: automatic dispatch, FINUFFT type3, or direct DFT
    buffer_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<std::array<target_mf_t, Rank>> target_mf_, int buf_size_,
             double tol_ = 1e-8, type_t type = type_t::automatic)
       : type_(type), fiw_vec(std::move(fiw_vec_)), tol(tol_) {

      state_.buf_size  = buf_size_;
      state_.n_targets = static_cast<int64_t>(target_mf_.size());
      state_.x_arr.resize(Rank, buf_size_);
      state_.fx_arr.resize(buf_size_);
      state_.fk_vec.resize(state_.n_targets);

      dispatch_tol_digits([&]<int TolDigits>() { init_nonuniform_impl<TolDigits>(target_mf_); });
    }

    /// Convenience constructor for Rank=1: accepts vector of matsubara_freq directly
    buffer_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<target_mf_t> const &target_mf_, int buf_size_, double tol_ = 1e-8,
             type_t type = type_t::automatic)
      requires(Rank == 1)
       : buffer_t(
            std::move(fiw_vec_),
            [&] {
              std::vector<std::array<target_mf_t, 1>> result;
              result.reserve(target_mf_.size());
              for (auto const &mf : target_mf_) result.push_back({mf});
              return result;
            }(),
            buf_size_, tol_, type) {}

    ~buffer_t() {
      if (state_.buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
    }

    buffer_t(buffer_t const &)            = delete;
    buffer_t &operator=(buffer_t const &) = delete;
    buffer_t(buffer_t &&)                 = default;
    buffer_t &operator=(buffer_t &&rhs) noexcept {
      // array_view::operator= does deep copy; we need rebind via destroy + placement new
      if (this != &rhs) {
        std::destroy_at(this);
        std::construct_at(this, std::move(rhs));
      }
      return *this;
    }

    void rebind(array_view<dcomplex, Rank> new_fiw_arr) {
      flush();
      TRIQS_ASSERT((new_fiw_arr.shape() == fiw_arr.shape() or fiw_arr.empty())
                   and " Nfft Buffer: Rebind to array of different shape not allowed ");
      fiw_arr.rebind(new_fiw_arr);
    }

    [[gnu::always_inline]] inline void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {
      if (state_.x_arr.empty()) [[unlikely]] throw_uninitialized_buffer_error();

      using x_arr_t = std::remove_reference_t<decltype(state_.x_arr)>;
      static_assert(x_arr_t::is_stride_order_C());

      int const idx = state_.buf_counter;

      if (type_ == type_t::type1 && direct_vs_finufft_threshold_ == 0) {
        double tau_sum = 0.0;
        double const pi_over_beta     = M_PI / state_.beta;
        double const two_pi_over_beta = 2.0 * pi_over_beta;
        poet::static_for<Rank>([&](auto r) {
          double const tau = tau_arr[r];
          tau_sum += tau;
          state_.x_arr(r, idx) = std::fma(two_pi_over_beta, tau, -M_PI);
        });
        state_.fx_arr[idx] = cis<12>(pi_over_beta * tau_sum) * ftau;
      } else {
        poet::static_for<Rank>([&](auto r) { state_.x_arr(r, idx) = tau_arr[r]; });
        state_.fx_arr[idx] = ftau;
      }

      if (++state_.buf_counter >= state_.buf_size) flush();
    }

    void flush() {
      if (state_.x_arr.empty()) [[unlikely]] throw_uninitialized_buffer_error();
      if (state_.buf_counter == 0) return;
      if (not do_nfft_fn_) [[unlikely]] throw_uninitialized_backend_error();
      do_nfft_fn_(*this);
      state_.buf_counter = 0;
    }

    private:
    static constexpr int64_t max_type1_dispatch_targets = 100'000;

    [[gnu::cold, noreturn]] static void throw_uninitialized_buffer_error() {
      NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";
    }

    [[gnu::cold, noreturn]] static void throw_uninitialized_backend_error() {
      NDA_RUNTIME_ERROR << " Nfft Buffer backend was not initialized\n";
    }

    type_t type_ = type_t::type1;
    // Automatic non-uniform dispatch first chooses a direct kernel family
    // (direct_type1 vs NAF/chain), then chooses the FINUFFT path
    // (type3 vs type1_gather) for buffer sizes above the direct crossover.
    bool selected_sparse_direct_is_chain_      = false;
    bool use_type3_below_switch_threshold_     = false;
    shared_state_t<Rank> state_;

    nda::array_view<dcomplex, Rank> fiw_arr;
    nda::array_view<dcomplex, 1> fiw_vec;
    std::array<int64_t, Rank> niws{};
    nda::array<dcomplex, Rank> fk_arr;
    int common_factor = 1;
    double tol        = 1e-8;

    finufft_variant_t finufft_kernel_;
    std::optional<kernel_direct_type1_t<Rank>> direct_type1_kernel_;
    std::optional<kernel_chain_t<Rank>> chain_kernel_;
    std::optional<kernel_naf_t<Rank>> naf_kernel_;

    int direct_vs_finufft_threshold_ = 0;
    int finufft_path_switch_threshold_ = 0;
    do_nfft_fn_t do_nfft_fn_          = nullptr;

    template <int TolDigits> auto &emplace_finufft() {
      return finufft_kernel_.template emplace<finufft_kernel_t<TolDigits>>();
    }

    template <int TolDigits> auto &get_finufft() { return std::get<finufft_kernel_t<TolDigits>>(finufft_kernel_); }

    template <typename Builder> void dispatch_tol_digits(Builder &&builder) {
      auto params = std::make_tuple(poet::DispatchParam<supported_tol_digits_t>{bucket_tol_digits(tol)});
      poet::dispatch(poet::throw_t, std::forward<Builder>(builder), params);
    }

    template <int TolDigits> void init_type1_impl(double beta_) {
      auto &finufft = emplace_finufft<TolDigits>();
      finufft.init_type1(niws, state_.buf_size, tol);

      int64_t n_targets = 1;
      for (auto n : niws) n_targets *= n;
      if (n_targets <= max_type1_dispatch_targets) {
        state_.n_targets = n_targets;
        state_.fk_vec.resize(n_targets);

        auto target_mf = build_uniform_target_mf(beta_);
        state_.init_direct_common(target_mf);
        direct_type1_kernel_.emplace(state_, target_mf);
        direct_vs_finufft_threshold_ = calibrate_dispatch_type1<12>(state_, *direct_type1_kernel_, finufft, fk_arr, common_factor);
      }

      do_nfft_fn_ = &buffer_t::template do_nfft_impl<TolDigits>;
    }

    template <int TolDigits> void init_nonuniform_impl(std::vector<std::array<target_mf_t, Rank>> const &target_mf) {
      if (type_ == type_t::type3) {
        emplace_finufft<TolDigits>().init_type3(target_mf, state_.n_targets, tol);

      } else if (type_ == type_t::direct_type1) {
        state_.init_direct_common(target_mf);
        direct_type1_kernel_.emplace(state_, target_mf);

      } else if (type_ == type_t::direct_chain) {
        state_.init_direct_common(target_mf);
        chain_kernel_.emplace(state_);

      } else if (type_ == type_t::direct_type3) {
        state_.init_direct_common(target_mf);
        naf_kernel_.emplace(state_, state_.buf_size);

      } else if (type_ == type_t::type1_gather) {
        state_.beta = target_mf[0][0].beta;
        emplace_finufft<TolDigits>().init_type1_gather(target_mf, state_.n_targets, tol);

      } else if (type_ == type_t::automatic) {
        auto &finufft = emplace_finufft<TolDigits>();
        finufft.init_type1_gather_and_type3(target_mf, state_.n_targets, tol);
        state_.init_direct_common(target_mf);
        naf_kernel_.emplace(state_, state_.buf_size);
        direct_type1_kernel_.emplace(state_, target_mf);
        chain_kernel_.emplace(state_);

        auto &naf_kernel   = *naf_kernel_;
        auto &chain_kernel = *chain_kernel_;

        // Automatic mode is chosen in two stages:
        // 1. Compare the two sparse-direct candidates (`chain` and `naf`) against `direct_type1`.
        //    Keep the sparse family with the larger direct-vs-FINUFFT crossover.
        // 2. For that winning family, honor the calibration result that says whether the best
        //    direct kernel is actually `direct_type1` or the sparse kernel itself.
        auto const chain_dispatch_plan = calibrate_dispatch_nonuniform<TolDigits>(state_, *direct_type1_kernel_, chain_kernel, finufft);
        auto const naf_dispatch_plan   = calibrate_dispatch_nonuniform<TolDigits>(state_, *direct_type1_kernel_, naf_kernel, finufft);

        auto const direct_path_cost = [](auto const &plan) { return plan.direct_path_lo_time + plan.direct_path_hi_time; };
        double const chain_direct_cost = direct_path_cost(chain_dispatch_plan);
        double const naf_direct_cost   = direct_path_cost(naf_dispatch_plan);

        if (!chain_dispatch_plan.select_direct_type1 && !naf_dispatch_plan.select_direct_type1) {
          selected_sparse_direct_is_chain_ =
             (chain_direct_cost < naf_direct_cost) ||
             ((chain_direct_cost == naf_direct_cost) &&
              (chain_dispatch_plan.direct_vs_finufft_threshold >= naf_dispatch_plan.direct_vs_finufft_threshold));
        } else {
          selected_sparse_direct_is_chain_ =
             (chain_dispatch_plan.direct_vs_finufft_threshold > naf_dispatch_plan.direct_vs_finufft_threshold) ||
             ((chain_dispatch_plan.direct_vs_finufft_threshold == naf_dispatch_plan.direct_vs_finufft_threshold) &&
              (chain_direct_cost <= naf_direct_cost));
        }
        auto const &selected_dispatch_plan = selected_sparse_direct_is_chain_ ? chain_dispatch_plan : naf_dispatch_plan;
        bool const select_direct_type1     = selected_dispatch_plan.select_direct_type1;
        direct_vs_finufft_threshold_       = selected_dispatch_plan.direct_vs_finufft_threshold;
        finufft_path_switch_threshold_     = selected_dispatch_plan.finufft_path_switch_threshold;
        use_type3_below_switch_threshold_  = selected_dispatch_plan.use_type3_below_switch_threshold;

        // Cleanup after the selected automatic dispatch policy is fully determined.
        prune_unselected_direct_kernels(select_direct_type1);

        // Release unused FINUFFT plans only if one path dominates over the full buffer range.
        if (finufft_path_switch_threshold_ <= 0) {
          if (use_type3_below_switch_threshold_)
            finufft.release_type3();
          else
            finufft.release_type1_gather();
        } else if (finufft_path_switch_threshold_ > state_.buf_size) {
          if (use_type3_below_switch_threshold_)
            finufft.release_type1_gather();
          else
            finufft.release_type3();
        }

      } else {
        NDA_RUNTIME_ERROR << "buffer_t: unsupported type_t for non-uniform target constructor\n";
      }

      do_nfft_fn_ = &buffer_t::template do_nfft_impl<TolDigits>;
    }

    int64_t estimate_naf_multiplies() const {
      int64_t naf_estimate = 0;
      for (int r = 0; r < Rank; ++r) {
        std::vector<unsigned long> unique_exp;
        unique_exp.reserve(state_.n_targets);
        for (int64_t d = 0; d < state_.n_targets; ++d) unique_exp.push_back(odd_exponent_abs(state_.target_n(r, d)));
        std::sort(unique_exp.begin(), unique_exp.end());
        unique_exp.erase(std::unique(unique_exp.begin(), unique_exp.end()), unique_exp.end());
        if (unique_exp.empty()) continue;

        naf_estimate += static_cast<int64_t>(std::bit_width(unique_exp.back()) - 1);
        for (unsigned long exp : unique_exp) naf_estimate += static_cast<int64_t>(compute_naf(exp).size()) - 1;
      }
      return naf_estimate;
    }

    void prune_unselected_direct_kernels(bool select_direct_type1) {
      if (select_direct_type1) {
        naf_kernel_.reset();
        chain_kernel_.reset();
        return;
      }

      direct_type1_kernel_.reset();
      if (selected_sparse_direct_is_chain_)
        naf_kernel_.reset();
      else
        chain_kernel_.reset();
    }

    bool should_use_type3_for_n(int n) const {
      if (finufft_path_switch_threshold_ <= 0) return !use_type3_below_switch_threshold_;
      if (finufft_path_switch_threshold_ > state_.buf_size) return use_type3_below_switch_threshold_;
      return use_type3_below_switch_threshold_ ? (n < finufft_path_switch_threshold_) : (n >= finufft_path_switch_threshold_);
    }

    template <int TolDigits> static void do_nfft_impl(buffer_t &self) {
      if (self.type_ == type_t::type1) {
        auto &finufft = self.template get_finufft<TolDigits>();
        if (self.state_.buf_counter < self.direct_vs_finufft_threshold_)
          self.template run_direct_type1<TolDigits>(*self.direct_type1_kernel_);
        else {
          if (self.direct_vs_finufft_threshold_ > 0) self.template prepare_type1_coords<12>();
          finufft.execute_type1(self.state_, self.fiw_arr, self.fk_arr, self.common_factor);
        }
      } else if (self.type_ == type_t::type1_gather) {
        auto &finufft = self.template get_finufft<TolDigits>();
        self.template prepare_type1_coords<TolDigits>();
        finufft.execute_type1_gather(self.state_, self.fiw_vec);
      } else if (self.type_ == type_t::automatic) {
        if (self.state_.buf_counter < self.direct_vs_finufft_threshold_) {
          if (self.direct_type1_kernel_)
            self.template run_direct<TolDigits>(*self.direct_type1_kernel_);
          else if (self.selected_sparse_direct_is_chain_ && self.chain_kernel_)
            self.template run_direct<TolDigits>(*self.chain_kernel_);
          else
            self.template run_direct<TolDigits>(*self.naf_kernel_);
        } else {
          auto &finufft = self.template get_finufft<TolDigits>();
          if (self.should_use_type3_for_n(self.state_.buf_counter)) {
            finufft.execute_type3(self.state_, self.fiw_vec);
          } else {
            self.template prepare_type1_coords<TolDigits>();
            finufft.execute_type1_gather(self.state_, self.fiw_vec);
          }
        }
      } else if (self.type_ == type_t::type3) {
        auto &finufft = self.template get_finufft<TolDigits>();
        finufft.execute_type3(self.state_, self.fiw_vec);
      } else if (self.type_ == type_t::direct_type1)
        self.template run_direct<TolDigits>(*self.direct_type1_kernel_);
      else if (self.type_ == type_t::direct_type3)
        self.template run_direct<TolDigits>(*self.naf_kernel_);
      else if (self.type_ == type_t::direct_chain)
        self.template run_direct<TolDigits>(*self.chain_kernel_);
      else
        self.template run_direct_type1<TolDigits>(*self.direct_type1_kernel_);
    }

    template <int TolDigits> void run_direct(kernel_direct_type1_t<Rank> &kernel) {
      state_.fk_vec = 0;
      kernel.template execute<12>(state_);
      fiw_vec += state_.fk_vec;
    }

    template <int TolDigits, typename Kernel> void run_direct(Kernel &kernel) {
      // Pad buffer counter to SIMD boundary to eliminate scalar tail loops
      int const buf_counter_padded = std::min(shared_state_t<Rank>::round_up_simd(state_.buf_counter), state_.buf_size);

      for (int j = state_.buf_counter; j < buf_counter_padded; ++j) {
        state_.fx_arr[j] = dcomplex{0.0, 0.0};
        for (int r = 0; r < Rank; ++r) state_.x_arr(r, j) = 0.0;
      }

      state_.fk_vec = 0;
      kernel.template execute<TolDigits>(state_);
      fiw_vec += state_.fk_vec;
    }

    // Direct type1 path: compute into flat fk_vec, then scatter to fiw_arr
    template <int TolDigits> void run_direct_type1(kernel_direct_type1_t<Rank> &kernel) {
      state_.fk_vec = 0;
      kernel.template execute<12>(state_);
      scatter_to_arr();
    }

    void scatter_to_arr() {
      if constexpr (Rank == 1) {
        fiw_arr += state_.fk_vec;
      } else {
        fiw_arr += nda::reshape(state_.fk_vec, fiw_arr.shape());
      }
    }

    template <int TolDigits> void prepare_type1_coords() { apply_type1_coord_transform<TolDigits>(state_, state_.buf_counter); }

    std::vector<std::array<target_mf_t, Rank>> build_uniform_target_mf(double beta) const {
      std::vector<std::array<target_mf_t, Rank>> target_mf;
      target_mf.reserve(state_.n_targets);

      auto recurse = [&](this auto &&self, std::array<target_mf_t, Rank> &mf, int r) -> void {
        if (r == Rank) {
          target_mf.push_back(mf);
          return;
        }
        for (int64_t k = 0; k < niws[r]; ++k) {
          mf[r] = target_mf_t(static_cast<int>(k - niws[r] / 2), beta, mesh::Fermion);
          self(mf, r + 1);
        }
      };
      std::array<target_mf_t, Rank> mf{};
      recurse(mf, 0);
      return target_mf;
    }
  };

  template <int Rank> buffer_t(nda::array_view<dcomplex, Rank>, int, double, double) -> buffer_t<Rank>;

  template <std::size_t N>
  buffer_t(nda::array_view<dcomplex, 1>, std::vector<std::array<target_mf_t, N>>, int, double, type_t) -> buffer_t<static_cast<int>(N)>;

  template <std::size_t N>
  buffer_t(nda::array_view<dcomplex, 1>, std::vector<std::array<target_mf_t, N>>, int, double) -> buffer_t<static_cast<int>(N)>;

  template <std::size_t N>
  buffer_t(nda::array_view<dcomplex, 1>, std::vector<std::array<target_mf_t, N>>, int) -> buffer_t<static_cast<int>(N)>;

  buffer_t(nda::array_view<dcomplex, 1>, std::vector<target_mf_t> const &, int, double, type_t) -> buffer_t<1>;

  buffer_t(nda::array_view<dcomplex, 1>, std::vector<target_mf_t> const &, int, double) -> buffer_t<1>;

  buffer_t(nda::array_view<dcomplex, 1>, std::vector<target_mf_t> const &, int) -> buffer_t<1>;

} // namespace triqs::utility::nfft
