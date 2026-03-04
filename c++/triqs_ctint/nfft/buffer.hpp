// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "common.hpp"
#include "kernels.hpp"
#include "calibration.hpp"
#include <optional>

namespace triqs::utility::nfft {

  template <int Rank> struct buffer_t {

    static_assert(Rank >= 1 and Rank <= 3, "buffer_t only supports Rank 1, 2, and 3");

    buffer_t() = default;

    /// Type 1: non-uniform tau -> uniform Matsubara grid
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

      finufft_kernel_.emplace();
      finufft_kernel_->init_type1(niws, buf_size_, tol);
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

      if (type == type_t::type3) {
        finufft_kernel_.emplace();
        finufft_kernel_->init_type3(target_mf_, state_.n_targets, tol);

      } else if (type == type_t::direct_type1) {
        state_.init_direct_common(target_mf_);
        direct_type1_kernel_.emplace(state_, target_mf_);

      } else if (type == type_t::direct_bitwise) {
        state_.init_direct_common(target_mf_);
        bitwise_kernel_.emplace(state_);

      } else if (type == type_t::direct_prime) {
        state_.init_direct_common(target_mf_);
        prime_kernel_.emplace(state_);

      } else if (type == type_t::direct_type3) {
        state_.init_direct_common(target_mf_);
        naf_kernel_.emplace(state_, buf_size_);

      } else if (type == type_t::automatic) {
        finufft_kernel_.emplace();
        finufft_kernel_->init_type3(target_mf_, state_.n_targets, tol);
        state_.init_direct_common(target_mf_);
        naf_kernel_.emplace(state_, buf_size_);
        dispatch_buf_threshold = calibrate_dispatch(state_, *naf_kernel_, *finufft_kernel_);

      } else {
        NDA_RUNTIME_ERROR << "buffer_t: unsupported type_t for non-uniform target constructor\n";
      }
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

    void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {
      if (state_.x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      if (type_ == type_t::type1) {
        double tau_sum = 0.0;
        for (int r = 0; r < Rank; ++r) {
          state_.x_arr(r, state_.buf_counter) = 2 * M_PI * (tau_arr[r] / state_.beta - 0.5);
          tau_sum += tau_arr[r];
        }
        state_.fx_arr[state_.buf_counter] = std::exp(dcomplex(0, M_PI * tau_sum / state_.beta)) * ftau;
      } else {
        for (int r = 0; r < Rank; ++r) state_.x_arr(r, state_.buf_counter) = tau_arr[r];
        state_.fx_arr[state_.buf_counter] = ftau;
      }

      if (++state_.buf_counter >= state_.buf_size) {
        do_nfft();
        state_.buf_counter = 0;
      }
    }

    void flush() {
      if (state_.x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";
      if (state_.buf_counter == 0) return;
      do_nfft();
      state_.buf_counter = 0;
    }

    private:
    type_t type_ = type_t::type1;
    shared_state_t<Rank> state_;

    // Type1-specific
    nda::array_view<dcomplex, Rank> fiw_arr;
    nda::array_view<dcomplex, 1> fiw_vec;
    std::array<int64_t, Rank> niws{};
    nda::array<dcomplex, Rank> fk_arr;
    int common_factor = 1;
    double tol        = 1e-8;

    // Kernels (only relevant ones initialized)
    std::optional<kernel_finufft_t<Rank>> finufft_kernel_;
    std::optional<kernel_direct_type1_t<Rank>> direct_type1_kernel_;
    std::optional<kernel_bitwise_t<Rank>> bitwise_kernel_;
    std::optional<kernel_prime_t<Rank>> prime_kernel_;
    std::optional<kernel_naf_t<Rank>> naf_kernel_;

    int dispatch_buf_threshold = 0;

    void do_nfft() {
      if (type_ == type_t::type1)
        finufft_kernel_->execute_type1(state_, fiw_arr, fk_arr, common_factor);
      else if (type_ == type_t::automatic) {
        if (state_.buf_counter < dispatch_buf_threshold)
          run_direct(*naf_kernel_);
        else
          finufft_kernel_->execute_type3(state_, fiw_vec);
      } else if (type_ == type_t::type3)
        finufft_kernel_->execute_type3(state_, fiw_vec);
      else if (type_ == type_t::direct_type3)
        run_direct(*naf_kernel_);
      else if (type_ == type_t::direct_bitwise)
        run_direct(*bitwise_kernel_);
      else if (type_ == type_t::direct_prime)
        run_direct(*prime_kernel_);
      else
        run_direct(*direct_type1_kernel_);
    }

    template <typename Kernel> void run_direct(Kernel &kernel) {
      state_.fk_vec = 0;
      kernel.execute(state_);
      fiw_vec += state_.fk_vec;
    }
  };

  // Deduction guides
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
