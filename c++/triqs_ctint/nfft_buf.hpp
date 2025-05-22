#pragma once
#include <numeric>
#include <nda/nda.hpp>
#include <array>
#include <cmath>

#include "finufft.h"

#define CHECK_ERROR(err)                                                                                                                             \
  if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n"

namespace triqs::utility {

  using nda::array_view;
  using nda::stdutil::sum;
  using dcomplex = std::complex<double>;

  template <int Rank> struct nfft_buf_t {

    // Possible future extensions:
    //  -Bosonic Matsubaras

    /// Default constructor, creates unusable buffer!
    nfft_buf_t() = default;

    /// Constructor
    nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_)
       : fiw_arr(std::move(fiw_arr_)),
         niws(nda::stdutil::make_std_array<int64_t>(fiw_arr.shape())),
         buf_size(buf_size_),
         beta(beta_),
         x_arr(Rank, buf_size),
         fx_arr(buf_size),
         fk_arr(fiw_arr.shape()) {

      // Capture frequency extents from fiw_arr and check that they are even ( i.e. fermionic matsubaras )
      for (int n : niws) {
        if (n % 2 != 0) NDA_RUNTIME_ERROR << " dimension with uneven frequency count not allowed in NFFT Buffer \n";
        common_factor *= (n / 2) % 2 ? -1 : 1; // Additional Minus sign for uneven Matsubara offset
      }

      // Init nfft_plan
      finufft_default_opts(&opts); // set default opts (must start with this)
      opts.nthreads = 1;           // enforce single-thread
      //opts.debug    = 1;                                       // print diagnostics or 2 prints some information about what finufft is doing
      auto Ns = std::vector(niws.rbegin(), niws.rend()); // Reverse order for FINUFFT
      CHECK_ERROR(finufft_makeplan(/*type =*/1, Rank, Ns.data(), /*iflag=*/1, /*ntrans =*/1, tol, &plan, &opts));
    }

    ~nfft_buf_t() {
      if (buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
      if (plan) finufft_destroy(plan);
    }

    // nfft_buffer needs to be uncopyable, because nfft_plan contains raw pointers
    nfft_buf_t(nfft_buf_t const &)            = delete;
    nfft_buf_t(nfft_buf_t &&)                 = default;
    nfft_buf_t &operator=(nfft_buf_t const &) = delete;
    nfft_buf_t &operator=(nfft_buf_t &&rhs) noexcept {
      fiw_arr.rebind(rhs.fiw_arr);
      niws = rhs.niws;
      std::swap(plan, rhs.plan);
      buf_size      = rhs.buf_size;
      beta          = rhs.beta;
      buf_counter   = rhs.buf_counter;
      common_factor = rhs.common_factor;
      opts          = std::move(rhs.opts);
      x_arr         = std::move(rhs.x_arr);
      fx_arr        = std::move(rhs.fx_arr);
      fk_arr        = std::move(rhs.fk_arr);
      return *this;
    }

    /// Rebind nfft buffer to new accumulation container of same shape
    void rebind(array_view<dcomplex, Rank> new_fiw_arr) {
      flush();
      TRIQS_ASSERT((get_shape(new_fiw_arr) == get_shape(fiw_arr) or get_shape(fiw_arr) == get_shape(array_view<dcomplex, Rank>{}))
                   and " Nfft Buffer: Rebind to array of different shape not allowed ");
      fiw_arr.rebind(new_fiw_arr);
    }

    /// Insert tau-vector {tau_1, tau_2, ... } \in [0,\beta)^Rank and corresponding f(tau) into the NFFT buffer
    void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {

      // Check if buffer has been properly initialized
      if (x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      // Write the set of shifted and normalized tau values (i. e. x values) to the NFFT buffer and sum taus
      double tau_sum = 0.0;
      for (int r = 0; r < Rank; ++r) {
        // Note: Nfft multi-arrays are stored in flattened arrays (c-order)
        x_arr(r, buf_counter) = 2 * M_PI * (tau_arr[r] / beta - 0.5); // \in [-PI, PI)
        tau_sum += tau_arr[r];                                        // Sum all tau values
      }

      // Write f(x), The prefactor accounts for the Pi/beta offset in fermionic Matsubaras
      fx_arr[buf_counter] = std::exp(dcomplex(0, M_PI * tau_sum / beta)) * ftau;

      ++buf_counter;

      // If buffer is full, perform transform
      if (is_full()) {
        do_nfft();
        buf_counter = 0;
      }
    }

    /// Flush contents of the nfft buffer
    void flush() {

      // Check if buffer has been properly initialized
      if (x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      // Don't do anything if buffer is empty
      if (is_empty()) return;

      // Execute the transform
      do_nfft();
      buf_counter = 0;
    }

    private:
    // Triqs array to contain the final NFFT output in matsubara frequencies
    nda::array_view<dcomplex, Rank> fiw_arr;

    // Dimensions of the output array
    std::array<int64_t, Rank> niws;

    // Finufft plan
    finufft_plan plan{nullptr};

    // Number of tau points for the nfft
    int buf_size;

    // Inverse temperature of fiw_arr
    double beta;

    // Counter for elements currently in the buffer
    int buf_counter = 0;

    // Common factor in container assignment
    int common_factor = 1;

    // FINUFFT options struct
    finufft_opts opts{};

    // Array containing x values for the NFFT transform
    nda::array<double, 2> x_arr;

    // Array containing f(x) values for the NFFT transform
    nda::vector<dcomplex> fx_arr;

    // Array containing the NFFT output h(k)
    nda::array<dcomplex, Rank> fk_arr;

    // Tolerance for the transformation
    double tol = 1e-14;

    // Function to check whether buffer is filled
    bool is_full() const { return buf_counter >= buf_size; }

    // Function to check whether buffer is filled
    bool is_empty() const { return buf_counter == 0; }

    // Perform NFFT transform and accumulate inside fiw_arr
    void do_nfft() {

      static_assert(Rank < 4, "NFFT Implemented only for Ranks 1, 2 and 3");

      // Execute transform
      auto _ = nda::range::all;
      if constexpr (Rank == 1) {
        CHECK_ERROR(finufft_setpts(plan, buf_counter, x_arr(0, _).data(), nullptr, nullptr, 0, nullptr, nullptr, nullptr));
      } else if constexpr (Rank == 2) {
        CHECK_ERROR(finufft_setpts(plan, buf_counter, x_arr(1, _).data(), x_arr(0, _).data(), nullptr, 0, nullptr, nullptr, nullptr));
      } else { // Rank == 3
        CHECK_ERROR(finufft_setpts(plan, buf_counter, x_arr(2, _).data(), x_arr(1, _).data(), x_arr(0, _).data(), 0, nullptr, nullptr, nullptr));
      }
      CHECK_ERROR(finufft_execute(plan, fx_arr.data(), fk_arr.data()));

      // Accumulate results in fiw_arr. Care to normalize results afterwards
      for (auto idx_tpl : fiw_arr.indices()) {
        auto idx_sum = std::apply([](auto... idx) { return (idx + ... + 0); }, idx_tpl);
        int factor   = common_factor * (idx_sum % 2 ? -1 : 1);
        std::apply(fiw_arr, idx_tpl) += std::apply(fk_arr, idx_tpl) * factor;
      }
    }
  };
} // namespace triqs::utility
