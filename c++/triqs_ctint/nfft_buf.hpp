#pragma once
#include <numeric>
#include <triqs/arrays.hpp>
#include <array>
#include <cmath>

//#define NFFT_PRECISION_DOUBLE
//#include "nfft3mp.h" // Multi Precision: Only available with nfft3.3.0+

#include "nfft3.h"

namespace triqs::utility {

  using triqs::arrays::array_view;
  using dcomplex = std::complex<double>;

  template <int Rank> struct nfft_buf_t {

    // Possible future extensions:
    //  -Bosonic Matsubaras
    //  -Move plan initialization to do_nfft for memory gain in case of large buf_size (performance penalty?)

    /// Default constructor, creates unusable buffer!
    nfft_buf_t() = default;

    /// Constructor
    nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, bool do_checks_ = false)
       : fiw_arr(std::move(fiw_arr_)), buf_size(buf_size_), beta(beta_), do_checks(do_checks_) {

      // Capture frequency extents from fiw_arr and check that they are even ( i.e. fermionic matsubaras )
      auto freq_extents = triqs::arrays::get_shape(fiw_arr).to_vector();
      std::vector<int> extents_int;
      for (int n : freq_extents) {
        if (n % 2 != 0) TRIQS_RUNTIME_ERROR << " dimension with uneven frequency count not allowed in NFFT Buffer \n";
        extents_int.push_back(n);
        common_factor *= (n / 2) % 2 ? -1 : 1; // Additional Minus sign for uneven Matsubara offset
      }

      // Init nfft_plan
      plan_ptr = std::make_unique<nfft_plan>();
      nfft_init(plan_ptr.get(), Rank, extents_int.data(), buf_size);
    }

    ~nfft_buf_t() {
      if (buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
      if (plan_ptr) nfft_finalize(plan_ptr.get());
    }

    // nfft_buffer needs to be uncopyable, because nfft_plan contains raw pointers
    nfft_buf_t(nfft_buf_t const &) = delete;
    nfft_buf_t(nfft_buf_t &&)      = default;
    nfft_buf_t &operator=(nfft_buf_t const &) = delete;
    nfft_buf_t &operator=(nfft_buf_t &&) = default;

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
      if (!plan_ptr) TRIQS_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      // Write the set of shifted and normalized tau values (i. e. x values) to the NFFT buffer and sum taus
      double tau_sum = 0.0;
      for (int r = 0; r < Rank; ++r) {
        // Note: Nfft multi-arrays are stored in flattened arrays (c-order)
        x_arr()[buf_counter * Rank + r] = tau_arr[r] / beta - 0.5; // \in [-0.5, 0.5) NOLINT
        tau_sum += tau_arr[r];                                     // Sum all tau values
      }

      // Write f(x) to nfft_plan-> The prefactor accounts for the Pi/beta offset in fermionic Matsubaras
      fx_arr()[buf_counter] = std::exp(1_j * M_PI * tau_sum / beta) * ftau; // NOLINT

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
      if (!plan_ptr) TRIQS_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      // Don't do anything if buffer is empty
      if (is_empty()) return;

      // Trivial initialization of the remaining points
      for (int i = buf_counter; i < buf_size; ++i) {
        fx_arr()[i] = 0.0; // NOLINT
        for (int r = 0; r < Rank; ++r) x_arr()[i * Rank + r] = -0.5 + double(i) / buf_size; // NOLINT
      }
      do_nfft();
      buf_counter = 0;
    }

    private:
    // Triqs array to contain the final NFFT output in matsubara frequencies
    array_view<dcomplex, Rank> fiw_arr;

    // Number of tau points for the nfft
    int buf_size;

    // Inverse temperature of fiw_arr
    double beta;

    // Switch for testing in nfft
    bool do_checks;

    // Nfft3 plan that allocates memory and performs NFFT transform
    std::unique_ptr<nfft_plan> plan_ptr;

    // Counter for elements currently in the buffer
    int buf_counter = 0;

    // Common factor in container assignment
    int common_factor = 1;

    // Get pointer to array containing x values for the NFFT transform
    double *x_arr() { return plan_ptr->x; }

    // Get pointer to array containing f(x) values for the NFFT transform
    dcomplex *fx_arr() { return reinterpret_cast<dcomplex *>(plan_ptr->f); } // NOLINT

    // Get pointer to array containing the NFFT output h(k)
    const dcomplex *fk_arr() const { return reinterpret_cast<dcomplex *>(plan_ptr->f_hat); } // NOLINT

    // Function to check whether buffer is filled
    bool is_full() const { return buf_counter >= buf_size; }

    // Function to check whether buffer is filled
    bool is_empty() const { return buf_counter == 0; }

    // Perform NFFT transform and accumulate inside fiw_arr
    void do_nfft() {

// --- NFFT Library precomputation and checks
#ifdef NFFT_OLD_API
      if (plan_ptr->nfft_flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan_ptr.get());
#else
      if (plan_ptr->flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan_ptr.get());
#endif
      if (do_checks) { // Check validity of NFFT parameters
        const char *error_str = nfft_check(plan_ptr.get());
        if (error_str != nullptr) TRIQS_RUNTIME_ERROR << "Error in NFFT module: " << error_str << "\n";
      }

      // Execute transform
      nfft_adjoint(plan_ptr.get());

      // Accumulate results in fiw_arr. Care to normalize results afterwards
      int count = 0;
      for (auto fiw_itr = fiw_arr.begin(); fiw_itr != fiw_arr.end(); ++fiw_itr) {
        int factor = common_factor * (sum(fiw_itr.indices()) % 2 ? -1 : 1);
        *fiw_itr += fk_arr()[count] * factor; // NOLINT
        ++count;
      }
    }
  };
} // namespace triqs_ctint
