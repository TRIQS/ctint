#pragma once
#include <numeric>
#include <triqs/arrays.hpp>
#include <array>
#include <cmath>

//#define NFFT_PRECISION_DOUBLE
//#include "nfft3mp.h" // Multi Precision: Only available with nfft3.3.0+

#include "nfft3.h"

using dcomplex = std::complex<double>;

namespace triqs::utility {

  using triqs::arrays::array_view;

  template <int Rank> struct nfft_buf_t {

    // Possible future extensions:
    //  -Bosonic Matsubaras
    //  -Possibly remove tau shift and exponential in push_back

    /// Constructor
    nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, bool do_checks_ = false)
       : fiw_arr(fiw_arr_), buf_size(buf_size_), beta(beta_), do_checks(do_checks_) {

      // Capture frequency extents from fiw_arr and check that they are even ( i.e. fermionic matsubaras )
      auto freq_extents = triqs::arrays::get_shape(fiw_arr).to_vector();
      std::vector<int> extents_int;
      for (int n : freq_extents) {
        if (n % 2 != 0) TRIQS_RUNTIME_ERROR << " dimension with uneven frequency count not allowed in nfft_buf_t \n";
        extents_int.push_back(n);
      }

      // Init nfft_plan
      plan_ptr = std::make_unique<nfft_plan>();
      int m    = 6; // Truncation order for the window functions
      int n[Rank];  // Size of fftw array. Each dimension should be 2-4 times larger than nfft freq extents (oversampling factor)
      for (int i = 0; i < Rank; i++) {
        int power = log2(extents_int.data()[i]);
        n[i]      = pow(2.0, 1.0 + ceil(power));
      }
      unsigned nfft_flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT /*|FFT_OUT_OF_PLACE*/ | NFFT_SORT_NODES;
      unsigned fftw_flags = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;
      nfft_init_guru(plan_ptr.get(), Rank, extents_int.data(), buf_size, n, m, nfft_flags, fftw_flags);
    }

    ~nfft_buf_t() {
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
      using triqs::arrays::get_shape;
      TRIQS_ASSERT(get_shape(new_fiw_arr) == get_shape(fiw_arr) && " Nfft Buffer: Rebind to array of different shape not allowed ");
      fiw_arr.rebind(new_fiw_arr);
    }

    /// Insert tau-vector {tau_1, tau_2, ... } \in [0,\beta)^Rank and corresponding f(tau) into the NFFT buffer
    void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {

      // Write the set of shifted and normalized tau values (i. e. x values) to the NFFT buffer and sum taus
      double tau_sum = 0.0;
      for (int r = 0; r < Rank; ++r) {
        // Note: Nfft multi-arrays are stored in flattened arrays (c-order)
        x_arr()[buf_counter * Rank + r] = tau_arr[r] / beta - 0.5; // \in [-0.5, 0.5)
        tau_sum += tau_arr[r];                                     // Sum all tau values
      }

      // Write f(x) to nfft_plan-> The prefactor accounts for the Pi/beta offset in fermionic Matsubaras
      fx_arr()[buf_counter] = std::exp(1_j * M_PI * tau_sum / beta) * ftau;

      ++buf_counter;

      // If buffer is full, perform transform
      if (is_full()) {
        do_nfft();
        buf_counter = 0;
      }
    }

    /// Flush contents of the nfft buffer
    void flush() {
      // Trivial initialization of the remaining points
      for (int i = buf_counter; i < buf_size; ++i) {
        fx_arr()[i] = 0.0;
        for (int r = 0; r < Rank; ++r) x_arr()[i * Rank + r] = -0.5 + double(i) / buf_size;
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

    // Counter for elements currently in buffer
    int buf_counter = 0;

    // Get pointer to array containing x values for the NFFT transform
    double *x_arr() { return plan_ptr->x; }

    // Get pointer to array containing f(x) values for the NFFT transform
    dcomplex *fx_arr() { return reinterpret_cast<dcomplex *>(plan_ptr->f); }

    // Get pointer to array containing the NFFT output h(k)
    const dcomplex *fk_arr() const { return reinterpret_cast<dcomplex *>(plan_ptr->f_hat); }

    // Function to check whether buffer is filled
    bool is_full() const { return (buf_counter >= buf_size); }

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
        if (error_str != 0) TRIQS_RUNTIME_ERROR << "Error in NFFT module: " << error_str << "\n";
      }

      // Execute transform
      nfft_adjoint(plan_ptr.get());

      // Accumulate results in fiw_arr. Care to normalize results afterwards
      int count = 0;
      for (auto fiw_itr = fiw_arr.begin(); fiw_itr != fiw_arr.end(); ++fiw_itr) {
        int factor = (sum(fiw_itr.indices()) % 2 ? -1 : 1);
        *fiw_itr += fk_arr()[count] * factor;
        ++count;
      }
    }
  };
} // namespace triqs_ctint
