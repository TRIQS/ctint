#pragma once
#include <nda/simd/simd.hpp>

namespace triqs_ctint::measures {

  namespace {
    template <long outer_loop, long inner_loop, bool equal_sizes>
    static void _impl_process_vectorized_loop(const mc_weight_t sign, const long bl1_size, const long bl2_size,
                                              const basic_array_view<const std::complex<double>, 2, C_layout> &M1_view_val,
                                              const basic_array_view<const std::complex<double>, 2, C_layout> &M1_view_if,
                                              const basic_array_view<const std::complex<double>, 2, C_layout> &M2_view,
                                              const basic_array_view<const std::complex<double>, 2, C_layout> &M2_view_if,
                                              basic_array_view<std::complex<double>, 4, C_stride_layout> &M4_view) {

      static_assert(outer_loop == 1, "outloop width 1 is supported.");
      static_assert(inner_loop == 1 or inner_loop == 2 or inner_loop == 4);
      static_assert(outer_loop == 1 or outer_loop == 2 or outer_loop == 4);
      using value_t                    = std::complex<double>;
      using simd_t                     = fixed_size_simd<value_t, inner_loop>;
      constexpr size_t inner_simd_size = simd_t::size();
      const long M1_view_if_stride     = M1_view_if.strides()[0]; // C_layout, stride length between columns.
      const long M4_view_stride        = M4_view.strides()[2];    // C_layout 2nd fastest dimension
      const long M2_view_stride        = M2_view.strides()[0];    // C_layout 2nd fastest dimension.
      std::array<simd_t, inner_simd_size> M4_block, M2_view_block;
      for (long i : range(bl1_size)) {
        for (long j : range(bl1_size)) {
          const auto M1val             = M1_view_val(j, i);
          const simd_t M1val_sign_simd = simd_t(sign * M1_view_val(j, i));

          long k = 0;
          for (; k + inner_simd_size <= bl2_size; k += inner_simd_size) {
            long l = 0;
            for (; l + inner_simd_size <= bl2_size; l += inner_simd_size) {
              // Here l and k indices are vectorized.

              // Load inner_simd_size*inner_simd_size chunk from M4_view and M2_view
              for (int m = 0; m < inner_simd_size; ++m) {
                M4_block[m].load_unaligned(&M4_view(i, j, k + m, l));
                M2_view_block[m].load_unaligned(&M2_view(l + m, k));
              }
              // Transpose the M2_view. Because 'l' index is the column index of M2_view.
              simd::kernel_transpose(M2_view_block);
              for (int m = 0; m < inner_simd_size; ++m) { M4_block[m] = simd::fma_add(M1val_sign_simd, M2_view_block[m], M4_block[m]); }

              if constexpr (equal_sizes) {
                //Since l is column. We need to call gather with correct stride. TODO: If you vectorize i as well we can call normal loads and then transpose.
                simd_t tmp = simd::gather<simd_t>(&M1_view_if(l, i), M1_view_if_stride);
                for (int m = 0; m < inner_simd_size; ++m) {
                  M4_block[m] = simd::fma_nadd(tmp * sign, simd_t(M2_view_if(j, k + m)), M4_block[m]);
                }
              }
              for (int m = 0; m < inner_simd_size; ++m) { M4_block[m].store_unaligned(&M4_view(i, j, k + m, l)); }
            }
            for (; l < bl2_size; ++l) {
              simd_t tmp = simd::gather<simd_t>(&M4_view(i, j, k, l), M4_view_stride);
              simd_t tmp2;
              tmp2.load_unaligned(&M2_view(l, k));
              tmp = simd::fma_add(M1val_sign_simd, tmp2, tmp);
              if constexpr (equal_sizes) {
                simd_t tmp3(M1_view_if(l, i));
                simd_t tmp4;
                tmp4.load_unaligned(&M2_view_if(j, k));
                tmp = simd::fma_nadd(tmp3 * sign, tmp4, tmp);
              }
              simd::scatter(tmp, &M4_view(i, j, k, l), M4_view_stride);
            }
          }
          for (; k < bl2_size; ++k) {
            long l = 0;
            for (; l + inner_simd_size <= bl2_size; l += inner_simd_size) {
              simd_t tmp;
              tmp.load_unaligned(&M4_view(i, j, k, l));
              tmp = simd::fma_add(M1val_sign_simd, simd::gather<simd_t>(&M2_view(l, k), M2_view_stride), tmp);
              if constexpr (equal_sizes) {
                tmp = simd::fma_nadd(simd::gather<simd_t>(&M1_view_if(l, i), M1_view_if_stride) * sign, simd_t(M2_view_if(j, k)), tmp);
              }
              tmp.store_unaligned(&M4_view(i, j, k, l));
            }
            for (; l < bl2_size; ++l) {
              M4_view(i, j, k, l) += sign * M1val * M2_view(l, k);
              if constexpr (equal_sizes) { M4_view(i, j, k, l) -= sign * M1_view_if(l, i) * M2_view_if(j, k); }
            }
          }
        }
      }
    }
  } // namespace
  template <bool equal_sizes>
  static void process_vectorized_loop(const mc_weight_t sign, const long bl1_size, const long bl2_size,
                                      const basic_array_view<const std::complex<double>, 2, C_layout> &M1_view_val,
                                      const basic_array_view<const std::complex<double>, 2, C_layout> &M1_view_if,
                                      const basic_array_view<const std::complex<double>, 2, C_layout> &M2_view,
                                      const basic_array_view<const std::complex<double>, 2, C_layout> &M2_view_if,
                                      basic_array_view<std::complex<double>, 4, C_stride_layout> &M4_view) {
    constexpr size_t max_width = native_simd<std::complex<double>>::size();
    if (bl2_size >= max_width or max_width == 1) { //NOLINT
      _impl_process_vectorized_loop<1, max_width, equal_sizes>(sign, bl1_size, bl2_size, M1_view_val, M1_view_if, M2_view, M2_view_if, M4_view);
    } else if (bl2_size >= max_width / 2 or max_width / 2 == 1) { //NOLINT
      _impl_process_vectorized_loop<1, max_width / 2, equal_sizes>(sign, bl1_size, bl2_size, M1_view_val, M1_view_if, M2_view, M2_view_if, M4_view);
    } else {
      _impl_process_vectorized_loop<1, 1, equal_sizes>(sign, bl1_size, bl2_size, M1_view_val, M1_view_if, M2_view, M2_view_if, M4_view);
    }
  }
} // namespace triqs_ctint::measures