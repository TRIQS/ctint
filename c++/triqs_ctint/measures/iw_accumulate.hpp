#pragma once
#include <xsimd/xsimd.hpp>

#if defined(_MSC_VER)
// For Microsoft Visual Studio
#define RESTRICT __restrict
#elif defined(__GNUC__)
// For GCC and Clang
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

namespace triqs_ctint::measures {
  namespace {
    template <class T, auto N> static constexpr auto min_simd_width() {
      // finds the smallest simd width that can handle N elements
      // simd size is batch size the SIMD width in xsimd terminology
      if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return min_simd_width<T, N * 2>();
      } else {
        return N;
      }
    };

    template <class T, auto N> static constexpr auto max_simd_width() {
      // finds the smallest simd width that can handle N elements
      // simd size is batch size the SIMD width in xsimd terminology
      if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return max_simd_width<T, N / 2>();
      } else {
        return N;
      }
    };

    static constexpr auto min_width = min_simd_width<double, 2>();
    static constexpr auto max_width = max_simd_width<double, 64>();

    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<double, min_width>>, "min_width is not valid");
    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<double, max_width>>, "max_width is not valid");

    template <auto N, class T = double> static constexpr auto get_arch_for_size() { return typename xsimd::make_sized_batch_t<T, N>::arch_type{}; }

    template <class T, auto N> using make_complex_sized_batch_t = xsimd::batch<std::complex<T>, decltype(get_arch_for_size<N, T>())>;

    static_assert(!std::is_void_v<make_complex_sized_batch_t<double, min_width>>, "failed to create min_width complex batch");
    static_assert(!std::is_void_v<make_complex_sized_batch_t<double, max_width>>, "failed to create max_width complex batch");

    template <auto bl1_batch, auto bl2_batch>
    void process_inner_loop(mc_weight_t sign, const auto &M1a, const auto &M2a, const auto &M1b, const auto &M2b, auto &M4a, const auto bl1_size,
                            const auto bl2_size, const auto bl1, const auto bl2) {
      using batch1_t = make_complex_sized_batch_t<double, std::max(std::min(bl1_batch, max_width), min_width)>;
      using batch2_t = make_complex_sized_batch_t<double, std::max(std::min(bl2_batch, max_width), min_width)>;

      for (auto i : range(bl1_size)) {
        for (auto j : range(bl1_size)) {
          uint64_t index = 0;

          const auto M1val     = M1a(j, i) * sign;
          const auto bl2square = bl2_size * bl2_size;
          if constexpr (bl1_batch >= min_width) {
            const auto M1_v           = batch1_t(M1val);
            const auto truncated_size = bl2square & (-batch1_t::size);
            for (; index < truncated_size; index += batch1_t::size) {
              auto *const RESTRICT m4_ptr = &M4a(i, j, 0, 0) + index;
              const auto batch            = batch1_t::load_unaligned(m4_ptr);
              const auto M2_batch         = batch1_t::load_unaligned(M2a.data() + index);
              const auto result           = xsimd::fma(M1_v, M2_batch, batch);
              result.store_unaligned(m4_ptr);
            }
          }
          for (; index < bl2square; index++) { (&M4a(i, j, 0, 0))[index] = xsimd::fma(M1val, M2a.data()[index], (&M4a(i, j, 0, 0))[index]); }

          if (bl1 == bl2) [[unlikely]] {
            for (const auto k : range(bl2_size)) {
              index             = 0;
              const auto M2sval = M2b(j, k) * sign;
              if constexpr (bl2_batch >= min_width) {
                const auto M2s_v          = batch2_t(M2sval);
                const auto truncated_size = bl2_size & (-batch2_t::size);
                for (; index < truncated_size; index += batch2_t::size) {
                  auto *const RESTRICT m4_ptr = &M4a(i, j, k, index);
                  const auto batch            = batch2_t::load_unaligned(m4_ptr);
                  const auto M1_batch         = batch2_t::load_unaligned(&M1b(index, i));
                  const auto result           = xsimd::fnma(M2s_v, M1_batch, batch);
                  result.store_unaligned(m4_ptr);
                }
              }
              for (; index < bl2_size; index++) { M4a(i, j, k, index) = xsimd::fnma(M2sval, M1b(index, i), M4a(i, j, k, index)); }
            }
          }
        }
      }
    }

    template <auto bl1_batch, auto bl2_batch>
    void iw4_accumulate_kernel(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2) {
      auto const &iw_mesh = std::get<0>(M4_iw(0, 0).mesh());
      auto const bl1_size = M[bl1].target_shape()[0];
      auto const bl2_size = M[bl2].target_shape()[0];
      auto const M1       = M[bl1];
      auto const M2       = M[bl2];
      auto &M4            = M4_iw(bl1, bl2);

      for (const auto &iw1 : iw_mesh) {
        for (const auto &iw2 : iw_mesh) {
          for (const auto &iw3 : iw_mesh) {
            const auto iw4 = iw1 + iw3 - iw2;
            const auto M1a = M1[iw2.value(), iw1];
            const auto M2a = M2[iw4, iw3];
            const auto M1b = M1[iw4, iw1];
            const auto M2b = M2[iw2.value(), iw3];
            auto M4a       = M4[iw1, iw2, iw3];
            process_inner_loop<bl1_batch, bl2_batch>(sign, M1a, M2a, M1b, M2b, M4a, bl1_size, bl2_size, bl1, bl2);
          }
        }
      }
    }

    template <auto bl1_batch, auto bl2_batch>
    void iw4ph_accumulate_kernel(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2) {
      auto const &iW_mesh = std::get<0>(M4_iw(0, 0).mesh());
      auto const &iw_mesh = std::get<1>(M4_iw(0, 0).mesh());
      auto const bl1_size = M[bl1].target_shape()[0];
      auto const bl2_size = M[bl2].target_shape()[0];
      auto const M1       = M[bl1];
      auto const M2       = M[bl2];
      auto &M4            = M4_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          for (auto iwp : iw_mesh) {
            const auto M1a = M1[iW + iw, iw.value()];
            const auto M2a = M2[iwp.value(), iW + iwp];
            const auto M1b = M1[iwp.value(), iw.value()];
            const auto M2b = M2[iW + iw, iW + iwp];
            auto M4a       = M4[iW, iw, iwp];
            process_inner_loop<bl1_batch, bl2_batch>(sign, M1a, M2a, M1b, M2b, M4a, bl1_size, bl2_size, bl1, bl2);
          }
        }
      }
    }

    template <auto bl1_batch, auto bl2_batch>
    void iw4pp_accumulate_kernel(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2) {
      auto const &iW_mesh = std::get<0>(M4_iw(0, 0).mesh());
      auto const &iw_mesh = std::get<1>(M4_iw(0, 0).mesh());
      auto const bl1_size = M[bl1].target_shape()[0];
      auto const bl2_size = M[bl2].target_shape()[0];
      auto const M1       = M[bl1];
      auto const M2       = M[bl2];
      auto &M4            = M4_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          for (auto iwp : iw_mesh) {
            const auto M1a = M1[iW - iwp, iw.value()];
            const auto M2a = M2[iwp.value(), iW - iw];
            const auto M1b = M1[iwp.value(), iw.value()];
            const auto M2b = M2[iW - iwp, iW - iw];
            auto M4a       = M4[iW, iw, iwp];
            process_inner_loop<bl1_batch, bl2_batch>(sign, M1a, M2a, M1b, M2b, M4a, bl1_size, bl2_size, bl1, bl2);
          }
        }
      }
    }

    template <auto bl1_batch, auto bl2_batch>
    void iw3pp_accumulate_kernel(mc_weight_t sign, const auto &GM, auto &M3_iw, const auto bl1, const auto bl2) {
      auto const [iW_mesh, iw_mesh] = M3_iw(0, 0).mesh();
      auto const bl1_size           = GM[bl1].target_shape()[0];
      auto const bl2_size           = GM[bl2].target_shape()[0];
      auto const GM1                = GM[bl1];
      auto const GM2                = GM[bl2];
      auto &M3                      = M3_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          const auto M1a = GM1[iw.value()];
          const auto M2a = GM2[iW - iw];
          auto M4a       = M3[iW, iw];
          process_inner_loop<bl1_batch, bl2_batch>(sign, M1a, M2a, M1a, M2a, M4a, bl1_size, bl2_size, bl1, bl2);
        }
      }
    }

    template <auto bl1_batch, auto bl2_batch>
    void iw3ph_accumulate_kernel(mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG, auto &M3_iw, const auto bl1,
                                 const auto bl2) {
      auto const [iW_mesh, iw_mesh] = M3_iw(0, 0).mesh();
      auto const bl1_size           = M[bl1].target_shape()[0];
      auto const bl2_size           = M[bl2].target_shape()[0];
      auto const M1                 = M[bl1];
      auto const GMG2               = GMG(bl2);
      auto const GM1                = GM[bl1];
      auto const MG2                = MG(bl2);
      auto &M3                      = M3_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          const auto M1a = M1[iW + iw, iw.value()];
          const auto M2a = GMG2;
          const auto M1b = GM1[iw.value()];
          const auto M2b = MG2[iW + iw];
          auto M4a       = M3[iW, iw];
          process_inner_loop<bl1_batch, bl2_batch>(sign, M1a, M2a, M1b, M2b, M4a, bl1_size, bl2_size, bl1, bl2);
        }
      }
    }
  } // namespace

  namespace simd {

    void iw4_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      // Dispatch to the correct SIMD instruction width based on the size of the blocks
      // It will try to use the widest SIMD instruction available for the given block sizes
      // TODO: fold expressions might be an option to simplify the code
      if (bl2_size >= 8) {
        return iw4_accumulate_kernel<8, 8>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 4) {
        return iw4_accumulate_kernel<8, 4>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 3) {
        return iw4_accumulate_kernel<8, 2>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 2) {
        return iw4_accumulate_kernel<4, 2>(sign, M, M4_iw, bl1, bl2);
      } else {
        return iw4_accumulate_kernel<1, 1>(sign, M, M4_iw, bl1, bl2);
      }
    }

    void iw4ph_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      // Dispatch to the correct SIMD instruction width based on the size of the blocks
      // It will try to use the widest SIMD instruction available for the given block sizes
      // TODO: fold expressions might be an option to simplify the code
      if (bl2_size >= 8) {
        return iw4ph_accumulate_kernel<8, 8>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 4) {
        return iw4ph_accumulate_kernel<8, 4>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 3) {
        return iw4ph_accumulate_kernel<8, 2>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 2) {
        return iw4ph_accumulate_kernel<4, 2>(sign, M, M4_iw, bl1, bl2);
      } else {
        return iw4ph_accumulate_kernel<1, 1>(sign, M, M4_iw, bl1, bl2);
      }
    }

    void iw4pp_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      // Dispatch to the correct SIMD instruction width based on the size of the blocks
      // It will try to use the widest SIMD instruction available for the given block sizes
      // TODO: fold expressions might be an option to simplify the code
      if (bl2_size >= 8) {
        return iw4pp_accumulate_kernel<8, 8>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 4) {
        return iw4pp_accumulate_kernel<8, 4>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 3) {
        return iw4pp_accumulate_kernel<8, 2>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 2) {
        return iw4pp_accumulate_kernel<4, 2>(sign, M, M4_iw, bl1, bl2);
      } else {
        return iw4pp_accumulate_kernel<1, 1>(sign, M, M4_iw, bl1, bl2);
      }
    }

    void iw3ph_accumulate(mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG, auto &M4_iw, const auto bl1,
                          const auto bl2, const auto bl2_size) {
      // Dispatch to the correct SIMD instruction width based on the size of the blocks
      // It will try to use the widest SIMD instruction available for the given block sizes
      // TODO: fold expressions might be an option to simplify the code
      if (bl2_size >= 8) {
        return iw3ph_accumulate_kernel<8, 8>(sign, M, GMG, GM, MG, M4_iw, bl1, bl2);
      } else if (bl2_size >= 4) {
        return iw3ph_accumulate_kernel<8, 4>(sign, M, GMG, GM, MG, M4_iw, bl1, bl2);
      } else if (bl2_size >= 3) {
        return iw3ph_accumulate_kernel<8, 2>(sign, M, GMG, GM, MG, M4_iw, bl1, bl2);
      } else if (bl2_size >= 2) {
        return iw3ph_accumulate_kernel<4, 2>(sign, M, GMG, GM, MG, M4_iw, bl1, bl2);
      } else {
        return iw3ph_accumulate_kernel<1, 1>(sign, M, GMG, GM, MG, M4_iw, bl1, bl2);
      }
    }

    void iw3pp_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      // Dispatch to the correct SIMD instruction width based on the size of the blocks
      // It will try to use the widest SIMD instruction available for the given block sizes
      // TODO: fold expressions might be an option to simplify the code
      if (bl2_size >= 8) {
        return iw3pp_accumulate_kernel<8, 8>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 4) {
        return iw3pp_accumulate_kernel<8, 4>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 3) {
        return iw3pp_accumulate_kernel<8, 2>(sign, M, M4_iw, bl1, bl2);
      } else if (bl2_size >= 2) {
        return iw3pp_accumulate_kernel<4, 2>(sign, M, M4_iw, bl1, bl2);
      } else {
        return iw3pp_accumulate_kernel<1, 1>(sign, M, M4_iw, bl1, bl2);
      }
    }
  } // namespace simd
} // namespace triqs_ctint::measures
