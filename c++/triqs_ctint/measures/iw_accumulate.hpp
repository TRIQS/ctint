#pragma once
#include <xsimd/xsimd.hpp>

#include "./../types.hpp"

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
    template <class T, auto N> constexpr auto min_simd_width() {
      // finds the smallest simd width that can handle N elements
      // simd size is batch size the SIMD width in xsimd terminology
      if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return min_simd_width<T, N * 2>();
      } else {
        return N;
      }
    };

    template <class T, auto N> constexpr auto max_simd_width() {
      // finds the smallest simd width that can handle N elements
      // simd size is batch size the SIMD width in xsimd terminology
      if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return max_simd_width<T, N / 2>();
      } else {
        return N;
      }
    };

    constexpr auto min_width = min_simd_width<std::complex<double>, 1>();
    constexpr auto max_width = max_simd_width<std::complex<double>, 128>();

    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<double, min_width>>, "min_width is not valid");
    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<double, max_width>>, "max_width is not valid");

    template <auto N, class T = double> constexpr auto get_arch_for_size() { return typename xsimd::make_sized_batch_t<T, N>::arch_type{}; }

    template <class T, auto N> using make_sized_batch_t = xsimd::batch<T, decltype(get_arch_for_size<N, T>())>;

    static_assert(!std::is_void_v<make_sized_batch_t<std::complex<double>, min_width>>, "failed to create min_width complex batch");
    static_assert(!std::is_void_v<make_sized_batch_t<std::complex<double>, max_width>>, "failed to create max_width complex batch");

    // === Fused, vectorized version for bl1_batch == bl2_batch ===
    template <int batch_size, typename M1a_t, typename M2a_t, typename M1b_t, typename M2b_t, typename M4a_t>
    void process_inner_loop_fused(mc_weight_t sign, const M1a_t &M1a, const M2a_t &M2a, const M1b_t &M1b, const M2b_t &M2b, M4a_t &M4a,
                                  auto bl1_size, auto bl2_size, auto bl1, auto bl2) {
      using batch_t = make_sized_batch_t<std::complex<double>, std::max(std::min(batch_size, max_width), min_width)>;

      for (auto i : range(bl1_size)) {
        for (auto j : range(bl1_size)) {
          const auto M1_v = batch_t{M1a(j, i) * sign};

          for (auto k : range(bl2_size)) {
            auto *M4_ptr  = &M4a(i, j, k, 0);
            auto *M2a_ptr = M2a.data() + k * bl2_size;

            const auto V   = batch_t::size;
            const auto end = bl2_size & ~(V - 1);
            auto idx       = 0;

            if (bl1 == bl2) {
              const auto M2s_v = batch_t{M2b(j, k) * sign};
              auto *M1b_ptr    = &M1b(0, i);

              // fused FMA + FNMA
              for (; idx < end; idx += V) {
                auto m4  = batch_t::load_unaligned(M4_ptr + idx);
                auto m2a = batch_t::load_unaligned(M2a_ptr + idx);
                auto m1b = batch_t::load_unaligned(M1b_ptr + idx);

                auto tmp    = xsimd::fma(M1_v, m2a, m4);
                auto result = xsimd::fnma(M2s_v, m1b, tmp);

                result.store_unaligned(M4_ptr + idx);
              }
              for (; idx < bl2_size; ++idx) {
                auto &m4 = M4_ptr[idx];
                m4       = xsimd::fnma(M2b(j, k) * sign, M1b(idx, i), xsimd::fma(M1a(j, i) * sign, M2a_ptr[idx], m4));
              }
            } else {
              // only FMA for M1a*M2a
              for (; idx < end; idx += V) {
                auto m4  = batch_t::load_unaligned(M4_ptr + idx);
                auto m2a = batch_t::load_unaligned(M2a_ptr + idx);

                auto result = xsimd::fma(M1_v, m2a, m4);
                result.store_unaligned(M4_ptr + idx);
              }
              for (; idx < bl2_size; ++idx) {
                auto &m4 = M4_ptr[idx];
                m4       = xsimd::fma(M1a(j, i) * sign, M2a_ptr[idx], m4);
              }
            }
          }
        }
      }
    }

    // === Fallback version for other batch sizes ===
    template <int bl1_batch, int bl2_batch, typename M1a_t, typename M2a_t, typename M1b_t, typename M2b_t, typename M4a_t>
    void process_inner_loop_fallback(mc_weight_t sign, const M1a_t &M1a, const M2a_t &M2a, const M1b_t &M1b, const M2b_t &M2b, M4a_t &M4a,
                                     auto bl1_size, auto bl2_size, auto bl1, auto bl2) {
      using batch1_t = make_sized_batch_t<std::complex<double>, std::max(std::min(bl1_batch, max_width), min_width)>;
      using batch2_t = make_sized_batch_t<std::complex<double>, std::max(std::min(bl2_batch, max_width), min_width)>;

      for (auto i : range(bl1_size)) {
        for (auto j : range(bl1_size)) {
          uint64_t index       = 0;
          const auto M1val     = M1a(j, i) * sign;
          const auto bl2square = bl2_size * bl2_size;

          // M1a * M2a update
          if constexpr (bl1_batch >= min_width) {
            const auto M1_v           = batch1_t{M1val};
            const auto truncated_size = bl2square & -(batch1_t::size);
            for (; index < truncated_size; index += batch1_t::size) {
              auto *const RESTRICT m4_ptr = &M4a(i, j, 0, 0) + index;
              auto batch                  = batch1_t::load_unaligned(m4_ptr);
              auto M2_batch               = batch1_t::load_unaligned(M2a.data() + index);
              auto result                 = xsimd::fma(M1_v, M2_batch, batch);
              result.store_unaligned(m4_ptr);
            }
          }
          for (; index < bl2square; ++index) { (&M4a(i, j, 0, 0))[index] = xsimd::fma(M1val, M2a.data()[index], (&M4a(i, j, 0, 0))[index]); }

          // M2b * M1b update only on diagonal blocks
          if (bl1 == bl2) {
            for (auto k : range(bl2_size)) {
              index             = 0;
              const auto M2sval = M2b(j, k) * sign;
              if constexpr (bl2_batch >= min_width) {
                const auto M2s_v          = batch2_t{M2sval};
                const auto truncated_size = bl2_size & -(batch2_t::size);
                for (; index < truncated_size; index += batch2_t::size) {
                  auto *const RESTRICT m4_ptr = &M4a(i, j, k, index);
                  auto batch                  = batch2_t::load_unaligned(m4_ptr);
                  auto M1_batch               = batch2_t::load_unaligned(&M1b(index, i));
                  auto result                 = xsimd::fnma(M2s_v, M1_batch, batch);
                  result.store_unaligned(m4_ptr);
                }
              }
              for (; index < bl2_size; ++index) { M4a(i, j, k, index) = xsimd::fnma(M2sval, M1b(index, i), M4a(i, j, k, index)); }
            }
          }
        }
      }
    }

    // === Wrapper that dispatches to fused or fallback ===
    template <int bl1_batch, int bl2_batch, typename M1a_t, typename M2a_t, typename M1b_t, typename M2b_t, typename M4a_t>
    void process_inner_loop(mc_weight_t sign, const M1a_t &M1a, const M2a_t &M2a, const M1b_t &M1b, const M2b_t &M2b, M4a_t &M4a,
                            auto bl1_size, auto bl2_size, auto bl1, auto bl2) {
      if constexpr (bl1_batch == bl2_batch && bl1_batch >= min_width) {
        process_inner_loop_fused<bl1_batch>(sign, M1a, M2a, M1b, M2b, M4a, bl1_size, bl2_size, bl1, bl2);
      } else {
        process_inner_loop_fallback<bl1_batch, bl2_batch>(sign, M1a, M2a, M1b, M2b, M4a, bl1_size, bl2_size, bl1, bl2);
      }
    }

    const auto iw4_accumulate_kernel = []<int bl1_batch, int bl2_batch>(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1,
                                                                        const auto bl2) {
      //auto const &iw_mesh = std::get<0>(M4_iw(0, 0).mesh());
      auto &[iw_mesh, _, _] = M4_iw(0, 0).mesh();
      auto const M1         = M[bl1];
      auto const M2         = M[bl2];
      auto const bl1_size   = M1.target_shape()[0];
      auto const bl2_size   = M2.target_shape()[0];
      auto &M4              = M4_iw(bl1, bl2);

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
    };

    const auto iw4ph_accumulate_kernel = []<int bl1_batch, int bl2_batch>(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1,
                                                                          const auto bl2) {
      //auto const &iW_mesh = std::get<0>(M4_iw(0, 0).mesh());
      //auto const &iw_mesh = std::get<1>(M4_iw(0, 0).mesh());
      auto const &[iW_mesh, iw_mesh, _] = M4_iw(0, 0).mesh();
      auto const M1                     = M[bl1];
      auto const M2                     = M[bl2];
      auto const bl1_size               = M1.target_shape()[0];
      auto const bl2_size               = M2.target_shape()[0];
      auto &M4                          = M4_iw(bl1, bl2);

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
    };

    const auto iw4pp_accumulate_kernel = []<int bl1_batch, int bl2_batch>(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1,
                                                                          const auto bl2) {
      //auto const &iW_mesh = std::get<0>(M4_iw(0, 0).mesh());
      //auto const &iw_mesh = std::get<1>(M4_iw(0, 0).mesh());
      auto const &[iW_mesh, iw_mesh, _] = M4_iw(0, 0).mesh();
      auto const bl1_size               = M[bl1].target_shape()[0];
      auto const bl2_size               = M[bl2].target_shape()[0];
      auto const M1                     = M[bl1];
      auto const M2                     = M[bl2];
      auto &M4                          = M4_iw(bl1, bl2);

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
    };

    const auto iw3pp_accumulate_kernel = []<int bl1_batch, int bl2_batch>(mc_weight_t sign, const auto &GM, auto &M3_iw, const auto bl1,
                                                                          const auto bl2) {
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
    };

    const auto iw3ph_accumulate_kernel = []<int bl1_batch, int bl2_batch>(mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM,
                                                                          const auto &MG, auto &M3_iw, const auto bl1, const auto bl2) {
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
    };
  } // namespace

  namespace simd {

    template <auto kernel> void kernel_dispatch(const auto bl2_size, auto &&...args) {
      // Dispatch to the correct SIMD instruction width based on the size of the blocks
      // It will try to use the widest SIMD instruction available for the given block sizes
      // TODO: fold expressions might be an option to simplify the code
      if (bl2_size >= 8) {
        return kernel.template operator()<8, 8>(args...);
      } else if (bl2_size >= 4) {
        return kernel.template operator()<8, 4>(args...);
      } else if (bl2_size >= 3) {
        return kernel.template operator()<8, 2>(args...);
      } else if (bl2_size >= 2) {
        return kernel.template operator()<4, 2>(args...);
      } else {
        return kernel.template operator()<1, 1>(args...);
      }
    }

    void iw4_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      kernel_dispatch<iw4_accumulate_kernel>(bl2_size, sign, M, M4_iw, bl1, bl2);
    }

    void iw4ph_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      kernel_dispatch<iw4ph_accumulate_kernel>(bl2_size, sign, M, M4_iw, bl1, bl2);
    }

    void iw4pp_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      kernel_dispatch<iw4pp_accumulate_kernel>(bl2_size, sign, M, M4_iw, bl1, bl2);
    }

    void iw3ph_accumulate(mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG, auto &M4_iw, const auto bl1,
                          const auto bl2, const auto bl2_size) {
      kernel_dispatch<iw3ph_accumulate_kernel>(bl2_size, sign, M, GMG, GM, MG, M4_iw, bl1, bl2);
    }

    void iw3pp_accumulate(mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) {
      kernel_dispatch<iw3pp_accumulate_kernel>(bl2_size, sign, M, M4_iw, bl1, bl2);
    }

  } // namespace simd
} // namespace triqs_ctint::measures
