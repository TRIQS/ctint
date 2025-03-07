#include "./M4_iw.hpp"
#include <nda/simd/simd.hpp>
#include <cmath>

namespace triqs_ctint::measures {

  M4_iw::M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M4};
    mesh::prod<imfreq, imfreq, imfreq> M4_iw_mesh{iw_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4_iw = make_block2_gf(M4_iw_mesh, params.gf_struct);
    M4_iw_.rebind(results->M4_iw.value());
    M4_iw_() = 0;

    // Construct Matsubara mesh for temporary Matrix
    mesh::imfreq iw_mesh_large{params.beta, Fermion, 3 * params.n_iw_M4};
    mesh::prod<imfreq, imfreq> M_mesh{iw_mesh_large, iw_mesh};

    // Initialize intermediate scattering matrix
    M = block_gf{M_mesh, params.gf_struct};

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), init_target_func};
    }
  }

  void M4_iw::accumulate(const mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int bl : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl],
               [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) { // Care for negative frequency in c transform (for M-objects)
                 buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), params.beta - double(c_i.tau)}, -Ginv_ji);
               })
        ;
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    auto const &iw_mesh = std::get<0>(M4_iw_(0, 0).mesh());

    for (int bl1 : range(params.n_blocks())) // FIXME c++17 Loops
      for (int bl2 : range(params.n_blocks())) {
        int bl1_size   = M[bl1].target_shape()[0];
        int bl2_size   = M[bl2].target_shape()[0];
        auto const &M1 = M[bl1];
        auto const &M2 = M[bl2];
        auto &M4       = M4_iw_(bl1, bl2);

        for (auto iw1 : iw_mesh)
          for (auto iw2 : iw_mesh)
            for (auto iw3 : iw_mesh)
              for (int i : range(bl1_size))
                for (int j : range(bl1_size)) {
                  auto M1val = M1[iw2.value(), iw1](j, i);

                  for (int k : range(bl2_size))
                    for (int l : range(bl2_size)) {
                      auto iw4 = iw1 + iw3 - iw2;
                      M4[iw1, iw2, iw3](i, j, k, l) += sign * M1val * M2[iw4, iw3](l, k);
                      if (bl1 == bl2) { M4[iw1, iw2, iw3](i, j, k, l) -= sign * M1[iw4, iw1](l, i) * M2[iw2.value(), iw3](j, k); }
                    }
                }
      }
  }

  void M4_iw::accumulate_v2(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int bl : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl],
               [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) { // Care for negative frequency in c transform (for M-objects)
                 buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), params.beta - double(c_i.tau)}, -Ginv_ji);
               })
        ;
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    auto const &iw_mesh        = std::get<0>(M4_iw_(0, 0).mesh());
    using value_t              = std::complex<double>;
    using simd_t               = native_simd<value_t>;
    constexpr size_t simd_size = simd_t::size();
    std::array<simd_t, simd_size> kernel_block, M4_block;
    simd_t sign_simd;
    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(sign)>, value_t>) {
      sign_simd = simd_t(sign);
    } else if constexpr (std::is_same_v<std::remove_cvref_t<decltype(sign)>, double>) {

      sign_simd = simd_t(value_t(sign, 0)); // TODO: Improve this by adding operator*(const Scalar& v) in complex simds.
    }
    for (int bl1 : range(params.n_blocks())) // FIXME c++17 Loops
      for (int bl2 : range(params.n_blocks())) {
        const bool equal_sizes = (bl1 == bl2);
        long bl1_size          = M[bl1].target_shape()[0];
        long bl2_size          = M[bl2].target_shape()[0];
        auto const &M1         = M[bl1];
        auto const &M2         = M[bl2];
        auto &M4               = M4_iw_(bl1, bl2);

        for (auto iw1 : iw_mesh) {
          for (auto iw2 : iw_mesh) {
            for (auto iw3 : iw_mesh) {
              const auto iw4         = iw1 + iw3 - iw2;
              auto M4_view           = M4[iw1, iw2, iw3];
              const auto M2_view     = M2[iw4, iw3];
              const auto M1_view     = M1[iw4, iw2];
              const auto M1_view_val = M1[iw2.value(), iw1];
              const auto M2_view_if  = M2[iw2.value(), iw3];
              const auto M1_view_if  = M1[iw4, iw1];
              //M4[iw1, iw2, iw3](i, j, k, l) += sign * M1val * M2[iw4, iw3](l, k);
              //                      if (bl1 == bl2) { M4[iw1, iw2, iw3](i, j, k, l) -= sign * M1[iw4, iw1](l, i) * M2[iw2.value(), iw3](j, k); }
              const long M1_view_if_stride = M1_view.strides()[0]; // C_layout, stride length between columns.
              const long M4_view_stride    = M4_view.strides()[2]; // C_layout 2nd fastest dimension
              const long M2_view_stride    = M2_view.strides()[0]; // C_layout 2nd fastest dimension.
              //TODO: What do I know about view sizes. Are they bl2_size * bl2_size matrices ?
              for (long i : range(bl1_size)) {
                for (long j : range(bl1_size)) {
                  const auto M1val             = M1_view_val(j, i);
                  const simd_t M1val_sign_simd = simd_t(sign * M1_view_val(j, i));

                  long k = 0;
                  for (; k + simd_size <= bl2_size; k += simd_size) {
                    long l = 0;
                    for (; l + simd_size <= bl2_size; l += simd_size) {
                      for (int m = 0; m < simd_size; ++m) {
                        M4_block[m]     = M4_view.load(i, j, k + m, l);
                        kernel_block[m] = M2_view.load(l + m, k);
                      }
                      kernel_block = simd::kernel_transpose(kernel_block);
                      for (int m = 0; m < simd_size; ++m) { M4_block[m] = simd::fma_add(M1val_sign_simd, kernel_block[m], M4_block[m]); }
                      if (equal_sizes) {
                        simd_t tmp  = simd::gather<simd_t>(&M1_view_if(l, i), M1_view_if_stride);
                        simd_t tmp2 = M2_view_if.load(j, k);
                        for (int m = 0; m < simd_size; ++m) { M4_block[m] = simd::fma_nadd(sign_simd * tmp, tmp2, M4_block[m]); }
                      }
                      for (int m = 0; m < simd_size; ++m) { M4_view.store(M4_block[m], i, j, k + m, l); }
                    }
                    for (; l < bl2_size; ++l) {
                      alignas(simd_t::alignment()) std::array<value_t, simd_size> tmp_array;
                      simd_t tmp = simd::gather<simd_t>(&M4_view(i, j, k, l), M4_view_stride);
                      tmp = simd::fma_add(M1val_sign_simd, M2_view.load(l, k), tmp);
                      if (equal_sizes) {
                        simd_t tmp2(M1_view_if(l, i));
                        tmp = simd::fma_nadd(sign_simd * tmp2, M2_view_if.load(j, k), tmp);
                      }
                      // TODO: Implement scatter with constant strides in simd classes.
                      tmp.store(tmp_array.data());
                      for (int m = 0; m < simd_size; ++m) { M4_view(i, j, k + m, l) = tmp_array[m]; }
                    }
                  }
                  for (; k < bl2_size; ++k) {
                    long l = 0;
                    for (; l + simd_size <= bl2_size; l += simd_size) {
                      simd_t tmp = M4_view.load(i, j, k, l);
                      tmp = simd::fma_add(M1val_sign_simd, simd::gather<simd_t>(&M2_view(l, k), M2_view_stride), tmp);
                      if (equal_sizes) {
                        tmp = simd::fma_nadd(sign_simd * simd::gather<simd_t>(&M1_view_if(l, i), M1_view_if_stride), simd_t(M2_view_if(j, k)), tmp);
                      }
                      M4_view.store(tmp, i, j, k, l);
                    }
                    for (; l < bl2_size; ++l) {
                      M4_view(i, j, k, l) += sign * M1val * M2_view(l, k);
                      if (equal_sizes) { M4_view(i, j, k, l) -= sign * M1_view_if(l, i) * M2_view_if(j, k); }
                    }
                  }
                }
              }
            }
          }
        }
      }
  }

  void M4_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z      = mpi::all_reduce(Z, comm);
    M4_iw_ = mpi::all_reduce(M4_iw_, comm);
    M4_iw_ = M4_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
