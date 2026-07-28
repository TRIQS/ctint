// Micro-benchmark for the M4 iw-accumulate kernels (measures/iw_accumulate.hpp).
// Sweeps n_orb (=block size) so the compile-time-length dispatch is exercised across
// the small-block regime where it matters and into the memory-bound large-block regime.
#include "../../c++/triqs_ctint/measures/iw_accumulate.hpp"

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <benchmark/benchmark.h>

#include <complex>

using namespace triqs_ctint;
using namespace triqs_ctint::measures;

namespace {

  // Accumulate cost scales as (2*n_iw)^3 * n_blocks^2 * n_orb^4; keep n_iw small so an
  // iteration is cheap while still filling the 3D frequency mesh.
  constexpr double beta = 100.0;
  constexpr int n_iw    = 6;

  struct fixture {
    gf_struct_t gf_struct;
    // Must match the real measures (M4_iw.hpp): the target matrix is transposed in memory, which
    // is what makes M2a.data()[k*N+l] == M2a(l,k) and &M1b(0,i)+l == M1b(l,i) in the kernels.
    using M_layout = nda::contiguous_layout_with_stride_order<nda::encode(std::array{0, 1, 3, 2})>;
    using M_t      = block_gf<mesh::prod<mesh::imfreq, mesh::imfreq>, matrix_valued, M_layout>;
    M_t M;
    block2_gf<mesh::prod<mesh::imfreq, mesh::imfreq, mesh::imfreq>, tensor_valued<4>> M4;
    int n_blocks;

    explicit fixture(int n_orb) : gf_struct{{"dn", n_orb}, {"up", n_orb}}, n_blocks(2) {
      mesh::imfreq iw_mesh{beta, Fermion, n_iw};
      mesh::prod<mesh::imfreq, mesh::imfreq, mesh::imfreq> M4_mesh{iw_mesh, iw_mesh, iw_mesh};
      M4   = make_block2_gf(M4_mesh, gf_struct);
      M4() = 0;

      mesh::imfreq iw_mesh_large{beta, Fermion, 3 * n_iw};
      mesh::prod<mesh::imfreq, mesh::imfreq> M_mesh{iw_mesh_large, iw_mesh};
      M = M_t{M_mesh, gf_struct};

      // Deterministic non-trivial fill so the arithmetic is real and reproducible.
      for (int bl = 0; bl < n_blocks; ++bl) {
        auto *p       = M[bl].data().data();
        std::size_t n = M[bl].data().size();
        for (std::size_t i = 0; i < n; ++i)
          p[i] = std::complex<double>(0.1 + 0.001 * double(i % 97), -0.2 + 0.002 * double(i % 53));
      }
    }

    template <auto accum> void accumulate(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) accum(sign, M, M4, bl1, bl2, M[bl2].target_shape()[0]);
    }
  };

  template <auto accum> void run(benchmark::State &st) {
    fixture fx(static_cast<int>(st.range(0)));
    const auto sign = mc_weight_t{1.0};
    for (auto _ : st) {
      fx.accumulate<accum>(sign);
      benchmark::DoNotOptimize(fx.M4(0, 0).data().data());
      benchmark::ClobberMemory();
    }
  }

  // Pass by forwarding reference: M/M4 must bind by reference, not be copied per call.
  void iw4(benchmark::State &st) { run<[](auto &&...a) { measures::simd::iw4_accumulate(decltype(a)(a)...); }>(st); }
  void iw4ph(benchmark::State &st) { run<[](auto &&...a) { measures::simd::iw4ph_accumulate(decltype(a)(a)...); }>(st); }
  void iw4pp(benchmark::State &st) { run<[](auto &&...a) { measures::simd::iw4pp_accumulate(decltype(a)(a)...); }>(st); }

  // Correctness probe: one deterministic pass on fresh data, report the summed M4 (re/im).
  void chk(benchmark::State &st) {
    const auto sign = mc_weight_t{1.0};
    for (auto _ : st) {
      fixture fx(static_cast<int>(st.range(0)));
      fx.accumulate<[](auto &&...a) { measures::simd::iw4_accumulate(decltype(a)(a)...); }>(sign);
      double re = 0, im = 0;
      for (int b1 = 0; b1 < fx.n_blocks; ++b1)
        for (int b2 = 0; b2 < fx.n_blocks; ++b2) {
          auto *p = fx.M4(b1, b2).data().data();
          for (std::size_t k = 0; k < fx.M4(b1, b2).data().size(); ++k) { re += p[k].real(); im += p[k].imag(); }
        }
      st.counters["chk_re"] = re;
      st.counters["chk_im"] = im;
    }
  }

} // namespace

BENCHMARK(chk)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(7)->Arg(8)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw4)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(7)->Arg(8)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw4ph)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(7)->Arg(8)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw4pp)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5)->Arg(6)->Arg(7)->Arg(8)->Unit(benchmark::kMicrosecond);
