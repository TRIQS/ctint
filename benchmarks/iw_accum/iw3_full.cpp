// Micro-benchmark for the M3 full-mesh iw-accumulate kernels (measures/iw_accumulate.hpp).
//
// Purpose: fef408b deleted the n_orb=1 scalar fast paths that 1a30154 had added to
// full_iw3ph/full_iw3pp ("~20% slower than the pre-SIMD scalar implementation for
// single-orbital models"). This measures the shipped general-N kernel against a hand-written
// n_orb=1 loop with identical arithmetic, in one binary so the two alternate under the same
// machine conditions.
#include "../../c++/triqs_ctint/measures/iw_accumulate.hpp"

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <benchmark/benchmark.h>

#include <complex>

using namespace triqs_ctint;
using namespace triqs_ctint::measures;

namespace {

  // Cost scales as (2*n_iw)^2 * n_blocks^2 * n_orb^4; M3 itself is (2*n_iw)^2 * n_orb^4 complex
  // per block pair, so n_iw stays small to keep the container in cache-plus-a-bit rather than GB.
  constexpr double beta = 100.0;
  constexpr int n_iw    = 32;

  struct fixture {
    gf_struct_t gf_struct;
    int n_blocks;

    // Layouts must match the real measures (M3{ph,pp}_iw_full.hpp): GM/MG store their target
    // transposed, which is what makes M2a.data()[k*n+l] == M2a(l,k) and &M1b(0,i)+l == M1b(l,i).
    using simd_layout = nda::contiguous_layout_with_stride_order<nda::encode(std::array{0, 2, 1})>;
    block_gf<mesh::prod<mesh::imfreq, mesh::imfreq>, matrix_valued> M;
    block_gf<mesh::imfreq, matrix_valued, simd_layout> GM, MG;
    array<array<dcomplex, 2, nda::F_layout>, 1> GMG;
    block2_gf<mesh::prod<mesh::imfreq, mesh::imfreq>, tensor_valued<4>> M3;

    explicit fixture(int n_orb) : gf_struct{{"dn", n_orb}, {"up", n_orb}}, n_blocks(2) {
      mesh::imfreq iw_mesh{beta, Fermion, n_iw};
      mesh::prod iw2_mesh{iw_mesh, iw_mesh};

      M3   = make_block2_gf(iw2_mesh, gf_struct);
      M3() = 0;
      M    = block_gf<mesh::prod<mesh::imfreq, mesh::imfreq>, matrix_valued>{iw2_mesh, gf_struct};
      GM   = decltype(GM){iw_mesh, gf_struct};
      MG   = decltype(MG){iw_mesh, gf_struct};
      GMG  = array_adapter{make_shape(n_blocks), [&](int) { return array<dcomplex, 2, nda::F_layout>(n_orb, n_orb); }};

      fill(M);
      fill(GM);
      fill(MG);
      for (int bl = 0; bl < n_blocks; ++bl) fill_raw(GMG(bl).data(), GMG(bl).size(), bl);

    }

    // Deterministic non-trivial fill so the arithmetic is real and reproducible.
    static void fill_raw(dcomplex *p, std::size_t n, std::size_t seed) {
      for (std::size_t i = 0; i < n; ++i) p[i] = dcomplex(0.1 + 0.001 * double((i + seed) % 97), -0.2 + 0.002 * double((i + seed) % 53));
    }
    void fill(auto &g) {
      for (int bl = 0; bl < n_blocks; ++bl) fill_raw(g[bl].data().data(), g[bl].data().size(), std::size_t(bl));
    }

    double checksum() {
      double s = 0;
      for (int b1 = 0; b1 < n_blocks; ++b1)
        for (int b2 = 0; b2 < n_blocks; ++b2) {
          auto *p = M3(b1, b2).data().data();
          for (std::size_t k = 0; k < M3(b1, b2).data().size(); ++k) s += p[k].real() + p[k].imag();
        }
      return s;
    }

    void ph_general(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) measures::simd::full_iw3ph_accumulate(sign, M, GMG, GM, MG, M3, bl1, bl2, GMG(bl2).shape()[0]);
    }
    void pp_general(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) measures::simd::full_iw3pp_accumulate(sign, GM, M3, bl1, bl2, GM[bl2].target_shape()[0]);
    }
    // Exactly the arithmetic of full_iw3ph_accumulate_kernel's N==1 branch, called without the
    // kernel plumbing: isolates the branch itself from the dispatch around it.
    void ph_rank1(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          auto &m3      = M3(bl1, bl2);
          auto *acc     = m3.data().data();
          const auto *m = M[bl1].data().data();
          const long n1 = m3.data().shape()[0], n2 = m3.data().shape()[1];
          const auto c1 = sign * GMG(bl2)(0, 0);
          if (bl1 == bl2) {
            const auto *gm1 = GM[bl1].data().data();
            const auto *mg2 = MG[bl2].data().data();
            for (long i1 = 0; i1 < n1; ++i1) add_scaled_minus_runtime(c1, m + i1 * n2, sign * gm1[i1], mg2, acc + i1 * n2, n2);
          } else {
            add_scaled_runtime(c1, m, acc, n1 * n2);
          }
        }
    }

    // n_orb=1 reference: the arithmetic accumulate_block<1, diagonal> performs, written out with
    // no dispatch, no view construction per mesh point, no batch machinery.
    void ph_scalar(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          auto const M1      = M[bl1];
          auto const GMG2_00 = GMG(bl2)(0, 0);
          auto const GM1     = GM[bl1];
          auto const MG2     = MG[bl2];
          auto &m3           = M3(bl1, bl2);
          if (bl1 == bl2) {
            for (auto mp : m3.mesh()) {
              auto [mp1, mp2] = mp;
              auto iw1 = mp1.value(), iw2 = mp2.value();
              m3[mp](0, 0, 0, 0) += sign * (M1[closest_mesh_pt(iw1, iw2)](0, 0) * GMG2_00 - GM1[iw1](0, 0) * MG2[iw2](0, 0));
            }
          } else {
            for (auto mp : m3.mesh()) {
              auto [mp1, mp2] = mp;
              auto iw1 = mp1.value(), iw2 = mp2.value();
              m3[mp](0, 0, 0, 0) += sign * M1[closest_mesh_pt(iw1, iw2)](0, 0) * GMG2_00;
            }
          }
        }
    }
    void pp_scalar(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          // On the diagonal both terms are M1a(0,0)*M2a(0,0) and cancel exactly (Pauli), so only
          // the off-diagonal pairs contribute anything at n_orb=1.
          if (bl1 == bl2) continue;
          auto const GM1 = GM[bl1];
          auto const GM2 = GM[bl2];
          auto &m3       = M3(bl1, bl2);
          for (auto mp : m3.mesh()) {
            auto [mp1, mp2] = mp;
            m3[mp](0, 0, 0, 0) += sign * GM1[mp1.value()](0, 0) * GM2[mp2.value()](0, 0);
          }
        }
    }
  };

  template <auto body> void run(benchmark::State &st) {
    fixture fx(static_cast<int>(st.range(0)));
    const auto sign = mc_weight_t{1.0};
    for (auto _ : st) {
      body(fx, sign);
      benchmark::DoNotOptimize(fx.M3(0, 0).data().data());
      benchmark::ClobberMemory();
    }
  }

  void iw3ph_full(benchmark::State &st) {
    run<[](fixture &f, mc_weight_t s) { f.ph_general(s); }>(st);
  }
  void iw3pp_full(benchmark::State &st) {
    run<[](fixture &f, mc_weight_t s) { f.pp_general(s); }>(st);
  }
  void iw3ph_full_rank1(benchmark::State &st) {
    run<[](fixture &f, mc_weight_t s) { f.ph_rank1(s); }>(st);
  }
  void iw3ph_full_scalar1(benchmark::State &st) {
    run<[](fixture &f, mc_weight_t s) { f.ph_scalar(s); }>(st);
  }
  void iw3pp_full_scalar1(benchmark::State &st) {
    run<[](fixture &f, mc_weight_t s) { f.pp_scalar(s); }>(st);
  }

  // Equality probe: general vs scalar reference must agree bit-for-bit at n_orb=1.
  void chk(benchmark::State &st) {
    for (auto _ : st) {
      fixture a(1), b(1);
      a.ph_general(1.0);
      b.ph_scalar(1.0);
      st.counters["ph_gen"] = a.checksum();
      st.counters["ph_ref"] = b.checksum();
      fixture e(1);
      e.ph_rank1(1.0);
      st.counters["ph_rk1"] = e.checksum();
      fixture c(1), d(1);
      c.pp_general(1.0);
      d.pp_scalar(1.0);
      st.counters["pp_gen"] = c.checksum();
      st.counters["pp_ref"] = d.checksum();
    }
  }

} // namespace

BENCHMARK(chk)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3ph_full)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3ph_full_rank1)->Arg(1)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3ph_full_scalar1)->Arg(1)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3pp_full)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3pp_full_scalar1)->Arg(1)->Unit(benchmark::kMicrosecond);
