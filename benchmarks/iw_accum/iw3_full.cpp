// Correctness probe and micro-benchmark for the M3 full-mesh iw-accumulate functions
// (M3ph_iw_full.cpp, M3pp_iw_full.cpp).
//
// chk_all_n checks both kernels element-wise against a general-N reference written straight from
// the update formula, for every block size the dispatch has a case for plus one beyond it: the
// SIMD branches only engage from n_orb >= 3 (prefer_simd) and the runtime-length path only at
// n_orb > max_block, so the sweep is what covers them.
//
// The timing side: fef408b deleted the n_orb=1 scalar fast paths that 1a30154 had added to
// full_iw3ph/full_iw3pp ("~20% slower than the pre-SIMD scalar implementation for
// single-orbital models"). It measures the shipped general-N kernel against a hand-written
// n_orb=1 loop with identical arithmetic, in one binary so the two alternate under the same
// machine conditions.
#include "../../c++/triqs_ctint/measures/M3ph_iw_full.hpp"
#include "../../c++/triqs_ctint/measures/M3pp_iw_full.hpp"
#include "../../c++/triqs_ctint/measures/iw_simd.hpp"

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <benchmark/benchmark.h>

#include <cmath>
#include <complex>
#include <string>
#include <utility>

using namespace triqs_ctint;
using namespace triqs_ctint::measures;

namespace {

  // Cost scales as (2*n_iw)^2 * n_blocks^2 * n_orb^4; M3 itself is (2*n_iw)^2 * n_orb^4 complex
  // per block pair, so n_iw stays small to keep the container in cache-plus-a-bit rather than GB.
  constexpr double beta = 100.0;
  constexpr int n_iw    = 32;
  // The general-N probe walks a scalar n_orb^4 loop per mesh point, so it uses a smaller mesh.
  constexpr int probe_n_iw = 6;

  struct fixture {
    gf_struct_t gf_struct;
    int n_blocks;

    // The same container types the real measures use, so the transposed GM/MG layout the
    // accumulation relies on cannot drift out of sync with them.
    M3_M_full_t M;
    M3_G_t GM, MG;
    M3_GMG_t GMG;
    chi3_iw_t M3;

    // n_iw_ shrinks the mesh for the general-N probe: point count scales its cost but not the
    // branch coverage.
    explicit fixture(int n_orb, int n_iw_ = n_iw) : gf_struct{{"dn", n_orb}, {"up", n_orb}}, n_blocks(2) {
      mesh::imfreq iw_mesh{beta, Fermion, n_iw_};
      mesh::prod iw2_mesh{iw_mesh, iw_mesh};

      M3   = make_block2_gf(iw2_mesh, gf_struct);
      M3() = 0;
      M    = M3_M_full_t{iw2_mesh, gf_struct};
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

    void ph_general(const mc_weight_t sign, chi3_iw_v_t &m3) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) measures::full_iw3ph_accumulate(sign, M, GMG, GM, MG, m3, bl1, bl2, GMG(bl2).shape()[0]);
    }
    void pp_general(const mc_weight_t sign, chi3_iw_v_t &m3) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) measures::full_iw3pp_accumulate(sign, GM, m3, bl1, bl2, GM[bl2].target_shape()[0]);
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
            for (long i1 = 0; i1 < n1; ++i1) {
              const auto c2 = sign * gm1[i1];
              for (long i2 = 0; i2 < n2; ++i2) acc[i1 * n2 + i2] += c1 * m[i1 * n2 + i2] - c2 * mg2[i2];
            }
          } else {
            for (long k = 0; k < n1 * n2; ++k) acc[k] += c1 * m[k];
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
    // acc(i,j,k,l) += sign * M1a(j,i) * M2a(l,k) - sign * M1b(l,i) * M2b(j,k), the second term on
    // the diagonal only: the update the kernels implement, indexed straight out of the operands,
    // for any block size. The n_orb=1 references above only reach the scalar branch.
    void ph_reference(const mc_weight_t sign) {
      const long n_orb = M3(0, 0).target_shape()[0];
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          auto const &M1   = M[bl1];
          auto const &GMG2 = GMG(bl2);
          auto const &GM1  = GM[bl1];
          auto const &MG2  = MG[bl2];
          auto &m3         = M3(bl1, bl2);
          const bool diag  = bl1 == bl2;
          for (auto mp : m3.mesh()) {
            auto [mp1, mp2] = mp;
            const auto iw1 = mp1.value(), iw2 = mp2.value();
            auto const M1a = M1[closest_mesh_pt(iw1, iw2)];
            auto const M1b = GM1[iw1];
            auto const M2b = MG2[iw2];
            auto acc       = m3[mp];
            for (long i = 0; i < n_orb; ++i)
              for (long j = 0; j < n_orb; ++j)
                for (long k = 0; k < n_orb; ++k)
                  for (long l = 0; l < n_orb; ++l) {
                    auto v = (M1a(j, i) * sign) * GMG2(l, k);
                    if (diag) v -= (M2b(j, k) * sign) * M1b(l, i);
                    acc(i, j, k, l) += v;
                  }
          }
        }
    }
    void pp_reference(const mc_weight_t sign) {
      const long n_orb = M3(0, 0).target_shape()[0];
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          auto const &GM1 = GM[bl1];
          auto const &GM2 = GM[bl2];
          auto &m3        = M3(bl1, bl2);
          const bool diag = bl1 == bl2;
          for (auto mp : m3.mesh()) {
            auto [mp1, mp2] = mp;
            auto const M1a  = GM1[mp1.value()];
            auto const M2a  = GM2[mp2.value()];
            auto acc        = m3[mp];
            for (long i = 0; i < n_orb; ++i)
              for (long j = 0; j < n_orb; ++j)
                for (long k = 0; k < n_orb; ++k)
                  for (long l = 0; l < n_orb; ++l) {
                    auto v = (M1a(j, i) * sign) * M2a(l, k);
                    if (diag) v -= (M2a(j, k) * sign) * M1a(l, i);
                    acc(i, j, k, l) += v;
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
    chi3_iw_v_t m3{fx.M3};
    const auto sign = mc_weight_t{1.0};
    for (auto _ : st) {
      body(fx, m3, sign);
      benchmark::DoNotOptimize(fx.M3(0, 0).data().data());
      benchmark::ClobberMemory();
    }
  }

  void iw3ph_full(benchmark::State &st) {
    run<[](fixture &f, chi3_iw_v_t &m3, mc_weight_t s) { f.ph_general(s, m3); }>(st);
  }
  void iw3pp_full(benchmark::State &st) {
    run<[](fixture &f, chi3_iw_v_t &m3, mc_weight_t s) { f.pp_general(s, m3); }>(st);
  }
  void iw3ph_full_rank1(benchmark::State &st) {
    run<[](fixture &f, chi3_iw_v_t &, mc_weight_t s) { f.ph_rank1(s); }>(st);
  }
  void iw3ph_full_scalar1(benchmark::State &st) {
    run<[](fixture &f, chi3_iw_v_t &, mc_weight_t s) { f.ph_scalar(s); }>(st);
  }
  void iw3pp_full_scalar1(benchmark::State &st) {
    run<[](fixture &f, chi3_iw_v_t &, mc_weight_t s) { f.pp_scalar(s); }>(st);
  }

  // Element-wise count, but the size of the disagreement is measured against the scale of the
  // result: the pp diagonal cancels to exact zero in places, where a per-element ratio against a
  // roundoff-level reference value reports 1 or 2 for a difference of one ulp.
  std::pair<long, double> compare(fixture const &a, fixture const &b) {
    long bad      = 0;
    double absmax = 0, scale = 0;
    for (int b1 = 0; b1 < a.n_blocks; ++b1)
      for (int b2 = 0; b2 < a.n_blocks; ++b2) {
        auto const &x = a.M3(b1, b2).data();
        auto const &y = b.M3(b1, b2).data();
        for (long k = 0; k < long(x.size()); ++k) {
          const dcomplex u = x.data()[k], v = y.data()[k];
          scale = std::max(scale, std::abs(v));
          if (u == v) continue;
          ++bad;
          absmax = std::max(absmax, std::abs(u - v));
        }
      }
    return {bad, absmax / std::max(scale, 1e-300)};
  }

  // General-N equality probe: every block size the dispatch has a case for, plus one beyond it.
  void chk_all_n(benchmark::State &st) {
    for (auto _ : st) {
      for (int n_orb = 1; n_orb <= max_block + 1; ++n_orb) {
        const auto tag = std::to_string(n_orb);
        fixture a(n_orb, probe_n_iw), b(n_orb, probe_n_iw);
        chi3_iw_v_t av{a.M3};
        a.ph_general(1.0, av);
        b.ph_reference(1.0);
        auto [ph_bad, ph_rel]            = compare(a, b);
        st.counters["ph" + tag + "_ne"]  = double(ph_bad);
        st.counters["ph" + tag + "_rel"] = ph_rel;

        fixture c(n_orb, probe_n_iw), d(n_orb, probe_n_iw);
        chi3_iw_v_t cv{c.M3};
        c.pp_general(1.0, cv);
        d.pp_reference(1.0);
        auto [pp_bad, pp_rel]            = compare(c, d);
        st.counters["pp" + tag + "_ne"]  = double(pp_bad);
        st.counters["pp" + tag + "_rel"] = pp_rel;
      }
    }
  }

  // Equality probe: general vs scalar reference must agree bit-for-bit at n_orb=1.
  void chk(benchmark::State &st) {
    for (auto _ : st) {
      fixture a(1), b(1);
      chi3_iw_v_t av{a.M3};
      a.ph_general(1.0, av);
      b.ph_scalar(1.0);
      st.counters["ph_gen"] = a.checksum();
      st.counters["ph_ref"] = b.checksum();
      fixture e(1);
      e.ph_rank1(1.0);
      st.counters["ph_rk1"] = e.checksum();
      fixture c(1), d(1);
      chi3_iw_v_t cv{c.M3};
      c.pp_general(1.0, cv);
      d.pp_scalar(1.0);
      st.counters["pp_gen"] = c.checksum();
      st.counters["pp_ref"] = d.checksum();
    }
  }

} // namespace

BENCHMARK(chk)->Unit(benchmark::kMicrosecond);
BENCHMARK(chk_all_n)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(iw3ph_full)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3ph_full_rank1)->Arg(1)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3ph_full_scalar1)->Arg(1)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3pp_full)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3pp_full_scalar1)->Arg(1)->Unit(benchmark::kMicrosecond);
