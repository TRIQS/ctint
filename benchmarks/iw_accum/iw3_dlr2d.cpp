// Correctness probe and micro-benchmark for the M3 DLR2D iw-accumulate functions
// (M3ph_iw.cpp, M3pp_iw.cpp).
//
// The probe checks both kernels element-wise against a general-N reference written straight from
// the update formula, for every block size the dispatch has a case for plus one beyond it. The
// SIMD branches only engage from n_orb >= 3 (prefer_simd) and the runtime-length path only at
// n_orb > max_block, so the sweep is what covers them.
#include "../../c++/triqs_ctint/measures/M3ph_iw.hpp"
#include "../../c++/triqs_ctint/measures/M3pp_iw.hpp"
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

  // Small mesh: the reference is a scalar n_orb^4 loop per mesh point, and the point count only
  // scales the runtime, not the branch coverage.
  constexpr double beta = 10.0;
  constexpr double wmax = 2.0;
  constexpr double eps  = 1e-6;

  struct fixture {
    gf_struct_t gf_struct;
    int n_blocks;
    int n_orb;

    // The same container types the measures use, so the transposed layouts the accumulation
    // relies on cannot drift out of sync with them.
    M3_M_t M;
    M3_G_t GM, MG;
    M3_GMG_t GMG;
    chi3_dlr2d_iw_t M3;

    fixture(int n_orb_, auto channel) : gf_struct{{"dn", n_orb_}, {"up", n_orb_}}, n_blocks(2), n_orb(n_orb_) {
      mesh::dlr2d_imfreq mesh2{beta, wmax, eps, channel, false};
      mesh::imfreq iw_mesh{beta, Fermion, mesh2.max_n() + 1};

      M3   = make_block2_gf(mesh2, gf_struct);
      M3() = 0;
      M    = M3_M_t{mesh2, gf_struct};
      GM   = M3_G_t{iw_mesh, gf_struct};
      MG   = M3_G_t{iw_mesh, gf_struct};
      GMG  = array_adapter{make_shape(n_blocks), [&](int) { return array<dcomplex, 2, nda::F_layout>(n_orb, n_orb); }};

      fill(M);
      fill(GM);
      fill(MG);
      for (int bl = 0; bl < n_blocks; ++bl) fill_raw(GMG(bl).data(), GMG(bl).size(), std::size_t(bl) + 7);
    }

    // Deterministic non-trivial fill so the arithmetic is real and reproducible.
    static void fill_raw(dcomplex *p, std::size_t n, std::size_t seed) {
      for (std::size_t i = 0; i < n; ++i) p[i] = dcomplex(0.1 + 0.001 * double((i + seed) % 97), -0.2 + 0.002 * double((i + seed) % 53));
    }
    void fill(auto &g) {
      for (int bl = 0; bl < n_blocks; ++bl) fill_raw(g[bl].data().data(), g[bl].data().size(), std::size_t(bl));
    }

    long mesh_size() const { return M3(0, 0).mesh().size(); }

    void ph_general(const mc_weight_t sign, chi3_dlr2d_iw_v_t &m3) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) dlr2d_iw3ph_accumulate(sign, M, GMG, GM, MG, m3, bl1, bl2, GMG(bl2).shape()[0]);
    }
    void pp_general(const mc_weight_t sign, chi3_dlr2d_iw_v_t &m3) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) dlr2d_iw3pp_accumulate(sign, GM, m3, bl1, bl2, GM[bl2].target_shape()[0]);
    }

    // acc(i,j,k,l) += sign * M1a(j,i) * M2a(l,k) - sign * M1b(l,i) * M2b(j,k), the second term on
    // the diagonal only: the update the kernels implement, indexed straight out of the operands.
    void ph_reference(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          auto const &M1   = M[bl1];
          auto const &GMG2 = GMG(bl2);
          auto const &GM1  = GM[bl1];
          auto const &MG2  = MG[bl2];
          auto &m3         = M3(bl1, bl2);
          const bool diag  = bl1 == bl2;
          for (auto mp : m3.mesh()) {
            auto [iw1, iw2] = mp.value();
            auto const M1a  = M1[mp];
            auto const M1b  = GM1[iw1];
            auto const M2b  = MG2[iw2];
            auto acc        = m3[mp];
            for (int i = 0; i < n_orb; ++i)
              for (int j = 0; j < n_orb; ++j)
                for (int k = 0; k < n_orb; ++k)
                  for (int l = 0; l < n_orb; ++l) {
                    auto v = (M1a(j, i) * sign) * GMG2(l, k);
                    if (diag) v -= (M2b(j, k) * sign) * M1b(l, i);
                    acc(i, j, k, l) += v;
                  }
          }
        }
    }
    void pp_reference(const mc_weight_t sign) {
      for (int bl1 = 0; bl1 < n_blocks; ++bl1)
        for (int bl2 = 0; bl2 < n_blocks; ++bl2) {
          auto const &GM1 = GM[bl1];
          auto const &GM2 = GM[bl2];
          auto &m3        = M3(bl1, bl2);
          const bool diag = bl1 == bl2;
          for (auto mp : m3.mesh()) {
            auto [iw1, iw2] = mp.value();
            auto const M1a  = GM1[iw1];
            auto const M2a  = GM2[iw2];
            auto acc        = m3[mp];
            for (int i = 0; i < n_orb; ++i)
              for (int j = 0; j < n_orb; ++j)
                for (int k = 0; k < n_orb; ++k)
                  for (int l = 0; l < n_orb; ++l) {
                    auto v = (M1a(j, i) * sign) * M2a(l, k);
                    if (diag) v -= (M2a(j, k) * sign) * M1a(l, i);
                    acc(i, j, k, l) += v;
                  }
          }
        }
    }
  };

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

  void chk(benchmark::State &st) {
    for (auto _ : st) {
      for (int n_orb = 1; n_orb <= max_block + 1; ++n_orb) {
        const auto tag = std::to_string(n_orb);
        fixture a(n_orb, mesh::PH), b(n_orb, mesh::PH);
        chi3_dlr2d_iw_v_t av{a.M3};
        a.ph_general(1.0, av);
        b.ph_reference(1.0);
        auto [ph_bad, ph_rel]            = compare(a, b);
        st.counters["ph" + tag + "_ne"]  = double(ph_bad);
        st.counters["ph" + tag + "_rel"] = ph_rel;

        fixture c(n_orb, mesh::PP), d(n_orb, mesh::PP);
        chi3_dlr2d_iw_v_t cv{c.M3};
        c.pp_general(1.0, cv);
        d.pp_reference(1.0);
        auto [pp_bad, pp_rel]            = compare(c, d);
        st.counters["pp" + tag + "_ne"]  = double(pp_bad);
        st.counters["pp" + tag + "_rel"] = pp_rel;
        st.counters["npt"]               = double(a.mesh_size());
      }
    }
  }

  template <auto body> void run(benchmark::State &st) {
    fixture fx(static_cast<int>(st.range(0)), mesh::PH);
    chi3_dlr2d_iw_v_t m3{fx.M3};
    const auto sign = mc_weight_t{1.0};
    for (auto _ : st) {
      body(fx, m3, sign);
      benchmark::DoNotOptimize(fx.M3(0, 0).data().data());
      benchmark::ClobberMemory();
    }
  }
  void iw3ph_dlr2d(benchmark::State &st) {
    run<[](fixture &f, chi3_dlr2d_iw_v_t &m3, mc_weight_t s) { f.ph_general(s, m3); }>(st);
  }
  void iw3pp_dlr2d(benchmark::State &st) {
    run<[](fixture &f, chi3_dlr2d_iw_v_t &m3, mc_weight_t s) { f.pp_general(s, m3); }>(st);
  }

} // namespace

BENCHMARK(chk)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(iw3ph_dlr2d)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Unit(benchmark::kMicrosecond);
BENCHMARK(iw3pp_dlr2d)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Unit(benchmark::kMicrosecond);
