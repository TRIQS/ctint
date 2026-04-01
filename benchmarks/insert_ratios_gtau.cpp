#include <benchmark/benchmark.h>
#include <triqs_ctint/solver_core.hpp>

// Benchmark the G0_shift_tau evaluation pattern inside insert_ratios/insert2_ratios.
// This isolates the bottleneck observed in the 8x8 Hubbard chiAB_tau measurement.
//
// The dominant cost is building the B(N,K) and C(K,N) matrices inside det_manip,
// where each element requires evaluating G0hat_t::operator() which does:
//   cyclic_difference → closest_mesh_pt → gf data lookup → matrix element access

using namespace triqs_ctint;
using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace nda;

// Build a random G0_shift_tau on a uniform imtime mesh
static auto make_g0_shift_tau(double beta, long n_tau, long n_orb) {
  auto g0 = gf<imtime, matrix_real_valued>{{beta, Fermion, n_tau}, {n_orb, n_orb}};
  // Fill with plausible decaying values
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (auto tau : g0.mesh())
    for (int a = 0; a < n_orb; ++a)
      for (int b = 0; b < n_orb; ++b) g0[tau](a, b) = dist(rng) * std::exp(-double(tau) / beta);
  return g0;
}

// Fill a det_manip to a given size by inserting random operators
static void fill_det(det_t &det, long target_size, long n_orb, std::mt19937_64 &rng) {
  std::uniform_int_distribution<std::uint64_t> tau_dist(10, tau_t::n_max - 10);
  std::uniform_int_distribution<int> orb_dist(0, n_orb - 1);

  for (long i = 0; i < target_size; ++i) {
    auto tau1 = tau_t{tau_dist(rng)};
    auto tau2 = tau_t{tau_dist(rng)};
    auto c    = c_t{tau1, orb_dist(rng)};
    auto cdag = cdag_t{tau2, orb_dist(rng)};
    auto ratio = det.try_insert(i, i, c, cdag);
    if (std::abs(ratio) < 1e-15) {
      det.reject_last_try();
    } else {
      det.complete_operation();
    }
  }
}

class InsertRatiosGtau : public benchmark::Fixture {
  public:
  static constexpr double beta  = 2.0;
  static constexpr long n_tau   = 4001; // typical for 8x8 case
  static constexpr long n_orb   = 64;   // 8x8 cluster
  static constexpr long det_N   = 500;  // typical perturbation order
  static constexpr long L       = 35;   // DLR tau points
  static constexpr long E       = 128;  // orbital pairs per group

  gf<imtime, matrix_real_valued> g0;
  std::unique_ptr<det_t> det;
  std::mt19937_64 rng{12345};

  void SetUp(benchmark::State const &) {
    tau_t::set_beta(beta);
    g0 = make_g0_shift_tau(beta, n_tau, n_orb);

    // Create alpha (trivial: zeros)
    alpha_t alpha(1, 2, 2, 1);
    alpha() = 0.0;

    det = std::make_unique<det_t>(G0hat_t{g0, alpha}, 1000);
    fill_det(*det, det_N, n_orb, rng);
  }

  void TearDown(benchmark::State const &) { det.reset(); }
};

// Benchmark insert_ratios with rank-2 array (L, E) — mimics AABB A-side
BENCHMARK_F(InsertRatiosGtau, insert_ratios_rank2)(benchmark::State &state) {
  std::uniform_int_distribution<std::uint64_t> tau_dist(10, tau_t::n_max - 10);
  std::uniform_int_distribution<int> orb_dist(0, n_orb - 1);

  // Build rank-2 arrays (L, E)
  nda::array<c_t, 2> cs(L, E);
  nda::array<cdag_t, 2> cdags(L, E);
  for (long l = 0; l < L; ++l) {
    auto tau_c    = tau_t{tau_dist(rng)};
    auto tau_cdag = tau_t{tau_dist(rng)};
    for (long e = 0; e < E; ++e) {
      cs(l, e)    = c_t{tau_c, orb_dist(rng)};
      cdags(l, e) = cdag_t{tau_cdag, orb_dist(rng)};
    }
  }

  for (auto _ : state) {
    auto ratios = det->insert_ratios(0, 0, cs, cdags);
    benchmark::DoNotOptimize(ratios.data());
  }
}

// Benchmark insert2_ratios with broadcast (L,E) x (E) — mimics AAAA case
BENCHMARK_F(InsertRatiosGtau, insert2_ratios_broadcast)(benchmark::State &state) {
  std::uniform_int_distribution<std::uint64_t> tau_dist(10, tau_t::n_max - 10);
  std::uniform_int_distribution<int> orb_dist(0, n_orb - 1);

  // A-side: rank-2 (L, E) — tau varies
  nda::array<c_t, 2> c_A(L, E);
  nda::array<cdag_t, 2> cdag_A(L, E);
  for (long l = 0; l < L; ++l) {
    auto tau_c    = tau_t{tau_dist(rng)};
    auto tau_cdag = tau_t{tau_dist(rng)};
    for (long e = 0; e < E; ++e) {
      c_A(l, e)    = c_t{tau_c, orb_dist(rng)};
      cdag_A(l, e) = cdag_t{tau_cdag, orb_dist(rng)};
    }
  }

  // B-side: rank-1 (E) — fixed tau
  nda::array<c_t, 1> c_B(E);
  nda::array<cdag_t, 1> cdag_B(E);
  for (long e = 0; e < E; ++e) {
    c_B(e)    = c_t{tau_t::zero(), orb_dist(rng)};
    cdag_B(e) = cdag_t{tau_t::epsilon(), orb_dist(rng)};
  }

  for (auto _ : state) {
    auto ratios = det->insert2_ratios(0, 1, 0, 1, c_A, c_B, cdag_A, cdag_B);
    benchmark::DoNotOptimize(ratios.data());
  }
}

// Benchmark just the f-evaluation loop (B and C matrix construction) in isolation
BENCHMARK_F(InsertRatiosGtau, f_evaluation_loop)(benchmark::State &state) {
  std::uniform_int_distribution<std::uint64_t> tau_dist(10, tau_t::n_max - 10);
  std::uniform_int_distribution<int> orb_dist(0, n_orb - 1);

  long K = L * E;
  std::vector<c_t> xs(K);
  std::vector<cdag_t> ys(K);
  for (long m = 0; m < K; ++m) {
    xs[m] = c_t{tau_t{tau_dist(rng)}, orb_dist(rng)};
    ys[m] = cdag_t{tau_t{tau_dist(rng)}, orb_dist(rng)};
  }

  // Get the f functor from the det_manip
  auto const &f = det->get_function();
  long N = det->size();

  // Simulate B matrix build: f(x_values[l], ys[m]) for l in [0,N), m in [0,K)
  // We don't have access to x_values directly, so we use the new points as proxies
  for (auto _ : state) {
    double sum = 0;
    for (long l = 0; l < N; ++l)
      for (long m = 0; m < K; ++m) sum += f(xs[l % K], ys[m]);
    benchmark::DoNotOptimize(sum);
  }
}
