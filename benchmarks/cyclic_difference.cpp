#include <benchmark/benchmark.h>
#include <triqs_ctint/solver_core.hpp>

// Stub implementation of the necessary interface components of the old tau_t
// and old cyclic_difference before the improvements.
namespace triqs_ctint_legacy {

  struct tau_t {
    static constexpr uint32_t const n_max = std::numeric_limits<uint32_t>::max();
    static double beta;
    uint32_t n = 0;
    explicit operator double() const { return beta * n / n_max; }
    bool operator>(const tau_t &tau) const { return n > tau.n; }
    bool operator<=(const tau_t &tau) const { return n <= tau.n; }
  };

  double triqs_ctint_legacy::tau_t::beta = 0.0;

  std::pair<double, double> cyclic_difference(tau_t const &tau1, tau_t const &tau2) {
    // Assume tau1 > tau2 for the equal time case
    double const sign  = tau2 > tau1 ? -1.0 : 1.0;
    double const value = static_cast<double>(tau_t{tau1.n - tau2.n});
    return std::make_pair(sign, value);
  }

} // namespace triqs_ctint_legacy

class CyclicDifference : public benchmark::Fixture {
  public:
  uint32_t n1;
  uint32_t n2;
  double beta;
  void SetUp(benchmark::State const &) {
    // Generate random parameters (with fixed seed)
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<uint32_t> dist_n(0, UINT32_MAX);
    std::uniform_real_distribution<double> dist_beta(0.0, 10.0);
    n1   = dist_n(rng);
    n2   = dist_n(rng);
    beta = dist_beta(rng);
  }
};

BENCHMARK_F(CyclicDifference, bench_old)(benchmark::State &state) {
  using namespace triqs_ctint_legacy;
  tau_t::beta = beta;
  tau_t tau1{n1};
  tau_t tau2{n2};
  for (auto _ : state) {
    auto diff = cyclic_difference(tau1, tau2);
    benchmark::DoNotOptimize(diff);
  }
}

BENCHMARK_F(CyclicDifference, bench_new)(benchmark::State &state) {
  using namespace triqs_ctint;
  tau_t::beta = beta;
  tau_t tau1{n1};
  tau_t tau2{n2};
  for (auto _ : state) {
    auto diff = cyclic_difference(tau1, tau2);
    benchmark::DoNotOptimize(diff);
  }
}
