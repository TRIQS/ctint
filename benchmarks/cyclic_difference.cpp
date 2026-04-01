#include <benchmark/benchmark.h>
#include <triqs_ctint/solver_core.hpp>

// Stub implementation of the old tau_t (uint32_t-based) and old cyclic_difference
// before the upgrade to uint64_t-based triqs::utility::tau_t.
namespace triqs_ctint_legacy {

  struct tau_t {
    static constexpr std::uint32_t n_max = std::numeric_limits<std::uint32_t>::max();
    inline static double _beta = 0.0;
    std::uint32_t n = 0;
    static void set_beta(double b) { _beta = b; }
    explicit operator double() const { return _beta * n / n_max; }
    bool operator>(const tau_t &tau) const { return n > tau.n; }
    bool operator<=(const tau_t &tau) const { return n <= tau.n; }
  };

  std::pair<double, double> cyclic_difference(tau_t const &tau1, tau_t const &tau2) {
    double const sign  = tau2 > tau1 ? -1.0 : 1.0;
    double const value = static_cast<double>(tau_t{tau1.n - tau2.n});
    return std::make_pair(sign, value);
  }

} // namespace triqs_ctint_legacy

class CyclicDifference : public benchmark::Fixture {
  public:
  std::uint64_t n1;
  std::uint64_t n2;
  double beta;
  void SetUp(benchmark::State const &) {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<std::uint64_t> dist_n(0, std::numeric_limits<std::uint64_t>::max());
    std::uniform_real_distribution<double> dist_beta(0.0, 10.0);
    n1   = dist_n(rng);
    n2   = dist_n(rng);
    beta = dist_beta(rng);
  }
};

BENCHMARK_F(CyclicDifference, bench_old)(benchmark::State &state) {
  using namespace triqs_ctint_legacy;
  tau_t::set_beta(beta);
  tau_t tau1{static_cast<std::uint32_t>(n1)};
  tau_t tau2{static_cast<std::uint32_t>(n2)};
  for (auto _ : state) {
    auto diff = cyclic_difference(tau1, tau2);
    benchmark::DoNotOptimize(diff);
  }
}

BENCHMARK_F(CyclicDifference, bench_new)(benchmark::State &state) {
  using namespace triqs_ctint;
  tau_t::set_beta(beta);
  tau_t tau1{n1};
  tau_t tau2{n2};
  for (auto _ : state) {
    auto diff = cyclic_difference(tau1, tau2);
    benchmark::DoNotOptimize(diff);
  }
}
