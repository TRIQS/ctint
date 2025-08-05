#include <triqs_ctint/solver_core.hpp>

#include <h5/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs_ctint;

int main(int argc, char *argv[]) {
  mpi::environment env(argc, argv);

  // System Parameters
  double U  = 2.0;
  double mu = U / 2;

  // Discrete bath energies and hoppings
  std::vector<double> energ = {-0.7, -0.15, 0.15, 0.7}; // Vertex Paper ( optimized for Beta=20 ) to fit box DOS
  std::vector<double> hopp  = {0.45, 0.34, 0.34, 0.45};

  // The alpha function
  alpha_t alpha(2);
  double delta = 0.35;
  alpha[0]     = {{0.5 + delta, 0.5 - delta}};
  alpha[1]     = {{0.5 - delta, 0.5 + delta}};

  // Construct Parameters
  constr_params_t pc;
  pc.beta      = 20.0;
  pc.gf_struct = {{"up", 1}, {"down", 1}};
  pc.n_tau     = 10000;
  pc.n_iw      = 500;
  pc.use_D     = false;
  pc.use_Jperp = false;

  // Solve Parameters
  solve_params_t ps;
  ps.h_int                = U * n("up", 0) * n("down", 0);
  ps.n_s                  = 2;
  ps.alpha                = alpha;
  ps.n_cycles             = 1000000;
  ps.length_cycle         = 10;
  ps.n_warmup_cycles      = 10000;
  ps.use_double_insertion = false;
  ps.max_time             = -1;
  ps.post_process         = true;
  ps.n_cheb_coeffs        = 50;

  ps.measure_histogram = true;

  ps.measure_M_tau      = false;
  ps.measure_M_tau_cheb = true;
  ps.measure_M_iw       = false;

  //ps.measure_M4_iw  = true;
  //ps.n_iw_M4        = 32;

  //ps.measure_M3pp_iw  = true;
  //ps.measure_M3ph_iw  = true;
  //ps.n_iw_M3          = 64;

  //ps.measure_M3pp_tau  = true;
  //ps.measure_M3ph_tau  = true;
  //ps.n_tau_M3           = 1000;

  //ps.measure_chi2pp_tau  = true;
  //ps.measure_chi2ph_tau  = true;
  //ps.n_tau_chi2          = 1000;
  //ps.n_iw_chi2           = 128;

  ps.nfft_buf_size = 50;

  solver_core S(pc);

  // set up hybridization delta(i\omega_n)
  auto delta_w = gf<imfreq, matrix_valued>{{pc.beta, Fermion, pc.n_iw}, make_shape(1, 1)};
  delta_w()    = 0;
  for (auto i : range(energ.size())) {
    auto term = gf<imfreq, matrix_valued>{{pc.beta, Fermion, pc.n_iw}, make_shape(1, 1)};
    term(iw_) << hopp[i] * hopp[i] / (iw_ - energ[i]);
    delta_w += term;
  }

  // set up non-interacting Green function g0(i\omega_n)
  S.G0_iw()[0](iw_) << 1.0 / (iw_ + mu + (-1.0) * delta_w(iw_));
  S.G0_iw()[1](iw_) << 1.0 / (iw_ + mu + (-1.0) * delta_w(iw_));

  // Solve!
  S.solve(ps);

  //if (mpi::communicator{}.rank() == 0) {
  //std::cout << "size= " << S.tau_samples[0](0, 0).size() << std::endl;
  //std::cout << "size= " << S.curlyG[0](0, 0).size() << std::endl;
  //for (auto [a, b] : zip(S.tau_samples[0](0, 0), S.weight_samples[0](0, 0))) { std::cout << "dtau= " << a << " weight= " << b << std::endl; }

  //auto archive = h5::file("anderson.out.h5", 'w');
  //h5_write(archive, "S", S);
  //}
}
