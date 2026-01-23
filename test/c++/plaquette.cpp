#include <triqs_ctint/solver_core.hpp>
#include <triqs_ctint/post_process.hpp>
#include <triqs/gfs/functions/dlr2d.hpp>

#include <h5/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs_ctint;

TEST(CtInt, Plaquette) { // NOLINT
  // --------- physical parameters ----------
  double const U    = 1.0;     // Density-density interaction
  double const t    = 1.0;     // Hopping
  double const mu   = U / 2.0; // Chemical Potential
  double const beta = 10;      // Inverse temperature

  // --------- simulation parameters ----------
  int const n_cyc = 1000;

  // --------- Define hopping matrix and interaction hamiltonian ----------

  long const Nx    = 2;
  long const Ny    = 2;
  long const n_orb = Nx * Ny;

  auto hloc0_mat = nda::matrix<double>::zeros(n_orb, n_orb);

  // diagonal terms
  diagonal(hloc0_mat) -= mu;

  auto hloc0_arr = reshape(hloc0_mat, Nx, Ny, Nx, Ny);
  auto _         = nda::range::all;
  for (long x : range(0, Nx - 1)) {
    diagonal(hloc0_arr(x, _, x + 1, _)) -= t;
    diagonal(hloc0_arr(x + 1, _, x, _)) -= t;
  }
  for (long y : range(0, Ny - 1)) {
    diagonal(hloc0_arr(_, y, _, y + 1)) -= t;
    diagonal(hloc0_arr(_, y + 1, _, y)) -= t;
  }

  many_body_operator h_int;
  for (long i = 0; i < n_orb; ++i) { h_int += U * n("up", i) * n("dn", i); }

  // --------- set up block structure ---------

  gf_struct_t gf_struct{{"dn", n_orb}, {"up", n_orb}};

  // --------- Construct the ctint solver ----------
  constr_params_t pc;
  pc.beta      = beta;
  pc.gf_struct = gf_struct;
  pc.n_iw      = 21;
  pc.n_tau     = 201;

  solver_core S(pc);

  // --------- Initialize the non-interacting Green's function ----------
  for (auto &&g_bl : S.G0_iw) {
    g_bl(iw_) << iw_ - hloc0_mat;
    g_bl = inverse(g_bl);
  }

  // Define the alpha tensor
  auto alpha_bl = nda::array<double, 2>(n_orb, 2);
  auto alpha    = alpha_t{alpha_bl, alpha_bl};
  double delta   = 0.1;
  alpha[0](_, 0) = 0.5 + delta;
  alpha[0](_, 1) = 0.5 - delta;
  alpha[1](_, 0) = 0.5 - delta;
  alpha[1](_, 1) = 0.5 + delta;

  // --------- Solve! ----------

  solve_params_t ps;
  ps.h_int              = h_int;
  ps.n_s                = 2;
  ps.alpha              = alpha;
  ps.n_cycles           = n_cyc;
  ps.length_cycle       = 50;
  ps.n_warmup_cycles    = 100;
  ps.random_seed        = 34788;
  ps.measure_histogram  = true;
  ps.measure_density    = true;
  ps.measure_M_tau      = true;
  ps.measure_M_iw       = true;
  ps.measure_M4_iw      = true;
  ps.measure_M4pp_iw    = true;
  ps.measure_M4ph_iw    = true;
  ps.n_iw_M4            = 2;
  ps.n_iW_M4            = 2;
  ps.nfft_buf_size      = 100000;
  ps.measure_M3pp_iw    = true;
  ps.measure_M3ph_iw    = true;
  ps.measure_M3pp_tau   = false;
  ps.measure_M3ph_tau   = false;
  ps.measure_M3xph_tau  = false;
  ps.dlr_wmax_M3        = 1.0;
  ps.n_iw_M3            = 4;
  ps.n_iW_M3            = 4;
  ps.n_tau_M3           = 4;
  ps.measure_chi2pp_tau = true;
  ps.measure_chi2ph_tau = true;
  ps.n_iw_chi2          = 10;
  ps.n_tau_chi2         = 21;
  ps.measure_chiAB_tau  = true;
  ps.chi_A_vec          = {n("up", 0) + n("dn", 0)};
  ps.chi_B_vec          = {n("up", 0) + n("dn", 0)};
  ps.post_process       = true;

  S.solve(ps);

  mpi::communicator world;
  int n_blocks = static_cast<int>(gf_struct.size());

  // -------- Compare DLR2D vs full-mesh M3pp ---------
  if (world.rank() == 0) std::cout << "Comparing DLR2D vs full-mesh M3pp ..." << std::endl;
  {
    auto const &M3pp_dlr2d = S.M3pp_iw_nfft.value();
    auto const &M3pp_full  = S.M3pp_iw_nfft_full.value();
    auto M3pp_ref          = M3pp_dlr2d;
    for (int bl1 = 0; bl1 < n_blocks; ++bl1)
      for (int bl2 = 0; bl2 < n_blocks; ++bl2)
        for (auto mp : M3pp_ref(bl1, bl2).mesh()) {
          auto [iw1, iw2] = mp.value();
          M3pp_ref(bl1, bl2)[mp] = M3pp_full(bl1, bl2)[iw1, iw2];
        }
    EXPECT_BLOCK2_GF_NEAR(M3pp_dlr2d, M3pp_ref, 1e-14);
  }

  // -------- Compare DLR2D vs full-mesh M3ph ---------
  if (world.rank() == 0) std::cout << "Comparing DLR2D vs full-mesh M3ph ..." << std::endl;
  {
    auto const &M3ph_dlr2d = S.M3ph_iw_nfft.value();
    auto const &M3ph_full  = S.M3ph_iw_nfft_full.value();
    auto M3ph_ref          = M3ph_dlr2d;
    for (int bl1 = 0; bl1 < n_blocks; ++bl1)
      for (int bl2 = 0; bl2 < n_blocks; ++bl2)
        for (auto mp : M3ph_ref(bl1, bl2).mesh()) {
          auto [iw1, iw2] = mp.value();
          M3ph_ref(bl1, bl2)[mp] = M3ph_full(bl1, bl2)[iw1, iw2];
        }
    EXPECT_BLOCK2_GF_NEAR(M3ph_dlr2d, M3ph_ref, 1e-14);
  }

  // -------- Compute full-mesh chi3 from full-mesh M3 ---------
  if (world.rank() == 0) std::cout << "Computing chi3pp from full-mesh M3pp ..." << std::endl;
  auto chi3pp_full = chi3_from_M3<Chan_t::PP>(S.M3pp_iw_nfft_full.value(), S.M_iw.value(), S.G0_shift_iw, S.density.value(), S.M_hartree.value());
  if (world.rank() == 0) std::cout << "Computing chi3ph from full-mesh M3ph ..." << std::endl;
  auto chi3ph_full = chi3_from_M3<Chan_t::PH>(S.M3ph_iw_nfft_full.value(), S.M_iw.value(), S.G0_shift_iw, S.density.value(), S.M_hartree.value());

  // -------- Compute DLR2D interpolated chi3 on full mesh ---------
  if (world.rank() == 0) std::cout << "Computing DLR2D interpolated chi3pp on full mesh ..." << std::endl;
  auto chi3pp_coefs = make_gf_dlr2d(S.chi3pp_iw_nfft.value());
  auto chi3pp_eval  = make_gf_imfreq(chi3pp_coefs);

  if (world.rank() == 0) std::cout << "Computing DLR2D interpolated chi3ph on full mesh ..." << std::endl;
  auto chi3ph_coefs = make_gf_dlr2d(S.chi3ph_iw_nfft.value());
  auto chi3ph_eval  = make_gf_imfreq(chi3ph_coefs);

  // -------- Compare DLR2D vs full-mesh chi3 at DLR2D points ---------
  if (world.rank() == 0) std::cout << "Comparing DLR2D vs full-mesh chi3 at DLR2D points ..." << std::endl;
  {
    auto chi3pp_ref = S.chi3pp_iw_nfft.value();
    for (int bl1 = 0; bl1 < n_blocks; ++bl1)
      for (int bl2 = 0; bl2 < n_blocks; ++bl2)
        for (auto mp : chi3pp_ref(bl1, bl2).mesh()) {
          auto [iw1, iw2] = mp.value();
          chi3pp_ref(bl1, bl2)[mp] = chi3pp_full(bl1, bl2)[iw1, iw2];
        }
    EXPECT_BLOCK2_GF_NEAR(S.chi3pp_iw_nfft.value(), chi3pp_ref, 1e-12);
  }
  {
    auto chi3ph_ref = S.chi3ph_iw_nfft.value();
    for (int bl1 = 0; bl1 < n_blocks; ++bl1)
      for (int bl2 = 0; bl2 < n_blocks; ++bl2)
        for (auto mp : chi3ph_ref(bl1, bl2).mesh()) {
          auto [iw1, iw2] = mp.value();
          chi3ph_ref(bl1, bl2)[mp] = chi3ph_full(bl1, bl2)[iw1, iw2];
        }
    EXPECT_BLOCK2_GF_NEAR(S.chi3ph_iw_nfft.value(), chi3ph_ref, 1e-12);
  }

  // -------- Compare DLR2D interpolated chi3 vs full-mesh chi3 ---------
  if (world.rank() == 0) std::cout << "Comparing DLR2D interpolated chi3 vs full-mesh chi3 ..." << std::endl;
  double expected_tol = 40 / std::sqrt(n_cyc * world.size());
  EXPECT_BLOCK2_GF_NEAR(chi3pp_eval, chi3pp_full, expected_tol);
  EXPECT_BLOCK2_GF_NEAR(chi3ph_eval, chi3ph_full, expected_tol);

  // -------- Save in archive ---------
  if (world.rank() == 0) std::cout << "Saving results to HDF5 ..." << std::endl;
  if (world.rank() == 0) {
    auto archive = h5::file("plaquette.out.h5", 'w');
    h5_write(archive, "histogram", S.histogram);
    h5_write(archive, "density", S.density);
    h5_write(archive, "M_tau", S.M_tau);
    h5_write(archive, "M_iw", S.M_iw_nfft);
    h5_write(archive, "M4_iw", S.M4_iw);
    h5_write(archive, "M4pp_iw", S.M4pp_iw);
    h5_write(archive, "M4ph_iw", S.M4ph_iw);
    h5_write(archive, "M3pp_iw_nfft", S.M3pp_iw_nfft);
    h5_write(archive, "M3ph_iw_nfft", S.M3ph_iw_nfft);
    h5_write(archive, "M3pp_iw_nfft_full", S.M3pp_iw_nfft_full);
    h5_write(archive, "M3ph_iw_nfft_full", S.M3ph_iw_nfft_full);
    h5_write(archive, "M3pp_tau", S.M3pp_tau);
    h5_write(archive, "M3ph_tau", S.M3ph_tau);
    h5_write(archive, "M3xph_tau", S.M3xph_tau);
    h5_write(archive, "chi3pp_iw_nfft", S.chi3pp_iw_nfft);
    h5_write(archive, "chi3ph_iw_nfft", S.chi3ph_iw_nfft);
    h5_write(archive, "chi3pp_full", chi3pp_full);
    h5_write(archive, "chi3ph_full", chi3ph_full);
    h5_write(archive, "chi3pp_eval", chi3pp_eval);
    h5_write(archive, "chi3ph_eval", chi3ph_eval);
    h5_write(archive, "chi2pp_tau", S.chi2pp_tau);
    h5_write(archive, "chi2ph_iw", S.chi2ph_tau);
    h5_write(archive, "chiAB_tau", S.chiAB_tau);
  }
}

MAKE_MAIN;
