// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include <benchmark/benchmark.h>
#include <triqs_ctint/nfft/buffer.hpp>
#include <triqs/mesh.hpp>
#include <triqs/mesh/matsubara_freq.hpp>
#include <random>
#include <set>

using namespace triqs::utility;
using namespace triqs::utility::nfft;
using namespace triqs::mesh;
using dcomplex = std::complex<double>;

// Parameters from plaquette test
static constexpr double beta     = 20.0;
static constexpr double dlr_wmax = 1.0;
static constexpr double dlr_eps  = 1e-6;
static constexpr int buf_size    = 100000;
static constexpr int bl_size     = 4;
static constexpr double tol      = 1e-8;

// Shared mesh data for M3ph_iw pattern (DLR2D), initialized once
struct MeshData {
  int64_t n_mesh_points;
  int64_t max_n;
  int64_t n_un1;
  std::vector<std::array<matsubara_freq, 2>> target_mf_2d; // for M (rank-2)
  std::vector<matsubara_freq> target_mf_n1;                // for GM (rank-1)

  MeshData() {
    dlr2d_imfreq mesh{beta, dlr_wmax, dlr_eps, PH};
    n_mesh_points = mesh.size();
    max_n         = mesh.max_n();

    // Collect unique n1 from DLR2D mesh
    std::set<long> unique_n1_set;
    for (auto mp : mesh) {
      auto [n1, n2] = mp.index();
      unique_n1_set.insert(n1);
    }
    std::vector<long> unique_n1(unique_n1_set.begin(), unique_n1_set.end());
    n_un1 = static_cast<int64_t>(unique_n1.size());

    // 2D target matsubara_freq for M (rank-2)
    target_mf_2d.reserve(n_mesh_points);
    for (long d = 0; d < n_mesh_points; ++d) {
      auto [n1, n2] = mesh.to_index(d);
      target_mf_2d.push_back({matsubara_freq(n2, beta, Fermion), matsubara_freq(n1, beta, Fermion)});
    }

    // 1D target matsubara_freq for GM (rank-1)
    target_mf_n1.reserve(n_un1);
    for (int64_t k = 0; k < n_un1; ++k) target_mf_n1.push_back(matsubara_freq(unique_n1[k], beta, Fermion));
  }
};

// Shared mesh data for M_iw pattern (pure DLR mesh), initialized once
struct DLRMeshData {
  int64_t n_dlr_pts;
  std::vector<matsubara_freq> target_mf;

  DLRMeshData() {
    dlr_imfreq mesh{beta, Fermion, dlr_wmax, dlr_eps};
    n_dlr_pts = mesh.size();

    target_mf.reserve(n_dlr_pts);
    for (auto w : mesh) target_mf.push_back(w);
  }
};

static MeshData &get_mesh_data() {
  static MeshData data;
  return data;
}

static DLRMeshData &get_dlr_mesh_data() {
  static DLRMeshData data;
  return data;
}

// Benchmark args: {k} where k is the perturbation order

// --- Rank-1 benchmarks (GM/MG pattern) ---
// n_points = k: after accumulating over columns/rows, each GM/MG buffer gets one push per row/column

static void BM_Nfft_Rank1_Type1(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  int64_t grid_size = 2 * (md.max_n + 1);
  nda::array<dcomplex, 1> output(grid_size);
  output = 0;

  buffer_t<1> buf{output, buf_size, beta, tol};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank1_Type3(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(md.n_un1);
  output = 0;

  buffer_t<1> buf{output, md.target_mf_n1, buf_size, tol, type_t::type3};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

// --- Rank-2 benchmarks (M pattern) ---
// n_points = k² / bl_size² because k² determinant entries are distributed across bl_size² buffers

static void BM_Nfft_Rank2_Type1(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = std::max<int64_t>(k * k / (bl_size * bl_size), 1);

  int64_t grid_size = 2 * (md.max_n + 1);
  nda::array<dcomplex, 2> output(grid_size, grid_size);
  output = 0;

  buffer_t<2> buf{output, buf_size, beta, tol};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = {tau_dist(rng), tau_dist(rng)};
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank2_Type3(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = std::max<int64_t>(k * k / (bl_size * bl_size), 1);

  nda::array<dcomplex, 1> output(md.n_mesh_points);
  output = 0;

  buffer_t<2> buf{output, md.target_mf_2d, buf_size, tol, type_t::type3};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = {tau_dist(rng), tau_dist(rng)};
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

// --- Direct DFT benchmarks (tolerance-independent) ---

static void BM_Nfft_Rank1_DirectType1(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(md.n_un1);
  output = 0;

  buffer_t<1> buf{output, md.target_mf_n1, buf_size, tol, type_t::direct_type1};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank1_DirectBitwise(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(md.n_un1);
  output = 0;

  buffer_t<1> buf{output, md.target_mf_n1, buf_size, tol, type_t::direct_bitwise};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank2_DirectType1(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = std::max<int64_t>(k * k / (bl_size * bl_size), 1);

  nda::array<dcomplex, 1> output(md.n_mesh_points);
  output = 0;

  buffer_t<2> buf{output, md.target_mf_2d, buf_size, tol, type_t::direct_type1};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = {tau_dist(rng), tau_dist(rng)};
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank2_DirectPrime(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = std::max<int64_t>(k * k / (bl_size * bl_size), 1);

  nda::array<dcomplex, 1> output(md.n_mesh_points);
  output = 0;

  buffer_t<2> buf{output, md.target_mf_2d, buf_size, tol, type_t::direct_prime};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = {tau_dist(rng), tau_dist(rng)};
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

// --- M_iw pattern benchmarks (pure DLR mesh, rank-1) ---
// Uses dlr_imfreq mesh instead of DLR2D unique_n1

static void BM_Nfft_M_iw_Type3(benchmark::State &state) {
  auto &dlr_md     = get_dlr_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(dlr_md.n_dlr_pts);
  output = 0;

  buffer_t<1> buf{output, dlr_md.target_mf, buf_size, tol, type_t::type3};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_M_iw_DirectType1(benchmark::State &state) {
  auto &dlr_md     = get_dlr_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(dlr_md.n_dlr_pts);
  output = 0;

  buffer_t<1> buf{output, dlr_md.target_mf, buf_size, tol, type_t::direct_type1};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_M_iw_DirectBitwise(benchmark::State &state) {
  auto &dlr_md     = get_dlr_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(dlr_md.n_dlr_pts);
  output = 0;

  buffer_t<1> buf{output, dlr_md.target_mf, buf_size, tol, type_t::direct_bitwise};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

// --- NAF direct DFT benchmarks ---

static void BM_Nfft_Rank1_DirectType3(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(md.n_un1);
  output = 0;

  buffer_t<1> buf{output, md.target_mf_n1, buf_size, tol, type_t::direct_type3};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank2_DirectType3(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = std::max<int64_t>(k * k / (bl_size * bl_size), 1);

  nda::array<dcomplex, 1> output(md.n_mesh_points);
  output = 0;

  buffer_t<2> buf{output, md.target_mf_2d, buf_size, tol, type_t::direct_type3};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = {tau_dist(rng), tau_dist(rng)};
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_M_iw_DirectType3(benchmark::State &state) {
  auto &dlr_md     = get_dlr_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(dlr_md.n_dlr_pts);
  output = 0;

  buffer_t<1> buf{output, dlr_md.target_mf, buf_size, tol, type_t::direct_type3};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

// --- Automatic dispatch benchmarks ---

static void BM_Nfft_Rank1_Automatic(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(md.n_un1);
  output = 0;

  buffer_t<1> buf{output, md.target_mf_n1, buf_size, tol, type_t::automatic};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_Rank2_Automatic(benchmark::State &state) {
  auto &md         = get_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = std::max<int64_t>(k * k / (bl_size * bl_size), 1);

  nda::array<dcomplex, 1> output(md.n_mesh_points);
  output = 0;

  buffer_t<2> buf{output, md.target_mf_2d, buf_size, tol, type_t::automatic};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = {tau_dist(rng), tau_dist(rng)};
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

static void BM_Nfft_M_iw_Automatic(benchmark::State &state) {
  auto &dlr_md     = get_dlr_mesh_data();
  int64_t k        = state.range(0);
  int64_t n_points = k;

  nda::array<dcomplex, 1> output(dlr_md.n_dlr_pts);
  output = 0;

  buffer_t<1> buf{output, dlr_md.target_mf, buf_size, tol, type_t::automatic};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  std::vector<double> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    taus[i] = tau_dist(rng);
    vals[i] = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back({taus[i]}, vals[i]);
    buf.flush();
  }
}

// --- M4_iw pattern benchmarks (rectangular uniform grid, rank-2, Type1) ---
// Grid shape: (6*n_iw_M4, 2*n_iw_M4) - rectangular, not square
// Push pattern: buf.push_back({tau_j, beta - tau_i}, -Ginv_ji)

static void BM_Nfft_M4_iw_Type1(benchmark::State &state) {
  int64_t k        = state.range(0);
  int64_t n_iw_M4  = state.range(1);
  int64_t n_points = std::max<int64_t>(static_cast<int64_t>(k) * k / (static_cast<int64_t>(bl_size) * bl_size), 1);

  // Rectangular grid: (6*n_iw_M4, 2*n_iw_M4)
  int64_t grid_size_0 = 6 * n_iw_M4;
  int64_t grid_size_1 = 2 * n_iw_M4;
  nda::array<dcomplex, 2> output(grid_size_0, grid_size_1);
  output = 0;

  buffer_t<2> buf{output, buf_size, beta, tol};

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  // M4_iw uses push_back({tau_j, beta - tau_i}, -Ginv_ji)
  std::vector<std::array<double, 2>> taus(n_points);
  std::vector<dcomplex> vals(n_points);
  for (int64_t i = 0; i < n_points; ++i) {
    double tau_j = tau_dist(rng);
    double tau_i = tau_dist(rng);
    taus[i]      = {tau_j, beta - tau_i};
    vals[i]      = dcomplex(val_dist(rng), val_dist(rng));
  }

  for (auto _ : state) {
    output = 0;
    for (int64_t i = 0; i < n_points; ++i) buf.push_back(taus[i], vals[i]);
    buf.flush();
  }
}

// clang-format off
// Register benchmarks with perturbation order k
BENCHMARK(BM_Nfft_Rank1_Type1)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank1_Type3)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank1_DirectType1)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank1_DirectBitwise)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank1_DirectType3)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank2_Type1)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank2_Type3)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank2_DirectType1)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank2_DirectPrime)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank2_DirectType3)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);

// M_iw pattern (pure DLR mesh)
BENCHMARK(BM_Nfft_M_iw_Type3)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_M_iw_DirectType1)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_M_iw_DirectBitwise)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_M_iw_DirectType3)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);

// Automatic dispatch
BENCHMARK(BM_Nfft_Rank1_Automatic)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_Rank2_Automatic)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_Nfft_M_iw_Automatic)->Arg(16)->Arg(64)->Arg(256)->Arg(1024);

// M4_iw pattern (rectangular uniform grid) - Args: {k, n_iw_M4}
BENCHMARK(BM_Nfft_M4_iw_Type1)
    ->Args({16, 16})->Args({64, 16})->Args({256, 16})->Args({1024, 16})
    ->Args({16, 32})->Args({64, 32})->Args({256, 32})->Args({1024, 32})
    ->Args({16, 64})->Args({64, 64})->Args({256, 64})->Args({1024, 64});
// clang-format on
