// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "./types.hpp"

namespace triqs_ctint {

  /// Parameters used for constructing the solver class.
  struct constr_params_t {

    /// Number of imaginary-time points for the single-particle quantities.
    int n_tau = 5001;

    /// DLR bandwidth cutoff \f$ w_{max} = \Lambda / \beta \f$ for the single-particle quantities.
    double dlr_wmax;

    /// DLR error tolerance \f$ \epsilon \f$ for the single-particle quantities.
    double dlr_eps = 1e-10;

    /// Inverse temperature \f$ \beta \f$.
    double beta;

    /// Structure of the Green's function (names and sizes of blocks).
    gf_struct_t gf_struct;

    /// Use a dynamic density-density interaction?
    bool use_D = false;

    /// Use a dynamic spin-spin interaction?
    bool use_Jperp = false;

    /// Number of imaginary-time points for \f$ D_0(\tau) \f$ and \f$ J_\perp(\tau) \f$.
    int n_tau_dynamical_interactions = this->n_tau;

    /// Number of blocks of the Green's function.
    int n_blocks() const { return gf_struct.size(); }

    /// Names of the blocks of the Green's function.
    auto block_names() const {
      std::vector<std::string> v;
      for (auto const &bl : gf_struct) v.push_back(bl.first);
      return v;
    }

    /// Write constr_params_t to hdf5.
    friend void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from hdf5.
    friend void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &cp);
  };

  /// Parameters passed to the solve method of the solver class.
  struct solve_params_t {

    // ----------- System Specific -----------

    /// Interacting part of the local Hamiltonian.
    many_body_operator h_int;

    // ----------- QMC Specific -----------

    /// Number of auxiliary spins.
    int n_s = 2;

    /// The \f$ \alpha \f$ tensor used in the determinantal expansion.
    alpha_t alpha;

    /// Number of QMC cycles.
    int n_cycles;

    /// Length of a single QMC cycle.
    int length_cycle = 100;

    /// Number of cycles for thermalization.
    int n_warmup_cycles = 5000;

    /// Seed for the random number generator (shared by all MPI ranks; the rank is used as the
    /// spawn key to derive an independent stream per rank, see triqs::mc_tools::random_generator).
    int random_seed = 34788;

    /// Name of the random number generator.
    std::string random_name = "";

    /// Use double insertion?
    bool use_double_insertion = true;

    /// Types of insertions to use.
    std::vector<int> insertion_types = {};

    /// Use auxiliary spin-flip insertion (requires \f$ n_s = 2 \f$)?
    bool use_auxiliary_spin_flip = false;

    /// Maximum runtime in seconds, use -1 to set infinite.
    int max_time = -1;

    /// Maximum perturbation order accepted during insertion and removal moves (use -1 for unlimited).
    int max_order = -1;

    /// Verbosity level.
    int verbosity = mpi::communicator().rank() == 0 ? 3 : 0;

    /// Catch exceptions on the nodes and rethrow them on rank 0?
    bool rethrow_exception = true;

    // ----------- Measurements -----------

    /// Measure the sign only?
    bool measure_sign_only = false;

    /// Measure the Monte-Carlo sign?
    bool measure_average_sign = true;

    /// Measure the average perturbation order?
    bool measure_average_k = true;

    /// Measure the auto-correlation time?
    bool measure_auto_corr_time = true;

    /// Measure the perturbation-order distribution?
    bool measure_histogram = false;

    /// Measure the density by operator insertion?
    bool measure_density = true;

    /// Measure \f$ M(\tau) \f$?
    bool measure_M_tau = true;

    /// Measure \f$ M(i\omega) \f$ using NFFT?
    bool measure_M_iw = false;

    /// Measure \f$ M^{(4)}(i\omega) \f$ using NFFT?
    bool measure_M4_iw = false;
    /// Measure \f$ M^{(4)}_{pp}(i\omega) \f$ using NFFT?
    bool measure_M4pp_iw = false;
    /// Measure \f$ M^{(4)}_{ph}(i\omega) \f$ using NFFT?
    bool measure_M4ph_iw = false;
    /// Number of positive bosonic Matsubara frequencies in \f$ M^{(4)} \f$.
    int n_iW_M4 = 32;
    /// Number of positive fermionic Matsubara frequencies in \f$ M^{(4)} \f$.
    int n_iw_M4 = 32;

    /// Measure \f$ M^{(3)}_{pp}(i\omega) \f$?
    bool measure_M3pp_iw = false;
    /// Measure \f$ M^{(3)}_{ph}(i\omega) \f$?
    bool measure_M3ph_iw = false;
    /// Number of positive fermionic Matsubara frequencies in \f$ M^{(3)} \f$.
    int n_iw_M3 = 64;
    /// Number of positive bosonic Matsubara frequencies in \f$ M^{(3)} \f$.
    int n_iW_M3 = 32;
    /// Measure \f$ M^{(3)}_{pp}(\tau) \f$?
    bool measure_M3pp_tau = false;
    /// Measure \f$ M^{(3)}_{ph}(\tau) \f$?
    bool measure_M3ph_tau = false;
    /// Measure \f$ M^{(3)}_{xph}(\tau) \f$?
    bool measure_M3xph_tau = false;
    /// Number of imaginary-time points in \f$ M^{(3)} \f$.
    int n_tau_M3 = 201;

    /// Measure \f$ \chi^{(2)}_{pp}(\tau) \f$ by insertion?
    bool measure_chi2pp_tau = false;
    /// Measure \f$ \chi^{(2)}_{ph}(\tau) \f$ by insertion?
    bool measure_chi2ph_tau = false;
    /// Number of imaginary-time points in \f$ \chi^{(2)} \f$.
    int n_tau_chi2 = 201;
    /// Number of positive Matsubara frequencies in \f$ \chi^{(2)} \f$.
    int n_iw_chi2 = 32;

    /// Measure \f$ \chi_{AB}(\tau) \f$ by insertion?
    bool measure_chiAB_tau = false;
    /// List of all operators \f$ A \f$.
    std::vector<many_body_operator> chi_A_vec = {};
    /// List of all operators \f$ B \f$.
    std::vector<many_body_operator> chi_B_vec = {};

    /// Size of the NFFT buffer.
    int nfft_buf_size = 100000;

    /// Tolerance for the NFFT transform.
    double nfft_tol = 1e-8;

    /// Perform post-processing?
    bool post_process = true;

    /// The maximum size of the determinant matrix before a resize.
    int det_init_size = 1000;

    /// Maximum number of operations before testing the accuracy of \f$ \det(M) \f$ and \f$ M^{-1} \f$.
    int det_n_operations_before_check = 100;

    /// Threshold for determinant precision warnings.
    double det_precision_warning = 1.e-8;

    /// Threshold for determinant precision errors.
    double det_precision_error = 1.e-5;

    /// Bound for the determinant matrix being singular (if \f$ < 0 \f$, checks for subnormal numbers instead).
    double det_singular_threshold = -1;

    /// Write solve_params_t to hdf5.
    friend void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read solve_params_t from hdf5.
    friend void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &sp);
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_params_t, solve_params_t {
    params_t() = default;
    params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
       : constr_params_t(constr_params_), solve_params_t(solve_params_) {}
  };

} // namespace triqs_ctint
