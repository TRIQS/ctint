// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/utility/macros.hpp>
#include "./params.hpp"
#include "./qmc_config.hpp"
#include "./container_set.hpp"

namespace triqs_ctint {

  /// The CT-INT solver.
  class solver_core : public container_set {

    public:
    /// Non-interacting Green's function \f$ G_0(i\omega) \f$ in Matsubara frequencies.
    g_iw_t G0_iw;

    /// Inverse of the non-interacting Green's function \f$ G_0^{-1}(i\omega) \f$.
    g_iw_t G0_iw_inv;

    /// Dynamic density-density interaction \f$ D_0(i\omega) \f$ in Matsubara frequencies (DLR mesh).
    std::optional<block2_gf<mesh::dlr_imfreq, matrix_valued>> D0_iw;

    /// Dynamic spin-spin interaction \f$ J_\perp(i\omega) \f$ in Matsubara frequencies (DLR mesh).
    std::optional<gf<mesh::dlr_imfreq, matrix_valued>> Jperp_iw;

    /**
     * @brief Construct a CT-INT solver.
     *
     * @param constr_params_ Set of parameters used to construct the solver.
     */
    solver_core(constr_params_t const &constr_params_);

    // Delete assignement operator because of const members
    solver_core(solver_core const &p)            = default;
    solver_core(solver_core &&p)                 = default;
    ~solver_core()                               = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p)      = default;

    /**
     * Solve the impurity problem with a CT-INT calculation.
     *
     * @param solve_params Set of parameters used for the solve.
     */
    void solve(solve_params_t const &solve_params);

    /// The shifted non-interacting Green's function in Matsubara frequencies.
    g_iw_t G0_shift_iw;

    /// The shifted non-interacting Green's function in imaginary time.
    g_tau_t G0_shift_tau;

    /// Calculate the shifted non-interacting Green's function given \f$ G_0(i\omega) \f$.
    C2PY_IGNORE void prepare_G0_shift_iw(params_t const &params);

    /// Parameters used to construct the solver.
    constr_params_t constr_params;

    /// Parameters used in the last solve (empty until the solver has been run).
    std::optional<solve_params_t> last_solve_params;

    private:
    // Mpi Communicator
    mpi::communicator world;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    // Function to perform the post-processing steps
    void post_process(params_t const &p);

    public:
    /// Retrigger post-processing with the last set of parameters.
    void post_process() {
      if (not last_solve_params) TRIQS_RUNTIME_ERROR << "You need to run the solver once before you post-process";
      post_process({constr_params, last_solve_params.value()});
    }

    static std::string hdf5_format() {
#if defined(INTERACTION_IS_COMPLEX)
      return "CTINT_SolverCore_complex_all";
#elif defined(GTAU_IS_COMPLEX)
      return "CTINT_SolverCore_complex_gtau";
#else
      return "CTINT_SolverCore";
#endif
    }

    // Function that writes the solver_core to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
      auto grp = h5group.create_group(subgroup_name);
      write_hdf5_format(grp, s);
      h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(STRINGIZE(TRIQS_GIT_HASH)));
      h5_write_attribute(grp, "CTINT_GIT_HASH", std::string(STRINGIZE(CTINT_GIT_HASH)));
      h5_write(grp, "", s.result_set());
      h5_write(grp, "constr_params", s.constr_params);
      h5_write(grp, "last_solve_params", s.last_solve_params);
      h5_write(grp, "G0_iw", s.G0_iw);
      h5_write(grp, "G0_iw_inv", s.G0_iw_inv);
      h5_write(grp, "G0_shift_iw", s.G0_shift_iw);
      h5_write(grp, "G0_shift_tau", s.G0_shift_tau);
      h5_write(grp, "D0_iw", s.D0_iw);
      h5_write(grp, "Jperp_iw", s.Jperp_iw);
    }

    // Function that read all containers to hdf5 file
    C2PY_IGNORE
    static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name) {
      auto grp           = h5group.open_group(subgroup_name);
      auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
      auto s             = solver_core{constr_params};
      h5_read(grp, "", s.result_set());
      h5_read(grp, "last_solve_params", s.last_solve_params);
      h5_read(grp, "G0_iw", s.G0_iw);
      h5::try_read(grp, "G0_iw_inv", s.G0_iw_inv);
      h5_read(grp, "G0_shift_iw", s.G0_shift_iw);
      h5_read(grp, "G0_shift_tau", s.G0_shift_tau);
      h5_read(grp, "D0_iw", s.D0_iw);
      h5_read(grp, "Jperp_iw", s.Jperp_iw);
      return s;
    }
  };
} // namespace triqs_ctint
