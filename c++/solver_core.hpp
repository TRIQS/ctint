#pragma once
#include "./params.hpp"
#include "./qmc_config.hpp"
#include "./container_set.hpp"

namespace triqs_ctint {

  /// The Solver class
  class solver_core : public container_set {

    public:
    /// Noninteracting Green Function in Matsubara frequencies
    g_iw_t G0_iw;

    /// Dynamic density-density interaction in Matsubara frequencies
    std::optional<g_iw_t> D0_iw;

    /// Dynamic spin-spin interaction in Matsubara frequencies
    std::optional<gf<imfreq, matrix_valued>> Jperp_iw;

    /**
     * Construct a CTINT solver
     *
     * @param construct_parameters Set of parameters specific to the CTINT solver
     */
    CPP2PY_ARG_AS_DICT
    solver_core(constr_params_t const &constr_params_);

    // Delete assignement operator because of const members
    solver_core(solver_core const &p) = default;
    solver_core(solver_core &&p)      = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p) = default;

    /**
     * Solve method that performs CTINT calculation
     *
     * @param solve_params_t Set of parameters specific to the CTINT run
     */
    CPP2PY_ARG_AS_DICT
    void solve(solve_params_t const &solve_params);

    void solve();

    private:
    // The shifted noninteracting Green Function in Matsubara frequencies
    g_iw_t G0_shift_iw;

    // The shifted noninteracting Green Function in imaginary time
    g_tau_t G0_shift_tau;

    // Struct containing the parameters relevant for the solver construction
    constr_params_t constr_params;

    // Struct containing the parameters relevant for the solve process
    solve_params_t solve_params;

    // Mpi Communicator
    triqs::mpi::communicator world;

    // Calculate G0_shift_tau given G0_iw
    void prepare_G0_shift_tau(params_t const &params);

    // Function to perform the post-processing steps
    void post_process(params_t const &p);

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    // Function that writes the solver_core to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
      h5_write(grp, "", s.result_set());
      h5_write(grp, "GIT_SHA1", std::string(STRINGIZE(GIT_SHA1)));
      h5_write(grp, "constr_params", s.constr_params);
      h5_write(grp, "solve_params", s.solve_params);
      h5_write(grp, "G0_iw", s.G0_iw);
      h5_write(grp, "D0_iw", s.D0_iw);
      h5_write(grp, "Jperp_iw", s.Jperp_iw);
    }

    // Function that read all containers to hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, solver_core &s) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      h5_read(grp, "", s.result_set());
      h5_read(grp, "constr_params", s.constr_params);
      h5_read(grp, "solve_params", s.solve_params);
      h5_read(grp, "G0_iw", s.G0_iw);
      h5_read(grp, "D0_iw", s.D0_iw);
      h5_read(grp, "Jperp_iw", s.Jperp_iw);
    }
  };

} // namespace triqs_ctint
