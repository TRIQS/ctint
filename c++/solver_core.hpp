#pragma once
#include "./params.hpp"
#include "./qmc_config.hpp"
#include "./container_set.hpp"

namespace triqs_ctint {

  /// The Solver class
  class solver_core : public container_set {

    public:
    /// Noninteracting Green Function in Matsubara frequencies
    block_gf<imfreq, matrix_valued> G0_iw;

    /// Dynamic density-density interaction in Matsubara frequencies
    std::optional<block_gf<imfreq, matrix_valued>> D0_iw;

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

    private:
    // The shifted noninteracting Green Function in imaginary time
    block_gf<imtime, matrix_valued> G0_shift_tau;

    // Struct containing the parameters relevant for construction
    const constr_params_t constr_params;

    // Mpi Communicator
    triqs::mpi::communicator world;

    // Calculate G0_shift_tau given G0_iw
    void prepare_G0_shift_tau(params_t const &params);

    // Perform post-processing
    void post_process(qmc_config_t const &qmc_config, container_set *results);

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    /// Function that writes the solver_core to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
      h5_write(grp, "", s.result_set());
      h5_write(grp, "G0_iw", s.G0_iw);
      h5_write(grp, "D0_iw", s.D0_iw);
      h5_write(grp, "Jperp_iw", s.Jperp_iw);
    }

    /// Function that read all containers to hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, solver_core &s) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      h5_read(grp, "", s.result_set());
      h5_read(grp, "G0_iw", s.G0_iw);
      h5_read(grp, "D0_iw", s.D0_iw);
      h5_read(grp, "Jperp_iw", s.Jperp_iw);
    }
  };

} // namespace triqs_ctint
