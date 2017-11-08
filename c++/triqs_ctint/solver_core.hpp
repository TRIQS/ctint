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

    /// The inverse of the noninteracting Green Function
    g_iw_t G0_iw_inv;

    /// Dynamic density-density interaction in Matsubara frequencies
    std::optional<block2_gf<imfreq, matrix_valued>> D0_iw;

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
    ~solver_core()                    = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p) = default;

    /**
     * Solve method that performs CTINT calculation
     *
     * @param solve_params_t Set of parameters specific to the CTINT run
     */
    CPP2PY_ARG_AS_DICT
    void solve(solve_params_t const &solve_params);

    /// The shifted noninteracting Green Function in Matsubara frequencies
    g_iw_t G0_shift_iw;

    /// The shifted noninteracting Green Function in imaginary time
    g_tau_t G0_shift_tau;

    // Calculate G0_shift_tau given G0_iw
    CPP2PY_ARG_AS_DICT
    void prepare_G0_shift_iw(params_t const &params);

    // Struct containing the parameters relevant for the solver construction
    constr_params_t constr_params;

    // Struct containing the parameters relevant for the solve process
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
    // Allow the user to retrigger post-processing with the last set of parameters
    void post_process() {
      if (not last_solve_params) TRIQS_RUNTIME_ERROR << "You need to run the solver once before you post-process";
      post_process({constr_params, last_solve_params.value()});
    }

    static std::string hdf5_format() { return "CTINT_SolverCore"; }

    // Function that writes the solver_core to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write_attribute(grp, "Format", solver_core::hdf5_format());
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
    CPP2PY_IGNORE
    static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name) {
      auto grp           = h5group.open_group(subgroup_name);
      auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
      auto s             = solver_core{constr_params};
      h5_read(grp, "", s.result_set());
      h5_read(grp, "last_solve_params", s.last_solve_params);
      h5_read(grp, "G0_iw", s.G0_iw);
      h5_try_read(grp, "G0_iw_inv", s.G0_iw_inv);
      h5_read(grp, "G0_shift_iw", s.G0_shift_iw);
      h5_read(grp, "G0_shift_tau", s.G0_shift_tau);
      h5_read(grp, "D0_iw", s.D0_iw);
      h5_read(grp, "Jperp_iw", s.Jperp_iw);
      return s;
    }
  };
} // namespace triqs_ctint
