#pragma once

#include "./types.hpp"

namespace triqs_ctint {

  /// The parameters for the solver construction
  struct constr_params_t {

    /// Number of tau points for gf<imtime, matrix_valued>
    int n_tau = 5001;

    /// Number of Matsubara frequencies for gf<imfreq, matrix_valued>
    int n_iw = 500;

    /// Inverse temperature
    double beta;

    /// block structure of the gf
    gf_struct_t gf_struct;

    /// Switch for dynamic density-density interaction
    bool use_D = false;

    /// Switch for dynamic spin-spin interaction
    bool use_Jperp = false;

    /// Number of tau pts for D0_tau and jperp_tau
    int n_tau_dynamical_interactions = this->n_tau;

    /// Number of matsubara freqs for D0_iw and jperp_iw
    int n_iw_dynamical_interactions = this->n_iw;

    /// Number of block indeces for the Green function
    int n_blocks() const { return gf_struct.size(); }

    /// Names of block indeces for the Green function
    auto block_names() const {
      std::vector<std::string> v;
      for (auto const &bl : gf_struct) v.push_back(bl.first);
      return v;
    }

    /// Write constr_params_t to hdf5
    friend void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from hdf5
    friend void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &cp);
  };

  /// The parameters for the solve function
  struct solve_params_t {

    // ----------- System Specific -----------

    /// Interaction Hamiltonian
    many_body_operator h_int;

    // ----------- QMC Specific -----------

    /// Number of auxiliary spins
    int n_s = 1;

    /// Alpha tensor
    alpha_t alpha;

    /// Number of MC cycles
    int n_cycles;

    /// Length of a MC cycles
    int length_cycle = 50;

    /// Number of warmup cycles
    int n_warmup_cycles = 5000;

    /// Random seed of the random generator
    int random_seed = 34788 + 928374 * mpi::communicator().rank();

    /// Name of the random generator
    std::string random_name = "";

    /// Use double insertion
    bool use_double_insertion = false;

    /// Maximum running time in seconds (-1 : no limit)
    int max_time = -1;

    /// Maximum pertubation order which is accepted in move::insert/remove
    /// (<0 : unlimited)
    int max_order = -1;

    /// Verbosity
    int verbosity = mpi::communicator().rank() == 0 ? 3 : 0;

    // ----------- Measurements -----------

    /// Measure Sign only mode
    bool measure_sign_only = false;

    /// Measure the MC sign
    bool measure_average_sign = true;

    /// Measure the average perturbation order
    bool measure_average_k = true;

    /// Measure the auto-correlation time
    bool measure_auto_corr_time = true;

    /// Measure the average perturbation order distribution
    bool measure_histogram = false;

    /// Measure the density matrix by operator insertion
    bool measure_density = false;

    /// Measure M(tau)
    bool measure_M_tau = true;

    /// Measure M(iomega) using nfft
    bool measure_M_iw = false;

    /// Measure M4(iw) NFFT
    bool measure_M4_iw = false;
    /// Measure M4pp(iw) NFFT
    bool measure_M4pp_iw = false;
    /// Measure M4ph(iw) NFFT
    bool measure_M4ph_iw = false;
    /// Number of positive bosonic Matsubara frequencies in M4
    int n_iW_M4 = 32;
    /// Number of positive fermionic Matsubara frequencies in M4
    int n_iw_M4 = 32;

    /// Measure M3pp(iw)
    bool measure_M3pp_iw = false;
    /// Measure M3ph(iw)
    bool measure_M3ph_iw = false;
    /// Number of positive fermionic Matsubara frequencies in M3
    int n_iw_M3 = 64;
    /// Number of positive bosonic Matsubara frequencies in M3
    int n_iW_M3 = 32;

    /// Measure M3pp(tau)
    bool measure_M3pp_tau = false;
    /// Measure M3ph(tau)
    bool measure_M3ph_tau = false;
    /// Measure M3xph(tau)
    bool measure_M3xph_tau = false;
    /// Number of imaginary time points in M3
    int n_tau_M3 = 201;

    /// Measure of chi2pp by insertion
    bool measure_chi2pp_tau = false;
    /// Measure of chi2ph by insertion
    bool measure_chi2ph_tau = false;
    /// Number of imaginary time points in chi2
    int n_tau_chi2 = 201;
    /// Number of positive Matsubara frequencies in chi2
    int n_iw_chi2 = 32;

    /// Measure of chiAB by insertion
    bool measure_chiAB_tau = false;
    /// The list of all operators A
    std::vector<many_body_operator> chi_A_vec = {};
    /// The list of all operators B
    std::vector<many_body_operator> chi_B_vec = {};

    /// Size of the Nfft buffer
    int nfft_buf_size = 500;

    /// Perform post processing
    bool post_process = true;

    /// The maximum size of the determinant matrix before a resize
    int det_init_size = 1000;

    /// Max number of ops before the test of deviation of the det, M^-1 is performed.
    int det_n_operations_before_check = 100;

    /// Threshold for determinant precision warnings
    double det_precision_warning = 1.e-8;

    /// Threshold for determinant precision error
    double det_precision_error = 1.e-5;

    /// Bound for the determinant matrix being singular: abs(det) < singular_threshold.
    /// For negative threshold check if !isnormal(abs(det)).
    double det_singular_threshold = -1;

    /// Write constr_params_t to hdf5
    friend void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read constr_params_t from hdf5
    friend void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &sp);
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_params_t, solve_params_t {
    params_t() = default;
    params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
       : constr_params_t(constr_params_), solve_params_t(solve_params_) {}
  };

} // namespace triqs_ctint
