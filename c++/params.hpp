#pragma once

#include "./types.hpp"

namespace triqs_ctint {

  struct constr_params_t {

    /// Number of tau points for gf<imtime, matrix_valued>
    int n_tau = 10000;

    /// Number of Matsubara frequencies for gf<imfreq, matrix_valued>
    int n_iw = 500;

    /// Inverse temperature
    double beta;

    ///block structure of the gf
    gf_struct_t gf_struct;

    /// Switch for dynamic density-density interaction
    bool use_D = false;

    /// Switch for dynamic spin-spin interaction
    bool use_Jperp = false;

    /// Number of tau pts for D0_tau and jperp_tau
    int n_tau_dynamical_interactions = 10001;

    /// Number of matsubara freqs for D0_iw and jperp_iw
    int n_iw_dynamical_interactions = 200;

    /// Number of block indeces for the Green function
    int n_blocks() const { return gf_struct.size(); }

    /// Names of block indeces for the Green function
    auto block_names() const {
      std::vector<std::string> v;
      for (auto const &bl : gf_struct) v.push_back(bl.first);
      return v;
    }

    /// Write constr_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp);
  };

  struct solve_params_t {

    // ----------- System Specific -----------

    /// Shift of the chemical potential mu_sigma --> mu_sigma + hartree_shift_sigma
    std::vector<double> hartree_shift = std::vector<double>{};

    /// Interaction Hamiltonian
    many_body_operator h_int;

    // ----------- QMC Specific -----------

    /// Switch for the use of the alpha function. Compare Sec. 1.3 in Notes.
    bool use_alpha = false;

    /// Number of auxiliary spins
    int n_s = 2;

    /// Alpha parameter
    alpha_t alpha;

    /// Number of MC cycles
    int n_cycles;

    /// Length of a MC cycles
    int length_cycle = 50;

    /// Number of warmup cycles
    int n_warmup_cycles = 5000;

    /// Random seed of the random generator
    int random_seed = 34788 + 928374 * triqs::mpi::communicator().rank();

    /// Name of the random generator
    std::string random_name = "";

    /// Use double insertion
    bool use_double_insertion = false;

    /// Maximum running time in seconds (-1 : no limit)
    int max_time = -1;

    /// Verbosity
    int verbosity = triqs::mpi::communicator().rank() == 0 ? 3 : 0;

    // ----------- Measurements -----------

    /// Measure the MC sign
    bool measure_average_sign = true;

    /// Measure M(tau)
    bool measure_M_tau = false;

    /// Measure M(iomega) using nfft
    bool measure_M_iw = false;

    /// Measure F(tau)
    bool measure_F_tau = false;

    /// Measure M4(iw) NFFT
    bool measure_M4_iw = false;
    int n_iw_M4        = 32;

    /// Measure M3(iw)
    bool measure_M3pp_iw  = false;
    bool measure_M3ph_iw  = false;
    int n_iw_M3           = 64;

    /// Measure M3(iw)
    bool measure_M3pp_tau  = false;
    bool measure_M3ph_tau  = false;
    int n_tau_M3           = 1000;

    /// Measure M2(tau)
    bool measure_M2pp_tau  = false;
    bool measure_M2ph_tau  = false;
    bool measure_M2xph_tau = false;
    int n_tau_M2           = 10000;
    int n_iw_M2            = 128;

    /// Size of the Nfft buffer
    int nfft_buf_size = 500;

    /// Perform post processing
    bool post_process = true;

    /// Write constr_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read constr_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp);
  };

  struct params_t : constr_params_t, solve_params_t {
    params_t(constr_params_t constr_params_, solve_params_t solve_params_) : constr_params_t(constr_params_), solve_params_t(solve_params_) {}

    /// Return matrix for static density-density interaction
    array<array<dcomplex, 2>, 2> get_U() const;
  };

} // namespace triqs_ctint
