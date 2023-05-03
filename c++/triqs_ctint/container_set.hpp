#pragma once
#include "./types.hpp"
#include <optional>

namespace triqs_ctint {

  /// The set of all (optional) measurement containers in solver_core
  struct container_set {

    //============ Containers for measurements

    /// Average sign of the CTINT
    mc_weight_t average_sign;

    /// Average perturbation order
    double average_k;

    /// Auto-correlation time
    double auto_corr_time;

    /// Average perturbation order distribution
    std::optional<std::vector<double>> histogram;

    /// The density matrix (measured by operator insertion)
    std::optional<block_matrix_t> density;

    /// Building block for the Green function in imaginary time (Eq. (23) in Notes)
    std::optional<block_gf<imtime, M_tau_target_t>> M_tau;

    /// Hartree-term of M_tau
    std::optional<block_matrix_t> M_hartree;

    /// Same as M_tau, but measured directly in Matsubara frequencies using NFFT
    std::optional<g_iw_t> M_iw_nfft;

    /// Building block for the full vertex function measured directly in Matsubara frequencies using NFFT
    std::optional<chi4_iw_t> M4_iw;

    /// Building block for the full vertex function (pp channel) measured directly in Matsubara frequencies using NFFT
    std::optional<chi4_iw_t> M4pp_iw;

    /// Building block for the full vertex function (ph channel) measured directly in Matsubara frequencies using NFFT
    std::optional<chi4_iw_t> M4ph_iw;

    /// Building block for the fermion boson vertex (pp channel) in Matsubara frequencies using NFFT
    std::optional<chi3_iw_t> M3pp_iw_nfft;

    /// Building block for the fermion boson vertex (ph channel) in Matsubara frequencies using NFFT
    std::optional<chi3_iw_t> M3ph_iw_nfft;

    /// Building block for the fermion boson vertex (pp channel) in imaginary time
    std::optional<chi3_tau_t> M3pp_tau;

    /// Building block for the fermion boson vertex (ph channel) in imaginary time
    std::optional<chi3_tau_t> M3ph_tau;

    /// Building block for the fermion boson vertex (xph channel) in imaginary time
    std::optional<chi3_tau_t> M3xph_tau;

    /// Equal-time peak in M3pp_tau
    std::optional<chi2_tau_t> M3pp_delta;

    /// Equal-time peak in M3ph_tau
    std::optional<chi2_tau_t> M3ph_delta;

    /// Equal-time peak in M3ph_tau
    std::optional<chi2_tau_t> M3xph_delta;

    /// The equal time correlator $\chi_2$ in the particle-particle channel in imaginary times as obtained by operator insertion
    std::optional<chi2_tau_t> chi2pp_tau;

    /// The equal time correlator $\chi_2$ in the particle-hole channel in imaginary times as obtained by operator insertion
    std::optional<chi2_tau_t> chi2ph_tau;

    /// The correlation function $\chi_AB$ in imaginary times
    std::optional<gf<imtime>> chiAB_tau;

    //============ Containers dependent on measured quantities

    /// The Fourier-transform of M_tau. Dependent on M_tau
    std::optional<g_iw_t> M_iw;

    /// Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_iw
    g_iw_t G_iw;

    /// Self-energy in Matsubara frequencies. Dependent on M_iw
    g_iw_t Sigma_iw;

    /// Building block for the fermion boson vertex (pp channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3pp_iw;

    /// Building block for the fermion boson vertex (ph channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3ph_iw;

    /// Building block for the fermion boson vertex (xph channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3xph_iw;

    /// The two-particle vertex function in purely fermionic notation (iw1, iw2, iw3)
    std::optional<chi4_iw_t> F_iw;

    /// The two-particle vertex function (pp channel)
    std::optional<chi4_iw_t> Fpp_iw;

    /// The two-particle vertex function (ph channel)
    std::optional<chi4_iw_t> Fph_iw;

    /// The connected part of the two-particle Green function
    std::optional<chi4_iw_t> G2c_iw;

    /// The connected part of the two-particle Green function (pp channel)
    std::optional<chi4_iw_t> G2ppc_iw;

    /// The connected part of the two-particle Green function (ph channel)
    std::optional<chi4_iw_t> G2phc_iw;

    /// The two-particle Green function
    std::optional<chi4_iw_t> G2_iw;

    /// The two-particle Green function (pp channel)
    std::optional<chi4_iw_t> G2pp_iw;

    /// The two-particle Green function (ph channel)
    std::optional<chi4_iw_t> G2ph_iw;

    /// The equal time correlator $\chi_2$ in the particle-particle channel in Matsubara frequencies
    std::optional<chi2_iw_t> chi2pp_iw;

    /// The equal time correlator $\chi_2$ in the particle-hole channel in Matsubara frequencies
    std::optional<chi2_iw_t> chi2ph_iw;

    /// M2 in the particle-particle channel in imaginary time as obtained from M3
    std::optional<chi2_tau_t> chi2pp_conn_tau_from_M3;

    /// M2 in the particle-hole channel in imaginary time as obtained from M3
    std::optional<chi2_tau_t> chi2ph_conn_tau_from_M3;

    /// M2 in the particle-hole-cross channel in imaginary time as obtained from M3
    std::optional<chi2_tau_t> chi2xph_conn_tau_from_M3;

    /// The equal time correlator $\chi_2$ in the particle-particle channel in imaginary times as obtained from M3pp_tau
    std::optional<chi2_tau_t> chi2pp_tau_from_M3;

    /// The equal time correlator $\chi_2$ in the particle-hole channel in imaginary times as obtained from M3ph_tau
    std::optional<chi2_tau_t> chi2ph_tau_from_M3;

    /// The equal time correlator $\chi_2$ in the particle-hole-cross channel in imaginary times as obtained from M3ph_tau
    std::optional<chi2_tau_t> chi2xph_tau_from_M3;

    /// The equal time correlator $\chi_2$ in the particle-particle channel in imaginary frequencies as obtained from M3pp_tau
    std::optional<chi2_iw_t> chi2pp_iw_from_M3;

    /// The equal time correlator $\chi_2$ in the particle-hole channel in imaginary frequencies as obtained from M3ph_tau
    std::optional<chi2_iw_t> chi2ph_iw_from_M3;

    /// The equal time correlator $\chi_2$ in the particle-hole-cross channel in imaginary frequencies as obtained from M3ph_tau
    std::optional<chi2_iw_t> chi2xph_iw_from_M3;

    /// The correlation function $\chi_AB$ in imaginary frequencies
    std::optional<gf<imfreq>> chiAB_iw;

    /// The equal time correlator $\chi_3$ in the particle-particle channel in Matsubara frequencies
    std::optional<chi3_iw_t> chi3pp_iw;

    /// The equal time correlator $\chi_3$ in the particle-hole channel in Matsubara frequencies
    std::optional<chi3_iw_t> chi3ph_iw;

    /// The equal time correlator $\chi_3$ in the particle-hole-cross channel in Matsubara frequencies
    std::optional<chi3_iw_t> chi3xph_iw;

    /// The equal time correlator $\chi_3$ in the particle-particle channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement
    std::optional<chi3_iw_t> chi3pp_iw_nfft;

    /// The equal time correlator $\chi_3$ in the particle-hole channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement
    std::optional<chi3_iw_t> chi3ph_iw_nfft;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write(grp, "average_sign", c.average_sign);
      h5_write(grp, "average_k", c.average_k);
      h5_write(grp, "auto_corr_time", c.auto_corr_time);
      h5_write(grp, "histogram", c.histogram);
      h5_write(grp, "density", c.density);
      h5_write(grp, "M_tau", c.M_tau);
      h5_write(grp, "M_hartree", c.M_hartree);
      h5_write(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_write(grp, "M4_iw", c.M4_iw);
      h5_write(grp, "M4pp_iw", c.M4pp_iw);
      h5_write(grp, "M4ph_iw", c.M4ph_iw);
      h5_write(grp, "M3pp_tau", c.M3pp_tau);
      h5_write(grp, "M3ph_tau", c.M3ph_tau);
      h5_write(grp, "M3xph_tau", c.M3xph_tau);
      h5_write(grp, "M3pp_delta", c.M3pp_delta);
      h5_write(grp, "M3ph_delta", c.M3ph_delta);
      h5_write(grp, "M3xph_delta", c.M3xph_delta);
      h5_write(grp, "M3pp_iw_nfft", c.M3pp_iw_nfft);
      h5_write(grp, "M3ph_iw_nfft", c.M3ph_iw_nfft);
      h5_write(grp, "chi2pp_tau", c.chi2pp_tau);
      h5_write(grp, "chi2ph_tau", c.chi2ph_tau);
      h5_write(grp, "chiAB_tau", c.chiAB_tau);
      h5_write(grp, "M_iw", c.M_iw);
      h5_write(grp, "G_iw", c.G_iw);
      h5_write(grp, "Sigma_iw", c.Sigma_iw);
      h5_write(grp, "M3pp_iw", c.M3pp_iw);
      h5_write(grp, "M3ph_iw", c.M3ph_iw);
      h5_write(grp, "M3xph_iw", c.M3xph_iw);
      h5_write(grp, "F_iw", c.F_iw);
      h5_write(grp, "Fpp_iw", c.Fpp_iw);
      h5_write(grp, "Fph_iw", c.Fph_iw);
      h5_write(grp, "G2_iw", c.G2_iw);
      h5_write(grp, "G2pp_iw", c.G2pp_iw);
      h5_write(grp, "G2ph_iw", c.G2ph_iw);
      h5_write(grp, "G2c_iw", c.G2c_iw);
      h5_write(grp, "G2ppc_iw", c.G2ppc_iw);
      h5_write(grp, "G2phc_iw", c.G2phc_iw);
      h5_write(grp, "chi2pp_iw", c.chi2pp_iw);
      h5_write(grp, "chi2ph_iw", c.chi2ph_iw);
      h5_write(grp, "chi2pp_conn_tau_from_M3", c.chi2pp_conn_tau_from_M3);
      h5_write(grp, "chi2ph_conn_tau_from_M3", c.chi2ph_conn_tau_from_M3);
      h5_write(grp, "chi2xph_conn_tau_from_M3", c.chi2xph_conn_tau_from_M3);
      h5_write(grp, "chi2pp_tau_from_M3", c.chi2pp_tau_from_M3);
      h5_write(grp, "chi2ph_tau_from_M3", c.chi2ph_tau_from_M3);
      h5_write(grp, "chi2xph_tau_from_M3", c.chi2xph_tau_from_M3);
      h5_write(grp, "chi2pp_iw_from_M3", c.chi2pp_iw_from_M3);
      h5_write(grp, "chi2ph_iw_from_M3", c.chi2ph_iw_from_M3);
      h5_write(grp, "chi2xph_iw_from_M3", c.chi2xph_iw_from_M3);
      h5_write(grp, "chiAB_iw", c.chiAB_iw);
      h5_write(grp, "chi3pp_iw", c.chi3pp_iw);
      h5_write(grp, "chi3ph_iw", c.chi3ph_iw);
      h5_write(grp, "chi3xph_iw", c.chi3xph_iw);
      h5_write(grp, "chi3pp_iw_nfft", c.chi3pp_iw_nfft);
      h5_write(grp, "chi3ph_iw_nfft", c.chi3ph_iw_nfft);
    }

    /// Function that reads all containers from hdf5 file
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set &c) {
      auto grp = h5group.open_group(subgroup_name);
      h5_read(grp, "average_sign", c.average_sign);
      h5::try_read(grp, "average_k", c.average_k);
      h5::try_read(grp, "auto_corr_time", c.auto_corr_time);
      h5_read(grp, "histogram", c.histogram);
      h5::try_read(grp, "density", c.density);
      h5_read(grp, "M_tau", c.M_tau);
      h5_read(grp, "M_hartree", c.M_hartree);
      h5_read(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_read(grp, "M4_iw", c.M4_iw);
      h5_read(grp, "M4pp_iw", c.M4pp_iw);
      h5_read(grp, "M4ph_iw", c.M4ph_iw);
      h5_read(grp, "M3pp_tau", c.M3pp_tau);
      h5_read(grp, "M3ph_tau", c.M3ph_tau);
      h5::try_read(grp, "M3xph_tau", c.M3xph_tau);
      h5_read(grp, "M3pp_delta", c.M3pp_delta);
      h5_read(grp, "M3ph_delta", c.M3ph_delta);
      h5::try_read(grp, "M3xph_delta", c.M3xph_delta);
      h5_read(grp, "M3pp_iw_nfft", c.M3pp_iw_nfft);
      h5_read(grp, "M3ph_iw_nfft", c.M3ph_iw_nfft);
      h5_read(grp, "chi2pp_tau", c.chi2pp_tau);
      h5_read(grp, "chi2ph_tau", c.chi2ph_tau);
      h5_read(grp, "chiAB_tau", c.chiAB_tau);
      h5_read(grp, "M_iw", c.M_iw);
      h5_read(grp, "G_iw", c.G_iw);
      h5_read(grp, "Sigma_iw", c.Sigma_iw);
      h5_read(grp, "M3pp_iw", c.M3pp_iw);
      h5_read(grp, "M3ph_iw", c.M3ph_iw);
      h5::try_read(grp, "M3xph_iw", c.M3xph_iw);
      h5_read(grp, "F_iw", c.F_iw);
      h5_read(grp, "Fpp_iw", c.Fpp_iw);
      h5_read(grp, "Fph_iw", c.Fph_iw);
      h5_read(grp, "G2_iw", c.G2_iw);
      h5_read(grp, "G2pp_iw", c.G2pp_iw);
      h5_read(grp, "G2ph_iw", c.G2ph_iw);
      h5_read(grp, "G2c_iw", c.G2c_iw);
      h5_read(grp, "G2ppc_iw", c.G2ppc_iw);
      h5_read(grp, "G2phc_iw", c.G2phc_iw);
      h5_read(grp, "chi2pp_iw", c.chi2pp_iw);
      h5_read(grp, "chi2ph_iw", c.chi2ph_iw);
      h5_read(grp, "chi2pp_conn_tau_from_M3", c.chi2pp_conn_tau_from_M3);
      h5_read(grp, "chi2ph_conn_tau_from_M3", c.chi2ph_conn_tau_from_M3);
      h5::try_read(grp, "chi2xph_conn_tau_from_M3", c.chi2xph_conn_tau_from_M3);
      h5_read(grp, "chi2pp_tau_from_M3", c.chi2pp_tau_from_M3);
      h5_read(grp, "chi2ph_tau_from_M3", c.chi2ph_tau_from_M3);
      h5::try_read(grp, "chi2xph_tau_from_M3", c.chi2xph_tau_from_M3);
      h5_read(grp, "chi2pp_iw_from_M3", c.chi2pp_iw_from_M3);
      h5_read(grp, "chi2ph_iw_from_M3", c.chi2ph_iw_from_M3);
      h5::try_read(grp, "chi2xph_iw_from_M3", c.chi2xph_iw_from_M3);
      // For backward compatibility we keep these additional reads
      if (!c.chi2pp_conn_tau_from_M3) h5_read(grp, "M2pp_tau", c.chi2pp_conn_tau_from_M3);
      if (!c.chi2ph_conn_tau_from_M3) h5_read(grp, "M2ph_tau", c.chi2ph_conn_tau_from_M3);
      if (!c.chi2pp_tau_from_M3) h5_read(grp, "chi2pp_new_tau", c.chi2pp_tau_from_M3);
      if (!c.chi2ph_tau_from_M3) h5_read(grp, "chi2ph_new_tau", c.chi2ph_tau_from_M3);
      if (!c.chi2pp_iw_from_M3) h5_read(grp, "chi2pp_new_iw", c.chi2pp_iw_from_M3);
      if (!c.chi2ph_iw_from_M3) h5_read(grp, "chi2ph_new_iw", c.chi2ph_iw_from_M3);
      h5_read(grp, "chiAB_iw", c.chiAB_iw);
      h5_read(grp, "chi3pp_iw", c.chi3pp_iw);
      h5_read(grp, "chi3ph_iw", c.chi3ph_iw);
      h5_read(grp, "chi3xph_iw", c.chi3xph_iw);
      h5_read(grp, "chi3pp_iw_nfft", c.chi3pp_iw_nfft);
      h5_read(grp, "chi3ph_iw_nfft", c.chi3ph_iw_nfft);
    }
  };

} // namespace triqs_ctint
