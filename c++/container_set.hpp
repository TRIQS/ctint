#pragma once
#include <optional>

namespace triqs_ctint {

  struct container_set {

    //============ Containers for measurements

    /// Average sign of the CTINT
    double average_sign = 0.0;

    /// Building block for the Green function in imaginary time (Eq. (23) in Notes)
    std::optional<g_tau_t> M_tau;

    /// Same as M_tau, but measured directly in Matsubara frequencies using NFFT
    std::optional<g_iw_t> M_iw_nfft;

    /// The improved estimator F_tau
    std::optional<g_tau_t> F_tau;

    /// Same as M4_tau, but measured directly in Matsubara frequencies using NFFT
    std::optional<chi4_iw_t> M4_iw;

    /// Building block for the fermion boson vertex (pp channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3pp_iw;

    /// Building block for the fermion boson vertex (ph channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3ph_iw;

    /// Building block for the susceptibility (pp channel) in imaginary time
    std::optional<chi2_tau_t> M2pp_tau;

    /// Building block for the susceptibility (ph channel) in imaginary time
    std::optional<chi2_tau_t> M2ph_tau;

    /// Building block for the susceptibility (xph channel) in imaginary time
    std::optional<chi2_tau_t> M2xph_tau;

    //============ Containers dependent on measured quantities

    /// The Fourier-transform of M_tau. Dependent on M_tau
    std::optional<g_iw_t> M_iw;

    /// Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_iw
    std::optional<g_iw_t> G_iw;

    /// Self-energy in Matsubara frequencies. Dependent on M_iw
    std::optional<g_iw_t> Sigma_iw;

    /// Building block for the susceptibility (pp channel) in Matsubara frequencies
    std::optional<chi2_iw_t> M2pp_iw;

    /// Building block for the susceptibility (ph channel) in Matsubara frequencies
    std::optional<chi2_iw_t> M2ph_iw;

    /// Building block for the susceptibility (xph channel) in Matsubara frequencies
    std::optional<chi2_iw_t> M2xph_iw;

    /// Building block for the susceptibility (xph channel) in Matsubara frequencies
    std::optional<chi4_iw_t> F_iw;

    /// Building block for the fermion boson vertex (pp channel) in Matsubara frequencies
    std::optional<chi3_iw_t> chi3pp_iw;

    /// Building block for the fermion boson vertex (ph channel) in Matsubara frequencies
    std::optional<chi3_iw_t> chi3ph_iw;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set const &c) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
      h5_write(grp, "average_sign", c.average_sign);
      h5_write(grp, "M_tau", c.M_tau);
      h5_write(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_write(grp, "F_tau", c.F_tau);
      h5_write(grp, "M4_iw", c.M4_iw);
      h5_write(grp, "M3pp_iw", c.M3pp_iw);
      h5_write(grp, "M3ph_iw", c.M3ph_iw);
      h5_write(grp, "M2pp_tau", c.M2pp_tau);
      h5_write(grp, "M2ph_tau", c.M2ph_tau);
      h5_write(grp, "M2xph_tau", c.M2xph_tau);
      h5_write(grp, "M_iw", c.M_iw);
      h5_write(grp, "G_iw", c.G_iw);
      h5_write(grp, "Sigma_iw", c.Sigma_iw);
      h5_write(grp, "M2pp_iw", c.M2pp_iw);
      h5_write(grp, "M2ph_iw", c.M2ph_iw);
      h5_write(grp, "M2xph_iw", c.M2xph_iw);
      h5_write(grp, "F_iw", c.F_iw);
      h5_write(grp, "chi3pp_iw", c.chi3pp_iw);
      h5_write(grp, "chi3ph_iw", c.chi3ph_iw);
    }

    /// Function that read all containers to hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set &c) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      h5_read(grp, "average_sign", c.average_sign);
      h5_read(grp, "M_tau", c.M_tau);
      h5_read(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_read(grp, "F_tau", c.F_tau);
      h5_read(grp, "M4_iw", c.M4_iw);
      h5_read(grp, "M3pp_iw", c.M3pp_iw);
      h5_read(grp, "M3ph_iw", c.M3ph_iw);
      h5_read(grp, "M2pp_tau", c.M2pp_tau);
      h5_read(grp, "M2ph_tau", c.M2ph_tau);
      h5_read(grp, "M2xph_tau", c.M2xph_tau);
      h5_read(grp, "M_iw", c.M_iw);
      h5_read(grp, "G_iw", c.G_iw);
      h5_read(grp, "Sigma_iw", c.Sigma_iw);
      h5_read(grp, "M2pp_iw", c.M2pp_iw);
      h5_read(grp, "M2ph_iw", c.M2ph_iw);
      h5_read(grp, "M2xph_iw", c.M2xph_iw);
      h5_read(grp, "F_iw", c.F_iw);
      h5_read(grp, "chi3pp_iw", c.chi3pp_iw);
      h5_read(grp, "chi3ph_iw", c.chi3ph_iw);
    }
  };

} // namespace triqs_ctint
