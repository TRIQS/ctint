#pragma once
#include <triqs/utility/optional_compat.hpp>

namespace triqs_ctint {

  struct container_set {

    //============ Containers for measurements

    /// Average sign of the CTINT
    double average_sign = 0.0;

    /// Building block for the Green function in imaginary time (Eq. (23) in Notes)
    std::optional<block_gf<imtime, matrix_valued>> M_tau;

    /// Same as M_tau, but measured directly in frequencies
    std::optional<block_gf<imfreq, matrix_valued>> M_iw_nfft;

    /// The improved estimator F_tau
    std::optional<block_gf<imtime, matrix_valued>> F_tau;

    //============ Containers dependent on measured quantities

    /// The Fourier-transform of M_tau. Dependent on M_tau
    std::optional<block_gf<imfreq, matrix_valued>> M_iw;

    /// Greens function in Matsubara frequencies (Eq. (18) in Notes). Dependent on M_tau
    std::optional<block_gf<imfreq, matrix_valued>> Giw;

    /// Self-energy in Matsubara frequencies. Dependent on M_tau
    std::optional<block_gf<imfreq, matrix_valued>> Sigma_iw;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set const &c) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
      h5_write(grp, "average_sign", c.average_sign);
      h5_write(grp, "M_tau", c.M_tau);
      h5_write(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_write(grp, "F_tau", c.F_tau);
      h5_write(grp, "M_iw", c.M_iw);
      h5_write(grp, "Giw", c.Giw);
      h5_write(grp, "Sigma_iw", c.Sigma_iw);
    }

    /// Function that read all containers to hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set &c) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      h5_read(grp, "average_sign", c.average_sign);
      h5_read(grp, "M_tau", c.M_tau);
      h5_read(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_read(grp, "F_tau", c.F_tau);
      h5_read(grp, "M_iw", c.M_iw);
      h5_read(grp, "Giw", c.Giw);
      h5_read(grp, "Sigma_iw", c.Sigma_iw);
    }
  };

} // namespace triqs_ctint
