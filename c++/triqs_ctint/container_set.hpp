// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "./types.hpp"
#include <optional>

namespace triqs_ctint {

  /// Container for all (optional) quantities measured and post-processed by the solver.
  struct container_set {

    //============ Containers for measurements

    /// Average Monte-Carlo sign.
    mc_weight_t average_sign;

    /// Total number of measurements.
    uint64_t nmeasures;

    /// Average perturbation order.
    double average_k;

    /// Error bar on the average sign.
    std::optional<double> average_sign_error;

    /// Error bar on the average perturbation order.
    std::optional<double> average_k_error;

    /// Auto-correlation time.
    double auto_corr_time;

    /// Warmup time in seconds.
    double warmup_time;

    /// Accumulation time in seconds.
    double accumulation_time;

    /// Perturbation-order distribution.
    std::optional<std::vector<double>> histogram;

    /// The density (measured by operator insertion).
    std::optional<block_matrix_t> density;

    /// Building block for the Green's function in imaginary time \f$ M(\tau) \f$.
    std::optional<block_gf<imtime, M_tau_target_t>> M_tau;

    /// Hartree term of \f$ M(\tau) \f$.
    std::optional<block_matrix_t> M_hartree;

    /// Same as \f$ M(\tau) \f$, but measured directly in Matsubara frequencies using NFFT on the DLR grid.
    std::optional<g_dlr_iw_t> M_iw_nfft;

    /// Building block \f$ M^{(4)}(i\omega) \f$ for the full vertex, measured in Matsubara frequencies using NFFT.
    std::optional<chi4_iw_t> M4_iw;

    /// Building block \f$ M^{(4)}_{pp}(i\omega) \f$ (particle-particle channel) for the full vertex, measured using NFFT.
    std::optional<chi4_iw_t> M4pp_iw;

    /// Building block \f$ M^{(4)}_{ph}(i\omega) \f$ (particle-hole channel) for the full vertex, measured using NFFT.
    std::optional<chi4_iw_t> M4ph_iw;

    /// Building block \f$ M^{(3)}_{pp}(i\omega) \f$ (particle-particle channel) for the fermion-boson vertex, measured using NFFT.
    std::optional<chi3_iw_t> M3pp_iw_nfft;

    /// Building block \f$ M^{(3)}_{ph}(i\omega) \f$ (particle-hole channel) for the fermion-boson vertex, measured using NFFT.
    std::optional<chi3_iw_t> M3ph_iw_nfft;

    /// Building block \f$ M^{(3)}_{pp}(\tau) \f$ (particle-particle channel) for the fermion-boson vertex in imaginary time.
    std::optional<chi3_tau_t> M3pp_tau;

    /// Building block \f$ M^{(3)}_{ph}(\tau) \f$ (particle-hole channel) for the fermion-boson vertex in imaginary time.
    std::optional<chi3_tau_t> M3ph_tau;

    /// Building block \f$ M^{(3)}_{xph}(\tau) \f$ (particle-hole-cross channel) for the fermion-boson vertex in imaginary time.
    std::optional<chi3_tau_t> M3xph_tau;

    /// Equal-time peak of \f$ M^{(3)}_{pp}(\tau) \f$.
    std::optional<chi2_tau_t> M3pp_delta;

    /// Equal-time peak of \f$ M^{(3)}_{ph}(\tau) \f$.
    std::optional<chi2_tau_t> M3ph_delta;

    /// Equal-time peak of \f$ M^{(3)}_{xph}(\tau) \f$.
    std::optional<chi2_tau_t> M3xph_delta;

    /// Correlator \f$ \chi^{(2)}_{pp}(\tau) \f$ (particle-particle channel) in imaginary time, obtained by operator insertion.
    std::optional<chi2_tau_t> chi2pp_tau;

    /// Correlator \f$ \chi^{(2)}_{ph}(\tau) \f$ (particle-hole channel) in imaginary time, obtained by operator insertion.
    std::optional<chi2_tau_t> chi2ph_tau;

    /// Correlation function \f$ \chi_{AB}(\tau) \f$ in imaginary time.
    std::optional<gf<imtime>> chiAB_tau;

    //============ Containers dependent on measured quantities

    /// Fourier transform of \f$ M(\tau) \f$.
    std::optional<g_iw_t> M_iw;

    /// Green's function \f$ G(i\omega) \f$ in Matsubara frequencies.
    g_iw_t G_iw;

    /// Dynamic self-energy \f$ \Sigma_{dyn}(i\omega) \f$ in Matsubara frequencies (DLR, decays to zero).
    g_iw_t Sigma_dyn_iw;

    /// Static (Hartree) part of the self-energy, \f$ \Sigma = \Sigma_{dyn} + \Sigma_{hartree} \f$.
    std::optional<block_matrix_t> Sigma_hartree;

    /// Building block \f$ M^{(3)}_{pp}(i\omega) \f$ (particle-particle channel) for the fermion-boson vertex in Matsubara frequencies.
    std::optional<chi3_iw_t> M3pp_iw;

    /// Building block \f$ M^{(3)}_{ph}(i\omega) \f$ (particle-hole channel) for the fermion-boson vertex in Matsubara frequencies.
    std::optional<chi3_iw_t> M3ph_iw;

    /// Building block \f$ M^{(3)}_{xph}(i\omega) \f$ (particle-hole-cross channel) for the fermion-boson vertex in Matsubara frequencies.
    std::optional<chi3_iw_t> M3xph_iw;

    /// The two-particle vertex function \f$ F \f$ in purely fermionic notation.
    std::optional<chi4_iw_t> F_iw;

    /// The two-particle vertex function \f$ F \f$ (particle-particle channel).
    std::optional<chi4_iw_t> Fpp_iw;

    /// The two-particle vertex function \f$ F \f$ (particle-hole channel).
    std::optional<chi4_iw_t> Fph_iw;

    /// Connected part of the two-particle Green's function \f$ G^{(2)} \f$.
    std::optional<chi4_iw_t> G2_conn_iw;

    /// Connected part of the two-particle Green's function \f$ G^{(2)} \f$ (particle-particle channel).
    std::optional<chi4_iw_t> G2pp_conn_iw;

    /// Connected part of the two-particle Green's function \f$ G^{(2)} \f$ (particle-hole channel).
    std::optional<chi4_iw_t> G2ph_conn_iw;

    /// The two-particle Green's function \f$ G^{(2)} \f$.
    std::optional<chi4_iw_t> G2_iw;

    /// The two-particle Green's function \f$ G^{(2)} \f$ (particle-particle channel).
    std::optional<chi4_iw_t> G2pp_iw;

    /// The two-particle Green's function \f$ G^{(2)} \f$ (particle-hole channel).
    std::optional<chi4_iw_t> G2ph_iw;

    /// Correlator \f$ \chi^{(2)}_{pp}(i\omega) \f$ (particle-particle channel) in Matsubara frequencies.
    std::optional<chi2_iw_t> chi2pp_iw;

    /// Correlator \f$ \chi^{(2)}_{ph}(i\omega) \f$ (particle-hole channel) in Matsubara frequencies.
    std::optional<chi2_iw_t> chi2ph_iw;

    /// Connected \f$ \chi^{(2)} \f$ (particle-particle channel) in imaginary time, obtained from \f$ M^{(3)} \f$.
    std::optional<chi2_tau_t> chi2pp_conn_tau_from_M3;

    /// Connected \f$ \chi^{(2)} \f$ (particle-hole channel) in imaginary time, obtained from \f$ M^{(3)} \f$.
    std::optional<chi2_tau_t> chi2ph_conn_tau_from_M3;

    /// Connected \f$ \chi^{(2)} \f$ (particle-hole-cross channel) in imaginary time, obtained from \f$ M^{(3)} \f$.
    std::optional<chi2_tau_t> chi2xph_conn_tau_from_M3;

    /// Correlator \f$ \chi^{(2)}_{pp}(\tau) \f$ (particle-particle channel), obtained from \f$ M^{(3)}_{pp}(\tau) \f$.
    std::optional<chi2_tau_t> chi2pp_tau_from_M3;

    /// Correlator \f$ \chi^{(2)}_{ph}(\tau) \f$ (particle-hole channel), obtained from \f$ M^{(3)}_{ph}(\tau) \f$.
    std::optional<chi2_tau_t> chi2ph_tau_from_M3;

    /// Correlator \f$ \chi^{(2)}_{xph}(\tau) \f$ (particle-hole-cross channel), obtained from \f$ M^{(3)}_{xph}(\tau) \f$.
    std::optional<chi2_tau_t> chi2xph_tau_from_M3;

    /// Correlator \f$ \chi^{(2)}_{pp}(i\omega) \f$ (particle-particle channel), obtained from \f$ M^{(3)}_{pp}(\tau) \f$.
    std::optional<chi2_iw_t> chi2pp_iw_from_M3;

    /// Correlator \f$ \chi^{(2)}_{ph}(i\omega) \f$ (particle-hole channel), obtained from \f$ M^{(3)}_{ph}(\tau) \f$.
    std::optional<chi2_iw_t> chi2ph_iw_from_M3;

    /// Correlator \f$ \chi^{(2)}_{xph}(i\omega) \f$ (particle-hole-cross channel), obtained from \f$ M^{(3)}_{xph}(\tau) \f$.
    std::optional<chi2_iw_t> chi2xph_iw_from_M3;

    /// Correlation function \f$ \chi_{AB}(i\omega) \f$ in Matsubara frequencies.
    std::optional<gf<imfreq>> chiAB_iw;

    /// Correlator \f$ \chi^{(3)}_{pp}(i\omega) \f$ (particle-particle channel) in Matsubara frequencies.
    std::optional<chi3_iw_t> chi3pp_iw;

    /// Correlator \f$ \chi^{(3)}_{ph}(i\omega) \f$ (particle-hole channel) in Matsubara frequencies.
    std::optional<chi3_iw_t> chi3ph_iw;

    /// Correlator \f$ \chi^{(3)}_{xph}(i\omega) \f$ (particle-hole-cross channel) in Matsubara frequencies.
    std::optional<chi3_iw_t> chi3xph_iw;

    /// Correlator \f$ \chi^{(3)}_{pp}(i\omega) \f$ (particle-particle channel), obtained from the NFFT \f$ M^{(3)} \f$ measurement.
    std::optional<chi3_iw_t> chi3pp_iw_nfft;

    /// Correlator \f$ \chi^{(3)}_{ph}(i\omega) \f$ (particle-hole channel), obtained from the NFFT \f$ M^{(3)} \f$ measurement.
    std::optional<chi3_iw_t> chi3ph_iw_nfft;

    /// Function that writes all containers to hdf5 file.
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write(grp, "average_sign", c.average_sign);
      h5_write(grp, "nmeasures", c.nmeasures);
      h5_write(grp, "average_k", c.average_k);
      h5_write(grp, "average_sign_error", c.average_sign_error);
      h5_write(grp, "average_k_error", c.average_k_error);
      h5_write(grp, "auto_corr_time", c.auto_corr_time);
      h5_write(grp, "warmup_time", c.warmup_time);
      h5_write(grp, "accumulation_time", c.accumulation_time);
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
      h5_write(grp, "Sigma_dyn_iw", c.Sigma_dyn_iw);
      h5_write(grp, "Sigma_hartree", c.Sigma_hartree);
      h5_write(grp, "M3pp_iw", c.M3pp_iw);
      h5_write(grp, "M3ph_iw", c.M3ph_iw);
      h5_write(grp, "M3xph_iw", c.M3xph_iw);
      h5_write(grp, "F_iw", c.F_iw);
      h5_write(grp, "Fpp_iw", c.Fpp_iw);
      h5_write(grp, "Fph_iw", c.Fph_iw);
      h5_write(grp, "G2_iw", c.G2_iw);
      h5_write(grp, "G2pp_iw", c.G2pp_iw);
      h5_write(grp, "G2ph_iw", c.G2ph_iw);
      h5_write(grp, "G2_conn_iw", c.G2_conn_iw);
      h5_write(grp, "G2pp_conn_iw", c.G2pp_conn_iw);
      h5_write(grp, "G2ph_conn_iw", c.G2ph_conn_iw);
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

    /// Function that reads all containers from hdf5 file.
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set &c) {
      auto grp = h5group.open_group(subgroup_name);
      h5_read(grp, "average_sign", c.average_sign);
      h5::try_read(grp, "nmeasures", c.nmeasures);
      h5::try_read(grp, "average_k", c.average_k);
      h5::try_read(grp, "average_sign_error", c.average_sign_error);
      h5::try_read(grp, "average_k_error", c.average_k_error);
      h5::try_read(grp, "auto_corr_time", c.auto_corr_time);
      h5::try_read(grp, "warmup_time", c.warmup_time);
      h5::try_read(grp, "accumulation_time", c.accumulation_time);
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
      h5::try_read(grp, "Sigma_dyn_iw", c.Sigma_dyn_iw);
      h5::try_read(grp, "Sigma_hartree", c.Sigma_hartree);
      h5_read(grp, "M3pp_iw", c.M3pp_iw);
      h5_read(grp, "M3ph_iw", c.M3ph_iw);
      h5::try_read(grp, "M3xph_iw", c.M3xph_iw);
      h5_read(grp, "F_iw", c.F_iw);
      h5_read(grp, "Fpp_iw", c.Fpp_iw);
      h5_read(grp, "Fph_iw", c.Fph_iw);
      h5_read(grp, "G2_iw", c.G2_iw);
      h5_read(grp, "G2pp_iw", c.G2pp_iw);
      h5_read(grp, "G2ph_iw", c.G2ph_iw);
      h5_read(grp, "G2_conn_iw", c.G2_conn_iw);
      h5_read(grp, "G2pp_conn_iw", c.G2pp_conn_iw);
      h5_read(grp, "G2ph_conn_iw", c.G2ph_conn_iw);
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
