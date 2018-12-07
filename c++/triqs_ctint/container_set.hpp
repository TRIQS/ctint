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

    /// Average perturbation order distribution
    std::optional<std::vector<double>> histogram;

    /// The density matrix (measured by operator insertion)
    std::optional<std::vector<matrix<dcomplex>>> density;

    /// Building block for the Green function in imaginary time (Eq. (23) in Notes)
    std::optional<block_gf<imtime, M_tau_target_t>> M_tau;

    /// Same as M_tau, but measured directly in Matsubara frequencies using NFFT
    std::optional<g_iw_t> M_iw_nfft;

    /// Same as M4_tau, but measured directly in Matsubara frequencies using NFFT
    std::optional<chi4_iw_t> M4_iw;

    /// Building block for the fermion boson vertex (pp channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3pp_iw_nfft;

    /// Building block for the fermion boson vertex (ph channel) in Matsubara frequencies
    std::optional<chi3_iw_t> M3ph_iw_nfft;

    /// Building block for the fermion boson vertex (pp channel) in imaginary time
    std::optional<chi3_tau_t> M3pp_tau;

    /// Building block for the fermion boson vertex (ph channel) in imaginary time
    std::optional<chi3_tau_t> M3ph_tau;

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

    /// Building block for the susceptibility (pp channel) in Matsubara frequencies
    std::optional<chi2_iw_t> M2pp_iw;

    /// Building block for the susceptibility (ph channel) in Matsubara frequencies
    std::optional<chi2_iw_t> M2ph_iw;

    /// The two-particle vertex function in purely fermionic notation (iw1, iw2, iw3)
    std::optional<chi4_iw_t> F_iw;

    /// The connected part of the two-particle Green function
    std::optional<chi4_iw_t> G2c_iw;

    /// The two-particle Green function
    std::optional<chi4_iw_t> G2_iw;

    /// The equal time correlator $\chi_2$ in the particle-particle channel in Matsubara frequencies
    std::optional<chi2_iw_t> chi2pp_iw;

    /// The equal time correlator $\chi_2$ in the particle-hole channel in Matsubara frequencies
    std::optional<chi2_iw_t> chi2ph_iw;

    /// The correlation function $\chi_AB$ in imaginary frequencies
    std::optional<gf<imfreq>> chiAB_iw;

    /// The equal time correlator $\chi_3$ in the particle-particle channel in Matsubara frequencies
    std::optional<chi3_iw_t> chi3pp_iw;

    /// The equal time correlator $\chi_3$ in the particle-hole channel in Matsubara frequencies
    std::optional<chi3_iw_t> chi3ph_iw;

    /// The equal time correlator $\chi_3$ in the particle-particle channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement
    std::optional<chi3_iw_t> chi3pp_iw_nfft;

    /// The equal time correlator $\chi_3$ in the particle-hole channel in Matsubara frequencies as obtained by the NFFT $M_3$ measurement
    std::optional<chi3_iw_t> chi3ph_iw_nfft;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set const &c) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write(grp, "average_sign", c.average_sign);
      h5_write(grp, "average_k", c.average_k);
      h5_write(grp, "histogram", c.histogram);
      h5_write(grp, "density", c.density);
      h5_write(grp, "M_tau", c.M_tau);
      h5_write(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_write(grp, "M4_iw", c.M4_iw);
      h5_write(grp, "M3pp_tau", c.M3pp_tau);
      h5_write(grp, "M3ph_tau", c.M3ph_tau);
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
      h5_write(grp, "M2pp_iw", c.M2pp_iw);
      h5_write(grp, "M2ph_iw", c.M2ph_iw);
      h5_write(grp, "F_iw", c.F_iw);
      h5_write(grp, "G2_iw", c.G2_iw);
      h5_write(grp, "G2c_iw", c.G2c_iw);
      h5_write(grp, "chi2pp_iw", c.chi2pp_iw);
      h5_write(grp, "chi2ph_iw", c.chi2ph_iw);
      h5_write(grp, "chiAB_iw", c.chiAB_iw);
      h5_write(grp, "chi3pp_iw", c.chi3pp_iw);
      h5_write(grp, "chi3ph_iw", c.chi3ph_iw);
      h5_write(grp, "chi3pp_iw_nfft", c.chi3pp_iw_nfft);
      h5_write(grp, "chi3ph_iw_nfft", c.chi3ph_iw_nfft);
    }

    /// Function that reads all containers from hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set &c) {
      auto grp = h5group.open_group(subgroup_name);
      // Useful to keep backward compatibility for older solvers
      auto h5_try_read = [](triqs::h5::group grp, std::string key_name, auto &obj){
	if(grp.has_key(key_name)) h5_read(grp, key_name, obj);
	else std::cout << "WARNING: Could not find key " << key_name << " during h5_read!\n";
      };
      h5_read(grp, "average_sign", c.average_sign);
      h5_try_read(grp, "average_k", c.average_k);
      h5_read(grp, "histogram", c.histogram);
      h5_read(grp, "density", c.density);
      h5_read(grp, "M_tau", c.M_tau);
      h5_read(grp, "M_iw_nfft", c.M_iw_nfft);
      h5_read(grp, "M4_iw", c.M4_iw);
      h5_read(grp, "M3pp_tau", c.M3pp_tau);
      h5_read(grp, "M3ph_tau", c.M3ph_tau);
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
      h5_read(grp, "M2pp_iw", c.M2pp_iw);
      h5_read(grp, "M2ph_iw", c.M2ph_iw);
      h5_read(grp, "F_iw", c.F_iw);
      h5_read(grp, "G2_iw", c.G2_iw);
      h5_read(grp, "G2c_iw", c.G2c_iw);
      h5_read(grp, "chi2pp_iw", c.chi2pp_iw);
      h5_read(grp, "chi2ph_iw", c.chi2ph_iw);
      h5_read(grp, "chiAB_iw", c.chiAB_iw);
      h5_read(grp, "chi3pp_iw", c.chi3pp_iw);
      h5_read(grp, "chi3ph_iw", c.chi3ph_iw);
      h5_read(grp, "chi3pp_iw_nfft", c.chi3pp_iw_nfft);
      h5_read(grp, "chi3ph_iw_nfft", c.chi3ph_iw_nfft);
    }
  };

} // namespace triqs_ctint
