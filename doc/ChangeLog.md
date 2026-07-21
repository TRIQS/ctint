# Changelog

## Version 4.0.1

CTINT version 4.0.1 is a patch release that pins the build dependencies for
the 4.0.x release branch and loosens a numerical test tolerance so the
dynamic-interaction tests pass with the reference BLAS used in the conda-forge
build environment.

We thank all contributors: Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Pin c2py (0.9.x), triqs_hartree_fock (4.0.x) and poet (v0.0.0) for the 4.0.x release
* Loosen the h5diff tolerance of the dynamic-interaction tests (all, all_continuousBoson, densdens, jperp) from 5e-5 to 1e-4


## Version 4.0.0

CTINT version 4.0.0 is a compatibility release for TRIQS version 4.0.0 that
* ports the Python bindings from cpp2py to the new c2py + clair framework
* migrates the single-particle and bosonic quantities to the DLR Matsubara mesh
* introduces a new per-term alpha tensor with support for dynamic interactions
* replaces nfft with finufft for the non-uniform Fourier transforms
* uses the latest [app4triqs](https://github.com/TRIQS/app4triqs) skeleton
* migrates the Jenkins CI to the new k8s-based system
* fixes a number of application issues

We provide a more detailed description of the changes below.

### Bindings (c2py + clair)
* Port the Python bindings from cpp2py to c2py + clair
* Generate the `solver_core` and `post_process` Python modules via `c2py_add_module`
* Disable C++20 module scanning for clair-c2py compatibility
* Update the doxygen documentation in the public headers (`params.hpp`, `container_set.hpp`, `post_process.hpp`, `solver_core.hpp`) for the c2py docstring pipeline

### DLR
* Migrate the single-particle and bosonic quantities to the DLR Matsubara mesh
* Adjust the ctint solver to use the DLR mesh version of the Hartree-Fock solver
* Use the mesh-based `triqs_hartree_fock` API for the alpha determination

### Alpha tensor
* Implement a new per-term alpha tensor and start the search from the last alpha
* Restore dynamic D0 interactions for the new per-term alpha tensor
* Allow `n_s = 2` for non-density-density terms
* Add an auxiliary spin-flip move

### NFFT / finufft
* Replace nfft with finufft (v2.5.1) for the non-uniform Fourier transforms; set `nfft_tol=1e-14`
* Extend `nfft_buf_t` with type-3 and direct DFT modes and refactor the direct NUDFT kernel with rank-specialized algorithms
* Build finufft with a portable architecture by default (no `-march=native`)
* Relax the default finufft type-3 tolerance to 1e-13

### Measurements
* Migrate `auto_corr_time` to the new TRIQS `log_binning` API and add density observables
* Add error-bar estimation for `average_sign` and `average_k`
* Enable `measure_density` by default

### General
* Port ctint to TRIQS 4.0 (`triqs.gf` -> `triqs.gfs`)
* Use `parity_sort`/`insertion_sort` for insert and remove moves
* Add `rethrow_exception` to `solve_params_t` and backwards compatibility for insertion parameters
* Install `triqs_hartree_fock` as part of the ctint installation
* Convert the implementation notes to reStructuredText and update them for the new alpha tensor

### cmake / deps
* Bump the `triqs_hartree_fock` dependency version to 4.0
* Increase the h5diff test tolerance to 5e-5 for BLAS portability

### jenkins / ghactions
* Migrate the Jenkins CI to the new k8s-based system (#26)
* Synchronize the GitHub Actions configuration with app4triqs

We thank all contributors: Marco Barbone, Jennifer Coulter, Thomas Hahn, Alexander Hampel, Marcel Klett, Henri Menke, Dylan Simon, Nils Wentzell


## Version 3.1.0

CTINT version 3.1.0 is a compatibility
release for TRIQS version 3.1.0 that
* uses the latest [app4triqs/3.1.x](https://github.com/TRIQS/app4triqs) skeleton
* switches documentation to read-the-docs theme
* improves ghactions and jenkins configuration
* fixes a number of application issues

We provide a more detailed description of the changes below.

### General
* Port ctint to triqs 3.1
* Fix runaway memory in G2c_from_M4
* Expose det_manip parameters to solver
* Fix compiler warnings
* Add triqs_ctint::measures::auto_corr_time and enable by default
* Change measure_density default from true to false
* Adjust ctint to changes in gf_struct object
* Add missing h5 write/read for measure_density parameter
* Remove rerun.py test as identical to solver_h5io
* Be sure to ignore build directories in .dockerignore
* Fix type of n_tau_M3_del from double to int in post_processing
* Fix move constructor of nfft_buf_t to properly call rebind

### doc
* Fix sphinx recursion error and two build warnings
* Update implementation notes to latest version

### cmake
* Avoid policy warning for CMP0115
* Add option to compile with MSAN

### gh-actions
* brew install nfft from wentzell/triqs/nfft

We thank all contributors: Michel Ferrero, Alexander Hampel, Henri Menke, Dylan Simon, Nils Wentzell


## Version 3.0.1

Consistently put bosonic frequency first for objects with mixed bosonic and fermionic frequencies.


## Version 3.0.0

ctint version 3.0.0 is a compatibility
release for TRIQS version 3.0.0 that
* introduces compatibility with Python 3 (Python 2 no longer supported)
* adds a cmake-based dependency management
* fixes several application issues

## Version 2.2.0

ctint version 2.2.0 is a compatibility
release for TRIQS version 2.2.0

This is the initial release for this project.
