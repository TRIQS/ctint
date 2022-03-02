(changelog)=

# Changelog

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
