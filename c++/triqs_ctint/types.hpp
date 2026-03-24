// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.


#pragma once

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/arrays/block_matrix.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/operators/util/extractors.hpp>

#include <itertools/itertools.hpp>
#include <mpi/mpi.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <variant>

namespace triqs_ctint {

  using namespace std::complex_literals; // Complex Unity 1i
  using namespace triqs;
  using namespace triqs::gfs;
  using namespace nda;
  using namespace triqs::operators;
  using namespace triqs::operators::utils;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace h5;

  using namespace itertools;

  /// The channel type
  enum class Chan_t { PP, PH, XPH };

  /// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  using triqs::gfs::gf_struct_t;

  /// Container type of one-particle Green and Vertex functions in imaginary times
#ifdef GTAU_IS_COMPLEX
  using g_tau_t = block_gf<imtime, matrix_valued>;
#else
  using g_tau_t = block_gf<imtime, matrix_real_valued>;
#endif

  /// A const_view to a g_tau_t
  using g_tau_cv_t = g_tau_t::const_view_type;

  /// A view to a g_tau_t
  using g_tau_v_t = g_tau_t::view_type;

  /// Scalar type of g_tau
  using g_tau_scalar_t = g_tau_t::g_t::scalar_t;

  /// Container type of the alpha function. alpha(vertex_index, aux_spin)
  using alpha_t = array<g_tau_scalar_t, 4>;

  /// Container type of one-particle Green and Vertex functions on DLR Matsubara frequencies
  using g_iw_t    = block_gf<mesh::dlr_imfreq, matrix_valued>;
  using g_dlr_iw_t = g_iw_t; // explicit alias used by M_iw measurement

  /// Container type on DLR Matsubara frequencies (view types)
  using g_dlr_iw_cv_t = g_iw_t::const_view_type;
  using g_dlr_iw_v_t  = g_iw_t::view_type;

  /// Container type of one-particle Green and Vertex functions on regular Matsubara frequencies (for post-processing)
  using g_reg_iw_t    = block_gf<imfreq, matrix_valued>;
  using g_reg_iw_cv_t = g_reg_iw_t::const_view_type;

#if defined(INTERACTION_IS_COMPLEX) && !defined(GTAU_IS_COMPLEX)
  static_assert(false, "INTERACTION_IS_COMPLEX requires GTAU_IS_COMPLEX");
#endif

  /// Scalar type of all interaction vertices
#ifdef INTERACTION_IS_COMPLEX
  using U_scalar_t = dcomplex;
#else
  using U_scalar_t = double;
#endif

  /// Type of the Monte-Carlo weight. Either double or dcomplex
  using mc_weight_t = decltype(U_scalar_t{} * g_tau_scalar_t{});

  /// The type of a block_matrix (e.g. density)
  using block_matrix_t = std::vector<matrix<g_tau_scalar_t>>;

  /// A view to a block_matrix_t
  using block_matrix_v_t = std::vector<matrix_view<g_tau_scalar_t>>;

  /// The type of a block_vector (e.g. diagonal densities)
  using block_vector_t = std::vector<nda::array<g_tau_scalar_t, 1>>;

  /// A view to a block_vector_t
  using block_vector_v_t = std::vector<nda::array_view<g_tau_scalar_t, 1>>;

  /// Container type of $\chi_2$ in DLR Matsubara frequencies
  using chi2_iw_t = block2_gf<mesh::dlr_imfreq, tensor_valued<4>>;

  /// Container type of $\chi_2$ in DLR imaginary time
  using chi2_tau_t = block2_gf<mesh::dlr_imtime, tensor_valued<4>>;

  /// A view to a chi2_tau_t
  using chi2_tau_v_t = chi2_tau_t::view_type;

  /// A const_view to a chi2_tau_t
  using chi2_tau_cv_t = chi2_tau_t::const_view_type;

  /// Container type of $\chi_3$ in Matsubara frequencies
  using chi3_iw_t = block2_gf<prod<imfreq, imfreq>, tensor_valued<4>>;

  /// A view to a chi3_iw_t
  using chi3_iw_v_t = chi3_iw_t::view_type;

  /// A const_view to a chi3_iw_t
  using chi3_iw_cv_t = chi3_iw_t::const_view_type;

  /// Container type of $\chi_3$ in Matsubara frequencies on 2D DLR grid
  using chi3_dlr2d_iw_t = block2_gf<mesh::dlr2d_imfreq, tensor_valued<4>>;

  /// A view to a chi3_dlr2d_iw_t
  using chi3_dlr2d_iw_v_t = chi3_dlr2d_iw_t::view_type;

  /// A const_view to a chi3_dlr2d_iw_t
  using chi3_dlr2d_iw_cv_t = chi3_dlr2d_iw_t::const_view_type;

  /// Container type of $\chi_3$ in imaginary time
  using chi3_tau_t = block2_gf<prod<imtime, imtime>, tensor_valued<4>>;

  /// A view to a chi3_tau_t
  using chi3_tau_v_t = chi3_tau_t::view_type;

  /// A const_view to a chi3_tau_t
  using chi3_tau_cv_t = chi3_tau_t::const_view_type;

  /// Container type of two-particle Green and Vertex functions in Matsubara frequencies
  using chi4_iw_t = block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>;

  /// Container type of two-particle Green and Vertex functions in imaginary time
  using chi4_tau_t = block2_gf<prod<imtime, imtime, imtime>, tensor_valued<4>>;

  // Declare some placeholders for the rest of the code. Use anonymous namespace for proper linkage
  // in this code, all variables with trailing _ are placeholders by convention.
  namespace {
    using nda::clef::placeholder;
    const placeholder<0> i_;
    const placeholder<1> j_;
    const placeholder<2> k_;
    const placeholder<3> l_;
    const placeholder<4> iw_;
    const placeholder<5> iw1_;
    const placeholder<6> iw2_;
    const placeholder<7> iw3_;
    const placeholder<8> iw4_;
    const placeholder<9> t_;
    const placeholder<10> t1_;
    const placeholder<11> t2_;
    const placeholder<12> t3_;
    const placeholder<13> bl_;
    const placeholder<14> bl1_;
    const placeholder<15> bl2_;
    const placeholder<16> iW_;
    const placeholder<17> iwp_;
  } // anonymous namespace

} // namespace triqs_ctint

namespace triqs::gfs {

  /// The maximum's norm of a triqs array. Takes element-wise absolute values of the data and returns the maximum of that.
  //FIXME Implement with is_array trait
  template <typename Value_t, int Rank> double max_norm(array_const_view<Value_t, Rank> const &arr) {
    auto max_itr = std::max_element(arr.begin(), arr.end(), [](auto a, auto b) { return std::abs(a) < std::abs(b); });
    return std::abs(*max_itr);
  }
  //FIXME array_const_view should be constructable from array_view
  template <typename Value_t, int Rank> double max_norm(array_view<Value_t, Rank> const &arr) { return max_norm(make_const_view(arr)); }

  /// The maximum's norm of a triqs Green function. Returns the max_norm of the data array.
  template <typename Gf> std::enable_if_t<is_gf_v<Gf>, double> max_norm(Gf const &G) { return max_norm(G.data()); }

  template <typename M, typename Target = tensor_valued<4>> block2_gf<M, Target> make_block2_gf(M const &m, gf_struct_t const &gf_struct) {

    std::vector<std::vector<gf<M, Target>>> gf_vecvec;
    std::vector<std::string> block_names;

    for (auto const &[bl1, bl1_size] : gf_struct) {
      block_names.push_back(bl1);
      std::vector<gf<M, Target>> gf_vec;
      for (auto const &[bl2, bl2_size] : gf_struct) {
        if constexpr (Target::rank == 4)
          gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
        else
          gf_vec.emplace_back(m, make_shape(bl1_size, bl2_size));
      }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
  }

  template <typename M1, typename M2, typename Target = tensor_valued<4>>
  block2_gf<M1, Target> make_block2_gf(M1 const &m, block2_gf_const_view<M2, Target> g_in) {

    std::vector<std::vector<gf<M1, Target>>> gf_vecvec;

    int n_blocks0 = g_in.block_names()[0].size();
    int n_blocks1 = g_in.block_names()[1].size();

    for (int i : range(n_blocks0)) {
      std::vector<gf<M1, Target>> gf_vec;
      for (int j : range(n_blocks1)) { gf_vec.emplace_back(m, g_in(i, j).target_shape()); }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(g_in.block_names()[0], g_in.block_names()[1], std::move(gf_vecvec));
  }

  template <typename Scalar_t> std::vector<matrix<Scalar_t>> make_block_vector(gf_struct_t const &gf_struct) {

    std::vector<matrix<Scalar_t>> res;
    for (auto const &[bl, bl_size] : gf_struct) { res.emplace_back(nda::zeros<Scalar_t>(bl_size, bl_size)); }
    return res;
  }
} // namespace triqs::gfs

namespace triqs::operators {
  /// Check if monomial is density-density interaction
  inline bool is_densdens_interact(monomial_t m) { return m.size() == 4 and m[0].indices == m[3].indices and m[1].indices == m[2].indices; }

  /// Function that returns a pair of integer indices (block, non_block), given the index of a c/c^dag operator
  inline std::pair<int, int> get_int_indices(canonical_ops_t const &op, hilbert_space::gf_struct_t const &gf_struct) {

    // The Fundamental operator-set allows for easy check of index validity
    hilbert_space::fundamental_operator_set fs(gf_struct);
    if (!fs.has_indices(op.indices)) TRIQS_RUNTIME_ERROR << " Index of c/c^+ operator not compatible with Green Function structure ";

    // Get block-name with apply visitor, lambda(0) is called to determine return type ...
    std::string op_bl_name = visit([](auto idx) { return std::to_string(idx); }, op.indices[0]);
    long nonbl_int_idx     = std::get<long>(op.indices[1]);

    // Capture positions in block and nonblock list
    for (auto [bl_int_idx, bl] : itertools::enumerate(gf_struct)) {
      auto const &[bl_name, bl_size] = bl;
      if (bl_name == op_bl_name and 0 <= nonbl_int_idx and nonbl_int_idx < bl_size) { return std::make_pair(bl_int_idx, nonbl_int_idx); }
    }
    TRIQS_RUNTIME_ERROR << "Error: Failed to retrieve integer indices for operator";
  }

  // Classification of bilinear operator monomials
  enum class bilinear_type { cdag_c, cdag_cdag, c_c };

  // Function that takes a bilinear operator Op = Sum_i a_i X_i Y_i (where X,Y are c or c†)
  // and returns a vector<tuple> with v[i] = (coef_i, type_i, (bl_X, bl_Y), (idx_X, idx_Y))
  inline auto get_terms(many_body_operator const &A, hilbert_space::gf_struct_t const &gf_struct) {
    std::vector<std::tuple<std::complex<double>, bilinear_type, std::pair<int, int>, std::pair<int, int>>> terms;
    for (auto const &term : A) {
      auto const &m = term.monomial;
      if (m.size() != 2) TRIQS_RUNTIME_ERROR << "chiAB operator monomial must have exactly 2 operators, got " << m.size();

      bilinear_type type;
      if (m[0].dagger && !m[1].dagger)
        type = bilinear_type::cdag_c;
      else if (m[0].dagger && m[1].dagger)
        type = bilinear_type::cdag_cdag;
      else if (!m[0].dagger && !m[1].dagger)
        type = bilinear_type::c_c;
      else
        TRIQS_RUNTIME_ERROR << "Unexpected c c† monomial (not in canonical form)";

      auto [bl1, i] = get_int_indices(m[0], gf_struct);
      auto [bl2, j] = get_int_indices(m[1], gf_struct);
      terms.emplace_back(term.coef, type, std::make_pair(bl1, bl2), std::make_pair(i, j));
    }
    return terms;
  }

} // namespace triqs::operators

// Useful macros

#define STR(x) #x
#define STRINGIZE(x) STR(x)

#ifdef DEBUG_CTINT
#define TRIQS_EXCEPTION_SHOW_CPP_TRACE
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define DEBUG(X) std::cerr << AS_STRING(X) << " = " << X << "      at " << __FILENAME__ << ':' << __LINE__ << std::endl
#define BREAK(X)                                                                                                                                     \
  std::cerr << X << " ... " << std::endl;                                                                                                            \
  getchar()
#else
#define DEBUG(X)
#define BREAK(X)
#endif
