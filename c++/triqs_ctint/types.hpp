/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, N. Wentzell
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <triqs/gfs.hpp>
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
  using namespace triqs::gfs;
  using namespace triqs::arrays;
  using namespace triqs::operators;
  using namespace triqs::operators::utils;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace h5;

  using namespace itertools;

  /// The channel type
  enum class Chan_t { PP, PH, XPH };

  /// Container type of the alpha function. alpha(vertex_index, aux_spin)
  using alpha_t = array<double, 4>;

  /// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  using triqs::hilbert_space::gf_struct_t;

  /// Container type of one-particle Green and Vertex functions in imaginary times
#ifdef GTAU_IS_COMPLEX
  using g_tau_t = block_gf<imtime, matrix_valued>;
#else
  using g_tau_t        = block_gf<imtime, matrix_real_valued>;
#endif

  /// A const_view to a g_tau_t
  using g_tau_cv_t = g_tau_t::const_view_type;

  /// A view to a g_tau_t
  using g_tau_v_t = g_tau_t::view_type;

  /// Scalar type of g_tau
  using g_tau_scalar_t = g_tau_t::g_t::scalar_t;

  /// Container type of one-particle Green and Vertex functions in Matsubara frequencies
  using g_iw_t = block_gf<imfreq, matrix_valued>;

  /// A const_view to a g_iw_t
  using g_iw_cv_t = g_iw_t::const_view_type;

  /// A view to a g_iw_t
  using g_iw_v_t = g_iw_t::view_type;

  /// The target_type of the intermediate scattering matrices
#if defined GTAU_IS_COMPLEX || defined INTERACTION_IS_COMPLEX
  using M_tau_target_t = matrix_valued;
#else
  using M_tau_target_t = matrix_real_valued;
#endif
  using M_tau_scalar_t = M_tau_target_t::scalar_t;

  /// Scalar type of all interaction vertices
#ifdef INTERACTION_IS_COMPLEX
  using U_scalar_t = dcomplex;
#else
  using U_scalar_t     = double;
#endif

  /// Type of the Monte-Carlo weight. Either double or dcomplex
  using mc_weight_t = decltype(U_scalar_t{} * g_tau_scalar_t{});

  /// The type of a block_matrix (e.g. density)
  using block_matrix_t = std::vector<matrix<M_tau_scalar_t>>;

  /// A view to a block_matrix_t
  using block_matrix_v_t = std::vector<matrix_view<M_tau_scalar_t>>;

  /// Container type of $\chi_3$ in Matsubara frequencies
  using chi2_iw_t = block2_gf<imfreq, tensor_valued<4>>;

  /// Container type of $\chi_3$ in imaginary time
  using chi2_tau_t = block2_gf<imtime, tensor_valued<4>>;

  /// A view to a chi2_tau_t
  using chi2_tau_v_t = chi2_tau_t::view_type;

  /// A const_view to a chi2_tau_t
  using chi2_tau_cv_t = chi2_tau_t::const_view_type;

  /// Container type of $\chi_3$ in Matsubara frequencies
  using chi3_iw_t = block2_gf<cartesian_product<imfreq, imfreq>, tensor_valued<4>>;

  /// A view to a chi3_iw_t
  using chi3_iw_v_t = chi3_iw_t::view_type;

  /// A const_view to a chi3_iw_t
  using chi3_iw_cv_t = chi3_iw_t::const_view_type;

  /// Container type of $\chi_3$ in imaginary time
  using chi3_tau_t = block2_gf<cartesian_product<imtime, imtime>, tensor_valued<4>>;

  /// A view to a chi3_tau_t
  using chi3_tau_v_t = chi3_tau_t::view_type;

  /// A const_view to a chi3_tau_t
  using chi3_tau_cv_t = chi3_tau_t::const_view_type;

  /// Container type of two-particle Green and Vertex functions in Matsubara frequencies
  using chi4_iw_t = block2_gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;

  /// Container type of two-particle Green and Vertex functions in imaginary time
  using chi4_tau_t = block2_gf<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>;

  // Declare some placeholders for the rest of the code. Use anonymous namespace for proper linkage
  // in this code, all variables with trailing _ are placeholders by convention.
  namespace {
    triqs::clef::placeholder<0> i_;
    triqs::clef::placeholder<1> j_;
    triqs::clef::placeholder<2> k_;
    triqs::clef::placeholder<3> l_;
    triqs::clef::placeholder<4> iw_;
    triqs::clef::placeholder<5> iw1_;
    triqs::clef::placeholder<6> iw2_;
    triqs::clef::placeholder<7> iw3_;
    triqs::clef::placeholder<8> iw4_;
    triqs::clef::placeholder<9> t_;
    triqs::clef::placeholder<10> t1_;
    triqs::clef::placeholder<11> t2_;
    triqs::clef::placeholder<12> t3_;
    triqs::clef::placeholder<13> bl_;
    triqs::clef::placeholder<14> bl1_;
    triqs::clef::placeholder<15> bl2_;
    triqs::clef::placeholder_prime<0> iW_;
    triqs::clef::placeholder_prime<1> iwp_;
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

  template <typename Var_t, typename Target = tensor_valued<4>>
  block2_gf<Var_t, Target> make_block2_gf(gf_mesh<Var_t> const &m, triqs::hilbert_space::gf_struct_t const &gf_struct) {

    std::vector<std::vector<gf<Var_t, Target>>> gf_vecvec;
    std::vector<std::string> block_names;

    for (auto const &bl1 : gf_struct) {
      auto &bname  = bl1.first;
      int bl1_size = bl1.second.size();
      block_names.push_back(bname);
      std::vector<std::string> indices1;
      for (auto const &var : bl1.second) visit([&indices1](auto &&arg) { indices1.push_back(std::to_string(arg)); }, var);

      std::vector<gf<Var_t, Target>> gf_vec;
      for (auto const &bl2 : gf_struct) {
        int bl2_size = bl2.second.size();
        std::vector<std::string> indices2;
        for (auto const &var : bl2.second) visit([&indices2](auto &&arg) { indices2.push_back(std::to_string(arg)); }, var);
        if constexpr (Target::rank == 4)
          gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size),
                              std::vector<std::vector<std::string>>{indices1, indices1, indices2, indices2});
        else
          gf_vec.emplace_back(m, make_shape(bl1_size, bl2_size), std::vector<std::vector<std::string>>{indices1, indices2});
      }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
  }

  template <typename Var1_t, typename Var2_t, typename Target = tensor_valued<4>>
  block2_gf<Var1_t, Target> make_block2_gf(gf_mesh<Var1_t> const &m, block2_gf_const_view<Var2_t, Target> g_in) {

    std::vector<std::vector<gf<Var1_t, Target>>> gf_vecvec;

    int n_blocks0 = g_in.block_names()[0].size();
    int n_blocks1 = g_in.block_names()[1].size();

    for (int i : range(n_blocks0)) {
      std::vector<gf<Var1_t, Target>> gf_vec;
      for (int j : range(n_blocks1)) { gf_vec.emplace_back(m, g_in(i, j).target_shape(), g_in(i, j).indices()); }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(g_in.block_names()[0], g_in.block_names()[1], std::move(gf_vecvec));
  }

  template <typename Var_t> auto bin_to_mesh(double val, gf_mesh<Var_t> const &m) {
    return gf_closest_point<Var_t, int>::invoke(m, closest_mesh_pt(val));
  }

  template <typename Scalar_t> std::vector<matrix<Scalar_t>> make_block_vector(hilbert_space::gf_struct_t const &gf_struct) {

    std::vector<matrix<Scalar_t>> res;
    for (auto const &bl : gf_struct) {
      int bl_size = bl.second.size();
      res.emplace_back(zeros<Scalar_t>(make_shape(bl_size, bl_size)));
    }
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

    // Capture positions in block and nonblock list
    for (auto [bl_int_idx, bl] : itertools::enumerate(gf_struct)) {
      auto const &[bl_name, idx_lst] = bl;
      if (bl_name == op_bl_name) {
        int nonbl_int_idx = std::distance(idx_lst.cbegin(), std::find(idx_lst.cbegin(), idx_lst.cend(), op.indices[1]));
        return std::make_pair(bl_int_idx, nonbl_int_idx);
      }
    }
    TRIQS_RUNTIME_ERROR << "Error: Failed to retrieve integer indices for operator";
  }

  // Function that takes a bosonic operator Op = Sum_i a_i c^+(bi, ui) c(bi, vi)
  // and returns a vector<tuple> with v[i] = (b_i, a_i, (c^+(bi, ui), c(bi, vi)))
  inline auto get_terms(many_body_operator const &A, hilbert_space::gf_struct_t const &gf_struct) {
    std::vector<std::tuple<std::complex<double>, std::pair<int, int>, std::pair<int, int>>> terms;
    for (auto const &term : A) {
      auto const &m = term.monomial;
      if (m.size() != 2 or !m[0].dagger or m[1].dagger)
        TRIQS_RUNTIME_ERROR << " Monomial in bosonic operator of chiAB measurement not of the proper form c^+ c \n";
      auto [bl1, i] = get_int_indices(m[0], gf_struct);
      auto [bl2, j] = get_int_indices(m[1], gf_struct);
      auto bl_pair  = std::make_pair(bl1, bl2);
      auto idx_pair = std::make_pair(i, j);
      terms.emplace_back(term.coef, bl_pair, idx_pair);
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
#define PRINT(X) std::cerr << "\n ========= " << X << " ========= " << std::endl
#else
#define DEBUG(X)
#define BREAK(X)
#define PRINT(X)
#endif
