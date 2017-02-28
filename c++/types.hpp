#pragma once

#ifdef DEBUG_CTINT
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#define TRIQS_EXCEPTION_SHOW_CPP_TRACE
#include <iostream>
#include <string.h>
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define DEBUG(X) std::cerr << AS_STRING(X) << " = " << X << "      at " << __FILENAME__ << ':' << __LINE__ << '\n'
#define BREAK(X)                                                                                                                                     \
  std::cerr << X << " ... \n";                                                                                                                       \
  getchar()
#define PRINT(X) std::cerr << "\n ========= " << X << " ========= \n"
#else
#define DEBUG(X)
#define BREAK(X)
#define PRINT(X)
#endif

#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/operators/util/extractors.hpp>

namespace triqs::gfs {

  /// The structure of the gf : block_name -> [...]= list of indices (int/string). FIXME Change to pair of vec<str> and vec<int> or vec<pair<str,int>>
  using block_gf_structure_t = std::map<std::string, std::vector<triqs::utility::variant_int_string>>;

  // Function template for block_gf initialization
  template <typename Var_t> block_gf<Var_t, matrix_valued> make_block_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {

    std::vector<gf<Var_t, matrix_valued>> gf_vec;
    std::vector<std::string> block_names;

    //for (auto const & [ bname, idx_lst ] : gf_struct) { // C++17
    for (auto const &bl : gf_struct) {
      auto &bname  = bl.first;
      auto bl_size = bl.second.size();
      block_names.push_back(bname);
      gf_vec.emplace_back(m, make_shape(bl_size, bl_size));
    }

    return make_block_gf(std::move(block_names), std::move(gf_vec));
  }

  template <typename Var_t> block2_gf<Var_t, tensor_valued<4>> make_block2_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {

    std::vector<std::vector<gf<Var_t, tensor_valued<4>>>> gf_vecvec;
    std::vector<std::string> block_names;

    for (auto const &bl1 : gf_struct) {
      auto &bname  = bl1.first;
      int bl1_size = bl1.second.size();
      block_names.push_back(bname);

      std::vector<gf<Var_t, tensor_valued<4>>> gf_vec;
      for (auto const &bl2 : gf_struct) {
        int bl2_size = bl2.second.size();
        gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
      }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
  }

} // namespace triqs::gfs

namespace triqs_ctint {

  using namespace triqs::gfs;
  using namespace triqs::arrays;
  using namespace triqs::operators;
  using namespace triqs::operators::util;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace triqs::h5;

  /// The channel type
  enum class Chan_t { PP, PH, XPH };

  /// Type of the Monte-Carlo weight. Either double or dcomplex
  using mc_weight_t = double;

  /// Container type of the alpha function. alpha[block](orbital,aux_spin)
  using alpha_t = std::vector<array<double, 2>>;

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
    triqs::clef::placeholder<8> iW_;
    triqs::clef::placeholder<9> t_;
    triqs::clef::placeholder<10> t1_;
    triqs::clef::placeholder<11> t2_;
    triqs::clef::placeholder<12> t3_;
  } // anonymous namespace

} // namespace triqs_ctint
