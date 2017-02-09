#pragma once

#ifdef DEBUG_CTINT
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#define TRIQS_EXCEPTION_SHOW_CPP_TRACE
#include <iostream>
#include <string.h>
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define DEBUG(X) std::cerr << AS_STRING(X) << " = " << X << "      at " << __FILENAME__ << ':' << __LINE__ << '\n'
#define BREAK(X) std::cerr << X << " ... \n"; getchar() 
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

  /// The structure of the gf : block_name -> [...]= list of indices (int/string).
  using block_gf_structure_t = std::map<std::string, std::vector<triqs::utility::variant_int_string>>;

  // Function template for block_gf initialization
  template <typename Mesh> block_gf<typename Mesh::var_t, matrix_valued> make_block_gf(Mesh const &m, block_gf_structure_t const &gf_struct) {

    std::vector<gf<typename Mesh::var_t, matrix_valued>> gf_vec;
    std::vector<std::string> block_names;

    //for (auto const & [ bname, idx_lst ] : gf_struct) { // C++17
    for (auto it = gf_struct.begin(); it != gf_struct.end(); ++it) {
      auto &bname   = it->first;
      auto &idx_lst = it->second;
      block_names.push_back(bname);
      gf_vec.emplace_back(m, make_shape(idx_lst.size(), idx_lst.size()));
    }

    return make_block_gf(std::move(block_names), std::move(gf_vec));
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
