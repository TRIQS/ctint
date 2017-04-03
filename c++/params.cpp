#include "./params.hpp"

namespace triqs_ctint {

  array<array<dcomplex, 2>, 2> params_t::get_U() const { // TODO Clean and Rework!

    // Calculate U-Matrix U(block1,block2)(i,j)
    array<array<dcomplex, 2>, 2> U(n_blocks(), n_blocks());

    // Extract flat U-Matrix, where indeces combine the block (non-leading) and non-block (leading) index
    array<dcomplex, 2> Uall = dict_to_matrix(extract_U_dict2(h_int), gf_struct);

    int ibl = 0, jbl = 0, offset_col = 0, offset_row = 0;
    for (auto const &bl1 : gf_struct) {
      jbl        = 0;
      offset_col = 0;
      for (auto const &bl2 : gf_struct) {
        U(ibl, jbl).resize(bl1.second.size(), bl2.second.size());
        for (int i = 0; i < bl1.second.size(); i++)
          for (int j = 0; j < bl2.second.size(); j++) { U(ibl, jbl)(i, j) = Uall(offset_row + i, offset_col + j); } // j
        offset_col += bl2.second.size();
        jbl++;
      } // bl2
      offset_row += bl1.second.size();
      ibl++;
    } // bl1

    return U;
  }

} // namespace triqs_ctint
