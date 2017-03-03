#pragma once

namespace triqs_ctint {

  template <template <typename, typename> class Gf, typename M_out, typename M_in, typename T>
  typename Gf<M_out, T>::regular_type _make_gf_impl(typename Gf<M_in, T>::const_view_type const &G_in, int n) {
    auto out_mesh = gf_mesh<M_out>{G_in.mesh().domain(), n};
    auto G_out    = gf<M_out, T>{out_mesh, G_in.target_shape()};
    _fourier_impl(G_out(), G_in);
    return G_out;
  }

  // Factory to create Matsubara frequency Green functions from inverse fourier transform of block[2]_gf[_const][_view]
  template <template <typename, typename> class Gf, typename M_out, typename M_in, typename T>
  typename Gf<M_out, T>::regular_type _make_block_gf_impl(typename Gf<M_in, T>::const_view_type const &G_in, int n) {
    // block_gf[_const][_view]
    if
      constexpr(Gf<M_in, T>::arity == 1) {
        std::vector<gf<M_out, T>> G_out_vec;
        //for (auto const &bl : G_in) G_out_vec.emplace_back(make_gf_from_fourier_new(bl, n));
        for (auto const &bl : G_in) G_out_vec.emplace_back(_make_gf_impl<gf_const_view, M_out, M_in, T>(bl(), n));
        typename Gf<M_out, T>::regular_type G_out{G_in.block_names(), G_out_vec};
        return G_out;
      }
    // block2_gf[_const][_view]
    else if (Gf<M_in, T>::arity == 2) {
      std::vector<std::vector<gf<M_out, T>>> G_out_vecvec;
      for (int i : range(G_in.size1())) {
        std::vector<gf<M_out, T>> temp_vec;
        //for (int j : range(G_in.size2())) temp_vec.emplace_back(make_gf_from_fourier_new(G_in(i,j), n));
        for (int j : range(G_in.size2())) temp_vec.emplace_back(_make_gf_impl<gf_const_view, M_out, M_in, T>(G_in(i, j), n));
        G_out_vecvec.emplace_back(std::move(temp_vec));
      }
      typename Gf<M_out, T>::regular_type G_out{G_in.block_names(), G_out_vecvec};
      return G_out;
    }
  }

  // Factory to create Matsubara frequency Green functions from inverse fourier transform of block[2]_gf[_const][_view]
  template <template <typename, typename> class Gf, typename T>
  std::enable_if_t<is_block_gf_or_view<Gf<imtime, T>>::value, typename Gf<imfreq, T>::regular_type>
  make_gf_from_fourier(Gf<imtime, T> const &G_tau, int n_iw = -1) {
    return _make_block_gf_impl<Gf, imfreq, imtime, T>(G_tau, n_iw);
  }

  // Factory to create imaginary time Green functions from inverse fourier transform of block[2]_gf[_const][_view]
  template <template <typename, typename> class Gf, typename T>
  std::enable_if_t<is_block_gf_or_view<Gf<imfreq, T>>::value, typename Gf<imtime, T>::regular_type>
  make_gf_from_inverse_fourier(Gf<imfreq, T> const &G_iw, int n_tau = -1) {
    return _make_block_gf_impl<Gf, imtime, imfreq, T>(G_iw, n_tau);
  }

} // namespace triqs_ctint
