#include "./vertex_factories.hpp"
#include "types.hpp"

namespace triqs_ctint {

  std::vector<vertex_factory_t> make_vertex_factories(params_t const &params, triqs::mc_tools::random_generator &rng,
                                                      std::optional<block2_gf_const_view<imfreq, matrix_valued>> D0_iw,
                                                      std::optional<gf_const_view<imfreq, matrix_valued>> Jperp_iw) {

    std::vector<vertex_factory_t> vertex_factories;

    {
      std::vector<vertex_idx_t> indices;
      std::vector<U_scalar_t> amplitudes;

      // Loop over the interaction hamiltonian and insert the corresponding indices/amplitiudes into the lists for the factory
      for (auto const &term : params.h_int) {

        amplitudes.push_back(-U_scalar_t(term.coef)); // n_s not needed here, cancels against prop_proba
        auto const &m = term.monomial;
        if (m.size() != 4 or !(m[0].dagger && m[1].dagger && !m[2].dagger && !m[3].dagger))
          TRIQS_RUNTIME_ERROR << " Monimial in h_int is not of the form c^+ c^+ c c ";
        std::vector<std::pair<int, int>> vec;
        for (auto op : m) { vec.push_back(get_int_indices(op, params.gf_struct)); }
        // Careful: h_int monomials automatically ordered as c^+_1 c^+_2 c_3 c_4
        indices.push_back({vec[0].first, vec[0].second, vec[3].first, vec[3].second, vec[1].first, vec[1].second, vec[2].first, vec[2].second});
      }

      if (indices.size() > 0) {
        // Construct the factory and insert into list
        auto l = [ beta = params.beta, n_s = params.n_s, indices = std::move(indices), amplitudes = std::move(amplitudes), &rng ] {
          int n                     = rng(indices.size());
          auto &idx                 = indices[n];
          bool is_densdens_interact = (idx.b1 == idx.b2) && (idx.b3 == idx.b4) && (idx.u1 == idx.u2) && (idx.u3 == idx.u4);
          tau_t t                   = tau_t::get_random(rng);
          int s                     = is_densdens_interact ? rng(n_s) : 0;
          double prop_proba         = 1.0 / (beta * indices.size()); // n_s here not needed, cancels against amplitude
          return vertex_t{indices[n], t, t, t, t, amplitudes[n], prop_proba, is_densdens_interact, s};
        };

        vertex_factories.emplace_back(l);
      }
    } // clean temporaries

    // ------------ Create Vertex Factory for Dynamic Density-Density Interactions --------------
    if (D0_iw) {

      std::vector<vertex_idx_t> indices;
#ifdef INTERACTION_IS_COMPLEX
      std::vector<gf<imtime, scalar_valued>> D0_tau_lst;
#else
      std::vector<gf<imtime, scalar_real_valued>> D0_tau_lst;
#endif

      // Loop over block indices
      for (int bl1 : range((*D0_iw).size1()))
        for (int bl2 : range((*D0_iw).size1())) {

          auto D = make_gf_from_inverse_fourier((*D0_iw)(bl1, bl2), params.n_tau_dynamical_interactions);

          // Loop over non-block indices
          for (int a = 0; a < D.target_shape()[0]; a++)
            for (int b = 0; b < D.target_shape()[1]; b++) {
              auto d = slice_target_to_scalar(D, a, b);

              // Add Vertex generator only if d is non-zero
              if (max_norm(d) > 1e-10) {
                indices.push_back({bl1, a, bl1, a, bl2, b, bl2, b});
#ifdef INTERACTION_IS_COMPLEX
                D0_tau_lst.emplace_back(d);
#else
                if (!is_gf_real(d)) TRIQS_RUNTIME_ERROR << " Assuming real interaction values, but Imag(D0(tau)) > 0 ";
                D0_tau_lst.emplace_back(get_real(d));
#endif
              }
            }
        }

      if (indices.size() > 0) {
        auto l = [ beta = params.beta, n_s = params.n_s, indices = std::move(indices), D0_tau_lst = std::move(D0_tau_lst), &rng ] {
          int n             = rng(indices.size());
          tau_t t           = tau_t::get_random(rng);
          tau_t tp          = tau_t::get_random(rng);
          double d_tau      = cyclic_difference(t, tp);
          int s             = rng(n_s);
          double prop_proba = 1.0 / (beta * beta * indices.size() * n_s);
          return vertex_t{indices[n], t, t, tp, tp, -D0_tau_lst[n](d_tau) / n_s, prop_proba, true, s};
        };

        vertex_factories.emplace_back(l);
      }
    }

    // ------------ Create Vertex Factory for Dynamic Spin-Spin Interactions --------------
    if (Jperp_iw) {

      std::vector<vertex_idx_t> indices;
#ifdef INTERACTION_IS_COMPLEX
      std::vector<gf<imtime, scalar_valued>> Jperp_tau_lst;
#else
      std::vector<gf<imtime, scalar_real_valued>> Jperp_tau_lst;
#endif

      // Jperp requires that blocks are of same size and that they represent spins which means there must be exactly two blocks WHAT?
      auto J = make_gf_from_inverse_fourier(*Jperp_iw, params.n_tau_dynamical_interactions);

      // Loop over non-block indices
      for (int a = 0; a < J.target_shape()[0]; a++)
        for (int b = 0; b < J.target_shape()[1]; b++) {
          auto d = slice_target_to_scalar(J, a, b);

          if (max_norm(d) > 1e-10)
            for (int bl = 0; bl < 2; bl++) {                           // FIXME Ugly ...
              indices.push_back({1 - bl, a, bl, a, bl, b, 1 - bl, b}); //S^+_a(tau) S^-_b(tau')
#ifdef INTERACTION_IS_COMPLEX
              Jperp_tau_lst.emplace_back(d);
#else
              if (!is_gf_real(d)) TRIQS_RUNTIME_ERROR << " Assuming real interaction values, but Imag(Jperp(tau)) > 0 ";
              Jperp_tau_lst.emplace_back(get_real(d));
#endif
            }
        }

      if (indices.size() > 0) {
        auto l = [ beta = params.beta, indices = std::move(indices), Jperp_tau_lst = std::move(Jperp_tau_lst), &rng ] {
          int n             = rng(indices.size());
          tau_t t           = tau_t::get_random(rng);
          tau_t tp          = tau_t::get_random(rng);
          double d_tau      = cyclic_difference(t, tp);
          double prop_proba = 1.0 / (beta * beta * indices.size());
          return vertex_t{indices[n], t, t, tp, tp, Jperp_tau_lst[n](d_tau) / 4.0, prop_proba};
        };

        vertex_factories.emplace_back(l);
      }
    }

    return vertex_factories;
  }

} // namespace triqs_ctint
