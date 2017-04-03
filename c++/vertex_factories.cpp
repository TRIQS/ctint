#include "./vertex_factories.hpp"

namespace triqs_ctint {

  std::vector<vertex_factory_t> make_vertex_factories(params_t const &params, triqs::mc_tools::random_generator &rng,
                                                      std::optional<block_gf_const_view<imfreq, matrix_valued>> D0_iw,
                                                      std::optional<gf_const_view<imfreq, matrix_valued>> Jperp_iw) {

    std::vector<vertex_factory_t> vertex_factories;

    // ----- TODO Rework interface

    {
      std::vector<vertex_idx_t> indices;
      std::vector<dcomplex> amplitudes;
      auto U = params.get_U();

      for (int a = 0; a < U.shape(0); ++a)
        for (int b = 0; b < U.shape(1); ++b)
          for (int i = 0; i < U(a, b).shape(0); ++i)
            for (int j = 0; j < U(a, b).shape(1); ++j) {
              if (std::abs(U(a, b)(i, j)) <= std::numeric_limits<double>::epsilon()) continue;
              indices.push_back({a, i, a, i, b, j, b, j});
              amplitudes.push_back(-U(a, b)(i, j) / 2.0 / params.n_s);
            }

      // move capture: moving the two vectors into the lambda and rename them
      auto l = [ beta = params.beta, n_s = params.n_s, indices = std::move(indices), amplitudes = std::move(amplitudes), &rng ] {
        int n             = rng(indices.size());
        tau_t t           = tau_t::get_random(rng);
        int s             = rng(n_s);
        double prop_proba = 1.0 / (beta * indices.size() * n_s);
        return vertex_t{indices[n], t, t, t, t, real(amplitudes[n]), prop_proba, true, s}; // TODO Real/Imag
      };

      vertex_factories.push_back(l);
    } // clean temporaries

    // ------------ Create Vertex Factory for Dynamic Density-Density Interactions --------------
    if (D0_iw) {

      std::vector<vertex_idx_t> indices;
      std::vector<gf_const_view<imtime, scalar_valued>> ktau_s;

      //--density-density
      for (int sig = 0; sig < (*D0_iw).size(); sig++) {                                            // TODO rename, not sig!!!
        auto D = make_gf_from_inverse_fourier((*D0_iw)[sig], params.n_tau_dynamical_interactions); // Use block_gf fourier
        for (int a = 0; a < D.target_shape()[0]; a++)
          for (int b = 0; b < D.target_shape()[1]; b++) {

            int bl1 = sig / params.n_blocks();
            int bl2 = sig % params.n_blocks();

            indices.push_back({bl1, a, bl1, a, bl2, b, bl2, b});

            auto d = slice_target_to_scalar(D, a, b);
            ktau_s.push_back(d);
          }
      }

      auto l = [ beta = params.beta, n_s = params.n_s, indices = std::move(indices), ktau_s = std::move(ktau_s), &rng ] {
        int n             = rng(indices.size());
        tau_t t           = tau_t::get_random(rng);
        tau_t tp          = tau_t::get_random(rng);
        double d_tau      = cyclic_difference(t, tp);
        int s             = rng(n_s);
        double prop_proba = 1.0 / (beta * beta * indices.size() * n_s);
        return vertex_t{indices[n], t, t, tp, tp, -real(ktau_s[n](d_tau)) / 2.0 / n_s, prop_proba, true, s};
      };

      vertex_factories.push_back(l);
    }

    // ------------ Create Vertex Factory for Dynamic Spin-Spin Interactions --------------
    if (Jperp_iw) {

      auto J = make_gf_from_inverse_fourier(*Jperp_iw, params.n_tau_dynamical_interactions);
      std::vector<vertex_idx_t> indices;
      std::vector<gf_const_view<imtime, scalar_valued>> ktau_s;

      // Jperp requires that blocks are of same size and that they represent spins which means there must be exactly two blocks
      for (int a = 0; a < J.target_shape()[0]; a++)
        for (int b = 0; b < J.target_shape()[1]; b++) {
          auto d = slice_target_to_scalar(J, a, b);

          for (int bl = 0; bl < 2; bl++) {                           // FIXME Ugly ...
            indices.push_back({1 - bl, a, bl, a, bl, b, 1 - bl, b}); //S^+_a(tau) S^-_b(tau')
            ktau_s.push_back(d);
          }
        }

      auto l = [ beta = params.beta, indices = std::move(indices), ktau_s = std::move(ktau_s), &rng ] {
        int n             = rng(indices.size());
        tau_t t           = tau_t::get_random(rng);
        tau_t tp          = tau_t::get_random(rng);
        double d_tau      = cyclic_difference(t, tp);
        double prop_proba = 1.0 / (beta * beta * indices.size());
        return vertex_t{indices[n], t, t, tp, tp, real(ktau_s[n](d_tau)) / 4.0, prop_proba}; // CHECK why 1/4 ?? Check sign
      };

      vertex_factories.push_back(l);
    }

    return vertex_factories;
  }

} // namespace triqs_ctint
