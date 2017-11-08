// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./vertex_factories.hpp"
#include "types.hpp"

namespace triqs_ctint {

  std::vector<vertex_factory_t> make_vertex_factories(params_t const &params, triqs::mc_tools::random_generator &rng,
                                                      std::optional<block2_gf_const_view<imfreq, matrix_valued>> D0_iw,
                                                      std::optional<gf_const_view<imfreq, matrix_valued>> Jperp_iw) {

    std::vector<vertex_factory_t> vertex_factories;

    // ------------ Create Vertex Factory for Static Interactions --------------
    {
      std::vector<vertex_idx_t> indices;
      std::vector<U_scalar_t> amplitudes;

      // Loop over the interaction hamiltonian and insert the corresponding indices/amplitiudes into the lists for the factory
      for (auto const &term : params.h_int) {
        amplitudes.push_back(-U_scalar_t(term.coef));
        auto const &m = term.monomial;
        if (m.size() != 4 or !(m[0].dagger and m[1].dagger and !m[2].dagger and !m[3].dagger))
          TRIQS_RUNTIME_ERROR << " Monimial in h_int is not of the form c^+ c^+ c c \n";
        // Careful: h_int monomials automatically ordered as c^+_0 c^+_1 c_2 c_3
        auto [bl_cdag_0, idx_cdag_0] = get_int_indices(m[0], params.gf_struct);
        auto [bl_cdag_1, idx_cdag_1] = get_int_indices(m[1], params.gf_struct);
        auto [bl_c_1, idx_c_1]       = get_int_indices(m[2], params.gf_struct);
        auto [bl_c_0, idx_c_0]       = get_int_indices(m[3], params.gf_struct);

        indices.push_back({bl_cdag_0, idx_cdag_0, bl_c_0, idx_c_0, bl_cdag_1, idx_cdag_1, bl_c_1, idx_c_1});
      }

      if (indices.size() > 0) {
        // Construct the factory and insert into list
        auto l = [beta = params.beta, n_s = params.n_s, indices = std::move(indices), amplitudes = std::move(amplitudes), &rng] {
          int n                     = rng(indices.size());
          auto &idx                 = indices[n];
          bool is_densdens_interact = (idx.b1 == idx.b2) and (idx.b3 == idx.b4) and (idx.u1 == idx.u2) and (idx.u3 == idx.u4);
          tau_t t                   = tau_t::get_random(rng);
          int s                     = is_densdens_interact ? rng(n_s) : 0;
          double prop_proba         = 1.0 / (beta * indices.size() * n_s);
          return vertex_t{indices[n], t, t, t, t, amplitudes[n] / n_s, prop_proba, n, s};
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

      auto tau_mesh = mesh::imtime{params.beta, Boson, params.n_tau_dynamical_interactions};
      auto D0_tau   = make_gf_from_fourier(*D0_iw, tau_mesh, make_zero_tail(*D0_iw, 2));

      // Loop over block indices
      for (int bl1 : range((*D0_iw).size1()))
        for (int bl2 : range((*D0_iw).size1())) {

          // Loop over non-block indices
          for (int a = 0; a < D0_tau(bl1, bl2).target_shape()[0]; a++)
            for (int b = 0; b < D0_tau(bl1, bl2).target_shape()[1]; b++) {
              auto d = slice_target_to_scalar(D0_tau(bl1, bl2), a, b);

              // Add Vertex generator only if d is non-zero
              if (max_norm(d) > 1e-10) {
                indices.push_back({bl1, a, bl1, a, bl2, b, bl2, b});
#ifdef INTERACTION_IS_COMPLEX
                D0_tau_lst.emplace_back(d);
#else
                if (!is_gf_real(d, 1e-8))
                  std::cerr << "WARNING: Assuming real interaction values, but found Imag(D0(tau)) > 1e-8. Casting to Real.\n";
                D0_tau_lst.emplace_back(real(d));
#endif
              }
            }
        }

      if (indices.size() > 0) {
        auto l = [beta = params.beta, n_s = params.n_s, indices = std::move(indices), D0_tau_lst = std::move(D0_tau_lst), &rng] {
          int n             = rng(indices.size());
          tau_t t           = tau_t::get_random(rng);
          tau_t tp          = tau_t::get_random(rng);
          auto [sig, dtau]  = cyclic_difference(t, tp);
          int s             = rng(n_s);
          double prop_proba = 1.0 / (beta * beta * indices.size() * n_s);
          return vertex_t{indices[n], t, t, tp, tp, -D0_tau_lst[n](dtau) / n_s, prop_proba, true, s};
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

      if (params.n_blocks() != 2) TRIQS_RUNTIME_ERROR << "Jperp requires exactly two blocks corresponding to the spins";

      auto tau_mesh  = mesh::imtime{params.beta, Boson, params.n_tau_dynamical_interactions};
      auto Jperp_tau = make_gf_from_fourier(*Jperp_iw, tau_mesh, make_zero_tail(*Jperp_iw, 2));

      // Loop over non-block indices
      for (int a = 0; a < Jperp_tau.target_shape()[0]; a++)
        for (int b = 0; b < Jperp_tau.target_shape()[1]; b++) {
          auto d = slice_target_to_scalar(Jperp_tau, a, b);

          if (max_norm(d) > 1e-10) {
            indices.push_back({1, a, 0, a, 0, b, 1, b}); //S^+_a(tau) S^-_b(tau')
            indices.push_back({0, a, 1, a, 1, b, 0, b}); //S^-_a(tau) S^+_b(tau')
#ifdef INTERACTION_IS_COMPLEX
            Jperp_tau_lst.emplace_back(d);
            Jperp_tau_lst.emplace_back(d);
#else
            if (!is_gf_real(d, 1e-8)) std::cerr << "WARNING: Assuming real interaction values, but found Imag(Jperp(tau)) > 1e-8. Casting to Real.\n";
            Jperp_tau_lst.emplace_back(real(d));
            Jperp_tau_lst.emplace_back(real(d));
#endif
          }
        }

      if (indices.size() > 0) {
        auto l = [beta = params.beta, indices = std::move(indices), Jperp_tau_lst = std::move(Jperp_tau_lst), &rng] {
          int n             = rng(indices.size());
          tau_t t           = tau_t::get_random(rng);
          tau_t tp          = tau_t::get_random(rng);
          auto [sig, dtau]  = cyclic_difference(t, tp);
          double prop_proba = 1.0 / (beta * beta * indices.size());
          return vertex_t{indices[n], t, t, tp, tp, Jperp_tau_lst[n](dtau) / 2.0, prop_proba}; // We add two identical terms above -> Divide by 2
        };

        vertex_factories.emplace_back(l);
      }
    }

    return vertex_factories;
  }

} // namespace triqs_ctint
