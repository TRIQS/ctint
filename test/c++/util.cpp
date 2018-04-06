#include <triqs_ctint/params.hpp>
#include <triqs/test_tools/many_body_operator.hpp>
#include <iostream>

using namespace triqs_ctint;

// Decompose an interaction hamiltonian into density and non-density term
TEST(util, hint_decomp) { // NOLINT

  double U = 1.0;
  double J = 1.0;

  auto sx_0 = c_dag("up", 0) * c("down", 0) + c_dag("down", 0) * c("up", 0);
  auto sx_1 = c_dag("up", 1) * c("down", 1) + c_dag("down", 1) * c("up", 1);

  auto h_int = U * n("up", 0) * n("down", 1) + J * sx_0 * sx_1;

  std::cout << " H_int reads : " << h_int << "\n\n";

  decltype(h_int) h_dens, h_nondens;

  for (auto const &term : h_int) {
    if (is_densdens_interact(term.monomial))
      h_dens += term;
    else
      h_nondens += term;
  }

  std::cout << " Density-Density terms:\n" << h_dens << "\n";
  std::cout << " Non-Density-Density terms\n" << h_nondens << "\n";

  EXPECT_OPERATOR_NEAR(h_dens, U * n("up", 0) * n("down", 1));
  EXPECT_OPERATOR_NEAR(h_nondens, J * sx_0 * sx_1);

}

// Given a gf_struct object, determine the integer indices of a given canonical operator
TEST(util, get_op_indices) { // NOLINT

  gf_struct_t gf_struct = {{"bl1", {"o1", "o2", "o3"}}, {"bl2", {"o4", "o5"}}};

  many_body_operator_generic<double> myop;

  std::cout << " gf block structure: { ";
  for (auto &bl : gf_struct) {
    std::cout << bl.first << ": { ";
    for (auto &idx : bl.second) {
      std::cout << to_string(idx) << " ";
      myop += c(bl.first, idx);
    }
    std::cout << "} ";
  }
  std::cout << "}\n\n";

  for (auto const &term : myop)
    for (auto const &op : term.monomial) {
      auto idx_pair = get_int_indices(op, gf_struct);
      std::cout << " Integer indices for " << op << " : " << idx_pair.first << ", " << idx_pair.second << "\n";
    }

  EXPECT_EQ(std::make_pair(1, 1), get_int_indices(canonical_ops_t{false, {"bl2", "o5"}}, gf_struct)); // NOLINT
  EXPECT_EQ(std::make_pair(0, 2), get_int_indices(canonical_ops_t{false, {"bl1", "o3"}}, gf_struct)); // NOLINT
}

MAKE_MAIN;
