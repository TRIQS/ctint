#pragma once
#include "vertex.hpp"
#include "dets.hpp"

namespace triqs_ctint {

  /**
   * A stack of lazy operations in the determinants.
   * The class keeps two lists of the c and c^\dagger operators to be inserted.
   * Additional verticies can be added with the << operator.
   * The order of insertion is of no importance, while the insertion in the
   * determinant is ordered using the binary comparison of arg_t.
   * A call to execute_try_insert/execute_try_remove performs the full
   * insert/remove operation in all the determinants involved.
   */
  class lazy_det_operation_t {

    private:
    // Vector of all determinants
    std::vector<det_t> *dets;

    // This type keeps track of lazy inserts for single determinant
    struct one_block {
      static constexpr int n_max = 2;

      // Array of $c$ (x) and $c^\dagger$ (y) operators
      std::array<arg_t, n_max> x, y;

      // Array of the operator positions
      std::array<int, n_max> px, py;

      // Counter for the $c$ and $c^\dagger$
      int x_count = 0, y_count = 0;

      // Lazy-add a $c$ operator with argument x_
      void lazy_add_x(arg_t const &x_) {
        TRIQS_ASSERT(x_count < n_max);
        x[x_count++] = x_;
      }

      // Lazy-add a $c^\dagger$ operator with argument y_
      void lazy_add_y(arg_t const &y_) {
        TRIQS_ASSERT(y_count < n_max);
        y[y_count++] = y_;
      }

      // Tries to perform all insertions into the det, returning the det ratio
      double execute_try_insert(det_t *d);

      // Tries to perform all removals from the det, returning the det ratio
      double execute_try_remove(det_t *d);
    };

    // The list of lazy operations, one for each determinant
    std::vector<one_block> lazy_op_lst;

    public:
    lazy_det_operation_t(std::vector<det_t> *dets) : dets(dets), lazy_op_lst(dets->size()) {}

    /// Clean all registered moves
    void reset() {
      for (auto &l : lazy_op_lst) {
        l.x_count = 0;
        l.y_count = 0;
      }
    }

    /// Register a vertex for insertion/removal into the configuration
    lazy_det_operation_t &operator<<(vertex_t const &v) {
      // Separately lazy-add all four operators associated with this vertex
      lazy_op_lst[v.idx.b1].lazy_add_y(arg_t{v.tau1, v.idx.i1, v.has_alpha_shift, v.s}); //c^\dagger
      lazy_op_lst[v.idx.b2].lazy_add_x(arg_t{v.tau2, v.idx.i2, v.has_alpha_shift, v.s}); //c
      lazy_op_lst[v.idx.b3].lazy_add_y(arg_t{v.tau3, v.idx.i3, v.has_alpha_shift, v.s}); //c^\dagger
      lazy_op_lst[v.idx.b4].lazy_add_x(arg_t{v.tau4, v.idx.i4, v.has_alpha_shift, v.s}); //c
      return *this;
    }

    /// Tries to perform the insertions into all dets, returning the det ratio
    mc_weight_t execute_try_insert() {
      mc_weight_t det_ratio = 1.0;
      for (auto i : range(dets->size())) det_ratio *= lazy_op_lst[i].execute_try_insert(&(*dets)[i]);
      reset();
      return det_ratio;
    }

    /// Tries to perform the removal from all dets, returning the det ratio
    mc_weight_t execute_try_remove() {
      mc_weight_t det_ratio = 1.0;
      for (auto i : range(dets->size())) det_ratio *= lazy_op_lst[i].execute_try_remove(&(*dets)[i]);
      reset();
      return det_ratio;
    }
  };
}
