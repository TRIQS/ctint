#pragma once
#include "vertex.hpp"
#include "dets.hpp"

namespace triqs_ctint {

  /**
   * A stack of lazy operations in the determinants.
   * The class keeps two lists of the c and c^\dagger operators to be inserted.
   * Additional verticies can be added with the << operator.
   * The order of insertion is of no importance, while the insertion in the
   * determinant is ordered using the binary comparison of c_t/cdag_t.
   * A call to execute_try_insert/execute_try_remove performs the full
   * insert/remove operation in all the determinants involved.
   */
  class lazy_det_operation_t {

    private:
    // Vector of all determinants
    std::vector<det_t> *dets;

    // This type keeps track of lazy inserts for single determinant
    struct one_block {
      // Vector of $c$ (x) and $c^\dagger$ (y) operators
      std::vector<c_t> c_lst = {};
      std::vector<cdag_t> cdag_lst = {};

      // Lazy-add a $c$ operator
      void lazy_add_c(c_t const &c_) {
        c_lst.emplace_back(c_);
      }

      // Lazy-add a $c^\dagger$ operator
      void lazy_add_cdag(cdag_t const &cdag_) {
        cdag_lst.emplace_back(cdag_);
      }

      // Tries to perform all insertions into the det, returning the det ratio
      g_tau_scalar_t execute_try_insert(det_t *d);

      // Tries to perform all removals from the det, returning the det ratio
      g_tau_scalar_t execute_try_remove(det_t *d);

      // Tries to perform all spinflips from the det, returning the det ratio
      g_tau_scalar_t execute_try_change_col_row(det_t *d);
    };

    // The list of lazy operations, one for each determinant
    std::vector<one_block> lazy_op_lst;

    public:
    lazy_det_operation_t(std::vector<det_t> *dets) : dets(dets), lazy_op_lst(dets->size()) {}

    /// Clean all registered moves
    void reset() {
      for (auto &l : lazy_op_lst) {
        l.c_lst.clear();
        l.cdag_lst.clear();
      }
    }

    /// Register a vertex for insertion/removal into the configuration
    lazy_det_operation_t &operator<<(vertex_t const &v) {
      // Separately lazy-add all four operators associated with this vertex
      lazy_op_lst[v.idx.b1].lazy_add_cdag(cdag_t{v.tau1, v.idx.u1, v.vertex_label, 0, v.s});  //c^\dagger
      lazy_op_lst[v.idx.b2].lazy_add_c(c_t{v.tau2, v.idx.u2, v.vertex_label, 0, v.s});        //c
      lazy_op_lst[v.idx.b3].lazy_add_cdag(cdag_t{v.tau3, v.idx.u3, v.vertex_label, 1, v.s}); //c^\dagger
      lazy_op_lst[v.idx.b4].lazy_add_c(c_t{v.tau4, v.idx.u4, v.vertex_label, 1, v.s});       //c
      return *this;
    }

    /// Tries to perform the insertions into all dets, returning the det ratio
    g_tau_scalar_t execute_try_insert() {
      g_tau_scalar_t det_ratio = 1.0;
      for (auto i : range(dets->size())) det_ratio *= lazy_op_lst[i].execute_try_insert(&(*dets)[i]);
      reset();
      return det_ratio;
    }

    /// Tries to perform the removal from all dets, returning the det ratio
    g_tau_scalar_t execute_try_remove() {
      g_tau_scalar_t det_ratio = 1.0;
      for (auto i : range(dets->size())) det_ratio *= lazy_op_lst[i].execute_try_remove(&(*dets)[i]);
      reset();
      return det_ratio;
    }

    /// Tries to perform the removal from all dets, returning the det ratio
    g_tau_scalar_t execute_try_change_col_row() {
      g_tau_scalar_t det_ratio = 1.0;
      for (auto i : range(dets->size())) det_ratio *= lazy_op_lst[i].execute_try_change_col_row(&(*dets)[i]);
      reset();
      return det_ratio;
    }
  };
} // namespace triqs_ctint
