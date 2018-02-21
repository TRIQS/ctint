#include "./lazy_det_operation.hpp"

namespace triqs_ctint {

  // lower_bound in a vector of length n, accessed by vec(i), where vec is a lambda
  // Complexity log.
  // Very similar to std::lower_bound, but computes a position and takes a lambda, not iterators.
  // Implementation adapted from cppreference.com
  template <typename Vec, typename Value> int lower_bound(int count, Vec const &vec, Value const &value) {
    int p = 0, first = 0;
    while (count > 0) {
      p        = first;
      int step = count / 2;
      p += step;
      if (vec(p) < value) {
        first = ++p;
        count -= step + 1;
      } else
        count = step;
    }
    return first;
  };

  // Find index of first row in determinant matrix that has an element not less than c.
  int get_c_lower_bound(det_t const *d, c_t const &c) {
    return lower_bound(d->size(), [d](int n) { return d->get_x(n); }, c);
  }

  // Find index of first column in determinant matrix that has an element not less than cdag.
  int get_cdag_lower_bound(det_t const *d, cdag_t const &cdag) {
    return lower_bound(d->size(), [d](int n) { return d->get_y(n); }, cdag);
  }

  g_tau_scalar_t lazy_det_operation_t::one_block::execute_try_insert(det_t *d) {

    // Equal number of creation and annihilation operators required
    if (c_count != cdag_count) TRIQS_RUNTIME_ERROR << " ERROR: Trying to insert unequal number of c and c_dag operators into block! ";

    // Trivial insert
    if (c_count == 0) return 1.0;

    // Prefactor to account for resorting of operators
    double prefactor = 1.0;

    // Sort, since operators have to be inserted in order
    if (c_count == 2) {
      if (c_lst[1] < c_lst[0]) {
        std::swap(c_lst[0], c_lst[1]);
        prefactor *= -1.0;
      }
      if (cdag_lst[1] < cdag_lst[0]) {
        std::swap(cdag_lst[0], cdag_lst[1]);
        prefactor *= -1.0;
      }
    }

    // Calculate the insertion positions
    for (int i = 0; i < c_count; ++i) {
      pos_c[i]    = i + get_c_lower_bound(d, c_lst[i]); // Shift by i to take into account of the insertion of previous ones.
      pos_cdag[i] = i + get_cdag_lower_bound(d, cdag_lst[i]);
      if ((pos_c[i] + pos_cdag[i]) % 2 != 0) prefactor *= -1.0;
    }

    // Perform single or double insertion
    switch (c_count) {
      case (1): return d->try_insert(pos_c[0], pos_cdag[0], c_lst[0], cdag_lst[0]) * prefactor;
      case (2): return d->try_insert2(pos_c[0], pos_c[1], pos_cdag[0], pos_cdag[1], c_lst[0], c_lst[1], cdag_lst[0], cdag_lst[1]) * prefactor;
      default:
        TRIQS_RUNTIME_ERROR << "Not implemented";
        return 0; // avoid compiler warning
    }
  }

  g_tau_scalar_t lazy_det_operation_t::one_block::execute_try_remove(det_t *d) {

    // Equal number of creation and annihilation operators required
    if (c_count != cdag_count) TRIQS_RUNTIME_ERROR << "Internal Error";

    // Trivial removal
    if (c_count == 0) return 1.0;

    // Prefactor to account for resorting of operators
    double prefactor = 1.0;

    // Sort, since operators have to be inserted in order
    if (c_count == 2) {
      if (c_lst[1] < c_lst[0]) {
        std::swap(c_lst[0], c_lst[1]);
        prefactor *= -1.0;
      }
      if (cdag_lst[1] < cdag_lst[0]) {
        std::swap(cdag_lst[0], cdag_lst[1]);
        prefactor *= -1.0;
      }
    }

    // Calculate the removal positions
    for (int i = 0; i < c_count; ++i) {
      pos_c[i]    = get_c_lower_bound(d, c_lst[i]);
      pos_cdag[i] = get_cdag_lower_bound(d, cdag_lst[i]);
      if ((pos_c[i] + pos_cdag[i]) % 2 != 0) prefactor *= -1.0;
    }

    // Perform single or double removal
    switch (c_count) {
      case (1): return d->try_remove(pos_c[0], pos_cdag[0]) * prefactor;
      case (2): return d->try_remove2(pos_c[0], pos_c[1], pos_cdag[0], pos_cdag[1]) * prefactor;
      default:
        TRIQS_RUNTIME_ERROR << "Not implemented";
        return 0; // avoid compiler warning
    }
  }

} // namespace triqs_ctint
