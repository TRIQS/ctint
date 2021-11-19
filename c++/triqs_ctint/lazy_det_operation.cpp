#include "./lazy_det_operation.hpp"

namespace triqs_ctint {

  // lower_bound in a vector of length n, accessed by vec(i), where vec is a lambda
  // Complexity log.
  // Very similar to std::lower_bound, but computes a position and takes a lambda, not iterators.
  // Implementation adapted from cppreference.com
  template <typename Vec, typename Value> size_t lower_bound(int count, Vec const &vec, Value const &value) {
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
    TRIQS_ASSERT(first >= 0);
    return first;
  }

  // Find index of first row in determinant matrix that has an element not less than c.
  size_t get_c_lower_bound(det_t const *d, c_t const &c) {
    return lower_bound(
       d->size(), [d](int n) { return d->get_x(n); }, c);
  }

  // Find index of first column in determinant matrix that has an element not less than cdag.
  size_t get_cdag_lower_bound(det_t const *d, cdag_t const &cdag) {
    return lower_bound(
       d->size(), [d](int n) { return d->get_y(n); }, cdag);
  }

  g_tau_scalar_t lazy_det_operation_t::one_block::execute_try_insert(det_t *d) {

    if (c_lst.size() != cdag_lst.size()) TRIQS_RUNTIME_ERROR << "Trying to insert unequal number of c and c_dag operators into block!";

    size_t const c_count = c_lst.size();

    // Prefactor to account for resorting of operators
    double prefactor = 1.0;

    if (c_count == 0) {
      // Trivial insert
      return 1.0;
    } else if (c_count == 1) {
      // Perform single insertion
      size_t const pos_c    = get_c_lower_bound(d, c_lst[0]);
      size_t const pos_cdag = get_cdag_lower_bound(d, cdag_lst[0]);
      if ((pos_c + pos_cdag) % 2) prefactor *= -1.0;
      return d->try_insert(pos_c, pos_cdag, c_lst[0], cdag_lst[0]) * prefactor;
    } else {
      // Sort, since operators have to be inserted in order
      auto const predicate = [&prefactor](auto const &lhs, auto const &rhs) {
        bool const result = std::less<>{}(lhs, rhs);
        if (!result) {
          // if not ordered, operators are swapped with minus sign
          prefactor *= -1.0;
        }
        return result;
      };
      std::stable_sort(c_lst.begin(), c_lst.end(), predicate);
      std::stable_sort(cdag_lst.begin(), cdag_lst.end(), predicate);

      // Calculate the insertion positions
      std::vector<size_t> pos_c(c_count);
      std::vector<size_t> pos_cdag(c_count);
      for (size_t i = 0; i < c_count; ++i) {
        pos_c[i]    = i + get_c_lower_bound(d, c_lst[i]); // Shift by i to take into account of the insertion of previous ones.
        pos_cdag[i] = i + get_cdag_lower_bound(d, cdag_lst[i]);
        if ((pos_c[i] + pos_cdag[i]) % 2) prefactor *= -1.0;
      }

      return d->try_insert_k(pos_c, pos_cdag, c_lst, cdag_lst) * prefactor;
    }
  }

  g_tau_scalar_t lazy_det_operation_t::one_block::execute_try_remove(det_t *d) {

    if (c_lst.size() != cdag_lst.size()) TRIQS_RUNTIME_ERROR << "Trying to remove unequal number of c and c_dag operators from block!";

    size_t const c_count = c_lst.size();

    // Prefactor to account for resorting of operators
    double prefactor = 1.0;

    if (c_count == 0) {
      // Trivial removal
      return 1.0;
    } else if (c_count == 1) {
      // Perform single removal
      size_t const pos_c    = get_c_lower_bound(d, c_lst[0]);
      size_t const pos_cdag = get_cdag_lower_bound(d, cdag_lst[0]);
      if ((pos_c + pos_cdag) % 2) prefactor *= -1.0;
      return d->try_remove(pos_c, pos_cdag) * prefactor;
    } else {
      // Sort, since operators have to be inserted in order
      auto const predicate = [&prefactor](auto const &lhs, auto const &rhs) {
        bool const result = std::less<>{}(lhs, rhs);
        if (!result) {
          // if not ordered, operators are swapped with minus sign
          prefactor *= -1.0;
        }
        return result;
      };
      std::stable_sort(c_lst.begin(), c_lst.end(), predicate);
      std::stable_sort(cdag_lst.begin(), cdag_lst.end(), predicate);

      // Calculate the removal positions
      std::vector<size_t> pos_c(c_count);
      std::vector<size_t> pos_cdag(c_count);
      for (size_t i = 0; i < c_count; ++i) {
        pos_c[i]    = get_c_lower_bound(d, c_lst[i]);
        pos_cdag[i] = get_cdag_lower_bound(d, cdag_lst[i]);
        if ((pos_c[i] + pos_cdag[i]) % 2) prefactor *= -1.0;
      }

      // Perform higher rank removal
      return d->try_remove_k(pos_c, pos_cdag) * prefactor;
    }
  }

} // namespace triqs_ctint
