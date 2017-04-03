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

  // Find index of first row in determinant matrix that has an arg_t not less than x.
  int get_x_lower_bound(det_t const *d, arg_t const &x) {
    return lower_bound(d->size(), [d](int n) { return d->get_x(n); }, x);
  }

  // Find index of first column in determinant matrix that has an arg_t not less than y.
  int get_y_lower_bound(det_t const *d, arg_t const &y) {
    return lower_bound(d->size(), [d](int n) { return d->get_y(n); }, y);
  }

  double lazy_det_operation_t::one_block::execute_try_insert(det_t *d) {

    // Equal number of creation and annihilation operators required
    if (x_count != y_count) TRIQS_RUNTIME_ERROR << "Internal Error";

    // Trivial insert
    if (x_count == 0) return 1.0;

    // Sort, since operators have to be inserted in order
    std::sort(x.begin(), x.begin() + x_count);
    std::sort(y.begin(), y.begin() + x_count);

    // Calculate the insertion positions
    for (int i = 0; i < x_count; ++i) {
      px[i] = i + get_x_lower_bound(d, x[i]); // Shift by i to take into account of the insertion of previous ones.
      py[i] = i + get_y_lower_bound(d, y[i]);
    }

    // Perform single or double insertion
    switch (x_count) {
      case (1): return d->try_insert(px[0], py[0], x[0], y[0]) * (((px[0] + py[0]) % 2) == 0 ? 1 : -1);
      case (2): return d->try_insert2(px[0], px[1], py[0], py[1], x[0], x[1], y[0], y[1]) * (((px[0] + px[1] + py[0] + py[1]) % 2) == 0 ? 1 : -1);
      default:
        TRIQS_RUNTIME_ERROR << "Not implemented";
        return 0; // avoid compiler warning
    }
  }

  double lazy_det_operation_t::one_block::execute_try_remove(det_t *d) {

    // Equal number of creation and annihilation operators required
    if (x_count != y_count) TRIQS_RUNTIME_ERROR << "Internal Error";

    // Trivial removal
    if (x_count == 0) return 1.0;

    // Calculate the removal positions
    for (int i = 0; i < x_count; ++i) {
      px[i] = get_x_lower_bound(d, x[i]);
      py[i] = get_y_lower_bound(d, y[i]);
    }

    // Perform single or double removal
    switch (x_count) {
      case (1): return d->try_remove(px[0], py[0]) * (((px[0] + py[0]) % 2) == 0 ? 1 : -1);
      case (2): return d->try_remove2(px[0], px[1], py[0], py[1]) * (((px[0] + px[1] + py[0] + py[1]) % 2) == 0 ? 1 : -1);
      default:
        TRIQS_RUNTIME_ERROR << "Not implemented";
        return 0; // avoid compiler warning
    }
  }

} // namespace triqs_ctint
