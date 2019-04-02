..
   Generated automatically by cpp2rst

.. highlight:: c


.. _G0hat_t:

G0hat_t
=========

Functor that evaluates the matrix elements, used by the det_manip.

 Cares for the possible alpha-shift along the diagonal of the matrix.

**Synopsis**:

.. code-block:: c

    class G0hat_t;
     Cares for the possible alpha-shift along the diagonal of the matrix.

Public members
--------------

+--------------+------------------------------------------------------+--------------------------------------------------+
| Member       | Type                                                 | Comment                                          |
+==============+======================================================+==================================================+
| G0_shift_tau | gf_const_view<triqs::gfs::imtime, g_tau_t::target_t> | /// The (shifted) non-interacting Green function |
+--------------+------------------------------------------------------+--------------------------------------------------+
| alpha        | array<double, 2>                                     | /// The alpha function                           |
+--------------+------------------------------------------------------+--------------------------------------------------+


Member functions
----------------

+----------------------------------------+---------+
| Member function                        | Comment |
+========================================+=========+
| :ref:`operator() <G0hat_t_operator()>` |         |
+----------------------------------------+---------+

.. toctree::
    :hidden:

    G0hat_t/operator()