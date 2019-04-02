..
   Generated automatically by cpp2rst

.. highlight:: c


.. _arg_t:

arg_t
=======

Type of row and column argument of the Green function matrix inside the determinant.

 G(x,y) is evaluated at the time-difference x.tau - y.tau.

**Synopsis**:

.. code-block:: c

    template<bool dag> class arg_t;
     G(x,y) is evaluated at the time-difference x.tau - y.tau.

Public members
--------------

+------------------+--------------------+-----------------------------------------------------------------------------------------------------+
| Member           | Type               | Comment                                                                                             |
+==================+====================+=====================================================================================================+
| tau              | triqs_ctint::tau_t | /// The imaginary time                                                                              |
+------------------+--------------------+-----------------------------------------------------------------------------------------------------+
| u                | int                | /// The orbital (or non-block) index                                                                |
+------------------+--------------------+-----------------------------------------------------------------------------------------------------+
| with_alpha_shift | bool               | /// Switch for alpha shift along the diagonal of the matrix. Relevant for density-type interactions |
+------------------+--------------------+-----------------------------------------------------------------------------------------------------+
| s                | int                | /// The auxiliary spin                                                                              |
+------------------+--------------------+-----------------------------------------------------------------------------------------------------+


Member functions
----------------

+--------------------------------------------------------+-------------------------------------------------------------------------------------------------+
| Member function                                        | Comment                                                                                         |
+========================================================+=================================================================================================+
| :ref:`operator comparison <arg_t_operator comparison>` | Lexicographical sorting of arg_t. This determines the order of row and columns inside the dets. |
+--------------------------------------------------------+-------------------------------------------------------------------------------------------------+

.. toctree::
    :hidden:

    arg_t/operator comparison