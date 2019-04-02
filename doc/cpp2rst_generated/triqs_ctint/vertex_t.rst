..
   Generated automatically by cpp2rst

.. highlight:: c


.. _vertex_t:

vertex_t
==========

Type representing an interaction vertex of the microscopic model at hand.

 Can be inserted in the Monte-Carlo move.

**Synopsis**:

.. code-block:: c

    class vertex_t;
     Can be inserted in the Monte-Carlo move.

Public members
--------------

+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| Member            | Type                      | Comment                                                                                          |
+===================+===========================+==================================================================================================+
| idx               | triqs_ctint::vertex_idx_t | /// Object containing discrete quantum numbers for external legs, i.e. block and non-block index |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| tau1              | triqs_ctint::tau_t        | /// Imaginary times for all four external legs (c^\dagger, c, c^\dagger, c)                      |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| tau2              | triqs_ctint::tau_t        | /// Imaginary times for all four external legs (c^\dagger, c, c^\dagger, c)                      |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| tau3              | triqs_ctint::tau_t        | /// Imaginary times for all four external legs (c^\dagger, c, c^\dagger, c)                      |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| tau4              | triqs_ctint::tau_t        | /// Imaginary times for all four external legs (c^\dagger, c, c^\dagger, c)                      |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| amplitude         | triqs_ctint::U_scalar_t   | /// Amplitude of the vertex, i.e. U, U(tau1-tau2), etc...                                        |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| proposition_proba | double                    | /// Probability of proposition for this vertex                                                   |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| has_alpha_shift   | bool                      | /// Switch for alpha shift                                                                       |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+
| s                 | int                       | /// Value of auxiliary spin                                                                      |
+-------------------+---------------------------+--------------------------------------------------------------------------------------------------+