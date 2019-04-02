..
   Generated automatically by cpp2rst

.. highlight:: c


.. _tau_t:

tau_t
=======

A point in imaginary time, i.e. $\tau \in [0,\beta]$, but defined on a very fine grid.

 The position in the segment is given by an uint32_t, i.e. a very long integer.
 This allows exact comparisons, which notoriously dangerous on floating point number.

**Synopsis**:

.. code-block:: c

    class tau_t;
     The position in the segment is given by an uint32_t, i.e. a very long integer.
 This allows exact comparisons, which notoriously dangerous on floating point number.

Public members
--------------

+--------+----------+-----------------------------------------------------------------------+
| Member | Type     | Comment                                                               |
+========+==========+=======================================================================+
| n      | uint32_t | /// :math:`\tau` value, represented as an integer on a very fine grid |
+--------+----------+-----------------------------------------------------------------------+


Member functions
----------------

+----------------------------------------------------------+--------------------------------------------+
| Member function                                          | Comment                                    |
+==========================================================+============================================+
| :ref:`get_random <tau_t_get_random>`                     | Get a random point in :math:`[0,\beta[`    |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`operator== != <tau_t_operator== !=>`               |                                            |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`operator comparison <tau_t_operator comparison>`   |                                            |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`get_zero <tau_t_get_zero>`                         | Return \tau = 0                            |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`get_zero_plus <tau_t_get_zero_plus>`               | Return \tau = 0^{+} = 0 + \delta           |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`get_zero_plus_plus <tau_t_get_zero_plus_plus>`     | Return \tau = 0^{++} = 0 + 2*\delta        |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`get_beta <tau_t_get_beta>`                         | Return \tau = \beta                        |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`get_beta_minus <tau_t_get_beta_minus>`             | Return \tau = beta^{-} = \beta - \delta    |
+----------------------------------------------------------+--------------------------------------------+
| :ref:`get_beta_minus_minus <tau_t_get_beta_minus_minus>` | Return \tau = beta^{--} = \beta - 2*\delta |
+----------------------------------------------------------+--------------------------------------------+

.. toctree::
    :hidden:

    tau_t/get_random
    tau_t/operator== !=
    tau_t/operator comparison
    tau_t/get_zero
    tau_t/get_zero_plus
    tau_t/get_zero_plus_plus
    tau_t/get_beta
    tau_t/get_beta_minus
    tau_t/get_beta_minus_minus


Non Member functions
--------------------

+----------------------------------------+------------------------------------------+
| Non member function                    | Comment                                  |
+========================================+==========================================+
| :ref:`operator<< <tau_t_operator\<\<>` | Operator allowing output to std::ostream |
+----------------------------------------+------------------------------------------+

.. toctree::
    :hidden:

    tau_t/operator<<