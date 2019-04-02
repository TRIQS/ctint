..
   Generated automatically by cpp2rst

.. highlight:: c


.. _vertex_idx_t:

vertex_idx_t
==============

Type representing the set of discrete quantum numbers for the vertices of

 the microscopic model at hand. We distinguish between block-indeces
 (in which the bare Green function is diagonal) and non-block indeces.

**Synopsis**:

.. code-block:: c

    class vertex_idx_t;
     the microscopic model at hand. We distinguish between block-indeces
 (in which the bare Green function is diagonal) and non-block indeces.

Public members
--------------

+--------+------+--------------------------------------------------------------------------------+
| Member | Type | Comment                                                                        |
+========+======+================================================================================+
| b1     | int  | /// First operator of the vertex (outgoing, c^\dagger): block index, non-block |
+--------+------+--------------------------------------------------------------------------------+
| u1     | int  | /// First operator of the vertex (outgoing, c^\dagger): block index, non-block |
+--------+------+--------------------------------------------------------------------------------+
| b2     | int  | /// Second operator of the vertex (ingoing, c): block index, non-block         |
+--------+------+--------------------------------------------------------------------------------+
| u2     | int  | /// Second operator of the vertex (ingoing, c): block index, non-block         |
+--------+------+--------------------------------------------------------------------------------+
| b3     | int  | /// Third operator of the vertex (outgoing, c^\dagger): block index, non-block |
+--------+------+--------------------------------------------------------------------------------+
| u3     | int  | /// Third operator of the vertex (outgoing, c^\dagger): block index, non-block |
+--------+------+--------------------------------------------------------------------------------+
| b4     | int  | /// Fourth operator of the vertex (ingoing, c): block index, non-block         |
+--------+------+--------------------------------------------------------------------------------+
| u4     | int  | /// Fourth operator of the vertex (ingoing, c): block index, non-block         |
+--------+------+--------------------------------------------------------------------------------+