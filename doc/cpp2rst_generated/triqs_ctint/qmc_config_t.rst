..
   Generated automatically by cpp2rst

.. highlight:: c


.. _qmc_config_t:

qmc_config_t
==============

Type of the Monte-Carlo configuration



**Synopsis**:

.. code-block:: c

    class qmc_config_t;


Public members
--------------

+------------+-----------------------+-------------------------------------------------------+
| Member     | Type                  | Comment                                               |
+============+=======================+=======================================================+
| vertex_lst | std::vector<vertex_t> | /// Unordered list of all vertices currently inserted |
+------------+-----------------------+-------------------------------------------------------+
| dets       | std::vector<det_t>    | /// List containing the determinant for each block    |
+------------+-----------------------+-------------------------------------------------------+


Member functions
----------------

+-------------------------------------------------------------+----------------------------+
| Member function                                             | Comment                    |
+=============================================================+============================+
| :ref:`perturbation_order <qmc_config_t_perturbation_order>` | Current perturbation order |
+-------------------------------------------------------------+----------------------------+
| :ref:`constructor <qmc_config_t_constructor>`               |                            |
+-------------------------------------------------------------+----------------------------+

.. toctree::
    :hidden:

    qmc_config_t/perturbation_order
    qmc_config_t/constructor