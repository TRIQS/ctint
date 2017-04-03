Reference documentation
==========================

The main class is **solver_core**. It constructs Monte-Carlo solver which owns a set of **moves** and **measures**. These in turn act on a **configuration**. The main data for the QMC are stored in **qmc_data**.


.. toctree::
   :maxdepth: 2

   reference/triqs/applications/impurity_solvers/ctint/solver_core
   moves
   measures
   reference/triqs/applications/impurity_solvers/ctint/configuration
   reference/triqs/applications/impurity_solvers/ctint/qmc_data


Python interface
~~~~~~~~~~~~~~~~~


.. autoclass:: ctint.SolverCore
  :members:
