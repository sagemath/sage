Numerical Optimization
======================

.. toctree::
   :maxdepth: 1

   sage/numerical/knapsack
   sage/numerical/mip
   sage/numerical/sdp
   sage/numerical/linear_functions
   sage/numerical/linear_tensor
   sage/numerical/linear_tensor_element
   sage/numerical/linear_tensor_constraints
   sage/numerical/optimize
   sage/numerical/interactive_simplex_method
   sage/numerical/gauss_legendre

Linear Optimization (LP) and Mixed Integer Linear Optimization (MIP) Solver backends
------------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 1

   sage/numerical/backends/generic_backend
   sage/numerical/backends/interactivelp_backend
   sage/numerical/backends/glpk_backend
   sage/numerical/backends/glpk_exact_backend
   sage/numerical/backends/glpk_graph_backend
   sage/numerical/backends/ppl_backend
   sage/numerical/backends/cvxopt_backend

.. backends depending on optional packages
..   sage/numerical/backends/coin_backend
..   sage/numerical/backends/cplex_backend
..   sage/numerical/backends/gurobi_backend

Sage also supports, via optional packages, CBC (COIN-OR), CPLEX (ILOG), and Gurobi. In order to find out
how to use them in Sage, please refer to the `Thematic Tutorial on Linear
Programming
<http://doc.sagemath.org/html/en/thematic_tutorials/linear_programming.html>`_.

The following backend is used for debugging and testing purposes.

.. toctree::
   :maxdepth: 1

   sage/numerical/backends/logging_backend

Semidefinite Optimization (SDP) Solver backends
-----------------------------------------------

.. toctree::
   :maxdepth: 1

   sage/numerical/backends/generic_sdp_backend
   sage/numerical/backends/cvxopt_sdp_backend

For more details on CVXOPT, see `CVXOPT documentation <http://cvxopt.org/documentation/index.html>`_.

.. include:: ../footer.txt
