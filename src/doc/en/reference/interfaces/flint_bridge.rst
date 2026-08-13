.. nodoctest
FLINT Bridge Interface
======================

This module provides functions to convert from SageMath's numerical types
to Python-FLINT's arbitrary precision types. It solves the issue of direct 
conversion failing with code like:

.. code-block:: python

    flint.arb(RR(1))
    flint.acb(CC(1+I))

The bridge uses string-based conversion through mpmath to preserve precision,
with fallback to direct conversion when necessary.

Real Number Conversion
---------------------

.. autofunction:: sage.interfaces.flint.bridge.sage_to_flint_arb

Complex Number Conversion
------------------------

.. autofunction:: sage.interfaces.flint.bridge.sage_to_flint_acb

Matrix Conversion
---------------

.. autofunction:: sage.interfaces.flint.bridge.sage_matrix_to_flint

Module Reference
--------------

.. automodule:: sage.interfaces.flint.bridge
   :members:
   :undoc-members:
   :show-inheritance: