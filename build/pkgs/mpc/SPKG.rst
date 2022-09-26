mpc: Arithmetic of complex numbers with arbitrarily high precision and correct rounding
=======================================================================================

Description
-----------

From https://www.multiprecision.org/mpc: GNU MPC is a C library for the
arithmetic of complex numbers with arbitrarily high precision and
correct rounding of the result. It extends the principles of the
IEEE-754 standard for fixed precision real floating point numbers to
complex numbers, providing well-defined semantics for every operation.
At the same time, speed of operation at high precision is a major design
goal.

License
-------

LGPLv3+ for the code and GFDLv1.3+ (with no invariant sections) for the
documentation.


Upstream Contact
----------------

The MPC website is located at https://www.multiprecision.org/mpc .

The MPC team can be contacted via the MPC mailing list: mpc-discuss@inria.fr

Special Update/Build Instructions
---------------------------------

-  mpc_mul_faster.patch: Patch from Paul Zimmermann to speed up MPC
   multiplication (for small precisions) by reducing overhead in MPFR
   operations.
