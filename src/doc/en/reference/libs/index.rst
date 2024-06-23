C/C++ Library Interfaces
========================

An underlying philosophy in the development of Sage is that it
should provide unified library-level access to the some of the best
GPL'd C/C++ libraries. Sage provides access to many libraries which
are included with Sage.

The interfaces are implemented via shared libraries and data is
moved between systems purely in memory. In particular, there is no
interprocess interpreter parsing (e.g., ``pexpect``),
since everything is linked together and run as a single process.
This is much more robust and efficient than using ``pexpect``.

Each of these interfaces is used by other parts of Sage. For
example, eclib is used by the elliptic curves module to compute
ranks of elliptic curves and PARI is used for computation of class
groups. It is thus probably not necessary for a casual user of Sage
to be aware of the modules described in this chapter.

ECL
---
.. toctree::
   :maxdepth: 1

   sage/libs/ecl

eclib
-----
.. toctree::
   :maxdepth: 1

   sage/libs/eclib/interface
   sage/libs/eclib/mwrank
   sage/libs/eclib/mat
   sage/libs/eclib/newforms
   sage/libs/eclib/homspace
   sage/libs/eclib/constructor

FLINT
-----
.. toctree::
   :maxdepth: 1

   sage/libs/flint/fmpz_poly_sage
   sage/libs/flint/fmpq_poly_sage
   sage/libs/flint/arith_sage
   sage/libs/flint/qsieve_sage
   sage/libs/flint/ulong_extras_sage

Giac
----
.. toctree::
   :maxdepth: 1

   sage/libs/giac

GMP-ECM
-------
.. toctree::
   :maxdepth: 1

   sage/libs/libecm

GSL
---
.. toctree::
   :maxdepth: 1

   sage/libs/gsl/array

lcalc
-----
.. toctree::
   :maxdepth: 1

   sage/libs/lcalc/lcalc_Lfunction

libSingular
-----------
.. toctree::
   :maxdepth: 1

   sage/libs/singular/function
   sage/libs/singular/function_factory
   sage/libs/singular/singular
   sage/libs/singular/polynomial
   sage/libs/singular/option
   sage/libs/singular/ring
   sage/libs/singular/groebner_strategy

GAP
---
.. toctree::
   :maxdepth: 1

   sage/libs/gap/context_managers
   sage/libs/gap/gap_functions
   sage/libs/gap/test_long
   sage/libs/gap/util
   sage/libs/gap/libgap
   sage/libs/gap/test
   sage/libs/gap/element
   sage/libs/gap/saved_workspace

LinBox
------
.. toctree::
   :maxdepth: 1

   sage/libs/linbox/linbox_flint_interface

lrcalc
------
.. toctree::
   :maxdepth: 1

   sage/libs/lrcalc/lrcalc

mpmath
------
.. toctree::
   :maxdepth: 1

   sage/libs/mpmath/sage_utils

NTL
---
.. toctree::
   :maxdepth: 1

   sage/libs/ntl/all

PARI
----
.. toctree::
   :maxdepth: 1

   sage/libs/pari
   sage/libs/pari/convert_sage
   sage/rings/pari_ring

Symmetrica
----------
.. toctree::
   :maxdepth: 1

   sage/libs/symmetrica/symmetrica

.. Cannot be imported independently of mpmath: sage/libs/mpmath/ext_main sage/libs/mpmath/ext_impl sage/libs/mpmath/ext_libmp

.. Modules depending on optional packages: sage/libs/coxeter3/coxeter sage/libs/coxeter3/coxeter_group  sage/libs/homfly sage/libs/braiding

.. include:: ../footer.txt
