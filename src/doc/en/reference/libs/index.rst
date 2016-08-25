C/C++ Library Interfaces
========================

An underlying philosophy in the development of Sage is that it
should provide unified library-level access to the some of the best
GPL'd C/C++ libraries. Currently Sage provides some access to
MWRANK, NTL, PARI, and Hanke, each of which are included with
Sage.

The interfaces are implemented via shared libraries and data is
moved between systems purely in memory. In particular, there is no
interprocess interpreter parsing (e.g., ``expect``),
since everything is linked together and run as a single process.
This is much more robust and efficient than using
``expect``.

Each of these interfaces is used by other parts of Sage. For
example, mwrank is used by the elliptic curves module to compute
ranks of elliptic curves, and PARI is used for computation of class
groups. It is thus probably not necessary for a casual user of Sage
to be aware of the modules described in this chapter.

.. toctree::
   :maxdepth: 2

   sage/libs/eclib/interface
   sage/libs/eclib/mwrank
   sage/libs/eclib/mat
   sage/libs/eclib/newforms
   sage/libs/eclib/homspace
   sage/libs/eclib/constructor
   sage/libs/gmp/rational_reconstruction
   sage/libs/lcalc/lcalc_Lfunction
   sage/libs/ratpoints
   sage/libs/singular/function
   sage/libs/singular/function_factory
   sage/libs/singular/singular
   sage/libs/singular/polynomial
   sage/libs/singular/option
   sage/libs/singular/ring
   sage/libs/singular/groebner_strategy
   sage/libs/ppl
   sage/libs/linbox/linbox
   sage/libs/flint/flint
   sage/libs/flint/fmpz_poly
   sage/libs/flint/arith
   sage/libs/symmetrica/symmetrica
   sage/libs/mpmath/utils
   sage/libs/ntl/all
   sage/libs/libecm
   sage/libs/lrcalc/lrcalc
   sage/libs/pari/gen
   sage/libs/pari/pari_instance
   sage/libs/pari/closure
   sage/libs/pari/convert
   sage/libs/pari/handle_error
   sage/rings/pari_ring
   sage/libs/fplll/fplll
   sage/libs/readline
   sage/libs/gap/context_managers
   sage/libs/gap/gap_functions
   sage/libs/gap/test_long
   sage/libs/gap/util
   sage/libs/gap/libgap
   sage/libs/gap/test
   sage/libs/gap/element
   sage/libs/gap/saved_workspace
   sage/libs/ecl

   sage/gsl/gsl_array

   sage/ext/pselect

.. Cannot be imported independently of mpmath: sage/libs/mpmath/ext_main sage/libs/mpmath/ext_impl sage/libs/mpmath/ext_libmp

.. Modules depending on optional packages: sage/libs/coxeter3/coxeter sage/libs/coxeter3/coxeter_group sage/libs/fes sage/libs/homfly sage/libs/braiding

.. include:: ../footer.txt
