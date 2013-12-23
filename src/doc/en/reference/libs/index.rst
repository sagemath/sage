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

   sage/libs/gap/libgap
   sage/libs/gap/element
   sage/libs/flint/fmpz_poly
   sage/libs/libecm
   sage/libs/lrcalc/lrcalc
   sage/libs/ntl/all
   sage/libs/mwrank/all
   sage/libs/mwrank/mwrank
   sage/libs/mwrank/interface
   sage/libs/pari/pari_instance
   sage/libs/pari/gen
   sage/libs/ppl
   sage/ext/pselect
   sage/libs/singular/function
   sage/libs/singular/option

.. include:: ../footer.txt
