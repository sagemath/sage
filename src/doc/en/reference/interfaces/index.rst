Interpreter Interfaces
======================

Sage provides a unified interface to the best computational
software. This is accomplished using both C-libraries (see
`C/C++ Library Interfaces <../libs/index.html>`_)
and interpreter interfaces, which are
implemented using pseudo-tty's, system files, etc. This chapter is
about these interpreter interfaces.

.. note::

    Each interface requires that the corresponding software is
    installed on your computer. Sage includes GAP, PARI, Singular, and
    Maxima, but does not include Octave (very easy to install), MAGMA
    (non-free), Maple (non-free), or Mathematica (non-free).


    There is overhead associated with each call to one of these
    systems. For example, computing ``2+2`` thousands of times using
    the GAP interface will be slower than doing it directly in
    Sage. In contrast, the C-library interfaces of
    `C/C++ Library Interfaces <../libs/index.html>`_
    incur less overhead.


In addition to the commands described for each of the interfaces
below, you can also type e.g., ``%gap``,
``%magma``, etc., to directly interact with a given
interface in its state. Alternatively, if ``X`` is an
interface object, typing ``X.interact()`` allows you to
interact with it. This is completely different than
``X.console()`` which starts a complete new copy of
whatever program ``X`` interacts with. Note that the
input for ``X.interact()`` is handled by Sage, so the
history buffer is the same as for Sage, tab completion is as for
Sage (unfortunately!), and input that spans multiple lines must be
indicated using a backslash at the end of each line. You can pull
data into an interactive session with ``X`` using
``sage(expression)``.

The console and interact methods of an interface do very different
things. For example, using gap as an example:


#. ``gap.console()``: You are completely using another
   program, e.g., gap/magma/gp Here Sage is serving as nothing more
   than a convenient program launcher, similar to bash.

#. ``gap.interact()``: This is a convenient way to
   interact with a running gap instance that may be "full of" Sage
   objects. You can import Sage objects into this gap (even from the
   interactive interface), etc.


The console function is very useful on occasion, since you get the
exact actual program available (especially useful for tab completion
and testing to make sure nothing funny is going on).

.. toctree::
   :maxdepth: 1

   sage/interfaces/interface
   sage/interfaces/expect
   sage/interfaces/sagespawn
   sage/interfaces/abc
   sage/interfaces/axiom
   sage/interfaces/asir
   sage/interfaces/ecm
   sage/interfaces/four_ti_2
   sage/interfaces/fricas
   sage/interfaces/frobby
   sage/interfaces/gap
   sage/interfaces/gap3
   sage/interfaces/gfan
   sage/interfaces/giac
   sage/interfaces/gnuplot
   sage/interfaces/gp
   sage/interfaces/jmoldata
   sage/interfaces/kash
   sage/interfaces/kenzo
   sage/interfaces/latte
   sage/interfaces/lie
   sage/interfaces/lisp
   sage/interfaces/macaulay2
   sage/interfaces/magma
   sage/interfaces/magma_free
   sage/interfaces/maple
   sage/interfaces/mathematica
   sage/interfaces/mathics
   sage/interfaces/matlab
   sage/interfaces/maxima
   sage/interfaces/maxima_abstract
   sage/interfaces/maxima_lib
   sage/interfaces/mupad
   sage/interfaces/mwrank
   sage/interfaces/octave
   sage/interfaces/phc
   sage/interfaces/polymake
   sage/interfaces/povray
   sage/interfaces/psage
   sage/interfaces/qepcad
   sage/interfaces/r
   sage/interfaces/rubik
   sage/interfaces/sage0
   sage/interfaces/scilab
   sage/interfaces/singular
   sage/interfaces/sympy
   sage/interfaces/tachyon
   sage/interfaces/tides

   sage/interfaces/cleaner
   sage/interfaces/quit
   sage/interfaces/read_data

.. include:: ../footer.txt
