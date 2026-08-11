# sage.doctest: optional snappy
r"""
Interface to SnapPy

`SnapPy <https://snappy.computop.org/>`__ is an open source software for
low-dimensional topology. From the home-page:

.. NOTE::

    SnapPy is a program for studying the topology and geometry of 3-manifolds,
    with a focus on hyperbolic structures. It runs on Mac OS X, Linux, and
    Windows, and combines a link editor and 3D-graphics for Dirichlet domains
    and cusp neighborhoods with a powerful command-line interface based on the
    Python programming language. You can see it in action, learn how to
    install it, and watch the tutorial.

The SnapPy interface will only work if the optional Sage package SnapPy
is installed. The interface lets you send certain Sage objects to SnapPy,
run SnapPy functions, import certain SnapPy expressions to Sage,
or any combination of the above.

To send a Sage object ``sobj`` to SnapPy, call ``snappy(sobj)``.
This exports the Sage object to SnapPy and returns a new Sage object
wrapping the SnapPy expression/variable, so that you can use the
SnapPy variable from within Sage. You can then call SnapPy
functions on the new object; for example::

    sage: A = AbelianGroup([5,15,0,0]); A
    Multiplicative Abelian group isomorphic to C5 x C15 x Z x Z
    sage: As = snappy(A); As
    Z/5 + Z/15 + Z + Z
    sage: As.order()
    'infinite'

In the above example the order of the group is obtained using SnapPy's
``order`` method.

To see SnapPy's output you can simply print the SnapPy wrapper
object. However if you want to import SnapPy's output back to Sage,
call the SnapPy wrapper object's ``sage()`` method. This method returns
a native Sage object::

    sage: K = Knots().from_table(8, 21); K
    Knot represented by 8 crossings
    sage: Ks = snappy(K); Ks
    <Link: 1 comp; 8 cross>
    sage: Ks.goeritz_matrix()
    [-2  1  0]
    [ 1 -4  2]
    [ 0  2  1]
    sage: Ks.sage() == K
    True

If you want to run a SnapPy function and don't already have the input
in the form of a Sage object, then it might be simpler to input a string
``expr`` to ``snappy(expr)``. This string will be evaluated as if you had
typed it into SnapPy::

    sage: M1 = snappy("Manifold('m125')"); M1
    m125(0,0)(0,0)

Alternatively, all constructors of SnapPy classes can be used directly as attributes
of the interface::

    sage: M2 = snappy.Manifold('m125'); M2
    m125(0,0)(0,0)
    sage: M1 == M2
    True

Finally, if you just want to use a SnapPy command line from within
Sage, the IPython magic function ``%snappy`` dumps you into an interactive
command-line SnapPy session. As long as you work in this environment the
prompt is ``snappy:``. To finish the environment type ``CTRL+D``::

    sage: %snappy                                 # not tested

    --> Switching to SnapPy <--

    snappy: M = Manifold('9_42')
    None
    snappy: M.volume()
    4.05686022423682
    snappy: M.cusp_info('shape')
    [-4.27893631592295 + 1.95728679749950*I]

    --> Exiting back to Sage <--

    sage:                                         # not tested


Complicated translations
------------------------

The ``sobj.sage()`` method tries to convert a SnapPy object to a Sage
object. In many cases, it will just work. In particular, it should be able to
convert expressions entirely consisting of:

- numbers, i.e. integers, floats, complex numbers;
- functions and named constants also present in Sage, where:

    - Sage knows how to translate the function or constant's name from
      SnapPy's, or
    - the Sage name for the function or constant is trivially related to
      SnapPy's;

- symbolic variables whose names don't pathologically overlap with
  objects already defined in Sage.

This method will not work when SnapPy's output includes:

- strings;
- functions unknown to Sage;
- SnapPy functions with different parameters/parameter order to
  the Sage equivalent.


AUTHORS:

- Sebastian Oehms (2026): first version.
"""

##############################################################################
#       Copyright (C) 2026 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.interfaces.python_internal import (
    PythonInternalElement,
    PythonInternalInterface,
)
from sage.misc.instancedoc import instancedoc


class SnapPy(PythonInternalInterface):
    r"""
    Interface to the SnapPy interpreter.

    EXAMPLES::

        sage: K = Knots().from_table(8, 21)
        sage: Ks = snappyhp(K); Ks
        <Link: 1 comp; 8 cross>
        sage: M = Ks.exterior()
        sage: success, rho = M.verify_hyperbolicity(); success
        True

    More examples can be found in the module header.
    """
    def __init__(self, high_precision=False):
        r"""
        Python constructor.

        TESTS::

            sage: TestSuite(snappy).run(skip=['_test_pickling', '_test_category'])
        """
        super().__init__('snappy')
        self._high_precision = high_precision

    def _repr_(self):
        r"""
        Return representation string.

        EXAMPLES::

            sage: snappy
            SnapPy
        """
        return 'SnapPy'

    def _start(self):
        """
        Start up the SnapPy interpreter and sets the initial prompt and options.

        This is called the first time the SnapPy interface is actually used.

        EXAMPLES::

            sage: snappy._start()
            sage: type(snappy._namespace.AbelianGroup())
            <class 'SnapPy.AbelianGroup'>
        """
        if not self._interface_globals:
            from sage.features.interfaces import SnapPy
            SnapPy().module.require()
            import snappy
            import spherogram
            if self._high_precision:
                d = snappy.SnapPyHP.__dict__
            else:
                d = snappy.SnapPy.__dict__
            e = spherogram.__dict__
            D = dict(d.items())
            D.update(e)
            self._namespace = type('snappy_names', (object,), D)
            # set up the globals
            D[self._namespace.__name__] = self._namespace
            D[self.name()] = snappy
            self._interface_globals = D
            self._interface_modules = [self.name(), 'spherogram', 'SnapPy', 'SnapPyHP']

    def _install_hints(self):
        """
        Hints for installing snappy on your computer.

        EXAMPLES::

            sage: len(snappy._install_hints())
            99
        """
        return """
In order to use the SnapPy interface you need to have the
optional Sage package SnapPy installed.
"""

    def _object_class(self):
        r"""
        Return the element class of this parent.
        This is used in the interface class.

        EXAMPLES::

            sage: snappy._object_class()
            <class 'sage.interfaces.snappy.SnapPyElement'>
        """
        return SnapPyElement


@instancedoc
class SnapPyElement(PythonInternalElement):
    r"""
    Element class of the SnapPy interface.

    Its instances are usually constructed via the instance call of its parent.
    It wrapes the SnapPy library for this object. In a session SnapPy methods
    can be obtained using tab completion.

    EXAMPLES::

        sage: K = Knots().from_table(8, 21)
        sage: Ks = snappy(K)
        sage: Ms = Ks.exterior(); Ms
        unnamed link(0,0)
        sage: type(Ms)
        <class 'sage.interfaces.snappy.SnapPyElement'>
        sage: Gs = Ms.fundamental_group(); Gs
        Generators:
           a,b
        Relators:
           aabABBAbABabbaBabbaBAbABBAb
        sage: type(Gs)
        <class 'sage.interfaces.snappy.SnapPyElement'>
        sage: G = Gs.sage(); G
        Finitely presented group < a, b | a^2*b*a^-1*b^-2*a^-1*b*a^-1*(b^-1*a*b^2*a)^2*b^-1*a^-1*b*a^-1*b^-2*a^-1*b >
        sage: Fs = Ks.faces(); Fs
        [[<CS 7, 3>, <CS 6, 3>, <CS 5, 1>, <CS 0, 0>], [<CS 7, 2>, <CS 0, 1>, <CS 1, 1>, <CS 2, 1>, <CS 3, 0>],
         [<CS 7, 1>, <CS 3, 1>, <CS 6, 1>], [<CS 7, 0>, <CS 6, 2>], [<CS 6, 0>, <CS 3, 2>, <CS 4, 0>, <CS 5, 0>],
         [<CS 5, 3>, <CS 4, 1>], [<CS 5, 2>, <CS 4, 2>, <CS 2, 3>, <CS 1, 3>, <CS 0, 3>],
         [<CS 4, 3>, <CS 3, 3>, <CS 2, 2>], [<CS 2, 0>, <CS 1, 2>], [<CS 1, 0>, <CS 0, 2>]]
        sage: type(Fs)
        <class 'sage.interfaces.snappy.SnapPyElement'>
        sage: F = Fs.sage(); F
        [[<CS 7, 3>, <CS 6, 3>, <CS 5, 1>, <CS 0, 0>],
         [<CS 7, 2>, <CS 0, 1>, <CS 1, 1>, <CS 2, 1>, <CS 3, 0>],
         [<CS 7, 1>, <CS 3, 1>, <CS 6, 1>],
         [<CS 7, 0>, <CS 6, 2>],
         [<CS 6, 0>, <CS 3, 2>, <CS 4, 0>, <CS 5, 0>],
         [<CS 5, 3>, <CS 4, 1>],
         [<CS 5, 2>, <CS 4, 2>, <CS 2, 3>, <CS 1, 3>, <CS 0, 3>],
         [<CS 4, 3>, <CS 3, 3>, <CS 2, 2>],
         [<CS 2, 0>, <CS 1, 2>],
         [<CS 1, 0>, <CS 0, 2>]]
        sage: Fs00 = Fs[0][0]; Fs00
        <CS 7, 3>
        sage: type(Fs00)
        <class 'sage.interfaces.snappy.SnapPyElement'>
        sage: F00 = Fs00.sage()
        sage: type(F00)
        <class 'spherogram.links.links_base.CrossingStrand'>
        sage: Fs.sage()[0][0] == F00
        True

    TESTS::

        sage: As = snappy.AbelianGroup([3, 0, 7, 0])
        sage: TestSuite(As).run(skip='_test_category')
    """
    def _sage_(self, locals={}):
        r"""
        Attempt to return a Sage version of this object.

        This method works successfully when SnapPy returns a result
        or list of results that consist only of:

        - numbers, i.e. integers, floats, complex numbers;
        - functions and named constants also present in Sage, where:

          * Sage knows how to translate the function or constant's name
            from SnapPy's naming scheme, or
          * you provide a translation dictionary `locals`, or
          * the Sage name for the function or constant is simply the
            SnapPy name in lower case;

        - symbolic variables whose names do not pathologically overlap with
          objects already defined in Sage.

        This method will not work when SnapPy's output includes:

        - strings;
        - functions unknown to Sage;
        - SnapPy functions with different parameters/parameter order to
          the Sage equivalent. In this case, define a function to do the
          parameter conversion, and pass it in via the locals dictionary.

        EXAMPLES::

            sage: Ds = snappy.RationalTangle(3,5).denominator_closure(); Ds
            <Link: 1 comp; 4 cross>
            sage: D = Ds.sage(); D
            Knot represented by 4 crossings
            sage: Ns = snappy.RationalTangle(3,5).numerator_closure(); Ns
            <Link: 1 comp; 4 cross>
            sage: N = Ns.sage(); N
            Knot represented by 4 crossings
            sage: D.is_isotopic(N)
            False
            sage: D.get_knotinfo()
            KnotInfo['K4_1']
            sage: N.get_knotinfo()
            KnotInfo['K3_1m']
        """
        inst = self._inst
        if hasattr(inst, 'sage'):
            return inst.sage()
        if hasattr(inst, 'sage_link'):
            return inst.sage_link()
        if locals:
            # if locals are given we use `_sage_repr`
            # surely this only covers simple cases
            from sage.misc.sage_eval import sage_eval
            return sage_eval(self._sage_repr(), locals=locals)
        return inst


# An instance
snappy = SnapPy()
snappyhp = SnapPy(high_precision=True)
