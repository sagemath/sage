# sage.doctest: optional regina
r"""
Interface to Regina

`Regina <https://regina-normal.github.io/>`__ is an open source software for
low-dimensional topology. From the home-page:

.. NOTE::

    Regina is a software package for 3-manifold and 4-manifold
    topologists, with a focus on triangulations, knots and links,
    normal surfaces, and angle structures.

    For 3-manifolds, it includes high-level tasks such as 3-sphere
    recognition, connected sum decomposition and Hakenness testing,
    comes with a rich database of census manifolds, and incorporates
    the SnapPea kernel for working with hyperbolic manifolds.
    For 4-manifolds, it offers a range of combinatorial and algebraic
    tools, plus support for normal hypersurfaces. For knots and links,
    Regina can perform combinatorial manipulation, compute knot
    polynomials, and work with several import/export formats.

The Regina interface will only work if the optional Sage package Regina
is installed. The interface lets you send certain Sage objects to Regina,
run Regina functions, import certain Regina expressions to Sage,
or any combination of the above.

To send a Sage object ``sobj`` to Regina, call ``regina(sobj)``.
This exports the Sage object to Regina and returns a new Sage object
wrapping the Regina expression/variable, so that you can use the
Regina variable from within Sage. You can then call Regina
functions on the new object; for example::

    sage: F3 = FreeGroup(3)
    sage: F3r = regina(F3); F3r
    <regina.GroupPresentation: < a b c >>
    sage: (F3r.parent(), type(F3r))
    (Regina, <class 'sage.interfaces.regina.ReginaElement'>)
    sage: f = F3.an_element(); f
    x0*x1*x2
    sage: fr = regina(f); fr
    <regina.GroupExpression: g0 g1 g2>
    sage: rel = fr^2; rel
    <regina.GroupExpression: g0 g1 g2 g0 g1 g2>
    sage: (type(fr), type(rel))
    (<class 'sage.interfaces.regina.ReginaElement'>,
     <class 'sage.interfaces.regina.ReginaElement'>)
    sage: F3r.addRelation(rel); F3r
    <regina.GroupPresentation: < a b c | a b c a b c >>

In the above example the relations are added using Regina's
``addRelations`` method.

To see Regina's output you can simply print the Regina wrapper
object. However if you want to import Regina's output back to Sage,
call the Regina wrapper object's ``sage()`` method. This method returns
a native Sage object::

    sage: G3 = F3r.sage(); G3
    Finitely presented group < x0, x1, x2 | (x0*x1*x2)^2 >
    sage: type(G3)
    <class 'sage.groups.finitely_presented.FinitelyPresentedGroup_with_category'>
    sage: regina(G3) == F3r
    True
    sage: fr.sage() == f
    True
    sage: rel.sage() == f^2
    True

If you want to run a Regina function and don't already have the input
in the form of a Sage object, then it might be simpler to input a string
``expr`` to ``regina(expr)``. This string will be evaluated as if you had
typed it into Regina::

    sage: rL = regina("Link('dabcabcv-')"); rL
    <regina.Link: 3-crossing knot: +++ ( ^0 _1 ^2 _0 ^1 _2 )>

Alternatively, all constructors of Regina classes can be used directly as attributes
of the interface::

    sage: rL == regina.Link('dabcabcv-')
    True

Finally, if you just want to use a Regina command line from within
Sage, the IPython magic function ``%regina`` dumps you into an interactive
command-line Regina session. As long as you work in this environment the
prompt is ``regina:``. To finish the environment type ``CTRL+D``::

    sage: %regina                                 # not tested

    --> Switching to Regina <--

    regina: u = Link()
    None
    regina: u
    <regina.Link: Empty link>
    regina: type(u)
    <class 'regina.engine.Link'>
    regina: v = u.fromKnotSig('iabcdbefcdghaefghRsgF+m')
    None
    regina: v
    <regina.Link: 8-crossing knot: ++++--++ ( ^0 ^1 _2 ^3 _1 ^4 _5 ^2 _3 _6 ^7 _0 _4 ^5 ^6 _7 )>
    regina: type(v)
    <class 'regina.engine.Link'>
    regina: v.homfly()
    <regina.Laurent2: 2 x^-2 y^2 + 3 x^-2 - x^-4 y^4 - 3 x^-4 y^2 - 3 x^-4 + x^-6 y^2 + x^-6>
    regina: u.homfly()
    <regina.Laurent2: 0>
    regina: type(u)
    <class 'regina.engine.Link'>

    --> Exiting back to Sage <--

    sage:                                         # not tested


Complicated translations
------------------------

The ``robj.sage()`` method tries to convert a Regina object to a Sage
object. In many cases, it will just work. In particular, it should be able to
convert expressions entirely consisting of:

- numbers, i.e. integers, floats, complex numbers;
- functions and named constants also present in Sage, where:

    - Sage knows how to translate the function or constant's name from
      Regina's, or
    - the Sage name for the function or constant is trivially related to
      Regina's;

- symbolic variables whose names don't pathologically overlap with
  objects already defined in Sage.

This method will not work when Regina's output includes:

- strings;
- functions unknown to Sage;
- Regina functions with different parameters/parameter order to
  the Sage equivalent.


AUTHORS:

- Sebastian Oehms (2025): first version.
"""

##############################################################################
#       Copyright (C) 2025 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
##############################################################################

from enum import Enum

from sage.interfaces.python_internal import (
    PythonInternalElement,
    PythonInternalInterface,
)
from sage.misc.instancedoc import instancedoc


class AlgorithmExt(Enum):
    r"""
    Enum to select algorithm choices.

    This extends the ``Algorithm`` class of Regina.
    """
    ALG_SIMPLIFY = 4
    ALG_WIRTINGER = 5
    ALG_USE_EXTERIOR = 6


class Regina(PythonInternalInterface):
    r"""
    Interface to the Regina interpreter.

    EXAMPLES::

        sage: K = Knots().from_table(8, 21)
        sage: Kr = regina(K); Kr
        <regina.Link: 8-crossing knot: ----++-- ( _5 _0 ^1 _2 _3 ^6 _7 ^3 _4 ^5 _6 ^7 ^0 _1 ^2 ^4 )>
        sage: Kr.knotSig()
        'iabcdbefcdghaefghRsgF+m'

    More examples can be found in the module header.
    """
    def __init__(self):
        r"""
        Python constructor.

        TESTS::

            sage: TestSuite(regina).run(skip=['_test_pickling', '_test_category'])
        """
        super().__init__('regina')

    def _start(self):
        """
        Start up the Regina interpreter and sets the initial prompt and options.

        This is called the first time the Regina interface is actually used.

        EXAMPLES::

            sage: regina._start()
            sage: type(regina._namespace.Cyclotomic())
            <class 'regina.engine.Cyclotomic'>
        """
        if not self._interface_globals:
            from sage.features.interfaces import Regina
            Regina().module.require()
            import regina
            self._namespace = regina.engine
            d = self._namespace.__dict__
            # add extras to the fixed namespace
            d['AlgorithmExt'] = AlgorithmExt
            for e in AlgorithmExt:
                d[e.name] = e
            # set up the globals
            D = dict(d)
            self._interface_modules = [self.name()]
            D[self._namespace.__name__] = self._namespace
            D[self.name()] = regina
            self._interface_globals = D

    def _install_hints(self):
        """
        Hints for installing regina on your computer.

        EXAMPLES::

            sage: len(regina._install_hints())
            99
        """
        return """
In order to use the Regina interface you need to have the
optional Sage package Regina installed.
"""

    def _object_class(self):
        r"""
        Return the element class of this parent.
        This is used in the interface class.

        EXAMPLES::

            sage: regina._object_class()
            <class 'sage.interfaces.regina.ReginaElement'>
        """
        return ReginaElement


@instancedoc
class ReginaElement(PythonInternalElement):
    r"""
    Element class of the Regina interface.

    Its instances are usually constructed via the instance call of its parent.
    It wrapes the Regina library for this object. In a session Regina methods
    can be obtained using tab completion.

    EXAMPLES::

        sage: b = BraidGroup(3)((1,2,-1))
        sage: re = regina(b); re
        <regina.GroupExpression: g0 g1 g0^-1>
        sage: type(re)
        <class 'sage.interfaces.regina.ReginaElement'>
        sage: P = re.parent(); P
        Regina
        sage: type(P)
        <class 'sage.interfaces.regina.Regina'>

    Access to the Regina expression objects::

        sage: res = re._inst
        sage: type(res)
         <class 'regina.engine.GroupExpression'>

    Applying Regina methods::

        sage: re.cycleLeft(); re
        <regina.GroupExpression: g0^-1 g0 g1>

    Conversion to Sage::

        sage: re.sage() == b
        False
        sage: re.cycleRight()
        sage: re.sage() == b
        True

    TESTS::

        sage: p = regina.Polynomial([-2, ~7])
        sage: TestSuite(p).run(skip='_test_category')
    """
    def __deepcopy__(self, memo=None):
        r"""
        EXAMPLES::

            sage: C = CubicBraidGroup(3)
            sage: c = C.an_element()
            sage: cr = regina(c); cr
            <regina.GroupExpression: g0 g1>
            sage: from copy import deepcopy, copy
            sage: crd = deepcopy(cr)
            sage: crc = copy(cr)
            sage: cr.cycleRight(); cr
            <regina.GroupExpression: g1 g0>
            sage: (crd, crd == cr)
            (<regina.GroupExpression: g0 g1>, False)
            sage: (crc, crc == cr)
            (<regina.GroupExpression: g1 g0>, True)
        """
        P = self._check_valid()
        inst = self._inst
        new = inst.__class__(inst)
        res = P(new)
        res._sage_parent = self._sage_parent
        return res

    def _sage_(self, locals={}):
        r"""
        Attempt to return a Sage version of this object.

        This method works successfully when Regina returns a result
        or list of results that consist only of:

        - numbers, i.e. integers, floats, complex numbers;
        - functions and named constants also present in Sage, where:

          * Sage knows how to translate the function or constant's name
            from Regina's naming scheme, or
          * you provide a translation dictionary `locals`, or
          * the Sage name for the function or constant is simply the
            Regina name in lower case;

        - symbolic variables whose names do not pathologically overlap with
          objects already defined in Sage.

        This method will not work when Regina's output includes:

        - strings;
        - functions unknown to Sage;
        - Regina functions with different parameters/parameter order to
          the Sage equivalent. In this case, define a function to do the
          parameter conversion, and pass it in via the locals dictionary.

        EXAMPLES::

            sage: p = regina("Laurent2([(1,2,3), (2,-1,5)])"); p
            <regina.Laurent2: 5 x^2 y^-1 + 3 x y^2>
            sage: p.sage()
            3*x*y^2 + 5*x^2*y^-1
            sage: R.<u> = PolynomialRing(ZZ)
            sage: p = u**3 -2*u + 9
            sage: rp = regina(p); rp
            <regina.PolynomialRational: x^3 - 2 x + 9>
            sage: rp.sage() == p
            True
            sage: F3 = FreeGroup(3)
            sage: a, b, c = F3.gens()
            sage: f = a**2*b*~c
            sage: fr = regina(f); fr
            <regina.GroupExpression: g0^2 g1 g2^-1>
            sage: fr.sage() == f
            True
        """
        def from_detail_str(lc):
            r"""
            Regina provides a detail method for many of its classes.
            Here we try to use it to convert back to Sage.
            """
            if locals:
                lc.update(locals)
            from sage.misc.sage_eval import sage_eval
            from sage.repl.preparse import implicit_mul
            s = self.detail().split('\n')[0]
            s = s.replace(' ', '')
            v = list(lc)
            for i in v + ['(']:
                for j in v + [')']:
                    s = s.replace('%s%s' % (i, j), '%s*%s' % (i, j))
            s = implicit_mul(s)
            return sage_eval(s, locals=lc)

        P = self._check_valid()
        inst = self._inst
        nspc = P._namespace
        if hasattr(inst, 'sage'):
            return inst.sage()
        if isinstance(inst, (nspc.Polynomial, nspc.Laurent, nspc.Laurent2)):
            if self._sage_parent:
                R = self._sage_parent
                old_var_names = ['x', 'y']
                new_var_names = R.variable_names()
                lc = {old_var_names[i]: R.gens_dict()[new_var_names[i]] for i in range(len(new_var_names))}
            elif isinstance(inst, nspc.Polynomial):
                from sage.rings.integer_ring import ZZ
                from sage.rings.polynomial.polynomial_ring_constructor import (
                    PolynomialRing,
                )
                R = PolynomialRing(ZZ, 'x')
                lc = R.gens_dict()
            else:
                from sage.rings.integer_ring import ZZ
                from sage.rings.polynomial.laurent_polynomial_ring import (
                    LaurentPolynomialRing,
                )
                if isinstance(inst, nspc.Laurent):
                    R = LaurentPolynomialRing(ZZ, 'x')
                else:
                    R = LaurentPolynomialRing(ZZ, 'x, y')
                lc = R.gens_dict()
            return from_detail_str(lc)
        if isinstance(inst, nspc.GroupExpression):
            num_gens = max(t.generator for t in inst.terms()) + 1
            if self._sage_parent:
                F = self._sage_parent
            else:
                from sage.groups.free_group import FreeGroup
                F = FreeGroup(num_gens)
            gens = F.gens()
            lc = {'g%s' % i: gens[i] for i in range(num_gens)}
            return from_detail_str(lc)
        if isinstance(inst, nspc.Link):
            from sage.knots.link import Link
            return Link(inst.pdData())
        if hasattr(self, 'detail'):
            return from_detail_str(locals)
        if locals:
            # if locals are given we use `_sage_repr`
            # surely this only covers simple cases
            from sage.misc.sage_eval import sage_eval
            return sage_eval(self._sage_repr(), locals=locals)
        return inst


# An instance
regina = Regina()
