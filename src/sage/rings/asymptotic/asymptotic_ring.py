r"""
Asymptotic Ring

This module provides a ring (called :class:`AsymptoticRing`) for
computations with :wikipedia:`asymptotic expansions <Asymptotic_expansion>`.


.. _asymptotic_ring_definition:

(Informal) Definition
=====================

An asymptotic expansion is a sum such as

.. MATH::

    5z^3 + 4z^2 + O(z)

as `z \to \infty` or

.. MATH::

    3x^{42}y^2 + 7x^3y^3 + O(x^2) + O(y)

as `x` and `y` tend to `\infty`. It is a truncated series (after a
finite number of terms), which approximates a function.

The summands of the asymptotic expansions are partially ordered. In
this module these summands are the following:

- Exact terms `c\cdot g` with a coefficient `c` and an element `g` of
  a growth group (:ref:`see below <asymptotic_ring_growth>`).

- `O`-terms `O(g)` (see :wikipedia:`Big O notation <Big_O_notation>`;
  also called *Bachmann--Landau notation*) for a growth group
  element `g` (:ref:`again see below <asymptotic_ring_growth>`).

See
:wikipedia:`the Wikipedia article on asymptotic expansions <Asymptotic_expansion>`
for more details.
Further examples of such elements can be found :ref:`here <asymptotic_ring_intro>`.


.. _asymptotic_ring_growth:

Growth Groups and Elements
--------------------------

The elements of a :doc:`growth group <growth_group>` are equipped with
a partial order and usually contain a variable. Examples---the order
is described below these examples---are

- elements of the form `z^q` for some integer or rational `q`
  (growth groups with :ref:`description strings <growth_group_description>`
  ``z^ZZ`` or ``z^QQ``),

- elements of the form `\log(z)^q` for some integer or rational `q`
  (growth groups ``log(z)^ZZ`` or ``log(z)^QQ``),

- elements of the form `a^z` for some
  rational `a` (growth group ``QQ^z``), or

- more sophisticated constructions like products
  `x^r \cdot \log(x)^s \cdot a^y \cdot y^q`
  (this corresponds to an element of the growth group
  ``x^QQ * log(x)^ZZ * QQ^y * y^QQ``).

The order in all these examples is induced by the magnitude of the
elements as `x`, `y`, or `z` (independently) tend to `\infty`. For
elements only using the variable `z` this means that `g_1 \leq g_2` if

.. MATH::

    \lim_{z\to\infty} \frac{g_1}{g_2} \leq 1.

.. NOTE::

    Asymptotic rings where the variable tend to some value distinct from
    `\infty` are not yet implemented.

To find out more about

- growth groups,

- on how they are created and

- about the above used *descriptions strings*

see the top of the module :doc:`growth group <growth_group>`.


.. _asymptotic_ring_intro:

Introductory Examples
=====================

We start this series of examples by defining two asymptotic rings.


Two Rings
---------

A Univariate Asymptotic Ring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we construct the following (very simple) asymptotic ring in the variable `z`::

    sage: A.<z> = AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ); A
    Asymptotic Ring <z^QQ> over Integer Ring

A typical element of this ring is
::

    sage: A.an_element()
    z^(3/2) + O(z^(1/2))

This element consists of two summands: the exact term with coefficient
`1` and growth `z^{3/2}` and the `O`-term `O(z^{1/2})`. Note that the
growth of `z^{3/2}` is larger than the growth of `z^{1/2}` as
`z\to\infty`, thus this expansion cannot be simplified (which would
be done automatically, see below).

Elements can be constructed via the generator `z` and the function
:func:`~sage.rings.big_oh.O`, for example

::

    sage: 4*z^2 + O(z)
    4*z^2 + O(z)

A Multivariate Asymptotic Ring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we construct a more sophisticated asymptotic ring in the
variables `x` and `y` by
::

    sage: B.<x, y> = AsymptoticRing(growth_group='x^QQ * log(x)^ZZ * (QQ_+)^y * y^QQ', coefficient_ring=QQ); B
    Asymptotic Ring <x^QQ * log(x)^ZZ * QQ^y * y^QQ> over Rational Field

Again, we can look at a typical (nontrivial) element::

    sage: B.an_element()
    1/8*x^(3/2)*log(x)^3*(1/8)^y*y^(3/2) + O(x^(1/2)*log(x)*(1/2)^y*y^(1/2))

Again, elements can be created using the generators `x` and `y`, as well as
the function :func:`~sage.rings.big_oh.O`::

    sage: log(x)*y/42 + O(1/2^y)
    1/42*log(x)*y + O((1/2)^y)

Arithmetical Operations
-----------------------

In this section we explain how to perform various arithmetical
operations with the elements of the asymptotic rings constructed
above.


The Ring Operations Plus and Times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We start our calculations in the ring
::

    sage: A
    Asymptotic Ring <z^QQ> over Integer Ring

Of course, we can perform the usual ring operations `+` and `*`::

    sage: z^2 + 3*z*(1-z)
    -2*z^2 + 3*z
    sage: (3*z + 2)^3
    27*z^3 + 54*z^2 + 36*z + 8

In addition to that, special powers---our growth group ``z^QQ`` allows
the exponents to be out of `\QQ`---can also be computed::

    sage: (z^(5/2)+z^(1/7)) * z^(-1/5)
    z^(23/10) + z^(-2/35)

The central concepts of computations with asymptotic expansions is
that the `O`-notation can be used. For example, we have
::

    sage: z^3 + z^2 + z + O(z^2)
    z^3 + O(z^2)

where the result is simplified automatically. A more sophisticated example is
::

    sage: (z+2*z^2+3*z^3+4*z^4) * (O(z)+z^2)
    4*z^6 + O(z^5)


Division
^^^^^^^^

The asymptotic expansions support division. For example, we can
expand `1/(z-1)` to a geometric series::

    sage: 1 / (z-1)
    z^(-1) + z^(-2) + z^(-3) + z^(-4) + ... + z^(-20) + O(z^(-21))

A default precision (parameter ``default_prec`` of
:class:`AsymptoticRing`) is predefined. Thus, only the first `20`
summands are calculated. However, if we only want the first `5` exact
terms, we cut of the rest by using
::

    sage: (1 / (z-1)).truncate(5)
    z^(-1) + z^(-2) + z^(-3) + z^(-4) + z^(-5) + O(z^(-6))

or

::

    sage: 1 / (z-1) + O(z^(-6))
    z^(-1) + z^(-2) + z^(-3) + z^(-4) + z^(-5) + O(z^(-6))

Of course, we can work with more complicated expansions as well::

    sage: (4*z+1) / (z^3+z^2+z+O(z^0))
    4*z^(-2) - 3*z^(-3) - z^(-4) + O(z^(-5))

Not all elements are invertible, for instance,

::

    sage: 1 / O(z)
    Traceback (most recent call last):
    ...
    ZeroDivisionError: Cannot invert O(z).

is not invertible, since it includes `0`.


Powers, Expontials and Logarithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It works as simple as it can be; just use the usual operators ``^``,
``exp`` and ``log``. For example, we obtain the usual series expansion
of the logarithm
::

    sage: -log(1-1/z)
    z^(-1) + 1/2*z^(-2) + 1/3*z^(-3) + ... + O(z^(-21))

as `z \to \infty`.

Similarly, we can apply the exponential function of an asymptotic expansion::

    sage: exp(1/z)
    1 + z^(-1) + 1/2*z^(-2) + 1/6*z^(-3) + 1/24*z^(-4) + ... + O(z^(-20))

Arbitrary powers work as well; for example, we have
::

    sage: (1 + 1/z + O(1/z^5))^(1 + 1/z)
    1 + z^(-1) + z^(-2) + 1/2*z^(-3) + 1/3*z^(-4) + O(z^(-5))


Multivariate Arithmetic
^^^^^^^^^^^^^^^^^^^^^^^

Now let us move on to arithmetic in the multivariate ring

::

    sage: B
    Asymptotic Ring <x^QQ * log(x)^ZZ * QQ^y * y^QQ> over Rational Field

.. TODO::

    write this part


More Examples
=============


The mathematical constant e as a limit
--------------------------------------

The base of the natural logarithm `e` satisfies the equation

.. MATH::

    e = \lim_{n\to\infty} \left(1+\frac{1}{n}\right)^n

By using asymptotic expansions, we obtain the more precise result
::

    sage: E.<n> = AsymptoticRing(growth_group='n^ZZ', coefficient_ring=SR, default_prec=5); E
    Asymptotic Ring <n^ZZ> over Symbolic Ring
    sage: (1 + 1/n)^n
    e - 1/2*e*n^(-1) + 11/24*e*n^(-2) - 7/16*e*n^(-3) + 2447/5760*e*n^(-4) + O(n^(-5))


Selected Technical Details
==========================


Coercions and Functorial Constructions
--------------------------------------

The :class:`AsymptoticRing` fully supports
`coercion <../../../../coercion/index.html>`_. For example, the coefficient ring is automatically extended when needed::

    sage: A
    Asymptotic Ring <z^QQ> over Integer Ring
    sage: (z + 1/2).parent()
    Asymptotic Ring <z^QQ> over Rational Field

Here, the coefficient ring was extended to allow `1/2` as a
coefficient. Another example is
::

    sage: C.<c> = AsymptoticRing(growth_group='c^ZZ', coefficient_ring=ZZ['e'])
    sage: C.an_element()
    e^3*c^3 + O(c)
    sage: C.an_element() / 7
    1/7*e^3*c^3 + O(c)

Here the result's coefficient ring is the newly found
::

    sage: (C.an_element() / 7).parent()
    Asymptotic Ring <c^ZZ> over
    Univariate Polynomial Ring in e over Rational Field

Not only the coefficient ring can be extended, but the growth group as
well. For example, we can add/multiply elements of the asymptotic
rings ``A`` and ``C`` to get an expansion of new asymptotic ring::

    sage: r = c*z + c/2 + O(z); r
    c*z + 1/2*c + O(z)
    sage: r.parent()
    Asymptotic Ring <c^ZZ * z^QQ> over
    Univariate Polynomial Ring in e over Rational Field


Data Structures
---------------

The summands of an
:class:`asymptotic expansion <AsymptoticExpansion>` are wrapped
:doc:`growth group elements <growth_group>`.
This wrapping is done by the
:doc:`term monoid module <term_monoid>`.
However, inside an
:class:`asymptotic expansion <AsymptoticExpansion>` these summands
(terms) are stored together with their growth-relationship, i.e., each
summand knows its direct predecessors and successors. As a data
structure a special poset (namely a
:mod:`mutable poset <sage.data_structures.mutable_poset>`)
is used. We can have a look at this::

    sage: b = x^3*y + x^2*y + x*y^2 + O(x) + O(y)
    sage: print(b.summands.repr_full(reverse=True))
    poset(x*y^2, x^3*y, x^2*y, O(x), O(y))
    +-- oo
    |   +-- no successors
    |   +-- predecessors:   x*y^2, x^3*y
    +-- x*y^2
    |   +-- successors:   oo
    |   +-- predecessors:   O(x), O(y)
    +-- x^3*y
    |   +-- successors:   oo
    |   +-- predecessors:   x^2*y
    +-- x^2*y
    |   +-- successors:   x^3*y
    |   +-- predecessors:   O(x), O(y)
    +-- O(x)
    |   +-- successors:   x*y^2, x^2*y
    |   +-- predecessors:   null
    +-- O(y)
    |   +-- successors:   x*y^2, x^2*y
    |   +-- predecessors:   null
    +-- null
    |   +-- successors:   O(x), O(y)
    |   +-- no predecessors


Various
=======

AUTHORS:

- Benjamin Hackl (2015)
- Daniel Krenn (2015)
- Clemens Heuberger (2016)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.


Classes and Methods
===================
"""

# *****************************************************************************
# Copyright (C) 2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.ring import Algebra
from sage.structure.element import CommutativeAlgebraElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.defaults import series_precision
import sage.rings.abc
from sage.rings.all import RIF
from .misc import WithLocals


class NoConvergenceError(RuntimeError):
    r"""
    A special :python:`RuntimeError<library/exceptions.html#exceptions.RuntimeError>`
    which is raised when an algorithm does not converge/stop.
    """
    pass


class AsymptoticExpansion(CommutativeAlgebraElement):
    r"""
    Class for asymptotic expansions, i.e., the elements of an
    :class:`AsymptoticRing`.

    INPUT:

    - ``parent`` -- the parent of the asymptotic expansion.

    - ``summands`` -- the summands as a
      :class:`~sage.data_structures.mutable_poset.MutablePoset`, which
      represents the underlying structure.

    - ``simplify`` -- a boolean (default: ``True``). It controls
      automatic simplification (absorption) of the asymptotic expansion.

    - ``convert`` -- a boolean (default: ``True``). If set, then the
      ``summands`` are converted to the asymptotic ring (the parent of this
      expansion). If not, then the summands are taken as they are. In
      that case, the caller must ensure that the parent of the terms is
      set correctly.

    EXAMPLES:

    There are several ways to create asymptotic expansions; usually
    this is done by using the corresponding :class:`asymptotic rings <AsymptoticRing>`::

        sage: R_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: R_y.<y> = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=ZZ); R_y
        Asymptotic Ring <y^ZZ> over Integer Ring

    At this point, `x` and `y` are already asymptotic expansions::

        sage: type(x)
        <class 'sage.rings.asymptotic.asymptotic_ring.AsymptoticRing_with_category.element_class'>

    The usual ring operations, but allowing rational exponents (growth
    group ``x^QQ``) can be performed::

        sage: x^2 + 3*(x - x^(2/5))
        x^2 + 3*x - 3*x^(2/5)
        sage: (3*x^(1/3) + 2)^3
        27*x + 54*x^(2/3) + 36*x^(1/3) + 8

    One of the central ideas behind computing with asymptotic
    expansions is that the `O`-notation (see
    :wikipedia:`Big_O_notation`) can be used. For example, we have::

        sage: (x+2*x^2+3*x^3+4*x^4) * (O(x)+x^2)
        4*x^6 + O(x^5)

    In particular, :func:`~sage.rings.big_oh.O` can be used to
    construct the asymptotic expansions. With the help of the
    :meth:`summands`, we can also have a look at the inner structure
    of an asymptotic expansion::

        sage: expr1 = x + 2*x^2 + 3*x^3 + 4*x^4; expr2 = O(x) + x^2
        sage: print(expr1.summands.repr_full())
        poset(x, 2*x^2, 3*x^3, 4*x^4)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   x
        +-- x
        |   +-- predecessors:   null
        |   +-- successors:   2*x^2
        +-- 2*x^2
        |   +-- predecessors:   x
        |   +-- successors:   3*x^3
        +-- 3*x^3
        |   +-- predecessors:   2*x^2
        |   +-- successors:   4*x^4
        +-- 4*x^4
        |   +-- predecessors:   3*x^3
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   4*x^4
        |   +-- no successors
        sage: print(expr2.summands.repr_full())
        poset(O(x), x^2)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   O(x)
        +-- O(x)
        |   +-- predecessors:   null
        |   +-- successors:   x^2
        +-- x^2
        |   +-- predecessors:   O(x)
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   x^2
        |   +-- no successors
        sage: print((expr1 * expr2).summands.repr_full())
        poset(O(x^5), 4*x^6)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   O(x^5)
        +-- O(x^5)
        |   +-- predecessors:   null
        |   +-- successors:   4*x^6
        +-- 4*x^6
        |   +-- predecessors:   O(x^5)
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   4*x^6
        |   +-- no successors

    In addition to the monomial growth elements from above, we can
    also compute with logarithmic terms (simply by constructing the
    appropriate growth group)::

        sage: R_log = AsymptoticRing(growth_group='log(x)^QQ', coefficient_ring=QQ)
        sage: lx = R_log(log(SR.var('x')))
        sage: (O(lx) + lx^3)^4
        log(x)^12 + O(log(x)^10)

    .. SEEALSO::

        :doc:`growth_group`,
        :doc:`term_monoid`,
        :mod:`~sage.data_structures.mutable_poset`.
    """
    def __init__(self, parent, summands, simplify=True, convert=True):
        r"""
        See :class:`AsymptoticExpansion` for more information.

        TESTS::

            sage: R_x.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: R_y.<y> = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=ZZ)
            sage: R_x is R_y
            False
            sage: ex1 = x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5
            sage: ex2 = x + O(R_x(1))
            sage: ex1 * ex2
            5*x^6 + O(x^5)

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = TermMonoid('O', G, ZZ); ET = TermMonoid('exact', G, ZZ)
            sage: R = AsymptoticRing(G, ZZ)
            sage: lst = [ET(x, coefficient=1), ET(x^2, coefficient=2), OT(x^3), ET(x^4, coefficient=4)]
            sage: expr = R(lst, simplify=False); expr  # indirect doctest
            4*x^4 + O(x^3) + 2*x^2 + x
            sage: print(expr.summands.repr_full())
            poset(x, 2*x^2, O(x^3), 4*x^4)
            +-- null
            |   +-- no predecessors
            |   +-- successors:   x
            +-- x
            |   +-- predecessors:   null
            |   +-- successors:   2*x^2
            +-- 2*x^2
            |   +-- predecessors:   x
            |   +-- successors:   O(x^3)
            +-- O(x^3)
            |   +-- predecessors:   2*x^2
            |   +-- successors:   4*x^4
            +-- 4*x^4
            |   +-- predecessors:   O(x^3)
            |   +-- successors:   oo
            +-- oo
            |   +-- predecessors:   4*x^4
            |   +-- no successors
            sage: expr._simplify_(); expr
            4*x^4 + O(x^3)
            sage: print(expr.summands.repr_full())
            poset(O(x^3), 4*x^4)
            +-- null
            |   +-- no predecessors
            |   +-- successors:   O(x^3)
            +-- O(x^3)
            |   +-- predecessors:   null
            |   +-- successors:   4*x^4
            +-- 4*x^4
            |   +-- predecessors:   O(x^3)
            |   +-- successors:   oo
            +-- oo
            |   +-- predecessors:   4*x^4
            |   +-- no successors
            sage: R(lst, simplify=True) # indirect doctest
            4*x^4 + O(x^3)

        ::

            sage: R.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ)
            sage: e = R(x^2 + O(x))
            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticExpansion
            sage: S = AsymptoticRing(growth_group='x^QQ', coefficient_ring=ZZ)
            sage: for s in AsymptoticExpansion(S, e.summands).summands.elements_topological():
            ....:     print(s.parent())
            O-Term Monoid x^QQ with implicit coefficients in Integer Ring
            Exact Term Monoid x^QQ with coefficients in Integer Ring
            sage: for s in AsymptoticExpansion(S, e.summands,
            ....:         convert=False).summands.elements_topological():
            ....:     print(s.parent())
            O-Term Monoid x^QQ with implicit coefficients in Rational Field
            Exact Term Monoid x^QQ with coefficients in Rational Field

        ::

            sage: AsymptoticExpansion(S, R(1/2).summands)
            Traceback (most recent call last):
            ...
            ValueError: Cannot include 1/2 with parent
            Exact Term Monoid x^QQ with coefficients in Rational Field in
            Asymptotic Ring <x^QQ> over Integer Ring
            > *previous* ValueError: Cannot create ExactTerm(1)
              since given coefficient 1/2 is not valid in
              Exact Term Monoid x^QQ with coefficients in Integer Ring.
            >> *previous* TypeError: no conversion of this rational to integer

        Check :trac:`19921`::

            sage: CR.<Z> = QQ['Z']
            sage: CR_mod = CR.quotient((Z^2 - 1)*CR)
            sage: R.<x> = AsymptoticRing(growth_group='x^NN', coefficient_ring=CR)
            sage: R_mod = R.change_parameter(coefficient_ring=CR_mod)
            sage: e = 1 + x*(Z^2-1)
            sage: R_mod(e)
            1

        Check that :trac:`19999` is resolved::

            sage: A.<x> = AsymptoticRing('(QQ_+)^x * x^QQ * UU^x', QQ)
            sage: 1 + (-1)^x + 2^x + (-2)^x
            2^x + 2^x*(-1)^x + (-1)^x + 1

            sage: A.<x> = AsymptoticRing('QQ^x * x^QQ', QQ)
            sage: 1 + (-1)^x + 2^x + (-2)^x
            2^x + 2^x*(-1)^x + (-1)^x + 1
        """
        super(AsymptoticExpansion, self).__init__(parent=parent)

        from sage.data_structures.mutable_poset import MutablePoset
        if not isinstance(summands, MutablePoset):
            raise TypeError('Summands %s are not in a mutable poset as expected '
                            'when creating an element of %s.' % (summands, parent))

        if convert:
            from .misc import combine_exceptions
            from .term_monoid import ZeroCoefficientError

            def convert_terms(element):
                T = self.parent().term_monoid(element.parent())
                try:
                    return T(element)
                except ZeroCoefficientError:
                    return None
                except (ValueError, TypeError) as e:
                    raise combine_exceptions(
                        ValueError('Cannot include %s with parent %s in %s' %
                                   (element, element.parent(), parent)), e)
            new_summands = summands.copy()
            new_summands.map(convert_terms, topological=True, reverse=True)
            self._summands_ = new_summands
        else:
            self._summands_ = summands

        if simplify:
            self._simplify_()

    @property
    def summands(self):
        r"""
        The summands of this asymptotic expansion stored in the
        underlying data structure (a
        :class:`~sage.data_structures.mutable_poset.MutablePoset`).

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: expr = 7*x^12 + x^5 + O(x^3)
            sage: expr.summands
            poset(O(x^3), x^5, 7*x^12)

        .. SEEALSO::

            :class:`sage.data_structures.mutable_poset.MutablePoset`
        """
        return self._summands_


    def __hash__(self):
        r"""
        A hash value for this element.

        .. WARNING::

            This hash value uses the string representation and might not be
            always right.

        TESTS::

            sage: R_log = AsymptoticRing(growth_group='log(x)^QQ', coefficient_ring=QQ)
            sage: lx = R_log(log(SR.var('x')))
            sage: elt = (O(lx) + lx^3)^4
            sage: hash(elt) # random
            -4395085054568712393
        """
        return hash(str(self))


    def __bool__(self):
        r"""
        Return whether this asymptotic expansion is not identically zero.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        TESTS::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: bool(R(0))  # indirect doctest
            False
            sage: bool(x)  # indirect doctest
            True
            sage: bool(7*x^12 + x^5 + O(x^3))  # indirect doctest
            True
        """
        return bool(self._summands_)

    __nonzero__ = __bool__

    def __eq__(self, other):
        r"""
        Return whether this asymptotic expansion is equal to ``other``.

        INPUT:

        - ``other`` -- an object.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: (1 + 2*x + 3*x^2) == (3*x^2 + 2*x + 1)  # indirect doctest
            True
            sage: O(x) == O(x)
            False

        TESTS::

            sage: x == None
            False

        ::

            sage: x == 'x'
            False
        """
        if other is None:
            return False
        try:
            return not bool(self - other)
        except (TypeError, ValueError):
            return False


    def __ne__(self, other):
        r"""
        Return whether this asymptotic expansion is not equal to ``other``.

        INPUT:

        - ``other`` -- an object.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: (1 + 2*x + 3*x^2) != (3*x^2 + 2*x + 1)  # indirect doctest
            False
            sage: O(x) != O(x)
            True

        TESTS::

            sage: x != None
            True
        """
        return not self == other


    def has_same_summands(self, other):
        r"""
        Return whether this asymptotic expansion and ``other`` have the
        same summands.

        INPUT:

        - ``other`` -- an asymptotic expansion.

        OUTPUT:

        A boolean.

        .. NOTE::

            While for example ``O(x) == O(x)`` yields ``False``,
            these expansions *do* have the same summands and this method
            returns ``True``.

            Moreover, this method uses the coercion model in order to
            find a common parent for this asymptotic expansion and
            ``other``.

        EXAMPLES::

            sage: R_ZZ.<x_ZZ> = AsymptoticRing('x^ZZ', ZZ)
            sage: R_QQ.<x_QQ> = AsymptoticRing('x^ZZ', QQ)
            sage: sum(x_ZZ^k for k in range(5)) == sum(x_QQ^k for k in range(5))  # indirect doctest
            True
            sage: O(x_ZZ) == O(x_QQ)
            False

        TESTS::

            sage: x_ZZ.has_same_summands(None)
            False
        """
        if other is None:
            return False
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._has_same_summands_(other)

        from sage.structure.element import get_coercion_model
        return get_coercion_model().bin_op(self, other,
                                           lambda self, other:
                                           self._has_same_summands_(other))


    def _has_same_summands_(self, other):
        r"""
        Return whether this :class:`AsymptoticExpansion` has the same
        summands as ``other``.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method compares two :class:`AsymptoticExpansion`
            with the same parent.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: O(x).has_same_summands(O(x))
            True
            sage: (1 + x + 2*x^2).has_same_summands(2*x^2 + O(x))  # indirect doctest
            False
        """
        if len(self.summands) != len(other.summands):
            return False
        from builtins import zip
        return all(s == o for s, o in
                   zip(self.summands.elements_topological(),
                       other.summands.elements_topological()))

    def _simplify_(self):
        r"""
        Simplify this asymptotic expansion.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing, but modifies this asymptotic expansion.

        .. NOTE::

            This method is usually called during initialization of
            this asymptotic expansion.

        .. NOTE::

            This asymptotic expansion is simplified by letting
            `O`-terms that are included in this expansion absorb all
            terms with smaller growth.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: OT = TermMonoid('O', G, ZZ); ET = TermMonoid('exact', G, ZZ)
            sage: R = AsymptoticRing(G, ZZ)
            sage: lst = [ET(x, coefficient=1), ET(x^2, coefficient=2), OT(x^3), ET(x^4, coefficient=4)]
            sage: expr = R(lst, simplify=False); expr  # indirect doctest
            4*x^4 + O(x^3) + 2*x^2 + x
            sage: expr._simplify_(); expr
            4*x^4 + O(x^3)
            sage: R(lst)  # indirect doctest
            4*x^4 + O(x^3)
        """
        self._summands_.merge(reverse=True)


    def _repr_(self, latex=False):
        r"""
        A representation string for this asymptotic expansion.

        INPUT:

        - ``latex`` -- (default: ``False``) a boolean. If set, then
          LaTeX-output is returned.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: (5*x^2+12*x) * (x^3+O(x))  # indirect doctest
            5*x^5 + 12*x^4 + O(x^3)
            sage: (5*x^2-12*x) * (x^3+O(x))  # indirect doctest
            5*x^5 - 12*x^4 + O(x^3)
        """
        if latex:
            from sage.misc.latex import latex as latex_repr
            f = latex_repr
        else:
            f = repr
        s = ' + '.join(f(elem) for elem in
                       self.summands.elements_topological(reverse=True,
                                                          key=repr))
        s = s.replace('+ -', '- ')
        if not s:
            return '0'
        return s


    def _latex_(self):
        r"""
        A LaTeX-representation string for this asymptotic expansion.

        OUTPUT:

        A string.

        TESTS::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: latex((5*x^2+12*x) * (x^3+O(x)))  # indirect doctest
            5 x^{5} + 12 x^{4} + O\!\left(x^{3}\right)
            sage: latex((5*x^2-12*x) * (x^3+O(x)))  # indirect doctest
            5 x^{5} - 12 x^{4} + O\!\left(x^{3}\right)
        """
        return self._repr_(latex=True)


    def show(self):
        r"""
        Pretty-print this asymptotic expansion.

        OUTPUT:

        Nothing, the representation is printed directly on the
        screen.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('QQ^x * x^QQ * log(x)^QQ', SR.subring(no_variables=True))
            sage: (pi/2 * 5^x * x^(42/17) - sqrt(euler_gamma) * log(x)^(-7/8)).show()
            1/2*pi*5^x*x^(42/17) - sqrt(euler_gamma)*log(x)^(-7/8)

        TESTS::

            sage: A.<x> = AsymptoticRing('(e^x)^QQ * x^QQ', SR.subring(no_variables=True))
            sage: (zeta(3) * (e^x)^(-1/2) * x^42).show()
            zeta(3)*(e^x)^(-1/2)*x^42
        """
        from sage.repl.rich_output.pretty_print import pretty_print
        pretty_print(self)

    def monomial_coefficient(self, monomial):
        r"""
        Return the coefficient in the base ring of the given monomial
        in this expansion.

        INPUT:

        - ``monomial`` -- a monomial element which can be converted
          into the asymptotic ring of this element

        OUTPUT:

        An element of the coefficient ring.

        EXAMPLES::

            sage: R.<m, n> = AsymptoticRing("m^QQ*n^QQ", QQ)
            sage: ae = 13 + 42/n + 2/n/m + O(n^-2)
            sage: ae.monomial_coefficient(1/n)
            42
            sage: ae.monomial_coefficient(1/n^3)
            0
            sage: R.<n> = AsymptoticRing("n^QQ", ZZ)
            sage: ae.monomial_coefficient(1/n)
            42
            sage: ae.monomial_coefficient(1)
            13

        TESTS:

        Conversion of ``monomial`` the parent of this element must be
        possible::

            sage: R.<m> = AsymptoticRing("m^QQ", QQ)
            sage: S.<n> = AsymptoticRing("n^QQ", QQ)
            sage: m.monomial_coefficient(n)
            Traceback (most recent call last):
            ...
            ValueError: Cannot include n with parent Exact Term Monoid
            n^QQ with coefficients in Rational Field in Asymptotic Ring
            <m^QQ> over Rational Field
            > *previous* ValueError: Growth n is not valid in
              Exact Term Monoid m^QQ with coefficients in Rational Field.
            >> *previous* ValueError: n is not in Growth Group m^QQ.

        Only monomials are allowed::

            sage: R.<n> = AsymptoticRing("n^QQ", QQ)
            sage: (n + 4).monomial_coefficient(n + 5)
            Traceback (most recent call last):
            ...
            ValueError: n + 5 not a monomial
            sage: n.monomial_coefficient(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 not a monomial

        Cannot extract the coefficient of an O term::

            sage: O(n).monomial_coefficient(n)
            Traceback (most recent call last):
            ...
            AttributeError: 'OTermMonoid_with_category.element_class'
            object has no attribute 'coefficient'

        The ``monomial`` must be exact::

            sage: n.monomial_coefficient(O(n))
            Traceback (most recent call last):
            ...
            ValueError: non-exact monomial O(n)

        """
        monomial = self.parent()(monomial)
        if not monomial.is_exact():
            raise ValueError("non-exact monomial {}".format(monomial))

        if len(monomial.summands) != 1:
            raise ValueError("{} not a monomial".format(monomial))

        monomial_growth = next(monomial.summands.elements()).growth
        try:
            return self.summands.element(monomial_growth).coefficient
        except KeyError:
            return self.parent().coefficient_ring(0)

    def _add_(self, other):
        r"""
        Add ``other`` to this asymptotic expansion.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        The sum as an :class:`AsymptoticExpansion`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1._add_(expr2)
            x^321 + x^123
            sage: expr1 + expr2  # indirect doctest
            x^321 + x^123

        If an `O`-term is added to an asymptotic expansion, then
        the `O`-term absorbs everything it can::

            sage: x^123 + x^321 + O(x^555)  # indirect doctest
            O(x^555)

        TESTS::

            sage: x + O(x)
            O(x)
            sage: O(x) + x
            O(x)
        """
        return self.parent()(self.summands.union(other.summands),
                             simplify=True, convert=False)


    def _sub_(self, other):
        r"""
        Subtract ``other`` from this asymptotic expansion.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        The difference as an :class:`AsymptoticExpansion`.

        .. NOTE::

            Subtraction of two asymptotic expansions is implemented
            by means of addition: `e_1 - e_2 = e_1 + (-1)\cdot e_2`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1 - expr2  # indirect doctest
            -x^321 + x^123
            sage: O(x) - O(x)
            O(x)
        """
        return self + self.parent().coefficient_ring(-1)*other


    def _mul_term_(self, term):
        r"""
        Helper method: multiply this asymptotic expansion by the
        asymptotic term ``term``.

        INPUT:

        - ``term`` -- an asymptotic term (see
          :doc:`term_monoid`).

        OUTPUT:

        The product as an :class:`AsymptoticExpansion`.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import OTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: T = OTermMonoid(TermMonoid, R.growth_group, ZZ)
            sage: expr = 10*x^2 + O(x)
            sage: t = T(R.growth_group.gen())
            sage: expr._mul_term_(t)
            O(x^3)
        """
        simplify = not term.is_exact()
        return self.parent()(self.summands.mapped(lambda element: term * element),
                             simplify=simplify, convert=False)


    def _mul_(self, other):
        r"""
        Multiply this asymptotic expansion by another asymptotic expansion ``other``.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        The product as an :class:`AsymptoticExpansion`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: ex1 = 5*x^12
            sage: ex2 = x^3 + O(x)
            sage: ex1 * ex2  # indirect doctest
            5*x^15 + O(x^13)

        .. TODO::

            The current implementation is the standard long
            multiplication. More efficient variants like Karatsuba
            multiplication, or methods that exploit the structure
            of the underlying poset shall be implemented at a later
            point.

        TESTS::

            sage: R(1) * R(0)
            0
            sage: _.parent()
            Asymptotic Ring <x^ZZ> over Integer Ring
        """
        return sum(iter(self._mul_term_(term_other) for
                        term_other in other.summands.elements()),
                   self.parent().zero())

    def _lmul_(self, other):
        r"""
        Multiply this asymptotic expansion by an element ``other`` of its
        coefficient ring.

        INPUT:

        - ``other`` -- an element of the coefficient ring.

        OUTPUT:

        An :class:`AsymptoticExpansion`.

        TESTS::

            sage: A.<a> = AsymptoticRing(growth_group='QQ^a * a^QQ * log(a)^QQ', coefficient_ring=ZZ)
            sage: 2*a # indirect doctest
            2*a
        """
        if other.is_zero():
            return self.parent().zero()

        E = self.parent().term_monoid('exact')
        e = E(self.parent().growth_group.one(), coefficient=other)
        return self._mul_term_(e)

    def _div_(self, other):
        r"""
        Divide this element through ``other``.

        INPUT:

        - ``other`` -- an asymptotic expansion.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ, default_prec=5)
            sage: 1/x^42
            x^(-42)
            sage: (1 + 4*x) / (x + 2*x^2)
            2*x^(-1) - 1/2*x^(-2) + 1/4*x^(-3) - 1/8*x^(-4) + 1/16*x^(-5) + O(x^(-6))
            sage: x / O(x)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert O(x).

        TESTS:

        See :trac:`19521`::

            sage: A.<n> = AsymptoticRing('n^ZZ', SR.subring(no_variables=True))
            sage: (A.one() / 1).parent()
            Asymptotic Ring <n^ZZ> over Symbolic Constants Subring
        """
        return self * ~other


    def __invert__(self, precision=None):
        r"""
        Return the multiplicative inverse of this element.

        INPUT:

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. WARNING::

            Due to truncation of infinite expansions, the element
            returned by this method might not fulfill
            ``el * ~el == 1``.

        .. TODO::

            As soon as `L`-terms are implemented, this
            implementation has to be adapted as well in order to
            yield correct results.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ, default_prec=4)
            sage: ~x
            x^(-1)
            sage: ~(x^42)
            x^(-42)
            sage: ex = ~(1 + x); ex
            x^(-1) - x^(-2) + x^(-3) - x^(-4) + O(x^(-5))
            sage: ex * (1+x)
            1 + O(x^(-4))
            sage: ~(1 + O(1/x))
            1 + O(x^(-1))

        TESTS::

            sage: A.<a> = AsymptoticRing(growth_group='a^ZZ', coefficient_ring=ZZ)
            sage: (1 / a).parent()
            Asymptotic Ring <a^ZZ> over Rational Field
            sage: (a / 2).parent()
            Asymptotic Ring <a^ZZ> over Rational Field

        ::

            sage: ~A(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert 0 in
            Asymptotic Ring <a^ZZ> over Integer Ring.

        ::

            sage: B.<s, t> = AsymptoticRing(growth_group='s^ZZ * t^ZZ', coefficient_ring=QQ)
            sage: ~(s + t)
            Traceback (most recent call last):
            ...
            ValueError: Cannot determine main term of s + t since there
            are several maximal elements s, t.
        """
        if not self.summands:
            raise ZeroDivisionError(
                'Cannot invert {} in {}.'.format(self, self.parent()))

        (imax_elem, x) = self._main_term_relative_error_(return_inverse_main_term=True)
        one = x.parent().one()

        if x:
            import itertools
            result = AsymptoticExpansion._power_series_(
                coefficients=itertools.repeat(one),
                start=one,
                ratio=-x,
                ratio_start=one,
                precision=precision)
        else:
            result = one

        return result._mul_term_(imax_elem)


    invert = __invert__


    def truncate(self, precision=None):
        r"""
        Truncate this asymptotic expansion.

        INPUT:

        - ``precision`` -- a positive integer or ``None``. Number of
          summands that are kept. If ``None`` (default value) is
          given, then ``default_prec`` from the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            For example, truncating an asymptotic expansion with
            ``precision=20`` does not yield an expansion with exactly 20
            summands! Rather than that, it keeps the 20 summands
            with the largest growth, and adds appropriate
            `O`-Terms.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: ex = sum(x^k for k in range(5)); ex
            x^4 + x^3 + x^2 + x + 1
            sage: ex.truncate(precision=2)
            x^4 + x^3 + O(x^2)
            sage: ex.truncate(precision=0)
            O(x^4)
            sage: ex.truncate()
            x^4 + x^3 + x^2 + x + 1
        """
        if precision is None:
            precision = self.parent().default_prec

        if len(self.summands) <= precision:
            return self

        summands = self.summands.copy()
        def convert_terms(element):
            if convert_terms.count < precision:
                convert_terms.count += 1
                return element
            T = self.parent().term_monoid('O')
            return T(element)
        convert_terms.count = 0
        summands.map(convert_terms, topological=True, reverse=True)
        return self.parent()(summands, simplify=True, convert=False)


    def exact_part(self):
        r"""
        Return the expansion consisting of all exact terms of this
        expansion.

        INPUT:

        Nothing

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^QQ * log(x)^QQ', QQ)
            sage: (x^2 + O(x)).exact_part()
            x^2
            sage: (x + log(x)/2 + O(log(x)/x)).exact_part()
            x + 1/2*log(x)

        TESTS::

            sage: R.<x, y> = AsymptoticRing('x^QQ * y^QQ', QQ)
            sage: (x + y + O(1/(x*y))).exact_part()
            x + y
            sage: O(x).exact_part()
            0
        """
        exact_terms = self.summands.copy()
        for term in self.summands.elements_topological():
            if not term.is_exact():
                exact_terms.remove(term.growth)

        return self.parent(exact_terms)


    def __pow__(self, exponent, precision=None):
        r"""
        Calculate the power of this asymptotic expansion to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: Q.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ)
            sage: x^(1/7)
            x^(1/7)
            sage: (x^(1/2) + O(x^0))^15
            x^(15/2) + O(x^7)

        ::

            sage: Z.<y> = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=ZZ)
            sage: y^(1/7)
            y^(1/7)
            sage: _.parent()
            Asymptotic Ring <y^QQ> over Rational Field
            sage: (y^2 + O(y))^(1/2)
            y + O(1)
            sage: (y^2 + O(y))^(-2)
            y^(-4) + O(y^(-5))
            sage: (1 + 1/y + O(1/y^3))^pi
            1 + pi*y^(-1) + (1/2*pi*(pi - 1))*y^(-2) + O(y^(-3))

        ::

            sage: B.<z> = AsymptoticRing(growth_group='z^QQ * log(z)^QQ', coefficient_ring=QQ)
            sage: (z^2 + O(z))^(1/2)
            z + O(1)

        ::

            sage: A.<x> = AsymptoticRing('QQ^x * x^SR * log(x)^ZZ', QQ)
            sage: x * 2^x
            2^x*x
            sage: 5^x * 2^x
            10^x
            sage: 2^log(x)
            x^(log(2))
            sage: 2^(x + 1/x)
            2^x + log(2)*2^x*x^(-1) + 1/2*log(2)^2*2^x*x^(-2) + ... + O(2^x*x^(-20))
            sage: _.parent()
            Asymptotic Ring <QQ^x * x^SR * log(x)^QQ * Signs^x> over Symbolic Ring

        ::

            sage: C.<c> = AsymptoticRing(growth_group='QQ^c * c^QQ', coefficient_ring=QQ, default_prec=5)
            sage: (3 + 1/c^2)^c
            3^c + 1/3*3^c*c^(-1) + 1/18*3^c*c^(-2) - 4/81*3^c*c^(-3)
            - 35/1944*3^c*c^(-4) + O(3^c*c^(-5))
            sage: _.parent()
            Asymptotic Ring <QQ^c * c^QQ * Signs^c> over Rational Field
            sage: (2 + (1/3)^c)^c
            2^c + 1/2*(2/3)^c*c + 1/8*(2/9)^c*c^2 - 1/8*(2/9)^c*c
            + 1/48*(2/27)^c*c^3 + O((2/27)^c*c^2)
            sage: _.parent()
            Asymptotic Ring <QQ^c * c^QQ * Signs^c> over Rational Field

        TESTS:

        See :trac:`19110`::

            sage: O(x)^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take O(x) to exponent -1.
            > *previous* ZeroDivisionError: rational division by zero

        ::

            sage: B.<z> = AsymptoticRing(growth_group='z^QQ * log(z)^QQ', coefficient_ring=QQ, default_prec=5)
            sage: z^(1+1/z)
            z + log(z) + 1/2*z^(-1)*log(z)^2 + 1/6*z^(-2)*log(z)^3 +
            1/24*z^(-3)*log(z)^4 + O(z^(-4)*log(z)^5)
            sage: _.parent()
            Asymptotic Ring <z^QQ * log(z)^QQ> over Rational Field

        ::

            sage: B(0)^(-7)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take 0 to the negative exponent -7.
            sage: B(0)^SR.var('a')
            Traceback (most recent call last):
            ...
            NotImplementedError: Taking 0 to the exponent a not implemented.

        ::

            sage: C.<s, t> = AsymptoticRing(growth_group='s^QQ * t^QQ', coefficient_ring=QQ)
            sage: (s + t)^s
            Traceback (most recent call last):
            ...
            ValueError: Cannot take s + t to the exponent s.
            > *previous* ValueError: Cannot determine main term of s + t
            since there are several maximal elements s, t.

        Check that :trac:`19945` is fixed::

            sage: A.<n> = AsymptoticRing('QQ^n * n^QQ', ZZ)
            sage: (1/2)^n
            (1/2)^n

        Check that :trac:`19946` is fixed::

            sage: assume(SR.an_element() > 0)
            sage: A.<n> = AsymptoticRing('QQ^n * n^QQ', SR)
            sage: e = 2^n; e
            2^n
            sage: e.parent()
            Asymptotic Ring <SR^n * n^QQ * Signs^n> over Symbolic Ring
            sage: e = A(e); e
            2^n
            sage: e.parent()
            Asymptotic Ring <QQ^n * n^QQ * Signs^n> over Symbolic Ring
            sage: forget()

        :trac:`22120`::

            sage: A.<w> = AsymptoticRing('w^QQbar', QQ)
            sage: w^QQbar(sqrt(2))
            w^(1.414213562373095?)
        """
        from .misc import strip_symbolic
        exponent = strip_symbolic(exponent)

        if not self.summands:
            if exponent == 0:
                return self.parent().one()
            elif exponent > 0:
                return self.parent().zero()
            elif exponent < 0:
                raise ZeroDivisionError('Cannot take %s to the negative exponent %s.' %
                                        (self, exponent))
            else:
                raise NotImplementedError('Taking %s to the exponent %s not implemented.' %
                                          (self, exponent))

        elif exponent == 0:
            return self.parent().one()

        elif exponent == 1:
            return self

        elif len(self.summands) == 1:
            element = next(self.summands.elements())
            if isinstance(exponent, AsymptoticExpansion) and element.is_constant():

                return exponent.rpow(base=element.coefficient, precision=precision)
            try:
                return self.parent()._create_element_in_extension_(
                    element ** exponent, element.parent())
            except (ArithmeticError, TypeError, ValueError):
                if not isinstance(exponent, AsymptoticExpansion):
                    raise

        from sage.rings.integer_ring import ZZ
        try:
            exponent = ZZ(exponent)
        except (TypeError, ValueError):
            pass
        else:
            return super(AsymptoticExpansion, self).__pow__(exponent)

        from sage.rings.rational_field import QQ
        try:
            exponent = QQ(exponent)
        except (TypeError, ValueError):
            pass
        else:
            return self.__pow_number__(exponent, precision=precision)

        from sage.structure.element import Expression
        if isinstance(exponent, Expression) and exponent.is_constant():
            return self.__pow_number__(exponent, precision=precision)

        if isinstance(exponent, AsymptoticExpansion) and len(self.summands) != 1:
            try:
                return self.__pow_number__(exponent, precision=precision,
                                           check_convergence=True)
            except NoConvergenceError:
                pass

        try:
            return (exponent * self.log(precision=precision)).exp(precision=precision)
        except (TypeError, ValueError, ZeroDivisionError) as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                ValueError('Cannot take %s to the exponent %s.' % (self, exponent)), e)


    pow = __pow__


    def __pow_number__(self, exponent, precision=None, check_convergence=False):
        r"""
        Return the power of this asymptotic expansion to some
        number (``exponent``).

        Let `m` be the maximal element of this asymptotic expansion
        and `r` the remaining summands. This method calculates

        .. MATH::

            (m + r)^{\mathit{exponent}}
            = m^{\mathit{exponent}} \sum_{k=0}^K
            \binom{\mathit{exponent}}{k} (r/m)^k

        where `K` is chosen such that adding an additional summand
        does not change the result.

        INPUT:

        - ``exponent`` -- a numerical value (e.g. integer, rational)
          or other constant.

        - ``precision`` -- a non-negative integer.

        - ``check_convergence`` -- (default: ``False``) a boolean. If set,
          then an additional check on the input is performed to ensure
          that the calculated sum converges.

        OUTPUT:

        An asymptotic expansion.

        .. SEEALSO::

            :meth:`pow`

        TESTS::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: (1 + x).__pow_number__(4)
            x^4 + 4*x^3 + 6*x^2 + 4*x + 1
            sage: _.parent()
            Asymptotic Ring <x^ZZ> over Rational Field
            sage: (x + 1).__pow_number__(1/2, precision=5)
            x^(1/2) + 1/2*x^(-1/2) - 1/8*x^(-3/2) + 1/16*x^(-5/2)
            - 5/128*x^(-7/2) + O(x^(-9/2))
            sage: _.parent()
            Asymptotic Ring <x^QQ> over Rational Field
            sage: (8 + 1/x).__pow_number__(1/3, precision=5)
            2 + 1/12*x^(-1) - 1/288*x^(-2) + 5/20736*x^(-3)
            - 5/248832*x^(-4) + O(x^(-5))
            sage: _.parent()
            Asymptotic Ring <x^QQ> over Rational Field

        ::

            sage: R(0).__pow_number__(-3/2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take 0 to the negative exponent -3/2.
            sage: R(0).__pow_number__(RIF(-1,1))
            Traceback (most recent call last):
            ...
            ValueError: Possible division by zero, since sign of
            the exponent 0.? cannot be determined.
            sage: R(0)^0
            1

        ::

            sage: A.<a, b> = AsymptoticRing(growth_group='a^ZZ * b^ZZ', coefficient_ring=QQ)
            sage: (a + b).__pow_number__(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Cannot determine main term of a + b since
            there are several maximal elements a, b.

        ::

            sage: S.<s> = AsymptoticRing(growth_group='QQ^s * s^ZZ', coefficient_ring=QQ)
            sage: (2 + 2/s^2).__pow_number__(s, precision=7)
            2^s + 2^s*s^(-1) + 1/2*2^s*s^(-2) - 1/3*2^s*s^(-3)
            - 11/24*2^s*s^(-4) + 11/120*2^s*s^(-5)
            + 271/720*2^s*s^(-6) + O(2^s*s^(-7))
            sage: _.parent()
            Asymptotic Ring <QQ^s * s^QQ * Signs^s> over Rational Field

            sage: S.<s> = AsymptoticRing(growth_group='(QQ_+)^s * s^ZZ', coefficient_ring=QQ)
            sage: (2 + 2/s^2).__pow_number__(s, precision=7)
            2^s + 2^s*s^(-1) + 1/2*2^s*s^(-2) - 1/3*2^s*s^(-3)
            - 11/24*2^s*s^(-4) + 11/120*2^s*s^(-5)
            + 271/720*2^s*s^(-6) + O(2^s*s^(-7))
            sage: _.parent()
            Asymptotic Ring <QQ^s * s^QQ> over Rational Field
        """
        if not self.summands:
            if exponent > 0:
                return self.parent().zero()
            elif exponent.is_zero():
                return self.parent().one()
            elif exponent < 0:
                raise ZeroDivisionError(
                    'Cannot take {} to the negative '
                    'exponent {}.'.format(self, exponent))
            else:
                raise ValueError(
                    'Possible division by zero, since sign of the exponent '
                    '{} cannot be determined.'.format(exponent))

        elif len(self.summands) == 1:
            element = next(self.summands.elements())
            return self.parent()._create_element_in_extension_(
                element**exponent, element.parent())

        try:
            (max_elem, x) = self._main_term_relative_error_()
        except ValueError:
            if check_convergence:
                raise NoConvergenceError
            raise

        if check_convergence:
            if not (x * exponent).is_little_o_of_one():
                raise NoConvergenceError

        pmax = self.parent()(max_elem)**exponent

        import itertools
        def binomials(a):
            P = a.parent()
            a = a + 1
            f = P(1)
            for k in itertools.count(1):
                k = P(k)
                b = a - k
                if b == 0:
                    return
                f *= b / k
                yield f

        one = x.parent().one()

        result = AsymptoticExpansion._power_series_(
            coefficients=binomials(exponent),
            start=one,
            ratio=x,
            ratio_start=one,
            precision=precision)

        return result * pmax


    def sqrt(self, precision=None):
        r"""
        Return the square root of this asymptotic expansion.

        INPUT:

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: A.<s> = AsymptoticRing(growth_group='s^QQ', coefficient_ring=QQ)
            sage: s.sqrt()
            s^(1/2)
            sage: a = (1 + 1/s).sqrt(precision=6); a
            1 + 1/2*s^(-1) - 1/8*s^(-2) + 1/16*s^(-3)
            - 5/128*s^(-4) + 7/256*s^(-5) + O(s^(-6))

        .. SEEALSO::

            :meth:`pow`, :meth:`rpow`, :meth:`exp`.

        TESTS::

            sage: P.<p> = PowerSeriesRing(QQ, default_prec=6)
            sage: bool(SR(a.exact_part()).subs(s=1/x) -
            ....:      SR((1+p).sqrt().polynomial()).subs(p=x) == 0)
            True
            """
        from sage.rings.rational_field import QQ
        return self.pow(QQ(1)/QQ(2), precision=precision)


    def O(self):
        r"""
        Convert all terms in this asymptotic expansion to `O`-terms.

        INPUT:

        Nothing.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: O(x)
            O(x)
            sage: type(O(x))
            <class 'sage.rings.asymptotic.asymptotic_ring.AsymptoticRing_with_category.element_class'>
            sage: expr = 42*x^42 + x^10 + O(x^2); expr
            42*x^42 + x^10 + O(x^2)
            sage: expr.O()
            O(x^42)
            sage: (2*x).O()
            O(x)

        .. SEEALSO::

            :func:`sage.rings.power_series_ring.PowerSeriesRing`,
            :func:`sage.rings.laurent_series_ring.LaurentSeriesRing`.

        TESTS::

            sage: AR(0).O()
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large x.
        """
        if not self.summands:
            from .misc import NotImplementedOZero
            raise NotImplementedOZero(self.parent(), exact_part=self.parent().zero())
        return sum(self.parent().create_summand('O', growth=element)
                   for element in self.summands.maximal_elements())


    def log(self, base=None, precision=None, locals=None):
        r"""
        The logarithm of this asymptotic expansion.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            Computing the logarithm of an asymptotic expansion
            is possible if and only if there is exactly one maximal
            summand in the expansion.

        ALGORITHM:

        If the expansion has more than one summand,
        the asymptotic expansion for `\log(1+t)` as `t` tends to `0`
        is used.

        .. TODO::

            As soon as `L`-terms are implemented, this
            implementation has to be adapted as well in order to
            yield correct results.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ * log(x)^ZZ', coefficient_ring=QQ)
            sage: log(x)
            log(x)
            sage: log(x^2)
            2*log(x)
            sage: log(x-1)
            log(x) - x^(-1) - 1/2*x^(-2) - 1/3*x^(-3) - ... + O(x^(-21))

        The coefficient ring is automatically extended if needed::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ * log(x)^ZZ', coefficient_ring=ZZ, default_prec=3)
            sage: (49*x^3-1).log()
            3*log(x) + 2*log(7) - 1/49*x^(-3) - 1/4802*x^(-6) ... + O(x^(-12))
            sage: _.parent()
            Asymptotic Ring <x^ZZ * log(x)^ZZ> over Symbolic Ring

        If one wants to avoid this extending to the Symbolic Ring, then
        the following helps::

            sage: L.<log7> = ZZ[]
            sage: def mylog(z, base=None):
            ....:     try:
            ....:         if ZZ(z).is_power_of(7):
            ....:             return log(ZZ(z), 7) * log7
            ....:     except (TypeError, ValueError):
            ....:         pass
            ....:     return log(z, base)
            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ * log(x)^ZZ', coefficient_ring=L, default_prec=3)
            sage: (49*x^3-1).log(locals={'log': mylog})
            3*log(x) + 2*log7 - 1/49*x^(-3) - 1/4802*x^(-6) ... + O(x^(-12))

        A ``log``-function can also be specified to always be used with the
        asymptotic ring::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ * log(x)^ZZ', coefficient_ring=L, default_prec=3, locals={'log': mylog})
            sage: log(49*x^3-1)
            3*log(x) + 2*log7 - 1/49*x^(-3) - 1/4802*x^(-6) - 1/352947*x^(-9) + O(x^(-12))

        TESTS::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ * log(x)^ZZ', coefficient_ring=QQ)
            sage: log(R(1))
            0
            sage: log(R(0))
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot compute log(0) in
            Asymptotic Ring <x^ZZ * log(x)^ZZ> over Rational Field.
            sage: C.<s, t> = AsymptoticRing(growth_group='s^ZZ * t^ZZ', coefficient_ring=QQ)
            sage: log(s + t)
            Traceback (most recent call last):
            ...
            ValueError: Cannot determine main term of s + t since
            there are several maximal elements s, t.
        """
        P = self.parent()
        locals = P.locals(locals)
        log = locals['log']

        if not self.summands:
            raise ArithmeticError('Cannot compute log(0) in %s.' % (self.parent(),))

        elif len(self.summands) == 1:
            if self.is_one():
                return P.zero()
            element = next(self.summands.elements())
            return sum(P._create_element_in_extension_(l, element.parent())
                       for l in element.log_term(base=base,
                                                 locals=locals))

        (max_elem, x) = self._main_term_relative_error_()
        geom = -x

        from sage.rings.integer_ring import ZZ
        import itertools

        result = - AsymptoticExpansion._power_series_(
            coefficients=iter(1 / ZZ(k)
                              for k in itertools.count(2)),
            start=geom,
            ratio=geom,
            ratio_start=geom,
            precision=precision)

        if base:
            result = result / log(base)

        result += x.parent()(max_elem).log(base=base, locals=locals)

        return result


    def is_exact(self):
        r"""
        Return whether all terms of this expansion are exact.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('x^QQ * log(x)^QQ', QQ)
            sage: (x^2 + O(x)).is_exact()
            False
            sage: (x^2 - x).is_exact()
            True

        TESTS::

            sage: A(0).is_exact()
            True
            sage: A.one().is_exact()
            True
        """
        return all(T.is_exact() for T in self.summands)


    def is_little_o_of_one(self):
        r"""
        Return whether this expansion is of order `o(1)`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('x^ZZ * log(x)^ZZ', QQ)
            sage: (x^4 * log(x)^(-2) + x^(-4) * log(x)^2).is_little_o_of_one()
            False
            sage: (x^(-1) * log(x)^1234 + x^(-2) + O(x^(-3))).is_little_o_of_one()
            True
            sage: (log(x) - log(x-1)).is_little_o_of_one()
            True

        ::

            sage: A.<x, y> = AsymptoticRing('x^QQ * y^QQ * log(y)^ZZ', QQ)
            sage: (x^(-1/16) * y^32 + x^32 * y^(-1/16)).is_little_o_of_one()
            False
            sage: (x^(-1) * y^(-3) + x^(-3) * y^(-1)).is_little_o_of_one()
            True
            sage: (x^(-1) * y / log(y)).is_little_o_of_one()
            False
            sage: (log(y-1)/log(y) - 1).is_little_o_of_one()
            True

        .. SEEALSO::

            :meth:`limit`
        """
        return all(term.is_little_o_of_one() for term in self.summands.maximal_elements())


    def rpow(self, base, precision=None, locals=None):
        r"""
        Return the power of ``base`` to this asymptotic expansion.

        INPUT:

        - ``base`` -- an element or ``'e'``.

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: (1/x).rpow('e', precision=5)
            1 + x^(-1) + 1/2*x^(-2) + 1/6*x^(-3) + 1/24*x^(-4) + O(x^(-5))

        TESTS::

            sage: assume(SR.an_element() > 0)
            sage: y = SR.var('y')
            sage: x.rpow(y)
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot construct y^x in Growth Group x^ZZ
            > *previous* TypeError: unsupported operand parent(s) for *:
              'Growth Group x^ZZ' and 'Growth Group SR^x * Arg_SR^x'
            sage: assume(y > 0)
            sage: x.rpow(y)
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot construct y^x in Growth Group x^ZZ
            > *previous* TypeError: unsupported operand parent(s) for *:
              'Growth Group x^ZZ' and 'Growth Group SR^x'
            sage: forget()

        Check that :trac:`19946` is fixed::

            sage: A.<n> = AsymptoticRing('(QQ_+)^n * n^QQ', SR)
            sage: n.rpow(2)
            2^n
            sage: _.parent()
            Asymptotic Ring <QQ^n * n^QQ> over Symbolic Ring
        """
        if isinstance(base, AsymptoticExpansion):
            return base.__pow__(self, precision=precision)

        P = self.parent()

        # first: remove terms from a copy of this term such that a
        # term in o(1) remains

        expr_o = self.summands.copy()
        large_terms = []
        for term in self.summands.elements_topological():
            if not term.is_little_o_of_one():
                large_terms.append(term)
                expr_o.remove(term.growth)

        expr_o = P(expr_o)

        # next: try to take the exponential function of the large elements

        try:
            large_result = P.prod(
                P._create_element_in_extension_(term.rpow(base),
                                                term.parent())
                for term in large_terms)
        except (TypeError, ValueError) as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                ValueError('Cannot construct the power of %s to the '
                           'exponent %s in %s.' %
                           (base, self, self.parent())), e)

        # then: expand expr_o

        if not expr_o:
            return large_result


        if base == 'e':
            geom = expr_o
        else:
            log = self.parent().locals(locals)['log']
            geom = expr_o * log(base)
        P = geom.parent()

        from sage.rings.integer_ring import ZZ
        import itertools

        def inverted_factorials():
            f = ZZ(1)
            for k in itertools.count(1):
                f /= ZZ(k)
                yield f

        result = AsymptoticExpansion._power_series_(
            coefficients=inverted_factorials(),
            start=P.one(),
            ratio=geom,
            ratio_start=P.one(),
            precision=precision)

        return result * large_result


    def _main_term_relative_error_(self, return_inverse_main_term=False):
        r"""
        Split this asymptotic expansion into `m(1+x)` with `x=o(1)`.

        INPUT:

        - ``return_inverse_main_term`` -- (default: ``False``) a boolean.
          If set, then the pair `(m^{-1},x)` is returned instead of `(m,x)`.

        OUTPUT:

        A pair (``m``, ``x``) consisting of
        a :mod:`term <sage.rings.asymptotic.term_monoid>` ``m`` and
        an :class:`asymptotic expansion <AsymptoticExpansion>` ``x``.

        EXAMPLES::

            sage: R.<n> = AsymptoticRing('n^ZZ', QQ)
            sage: ex = 2*n^2 + n + O(1/n)
            sage: (m, x) = ex._main_term_relative_error_()
            sage: m
            2*n^2
            sage: x
            1/2*n^(-1) + O(n^(-3))
            sage: ex = 2*n^2 + n
            sage: (m, x) = ex._main_term_relative_error_()
            sage: m
            2*n^2
            sage: x
            1/2*n^(-1)
            sage: ex._main_term_relative_error_(return_inverse_main_term=True)
            (1/2*n^(-2), 1/2*n^(-1))
            sage: R(0)._main_term_relative_error_()
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot determine main term of 0.

        TESTS::

            sage: R.<m, n> = AsymptoticRing('n^ZZ*m^ZZ', QQ)
            sage: (m + n)._main_term_relative_error_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot determine main term of m + n since
            there are several maximal elements m, n.
        """
        if not self.summands:
            raise ArithmeticError("Cannot determine main term of 0.")

        max_elem = tuple(self.summands.maximal_elements())
        if len(max_elem) != 1:
            raise ValueError('Cannot determine main term of {} since there '
                             'are several maximal elements {}.'.format(
                             self, ', '.join(str(e) for e in
                                              sorted(max_elem, key=str))))
        max_elem = max_elem[0]

        imax_elem = ~max_elem
        if imax_elem.parent() is max_elem.parent():
            new_self = self
        else:
            new_self = self.parent()._create_element_in_extension_(
                imax_elem, max_elem.parent()).parent()(self)

        one = new_self.parent().one()
        x = - one + new_self._mul_term_(imax_elem)

        if return_inverse_main_term:
            return (imax_elem, x)
        else:
            return (max_elem, x)


    @staticmethod
    def _power_series_(coefficients, start, ratio, ratio_start, precision):
        r"""
        Return a taylor series.

        Let `c_k` be determined by the ``coefficients`` and set

        .. MATH::

            s_k = c_k \cdot \mathit{ratio\_start} \cdot \mathit{ratio}^k.

        The result is

        .. MATH::

            \mathit{start} + \sum_{k=1}^K s_k

        where `K` is chosen such that adding `s_{K+1}` does not change
        the result.

        INPUT:

        - ``coefficients`` -- an iterator.

        - ``start`` -- an asymptotic expansion.

        - ``ratio`` -- an asymptotic expansion.

        - ``ratio_start`` -- an asymptotic expansion.

        - ``precision`` -- a non-negative integer. All intermediate
          results are truncated to this precision.

        OUTPUT:

        An asymptotic expansion.

        TESTS::

            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticExpansion
            sage: from itertools import count
            sage: A.<g> = AsymptoticRing('g^ZZ', QQ)
            sage: AsymptoticExpansion._power_series_(
            ....:     coefficients=iter(ZZ(k) for k in count(1)),
            ....:     start=A(42),
            ....:     ratio=1/g,
            ....:     ratio_start=A(5),
            ....:     precision=4)
            42 + 5*g^(-1) + 10*g^(-2) + 15*g^(-3) + O(g^(-4))
            sage: AsymptoticExpansion._power_series_(
            ....:     coefficients=iter(ZZ(k) for k in count(1)),
            ....:     start=A(42),
            ....:     ratio=1/g+O(1/g^2),
            ....:     ratio_start=A(5),
            ....:     precision=4)
            42 + 5*g^(-1) + O(g^(-2))
            sage: AsymptoticExpansion._power_series_(
            ....:     coefficients=iter(ZZ(k) for k in count(1)),
            ....:     start=A(42),
            ....:     ratio=1/g+O(1/g^2),
            ....:     ratio_start=A(5),
            ....:     precision=1000000)
            42 + 5*g^(-1) + O(g^(-2))
        """
        result = start
        g = ratio_start
        for c in coefficients:
            g *= ratio
            new_result = (result + c*g).truncate(precision=precision)
            if new_result.has_same_summands(result):
                break
            result = new_result
        return result


    def exp(self, precision=None):
        r"""
        Return the exponential of (i.e., the power of `e` to) this asymptotic expansion.

        INPUT:

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            The exponential function of this expansion can only be
            computed exactly if the respective growth element can be
            constructed in the underlying growth group.

        ALGORITHM:

        If the corresponding growth can be constructed, return
        the exact exponential function. Otherwise, if this term
        is `o(1)`, try to expand the series and truncate
        according to the given precision.

        .. TODO::

            As soon as `L`-terms are implemented, this
            implementation has to be adapted as well in order to
            yield correct results.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('(e^x)^ZZ * x^ZZ * log(x)^ZZ', SR)
            sage: exp(x)
            e^x
            sage: exp(2*x)
            (e^x)^2
            sage: exp(x + log(x))
            e^x*x

        ::

            sage: (x^(-1)).exp(precision=7)
            1 + x^(-1) + 1/2*x^(-2) + 1/6*x^(-3) + ... + O(x^(-7))

        TESTS::

            sage: A.<x> = AsymptoticRing('(e^x)^ZZ * x^QQ * log(x)^QQ', SR)
            sage: exp(log(x))
            x
            sage: log(exp(x))
            x

        ::

            sage: exp(x+1)
            e*e^x

        See :trac:`19521`::

            sage: A.<n> = AsymptoticRing('n^ZZ', SR.subring(no_variables=True))
            sage: exp(O(n^(-3))).parent()
            Asymptotic Ring <n^ZZ> over Symbolic Constants Subring
        """
        return self.rpow('e', precision=precision)


    def substitute(self, rules=None, domain=None, **kwds):
        r"""
        Substitute the given ``rules`` in this asymptotic expansion.

        INPUT:

        - ``rules`` -- a dictionary.

        - ``kwds`` -- keyword arguments will be added to the
          substitution ``rules``.

        - ``domain`` -- (default: ``None``) a parent. The neutral
          elements `0` and `1` (rules for the keys ``'_zero_'`` and
          ``'_one_'``, see note box below) are taken out of this
          domain. If ``None``, then this is determined automatically.

        OUTPUT:

        An object.

        .. NOTE::

          The neutral element of the asymptotic ring is replaced by
          the value to the key ``'_zero_'``; the neutral element of
          the growth group is replaced by the value to the key
          ``'_one_'``.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing(growth_group='(e^x)^QQ * x^ZZ * log(x)^ZZ', coefficient_ring=QQ, default_prec=5)

        ::

            sage: (e^x * x^2 + log(x)).subs(x=SR('s'))
            s^2*e^s + log(s)
            sage: _.parent()
            Symbolic Ring

        ::

            sage: (x^3 + x + log(x)).subs(x=x+5).truncate(5)
            x^3 + 15*x^2 + 76*x + log(x) + 130 + O(x^(-1))
            sage: _.parent()
            Asymptotic Ring <(e^x)^QQ * x^ZZ * log(x)^ZZ> over Rational Field

        ::

            sage: (e^x * x^2 + log(x)).subs(x=2*x)
            4*(e^x)^2*x^2 + log(x) + log(2)
            sage: _.parent()
            Asymptotic Ring <(e^x)^QQ * x^QQ * log(x)^QQ> over Symbolic Ring

        ::

            sage: (x^2 + log(x)).subs(x=4*x+2).truncate(5)
            16*x^2 + 16*x + log(x) + 2*log(2) + 4 + 1/2*x^(-1) + O(x^(-2))
            sage: _.parent()
            Asymptotic Ring <(e^x)^QQ * x^ZZ * log(x)^ZZ> over Symbolic Ring

        ::

            sage: (e^x * x^2 + log(x)).subs(x=RIF(pi))
            229.534211738584?
            sage: _.parent()
            Real Interval Field with 53 bits of precision

        .. SEEALSO::

            :meth:`sage.symbolic.expression.Expression.subs`

        TESTS::

            sage: x.subs({'y': -1})
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute y in x since it is not a generator of
            Asymptotic Ring <(e^x)^QQ * x^ZZ * log(x)^ZZ> over Rational Field.
            sage: B.<u, v, w> = AsymptoticRing(growth_group='u^QQ * v^QQ * w^QQ', coefficient_ring=QQ)
            sage: (1/u).subs({'u': 0})
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply the substitution rules {u: 0} on u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            > *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            >> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Exact Term Monoid u^QQ * v^QQ * w^QQ with coefficients in Rational Field.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ * v^QQ * w^QQ.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ.
            >...> *previous* ZeroDivisionError: rational division by zero
            sage: (1/u).subs({'u': 0, 'v': SR.var('v')})
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply the substitution rules {u: 0, v: v} on u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            > *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            >> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Exact Term Monoid u^QQ * v^QQ * w^QQ with coefficients in Rational Field.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ * v^QQ * w^QQ.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ.
            >...> *previous* ZeroDivisionError: rational division by zero

        ::

            sage: u.subs({u: 0, 'v': SR.var('v')})
            0
            sage: v.subs({u: 0, 'v': SR.var('v')})
            v
            sage: _.parent()
            Symbolic Ring

        ::

            sage: u.subs({SR.var('u'): -1})
            Traceback (most recent call last):
            ...
            TypeError: Cannot substitute u in u since it is neither an
            asymptotic expansion nor a string
            (but a <class 'sage.symbolic.expression.Expression'>).

        ::

            sage: u.subs({u: 1, 'u': 1})
            1
            sage: u.subs({u: 1}, u=1)
            1
            sage: u.subs({u: 1, 'u': 2})
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute in u: duplicate key u.
            sage: u.subs({u: 1}, u=3)
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute in u: duplicate key u.

        ::

            sage: B(0).subs({'_zero_': None}) is None
            True
            sage: B(1).subs({'_one_': AA(1)}).parent() is AA
            True
        """
        # check if nothing to do
        if not rules and not kwds:
            return self

        # init and process keyword arguments
        gens = self.parent().gens()
        locals = kwds or dict()

        # update with rules
        if isinstance(rules, dict):
            for k, v in rules.items():
                if not isinstance(k, str) and k not in gens:
                    raise TypeError('Cannot substitute %s in %s '
                                    'since it is neither an '
                                    'asymptotic expansion '
                                    'nor a string (but a %s).' %
                                    (k, self, type(k)))
                k = str(k)
                if k in locals and locals[k] != v:
                    raise ValueError('Cannot substitute in %s: '
                                     'duplicate key %s.' % (self, k))
                locals[k] = v
        elif rules is not None:
            raise TypeError('Substitution rules %s have to be a dictionary.' %
                            (rules,))

        # fill up missing rules
        for g in gens:
            locals.setdefault(str(g), g)

        # check if all keys are generators
        gens_str = tuple(str(g) for g in gens)
        for k in locals:
            sk = str(k)
            if sk not in gens_str and not sk.startswith('_'):
                raise ValueError('Cannot substitute %s in %s '
                                 'since it is not a generator of %s.' %
                                 (k, self, self.parent()))

        # determine 0 and 1
        if domain is None and \
               ('_zero_' not in locals or '_one_' not in locals):
            P = self.parent()
            for g in gens:
                G = locals[str(g)].parent()
                if G is not P:
                    domain = G
                    break
            else:
                domain = P
        locals.setdefault('_zero_', domain.zero())
        locals.setdefault('_one_', domain.one())

        # do the actual substitution
        try:
            return self._substitute_(locals)
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import combine_exceptions
            rules = '{' + ', '.join(
                '%s: %s' % (k, v)
                for k, v in sorted(locals.items(),
                                   key=lambda k: str(k[0]))
                if not k.startswith('_') and
                not any(k == str(g) and v is g for g in gens)) + '}'
            raise combine_exceptions(
                TypeError('Cannot apply the substitution rules %s on %s '
                          'in %s.' % (rules, self, self.parent())), e)

    subs = substitute

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this asymptotic expansion.

        INPUT:

        - ``rules`` -- a dictionary.
          The neutral element of the asymptotic ring is replaced by the value
          to key ``'_zero_'``.

        OUTPUT:

        An object.

        TESTS::

            sage: A.<z> = AsymptoticRing(growth_group='z^QQ', coefficient_ring=QQ)
            sage: z._substitute_({'z': SR.var('a')})
            a
            sage: _.parent()
            Symbolic Ring
            sage: A(0)._substitute_({'_zero_': 'zero'})
            'zero'
            sage: (1/z)._substitute_({'z': 4})
            1/4
            sage: _.parent()
            Rational Field
            sage: (1/z)._substitute_({'z': 0})
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot substitute in z^(-1) in
            Asymptotic Ring <z^QQ> over Rational Field.
            > *previous* ZeroDivisionError: Cannot substitute in z^(-1) in
            Exact Term Monoid z^QQ with coefficients in Rational Field.
            >> *previous* ZeroDivisionError: Cannot substitute in z^(-1) in
            Growth Group z^QQ.
            >...> *previous* ZeroDivisionError: rational division by zero
        """
        if not self.summands:
            return rules['_zero_']
        from sage.symbolic.operators import add_vararg
        try:
            return add_vararg(
                *tuple(s._substitute_(rules)
                       for s in self.summands.elements_topological()))
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)


    def compare_with_values(self, variable, function, values,
                            rescaled=True, ring=RIF):
        """
        Compute the (rescaled) difference between this asymptotic
        expansion and the given values.

        INPUT:

        - ``variable`` -- an asymptotic expansion or a string.

        - ``function`` -- a callable or symbolic expression giving the
          comparison values.

        - ``values`` -- a list or iterable of values where the comparison
          shall be carried out.

        - ``rescaled`` -- (default: ``True``) determines whether
          the difference is divided by the error term of the asymptotic
          expansion.

        - ``ring`` -- (default: ``RIF``) the parent into which the
          difference is converted.

        OUTPUT:

        A list of pairs containing comparison points and (rescaled)
        difference values.

        EXAMPLES::

            sage: assume(SR.an_element() > 0)
            sage: A.<n> = AsymptoticRing('QQ^n * n^ZZ', SR)
            sage: catalan = binomial(2*x, x)/(x+1)
            sage: expansion = 4^n*(1/sqrt(pi)*n^(-3/2)
            ....:     - 9/8/sqrt(pi)*n^(-5/2)
            ....:     + 145/128/sqrt(pi)*n^(-7/2) + O(n^(-9/2)))
            sage: expansion.compare_with_values(n, catalan, srange(5, 10))
            [(5, 0.5303924444775?),
             (6, 0.5455279498787?),
             (7, 0.556880411050?),
             (8, 0.565710587724?),
             (9, 0.572775029098?)]
            sage: expansion.compare_with_values(n, catalan, [5, 10, 20], rescaled=False)
            [(5, 0.3886263699387?), (10, 19.1842458318?), (20, 931314.63637?)]
            sage: expansion.compare_with_values(n, catalan, [5, 10, 20], rescaled=False, ring=SR)
            [(5, 168/5*sqrt(5)/sqrt(pi) - 42),
             (10, 1178112/125*sqrt(10)/sqrt(pi) - 16796),
             (20, 650486218752/125*sqrt(5)/sqrt(pi) - 6564120420)]

        Instead of a symbolic expression, a callable function can
        be specified as well::

            sage: A.<n> = AsymptoticRing('n^ZZ * log(n)^ZZ', SR)
            sage: def H(n):
            ....:     return sum(1/k for k in srange(1, n+1))
            sage: H_expansion = (log(n) + euler_gamma + 1/(2*n)
            ....:                - 1/(12*n^2) + O(n^-4))
            sage: H_expansion.compare_with_values(n, H, srange(25, 30)) # rel tol 1e-6
            [(25, -0.008326995?),
             (26, -0.008327472?),
             (27, -0.008327898?),
             (28, -0.00832828?),
             (29, -0.00832862?)]
            sage: forget()

        .. SEEALSO::

            :meth:`plot_comparison`

        TESTS::

            sage: A.<x, y> = AsymptoticRing('x^ZZ*y^ZZ', QQ)
            sage: expansion = x^2 + O(x) + O(y)
            sage: expansion.compare_with_values(y, lambda z: z^2, srange(20, 30))
            Traceback (most recent call last):
            ....
            NotImplementedError: exactly one error term required
            sage: expansion = x^2
            sage: expansion.compare_with_values(y, lambda z: z^2, srange(20, 30))
            Traceback (most recent call last):
            ....
            NotImplementedError: exactly one error term required
            sage: expansion = x^2 + O(x)
            sage: expansion.compare_with_values(y, lambda z: z^2, srange(20, 30))
            Traceback (most recent call last):
            ....
            NameError: name 'x' is not defined
            sage: expansion.compare_with_values(x, lambda z: z^2, srange(20, 30))
            [(20, 0), (21, 0), ..., (29, 0)]
            sage: expansion.compare_with_values(x, SR('x*y'), srange(20, 30))
            Traceback (most recent call last):
            ....
            NotImplementedError: expression x*y has more than one variable
        """
        from .term_monoid import OTerm
        from sage.rings.integer_ring import ZZ

        main = self.exact_part()
        error = self - main
        error_terms = list(error.summands)
        if len(error_terms) != 1:
            raise NotImplementedError("exactly one error term required")
        if not isinstance(error_terms[0], OTerm):
            raise NotImplementedError("{} is not an O term".format(error))
        error_growth = error_terms[0].growth

        if hasattr(function, 'variables'):
            expr = function
            vars = expr.variables()
            if len(vars) > 1:
                raise NotImplementedError("expression {} has more than one "
                                          "variable".format(expr))
            elif len(vars) == 1:
                v = vars[0]
                def function(arg):
                    return expr.subs({v: arg})
            else:
                def function(arg):
                    return expr

        if rescaled:
            points = list(
                (k, ring((main.subs({variable: k}) - function(k)) /
                         error_growth._substitute_({str(variable): k,
                                                    '_one_': ZZ(1)})))
                for k in values)
        else:
            points = list(
                (k, ring(main.subs({variable: k}) - function(k)))
                for k in values)

        return points


    def plot_comparison(self, variable, function, values, rescaled=True,
                        ring=RIF, relative_tolerance=0.025, **kwargs):
        r"""
        Plot the (rescaled) difference between this asymptotic
        expansion and the given values.

        INPUT:

        - ``variable`` -- an asymptotic expansion or a string.

        - ``function`` -- a callable or symbolic expression giving the
          comparison values.

        - ``values`` -- a list or iterable of values where the comparison
          shall be carried out.

        - ``rescaled`` -- (default: ``True``) determines whether
          the difference is divided by the error term of the asymptotic
          expansion.

        - ``ring`` -- (default: ``RIF``) the parent into which the
          difference is converted.

        - ``relative_tolerance`` -- (default: ``0.025``). Raise error
          when relative error exceeds this tolerance.

        Other keyword arguments are passed to :func:`list_plot`.

        OUTPUT:

        A graphics object.

        .. NOTE::

            If rescaled (i.e. divided by the error term), the output
            should be bounded.

            This method is mainly meant to have an easily usable
            plausibility check for asymptotic expansion created in
            some way.

        EXAMPLES:

        We want to check the quality of the asymptotic expansion of
        the harmonic numbers::

            sage: A.<n> = AsymptoticRing('n^ZZ * log(n)^ZZ', SR)
            sage: def H(n):
            ....:     return sum(1/k for k in srange(1, n+1))
            sage: H_expansion = (log(n) + euler_gamma + 1/(2*n)
            ....:                - 1/(12*n^2) + O(n^-4))
            sage: H_expansion.plot_comparison(n, H, srange(1, 30))
            Graphics object consisting of 1 graphics primitive

        Alternatively, the unscaled (absolute) difference can be
        plotted as well::

            sage: H_expansion.plot_comparison(n, H, srange(1, 30),
            ....:                             rescaled=False)
            Graphics object consisting of 1 graphics primitive

        Additional keywords are passed to :func:`list_plot`::

            sage: H_expansion.plot_comparison(n, H, srange(1, 30),
            ....:                             plotjoined=True, marker='o',
            ....:                             color='green')
            Graphics object consisting of 1 graphics primitive

        .. SEEALSO::

            :meth:`compare_with_values`

        TESTS::

            sage: H_expansion.plot_comparison(n, H, [600])
            Traceback (most recent call last):
            ...
            ValueError: Numerical noise is too high, the comparison is inaccurate
            sage: H_expansion.plot_comparison(n, H, [600], relative_tolerance=2)
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.plot import list_plot
        points = self.compare_with_values(variable, function,
                                          values, rescaled=rescaled, ring=ring)

        if isinstance(ring, sage.rings.abc.RealIntervalField):
            if not all(p[1].relative_diameter() <= relative_tolerance for p in points):
                raise ValueError('Numerical noise is too high, the '
                                 'comparison is inaccurate')

            # RIFs cannot be plotted, they need to be converted to RR
            # (see #15011).
            points = [(p[0], p[1].center()) for p in points]

        return list_plot(points, **kwargs)


    def symbolic_expression(self, R=None):
        r"""
        Return this asymptotic expansion as a symbolic expression.

        INPUT:

        - ``R`` -- (a subring of) the symbolic ring or ``None``.
          The output will be an element of ``R``. If ``None``,
          then the symbolic ring is used.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: A.<x, y, z> = AsymptoticRing(growth_group='x^ZZ * y^QQ * log(y)^QQ * (QQ_+)^z * z^QQ', coefficient_ring=QQ)
            sage: SR(A.an_element())  # indirect doctest
            1/8*(1/8)^z*x^3*y^(3/2)*z^(3/2)*log(y)^(3/2) +
            Order((1/2)^z*x*sqrt(y)*sqrt(z)*sqrt(log(y)))

            sage: A.<x, y, z> = AsymptoticRing(growth_group='x^ZZ * y^QQ * log(y)^QQ * (QQ_+)^z * z^QQ', coefficient_ring=QQ)
            sage: SR(A.an_element())  # indirect doctest
            1/8*(1/8)^z*x^3*y^(3/2)*z^(3/2)*log(y)^(3/2) +
            Order((1/2)^z*x*sqrt(y)*sqrt(z)*sqrt(log(y)))

        TESTS::

            sage: a = A.an_element(); a
            1/8*x^3*y^(3/2)*log(y)^(3/2)*(1/8)^z*z^(3/2) +
            O(x*y^(1/2)*log(y)^(1/2)*(1/2)^z*z^(1/2))
            sage: a.symbolic_expression()
            1/8*(1/8)^z*x^3*y^(3/2)*z^(3/2)*log(y)^(3/2) +
            Order((1/2)^z*x*sqrt(y)*sqrt(z)*sqrt(log(y)))
            sage: _.parent()
            Symbolic Ring

        ::

            sage: from sage.symbolic.ring import SymbolicRing
            sage: class MySymbolicRing(SymbolicRing):
            ....:     pass
            sage: mySR = MySymbolicRing()
            sage: a.symbolic_expression(mySR).parent() is mySR
            True
        """
        if R is None:
            from sage.symbolic.ring import SR
            R = SR

        return self.substitute(dict((g, R(R.var(str(g))))
                                    for g in self.parent().gens()),
                               domain=R)


    _symbolic_ = symbolic_expression  # will be used by SR._element_constructor_


    def map_coefficients(self, f, new_coefficient_ring=None):
        r"""
        Return the asymptotic expansion obtained by applying ``f`` to
        each coefficient of this asymptotic expansion.

        INPUT:

        - ``f`` -- a callable. A coefficient `c` will be mapped to `f(c)`.

        - ``new_coefficient_ring`` -- (default: ``None``) a ring.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: A.<n> = AsymptoticRing(growth_group='n^ZZ', coefficient_ring=ZZ)
            sage: a = n^4 + 2*n^3 + 3*n^2 + O(n)
            sage: a.map_coefficients(lambda c: c+1)
            2*n^4 + 3*n^3 + 4*n^2 + O(n)
            sage: a.map_coefficients(lambda c: c-2)
            -n^4 + n^2 + O(n)

        TESTS::

            sage: a.map_coefficients(lambda c: 1/c, new_coefficient_ring=QQ)
            n^4 + 1/2*n^3 + 1/3*n^2 + O(n)
            sage: _.parent()
            Asymptotic Ring <n^ZZ> over Rational Field
            sage: a.map_coefficients(lambda c: 1/c)
            Traceback (most recent call last):
            ...
            ValueError: Cannot create ExactTerm(n^3) since
            given coefficient 1/2 is not valid in
            Exact Term Monoid n^ZZ with coefficients in Integer Ring.
            > *previous* TypeError: no conversion of this rational to integer
        """
        def mapping(term):
            T = term.parent().change_parameter(
                coefficient_ring=new_coefficient_ring)
            if hasattr(term, 'coefficient'):
                c = f(term.coefficient)
                if c.is_zero():
                    return None
                return T(term.growth, coefficient=c)
            else:
                return T(term.growth)

        P = self.parent().change_parameter(coefficient_ring=new_coefficient_ring)
        S = self.summands.copy()
        S.map(mapping)
        return P(S, simplify=False, convert=False)


    def factorial(self):
        r"""
        Return the factorial of this asymptotic expansion.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: A.<n> = AsymptoticRing(growth_group='n^ZZ * log(n)^ZZ', coefficient_ring=ZZ, default_prec=5)
            sage: n.factorial()
            sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(1/2)
            + 1/12*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-1/2)
            + 1/288*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-3/2)
            + O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: _.parent()
            Asymptotic Ring <(e^(n*log(n)))^QQ * (e^n)^QQ * n^QQ * log(n)^QQ>
            over Symbolic Constants Subring

        :wikipedia:`Catalan numbers <Catalan_number>`
        `\frac{1}{n+1}\binom{2n}{n}`::

            sage: (2*n).factorial() / n.factorial()^2 / (n+1)  # long time
            1/sqrt(pi)*(e^n)^(2*log(2))*n^(-3/2)
            - 9/8/sqrt(pi)*(e^n)^(2*log(2))*n^(-5/2)
            + 145/128/sqrt(pi)*(e^n)^(2*log(2))*n^(-7/2)
            + O((e^n)^(2*log(2))*n^(-9/2))

        Note that this method substitutes the asymptotic expansion into
        Stirling's formula. This substitution has to be possible which is
        not always guaranteed::

            sage: S.<s> = AsymptoticRing(growth_group='s^QQ * log(s)^QQ', coefficient_ring=QQ, default_prec=4)
            sage: log(s).factorial()
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply the substitution rules {s: log(s)} on
            sqrt(2)*sqrt(pi)*e^(s*log(s))*(e^s)^(-1)*s^(1/2)
            + O(e^(s*log(s))*(e^s)^(-1)*s^(-1/2)) in
            Asymptotic Ring <(e^(s*log(s)))^QQ * (e^s)^QQ * s^QQ * log(s)^QQ>
            over Symbolic Constants Subring.
            ...

        .. SEEALSO::

            :meth:`~sage.rings.asymptotic.asymptotic_expansion_generators.AsymptoticExpansionGenerators.Stirling`

        TESTS::

            sage: A.<m> = AsymptoticRing(growth_group='m^ZZ * log(m)^ZZ', coefficient_ring=QQ, default_prec=5)
            sage: m.factorial()
            sqrt(2)*sqrt(pi)*e^(m*log(m))*(e^m)^(-1)*m^(1/2)
            + 1/12*sqrt(2)*sqrt(pi)*e^(m*log(m))*(e^m)^(-1)*m^(-1/2)
            + 1/288*sqrt(2)*sqrt(pi)*e^(m*log(m))*(e^m)^(-1)*m^(-3/2)
            + O(e^(m*log(m))*(e^m)^(-1)*m^(-5/2))

        ::

            sage: A(1/2).factorial()
            1/2*sqrt(pi)
            sage: _.parent()
            Asymptotic Ring <m^ZZ * log(m)^ZZ> over Symbolic Ring

        ::

            sage: B.<a, b> = AsymptoticRing('a^ZZ * b^ZZ', QQ, default_prec=3)
            sage: b.factorial()
            O(e^(b*log(b))*(e^b)^(-1)*b^(1/2))
            sage: (a*b).factorial()
            Traceback (most recent call last):
            ...
            ValueError: Cannot build the factorial of a*b
            since it is not univariate.
        """
        vars = self.variable_names()

        if len(vars) == 0:
            if self.is_zero():
                return self.parent().one()
            assert len(self.summands) == 1
            element = next(self.summands.elements())
            return self.parent()._create_element_in_extension_(
                element._factorial_(), element.parent())

        if len(vars) == 1:
            from .asymptotic_expansion_generators import \
                asymptotic_expansions
            var = vars[0]
            S = asymptotic_expansions.Stirling(
                var, precision=self.parent().default_prec)
            from sage.structure.element import get_coercion_model
            cm = get_coercion_model()
            P = cm.common_parent(self, S)
            return S.subs({var: P.coerce(self)})

        else:
            raise ValueError(
                'Cannot build the factorial of {} since it is not '
                'univariate.'.format(self))


    def variable_names(self):
        r"""
        Return the names of the variables of this asymptotic expansion.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: A.<m, n> = AsymptoticRing('QQ^m * m^QQ * n^ZZ * log(n)^ZZ', QQ)
            sage: (4*2^m*m^4*log(n)).variable_names()
            ('m', 'n')
            sage: (4*2^m*m^4).variable_names()
            ('m',)
            sage: (4*log(n)).variable_names()
            ('n',)
            sage: (4*m^3).variable_names()
            ('m',)
            sage: (4*m^0).variable_names()
            ()
            sage: (4*2^m*m^4 + log(n)).variable_names()
            ('m', 'n')
            sage: (2^m + m^4 + log(n)).variable_names()
            ('m', 'n')
            sage: (2^m + m^4).variable_names()
            ('m',)
        """
        vars = sorted(sum(iter(s.variable_names()
                               for s in self.summands),
                          tuple()))
        from itertools import groupby
        return tuple(v for v, _ in groupby(vars))


    def _singularity_analysis_(self, var, zeta, precision=None):
        r"""
        Return the asymptotic growth of the coefficients of some
        generating function having this singular expansion around `\zeta`.

        INPUT:

        - ``var`` -- a string, the variable for the growth of the coefficients,
          or the generator of an asymptotic ring.

        - ``zeta`` -- location of the singularity

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the parent of this expansion is used.

        OUTPUT:

        An asymptotic expansion in ``var``.

        EXAMPLES::

            sage: C.<T> = AsymptoticRing('T^QQ', QQ)
            sage: ex = 2 - 2*T^(-1/2) + 2*T^(-1) - 2*T^(-3/2) + O(T^(-2))
            sage: ex._singularity_analysis_('n', 1/4, precision=2)
            1/sqrt(pi)*4^n*n^(-3/2) - 9/8/sqrt(pi)*4^n*n^(-5/2) + O(4^n*n^(-3))

        The parameter ``var`` can also be the generator of an asymptotic
        ring::

            sage: A.<n> = AsymptoticRing('n^QQ', QQ)
            sage: ex._singularity_analysis_(n, 1/4, precision=2)
            1/sqrt(pi)*4^n*n^(-3/2) - 9/8/sqrt(pi)*4^n*n^(-5/2) + O(4^n*n^(-3))

        If the parameter ``precision`` is omitted, the default precision
        of the parent of this expansion is used. ::

            sage: C.<T> = AsymptoticRing('T^QQ', QQ, default_prec=1)
            sage: ex = 2 - 2*T^(-1/2) + 2*T^(-1) - 2*T^(-3/2) + O(T^(-2))
            sage: ex._singularity_analysis_('n', 1/4)
            1/sqrt(pi)*4^n*n^(-3/2)  + O(4^n*n^(-5/2))

        .. SEEALSO::

            :meth:`AsymptoticRing.coefficients_of_generating_function`

        .. WARNING::

            Once singular expansions around points other than infinity
            are implemented (:trac:`20050`), this method will be
            renamed to ``singularity_analysis``, the parameter
            ``zeta`` will be dropped (as it will be part of the
            singular expansion) and expansions around infinity will no
            longer be accepted.

        TESTS::

            sage: C.<T> = AsymptoticRing('T^QQ', QQ)
            sage: (1/T)._singularity_analysis_('n', 1)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large n.
        """
        from .misc import NotImplementedOZero
        OZeroEncountered = False

        if precision is None:
            precision = self.parent().default_prec

        result = 0
        for s in self.summands:
            try:
                contribution = s._singularity_analysis_(
                    var=var, zeta=zeta,
                    precision=precision)
            except NotImplementedOZero as ozero:
                OZeroEncountered = True
                result += ozero.exact_part
            else:
                result += contribution

        if OZeroEncountered and (isinstance(result, int) and result == 0
                                 or result.is_exact()):
            raise NotImplementedOZero(var=var, exact_part=result)
        return result

    def limit(self):
        """
        Compute the limit of this asymptotic expansion.

        OUTPUT:

        An element of the coefficient ring.

        EXAMPLES::

            sage: A.<s> = AsymptoticRing("s^ZZ", SR, default_prec=3)
            sage: (3 + 1/s + O(1/s^2)).limit()
            3
            sage: ((1+1/s)^s).limit()
            e
            sage: (1/s).limit()
            0
            sage: (s + 3 + 1/s + O(1/s^2)).limit()
            Traceback (most recent call last):
            ...
            ValueError: Cannot determine limit of s + 3 + s^(-1) + O(s^(-2))
            sage: (O(s^0)).limit()
            Traceback (most recent call last):
            ...
            ValueError: Cannot determine limit of O(1)

        .. SEEALSO::

            :meth:`is_little_o_of_one`
        """
        non_o_one_terms = list(
            term for term in self.summands
            if not term.is_little_o_of_one()
        )
        if not non_o_one_terms:
            return self.parent().base_ring()(0)
        elif (
            len(non_o_one_terms) == 1
            and non_o_one_terms[0].growth.is_one()
            and non_o_one_terms[0].is_exact()
            ):
            return non_o_one_terms[0].coefficient
        else:
            raise ValueError("Cannot determine limit of {}".format(self))

    def B(self, valid_from=0):
        r"""
        Convert all terms in this asymptotic expansion to `B`-terms.

        INPUT:

        - ``valid_from`` -- dictionary mapping variable names to lower bounds
          for the corresponding variable. The bound implied by this term is valid when
          all variables are at least their corresponding lower bound. If a number
          is passed to ``valid_from``, then the lower bounds for all variables of
          the asymptotic expansion are set to this number

        OUTPUT:

        An asymptotic expansion

        EXAMPLES::

            sage: AR.<x, z> = AsymptoticRing(growth_group='x^ZZ * z^ZZ', coefficient_ring=ZZ)
            sage: AR.B(2*x^2, {x: 10}) # indirect doctest
            doctest:warning
            ...
            FutureWarning: This class/method/function is marked as experimental.
            It, its functionality or its interface might change without a formal deprecation.
            See https://trac.sagemath.org/31922 for details.
            B(2*x^2, x >= 10)
            sage: expr = 42*x^42 + x^10 + AR.B(x^2, 20); expr # indirect doctest
            42*x^42 + x^10 + B(x^2, x >= 20, z >= 20)
            sage: type(AR.B(x, 10)) # indirect doctest
            <class 'sage.rings.asymptotic.asymptotic_ring.AsymptoticRing_with_category.element_class'>
            sage: 2*z^3 + AR.B(5*z^2, {z: 20}) # indirect doctest
            2*z^3 + B(5*z^2, z >= 20)
            sage: (2*x).B({x: 20})
            B(2*x, x >= 20)
            sage: AR.B(4*x^2*z^3, valid_from=10) # indirect doctest
            B(4*x^2*z^3, x >= 10, z >= 10)
            sage: AR.B(42*x^2) # indirect doctest
            B(42*x^2, x >= 0, z >= 0)

        TESTS::
            sage: AR(0).B(20) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedBZero: got B(0)
            The error term B(0) means 0 for sufficiently large x, z.
        """
        if not self.summands:
            from .misc import NotImplementedBZero
            raise NotImplementedBZero(self.parent(), exact_part=self.parent().zero())
        return sum(self.parent().create_summand('B', growth=element, valid_from=valid_from)
                   for element in self.summands.elements())


class AsymptoticRing(Algebra, UniqueRepresentation, WithLocals):
    r"""
    A ring consisting of :class:`asymptotic expansions <AsymptoticExpansion>`.

    INPUT:

    - ``growth_group`` -- either a partially ordered group (see
      :doc:`growth_group`) or a string
      describing such a growth group (see
      :class:`~sage.rings.asymptotic.growth_group.GrowthGroupFactory`).

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the expansions.

    - ``default_prec`` -- a positive integer. This is the number of
      summands that are kept before truncating an infinite series.

    - ``category`` -- the category of the parent can be specified
      in order to broaden the base structure. It has to be a
      subcategory of ``Category of rings``. This is also the default
      category if ``None`` is specified.

    - ``term_monoid_factory`` -- a :class:`~sage.rings.asymptotic.term_monoid.TermMonoidFactory`.
      If ``None``, then :class:`~sage.rings.asymptotic.term_monoid.DefaultTermMonoidFactory`
      is used.

    - ``locals`` -- a dictionary which may contain the following keys and values:

      - ``'log'`` -- value: a function. If not given, then the usual
        :class:`log <sage.functions.log.Function_log>` is taken.
        (See also :meth:`AsymptoticExpansion.log`.)

    EXAMPLES:

    We begin with the construction of an asymptotic ring in various
    ways. First, we simply pass a string specifying the underlying
    growth group::

        sage: R1_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R1_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: x
        x

    This is equivalent to the following code, which explicitly
    specifies the underlying growth group::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G_QQ = GrowthGroup('x^QQ')
        sage: R2_x.<x> = AsymptoticRing(growth_group=G_QQ, coefficient_ring=QQ); R2_x
        Asymptotic Ring <x^QQ> over Rational Field

    Of course, the coefficient ring of the asymptotic ring and the
    base ring of the underlying growth group do not need to
    coincide::

        sage: R_ZZ_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=ZZ); R_ZZ_x
        Asymptotic Ring <x^QQ> over Integer Ring

    Note, we can also create and use logarithmic growth groups::

        sage: R_log = AsymptoticRing(growth_group='log(x)^ZZ', coefficient_ring=QQ); R_log
        Asymptotic Ring <log(x)^ZZ> over Rational Field

    Other growth groups are available. See :doc:`asymptotic_ring` for
    more examples.

    Below there are some technical details.

    According to the conventions for parents, uniqueness is ensured::

        sage: R1_x is R2_x
        True

    Furthermore, the coercion framework is also involved. Coercion
    between two asymptotic rings is possible (given that the
    underlying growth groups and coefficient rings are chosen
    appropriately)::

        sage: R1_x.has_coerce_map_from(R_ZZ_x)
        True

    Additionally, for the sake of convenience, the coefficient ring
    also coerces into the asymptotic ring (representing constant
    quantities)::

        sage: R1_x.has_coerce_map_from(QQ)
        True

    It is possible to customize the terms in an asymptotic expansion::

        sage: from sage.rings.asymptotic.term_monoid import ExactTermMonoid, OTermMonoid
        sage: from sage.rings.asymptotic.term_monoid import TermMonoidFactory
        sage: class MyExactTermMonoid(ExactTermMonoid):
        ....:     pass
        sage: class MyOTermMonoid(OTermMonoid):
        ....:     pass
        sage: MyTermMonoid = TermMonoidFactory('MyTermMonoid',
        ....:                                  exact_term_monoid_class=MyExactTermMonoid,
        ....:                                  O_term_monoid_class=MyOTermMonoid)
        sage: G = GrowthGroup('x^ZZ')
        sage: A.<n> = AsymptoticRing(growth_group=G, coefficient_ring=QQ, term_monoid_factory=MyTermMonoid)
        sage: a = A.an_element(); a
        1/8*x^3 + O(x)
        sage: for t in a.summands.elements_topological(reverse=True):
        ....:     print(t, type(t))
        1/8*x^3 <class '__main__.MyExactTermMonoid_with_category.element_class'>
        O(x) <class '__main__.MyOTermMonoid_with_category.element_class'>

    TESTS::

        sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRing as AR_class
        sage: class AR(AR_class):
        ....:     class Element(AR_class.Element):
        ....:         __eq__ = AR_class.Element.has_same_summands
        sage: A = AR(growth_group='z^QQ', coefficient_ring=QQ)
        sage: from itertools import islice
        sage: TestSuite(A).run(  # not tested  # long time  # see #19424
        ....:     verbose=True,
        ....:     elements=tuple(islice(A.some_elements(), int(10))),
        ....:     skip=('_test_some_elements',  # to many elements
        ....:           '_test_distributivity'))  # due to cancellations: O(z) != O(z^2)
    """

    # enable the category framework for elements
    Element = AsymptoticExpansion


    @staticmethod
    def __classcall__(cls, growth_group=None, coefficient_ring=None,
                      names=None, category=None, default_prec=None,
                      term_monoid_factory=None,
                      locals=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation of the parent.

        For more information see :class:`AsymptoticRing`.

        EXAMPLES:

        ``__classcall__`` unifies the input to the constructor of
        :class:`AsymptoticRing` such that the instances generated
        are unique. Also, this enables the use of the generation
        framework::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: MG = GrowthGroup('x^ZZ')
            sage: AR1 = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR2.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR1 is AR2
            True

        The bracket notation can only be used if the growth group
        has a generator::

            sage: AR.<lx> = AsymptoticRing(growth_group='log(x)^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError:  Growth Group log(x)^ZZ does not provide any
            generators but name 'lx' given.

        The names of the generators have to agree with the names used in
        the growth group except for univariate rings::

            sage: A.<icecream> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ); A
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: icecream
            x
            sage: A.<x, y> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ); A
            Asymptotic Ring <x^ZZ * y^ZZ> over Integer Ring
            sage: A.<y, x> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'y', 'x' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<a, b> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'a', 'b' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<x, b> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'x', 'b' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<x> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Name 'x' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<x, y, z> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'x', 'y', 'z' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.

        TESTS::

            sage: AsymptoticRing(growth_group=None, coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Growth group not specified. Cannot continue.
            sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring=None)
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring not specified. Cannot continue.
            sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring='icecream')
            Traceback (most recent call last):
            ...
            ValueError: icecream is not a ring. Cannot continue.

        ::

            sage: A.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=ZZ)
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory
            sage: A.term_monoid_factory is DefaultTermMonoidFactory
            True
        """
        from sage.categories.sets_cat import Sets
        from sage.categories.rings import Rings

        Sets_parent_class = Sets().parent_class
        while issubclass(cls, Sets_parent_class):
            cls = cls.__base__

        if isinstance(growth_group, str):
            from .growth_group import GrowthGroup
            growth_group = GrowthGroup(growth_group)

        if growth_group is None:
            raise ValueError('Growth group not specified. Cannot continue.')

        if coefficient_ring is None:
            raise ValueError('Coefficient ring not specified. Cannot continue.')
        if coefficient_ring not in Rings():
            raise ValueError('%s is not a ring. Cannot continue.' % (coefficient_ring,))

        strgens = tuple(str(g) for g in growth_group.gens_monomial())
        def format_names(N):
            return ('s ' if len(N) != 1 else ' ') + ', '.join("'%s'" % n for n in N)
        if names and not strgens:
            raise ValueError('%s does not provide any generators but name%s given.' %
                             (growth_group, format_names(names)))
        elif names is not None and len(names) == 1 and len(strgens) == 1:
            pass
        elif names is not None and names != strgens:
            raise ValueError('Name%s do not coincide with generator%s of %s.' %
                             (format_names(names), format_names(strgens), growth_group))

        if category is None:
            from sage.categories.commutative_algebras import CommutativeAlgebras
            from sage.categories.rings import Rings
            category = CommutativeAlgebras(Rings())

        if default_prec is None:
            default_prec = series_precision()

        if term_monoid_factory is None:
            from .term_monoid import DefaultTermMonoidFactory
            term_monoid_factory = DefaultTermMonoidFactory

        if locals is not None:
            locals = cls._convert_locals_(locals)

        return super(AsymptoticRing,
                     cls).__classcall__(cls, growth_group, coefficient_ring,
                                        category=category,
                                        default_prec=default_prec,
                                        term_monoid_factory=term_monoid_factory,
                                        locals=locals)

    def __init__(self, growth_group, coefficient_ring,
                 category, default_prec,
                 term_monoid_factory, locals):
        r"""
        See :class:`AsymptoticRing` for more information.

        TESTS::

            sage: R1 = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ); R1
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: R2.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R2
            Asymptotic Ring <x^QQ> over Rational Field
            sage: R1 is R2
            False

        ::

            sage: R3 = AsymptoticRing('x^ZZ')
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring not specified. Cannot continue.
        """
        self._coefficient_ring_ = coefficient_ring
        self._growth_group_ = growth_group
        self._default_prec_ = default_prec
        self._term_monoid_factory_ = term_monoid_factory
        self._locals_ = locals
        super(AsymptoticRing, self).__init__(base_ring=coefficient_ring,
                                             category=category)


    @property
    def growth_group(self):
        r"""
        The growth group of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.growth_group
            Growth Group x^ZZ

        .. SEEALSO::

            :doc:`growth_group`
        """
        return self._growth_group_


    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.coefficient_ring
            Integer Ring
        """
        return self._coefficient_ring_


    @property
    def default_prec(self):
        r"""
        The default precision of this asymptotic ring.

        This is the parameter used to determine how many summands
        are kept before truncating an infinite series (which occur
        when inverting asymptotic expansions).

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.default_prec
            20
            sage: AR = AsymptoticRing('x^ZZ', ZZ, default_prec=123)
            sage: AR.default_prec
            123
        """
        return self._default_prec_


    @property
    def term_monoid_factory(self):
        r"""
        The term monoid factory of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.term_monoid_factory
            Term Monoid Factory 'sage.rings.asymptotic.term_monoid.DefaultTermMonoidFactory'

        .. SEEALSO::

            :doc:`term_monoid`
        """
        return self._term_monoid_factory_


    def term_monoid(self, type):
        r"""
        Return the term monoid of this asymptotic ring of specified ``type``.

        INPUT:

        - ``type`` -- 'O' or 'exact', or an instance of an existing
          term monoid.
          See :class:`~sage.rings.asymptotic.term_monoid.TermMonoidFactory`
          for more details.

        OUTPUT:

        A term monoid object derived from
        :class:`~sage.rings.asymptotic.term_monoid.GenericTermMonoid`.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.term_monoid('exact')
            Exact Term Monoid x^ZZ with coefficients in Integer Ring
            sage: AR.term_monoid('O')
            O-Term Monoid x^ZZ with implicit coefficients in Integer Ring
            sage: AR.term_monoid(AR.term_monoid('exact'))
            Exact Term Monoid x^ZZ with coefficients in Integer Ring
        """
        TermMonoid = self.term_monoid_factory
        return TermMonoid(type, asymptotic_ring=self)


    def change_parameter(self, **kwds):
        r"""
        Return an asymptotic ring with a change in one or more of the given parameters.

        INPUT:

        - ``growth_group`` -- (default: ``None``) the new growth group.

        - ``coefficient_ring`` -- (default: ``None``) the new coefficient ring.

        - ``category`` -- (default: ``None``) the new category.

        - ``default_prec`` -- (default: ``None``) the new default precision.

        OUTPUT:

        An asymptotic ring.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: A.change_parameter(coefficient_ring=QQ)
            Asymptotic Ring <x^ZZ> over Rational Field

        TESTS::

            sage: A.change_parameter(coefficient_ring=ZZ) is A
            True
            sage: A.change_parameter(coefficient_ring=None) is A
            True
        """
        parameters = ('growth_group', 'coefficient_ring', 'default_prec',
                      'term_monoid_factory')
        values = dict()
        category = self.category()
        values['category'] = category
        locals = self._locals_
        values['locals'] = locals
        for parameter in parameters:
            default = getattr(self, parameter)
            values[parameter] = kwds.get(parameter, default)
            if values[parameter] is None:
                values[parameter] = default
        if isinstance(values['growth_group'], str):
            from .growth_group import GrowthGroup
            values['growth_group'] = GrowthGroup(values['growth_group'])
        if all(values[parameter] is getattr(self, parameter)
               for parameter in parameters) and values['category'] is category and values['locals'] is locals:
            return self
        return self._underlying_class()(**values)

    @staticmethod
    def _create_empty_summands_():
        r"""
        Create an empty data structure suitable for storing and working
        with summands.

        INPUT:

        Nothing.

        OUTPUT:

        A :class:`~sage.data_structures.mutable_poset.MutablePoset`.

        TESTS::

            sage: AsymptoticRing._create_empty_summands_()
            poset()
        """
        from sage.data_structures.mutable_poset import MutablePoset
        from .term_monoid import can_absorb, absorption
        return MutablePoset(key=lambda element: element.growth,
                            can_merge=can_absorb,
                            merge=absorption)


    def _create_element_in_extension_(self, term, old_term_parent=None):
        r"""
        Create an element in an extension of this asymptotic ring which
        is chosen according to the input.

        INPUT:

        - ``term`` -- the element data.

        - ``old_term_parent`` -- the parent of ``term`` is compared to this
          parent. If both are the same or ``old_parent`` is ``None``,
          then the result is an expansion in this (``self``) asymptotic ring.

        OUTPUT:

        An element.

        EXAMPLES::

            sage: A = AsymptoticRing('z^ZZ', ZZ)
            sage: a = next(A.an_element().summands.elements_topological())
            sage: B = AsymptoticRing('z^QQ', QQ)
            sage: b = next(B.an_element().summands.elements_topological())
            sage: c = A._create_element_in_extension_(a, a.parent())
            sage: next(c.summands.elements_topological()).parent()
            O-Term Monoid z^ZZ with implicit coefficients in Integer Ring
            sage: c = A._create_element_in_extension_(b, a.parent())
            sage: next(c.summands.elements_topological()).parent()
            O-Term Monoid z^QQ with implicit coefficients in Rational Field

        TESTS::

            sage: c = A._create_element_in_extension_(b, None)
            sage: next(c.summands.elements_topological()).parent()
            O-Term Monoid z^QQ with implicit coefficients in Rational Field
        """
        if old_term_parent is None or term.parent() is old_term_parent:
            parent = self
        else:
            # Insert an 'if' here once terms can have different
            # coefficient rings, as this will be for L-terms.
            parent = self.change_parameter(
                growth_group=term.parent().growth_group,
                coefficient_ring=term.parent().coefficient_ring)
        return parent(term, simplify=False, convert=False)


    def _element_constructor_(self, data, simplify=True, convert=True):
        r"""
        Convert a given object to this asymptotic ring.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``simplify`` -- (default: ``True``) if set, then the constructed
          element is simplified (terms are absorbed) automatically.

        - ``convert`` -- (default: ``True``) passed on to the element
          constructor.  If set, then the ``summands`` are converted to
          the asymptotic ring (the parent of this expansion). If not,
          then the summands are taken as they are. In that case, the
          caller must ensure that the parent of the terms is set
          correctly.

        OUTPUT:

        An element of this asymptotic ring.

        TESTS::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR(5)  # indirect doctest
            5
            sage: AR(3*x^2)  # indirect doctest
            3*x^2
            sage: x = ZZ['x'].gen(); x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: AR(x)
            x
            sage: y = ZZ['y'].gen(); AR(y)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Polynomial y is not in
            Asymptotic Ring <x^ZZ> over Integer Ring
            > *previous* ValueError: Growth y is not valid in
            Exact Term Monoid x^ZZ with coefficients in Integer Ring.
            >> *previous* ValueError: y is not in Growth Group x^ZZ.

        ::

            sage: A = AsymptoticRing(growth_group='p^ZZ', coefficient_ring=QQ)
            sage: P.<p> = QQ[]
            sage: A(p)  # indirect doctest
            p
            sage: A(p^11)  # indirect doctest
            p^11
            sage: A(2*p^11)  # indirect doctest
            2*p^11
            sage: A(3*p^4 + 7/3*p - 8)  # indirect doctest
            3*p^4 + 7/3*p - 8

        ::

            sage: S = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=QQ)
            sage: var('x, y')
            (x, y)
            sage: S(x + y)  # indirect doctest
            x + y
            sage: S(2*x - 4*x*y^6)  # indirect doctest
            -4*x*y^6 + 2*x

        ::

            sage: A.<a,b> = AsymptoticRing('a^ZZ * b^ZZ', QQ)
            sage: 1/a
            a^(-1)

        ::

            sage: P.<a, b, c> = ZZ[]
            sage: A(a + b)
            a + b
            sage: A(a + c)
            Traceback (most recent call last):
            ...
            ValueError: Polynomial a + c is not in
            Asymptotic Ring <a^ZZ * b^ZZ> over Rational Field
            > *previous* ValueError: Growth c is not valid in
            Exact Term Monoid a^ZZ * b^ZZ with coefficients in Rational Field.
            >> *previous* ValueError: c is not in Growth Group a^ZZ * b^ZZ.
            >...> *previous* ValueError: c is not in any of the factors
            of Growth Group a^ZZ * b^ZZ
            >...> *previous* ValueError: c is not in Growth Group a^ZZ.
            >...> *and* ValueError: c is not in Growth Group b^ZZ.

        ::

            sage: M = AsymptoticRing('m^ZZ', ZZ)
            sage: N = AsymptoticRing('n^ZZ', QQ)
            sage: N(M.an_element())  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot include m^3 with parent
            Exact Term Monoid m^ZZ with coefficients in Integer Ring
            in Asymptotic Ring <n^ZZ> over Rational Field
            > *previous* ValueError: Growth m^3 is not valid in
              Exact Term Monoid n^ZZ with coefficients in Rational Field.
            >> *previous* ValueError: m^3 is not in Growth Group n^ZZ.

        ::

            sage: M([1])  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: Not all list entries of [1] are asymptotic terms,
            so cannot create an asymptotic expansion in
            Asymptotic Ring <m^ZZ> over Integer Ring.
            sage: M(SR.var('a') + 1)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Symbolic expression a + 1 is not in
            Asymptotic Ring <m^ZZ> over Integer Ring.
            > *previous* ValueError: a is not in
            Exact Term Monoid m^ZZ with coefficients in Integer Ring.
            >> *previous* ValueError: Factor a of a is neither a coefficient
            (in Integer Ring) nor growth (in Growth Group m^ZZ).
        """
        from sage.data_structures.mutable_poset import MutablePoset
        if isinstance(data, MutablePoset):
            return self.element_class(self, data, simplify=simplify, convert=convert)

        if type(data) == self.element_class and data.parent() == self:
            return data

        if isinstance(data, AsymptoticExpansion):
            return self.element_class(self, data.summands,
                                      simplify=simplify, convert=convert)

        from .term_monoid import GenericTerm
        if isinstance(data, GenericTerm):
            data = (data,)

        if isinstance(data, (list, tuple)):
            if not all(isinstance(elem, GenericTerm) for elem in data):
                raise TypeError('Not all list entries of %s '
                                'are asymptotic terms, so cannot create an '
                                'asymptotic expansion in %s.' % (data, self))
            summands = AsymptoticRing._create_empty_summands_()
            summands.union_update(data)
            return self.element_class(self, summands,
                                      simplify=simplify, convert=convert)

        if not data:
            summands = AsymptoticRing._create_empty_summands_()
            return self.element_class(self, summands,
                                      simplify=simplify, convert=False)

        try:
            P = data.parent()
        except AttributeError:
            return self.create_summand('exact', data)

        from .misc import combine_exceptions
        from sage.symbolic.ring import SymbolicRing
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
        from sage.rings.power_series_ring import is_PowerSeriesRing

        if isinstance(P, SymbolicRing):
            from sage.symbolic.operators import add_vararg
            if data.operator() == add_vararg:
                summands = []
                for summand in data.operands():
                    # TODO: check if summand is an O-Term here
                    # (see #19425, #19426)
                    try:
                        summands.append(self.create_summand('exact', summand))
                    except ValueError as e:
                        raise combine_exceptions(
                            ValueError('Symbolic expression %s is not in %s.' %
                                       (data, self)), e)
                return sum(summands, self.zero())

        elif is_PolynomialRing(P):
            p = P.gen()
            try:
                return sum(iter(self.create_summand('exact', growth=p**i,
                                                    coefficient=c)
                                for i, c in enumerate(data)),
                           self.zero())
            except ValueError as e:
                raise combine_exceptions(
                    ValueError('Polynomial %s is not in %s' % (data, self)), e)

        elif is_MPolynomialRing(P):
            try:
                return sum(iter(self.create_summand('exact', growth=g, coefficient=c)
                                for c, g in iter(data)),
                           self.zero())
            except ValueError as e:
                raise combine_exceptions(
                    ValueError('Polynomial %s is not in %s' % (data, self)), e)

        elif is_PowerSeriesRing(P):
            raise NotImplementedError(
                'Cannot convert %s from the %s to an asymptotic expansion '
                'in %s, since growths at other points than +oo are not yet '
                'supported.' % (data, P, self))
            # Delete lines above as soon as we can deal with growths
            # other than the that at going to +oo.
            from sage.rings.infinity import PlusInfinity
            p = P.gen()
            try:
                result = self(data.polynomial())
            except ValueError as e:
                raise combine_exceptions(
                    ValueError('Powerseries %s is not in %s' % (data, self)), e)
            prec = data.precision_absolute()
            if prec < PlusInfinity():
                try:
                    result += self.create_summand('O', growth=p**prec)
                except ValueError as e:
                    raise combine_exceptions(
                        ValueError('Powerseries %s is not in %s' %
                                   (data, self)), e)
            return result

        return self.create_summand('exact', data)


    def _coerce_map_from_(self, R):
        r"""
        Return whether ``R`` coerces into this asymptotic ring.

        INPUT:

        - ``R`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            There are two possible cases: either ``R`` coerces in the
            :meth:`coefficient_ring` of this asymptotic ring, or ``R``
            itself is an asymptotic ring, where both the
            :meth:`growth_group` and the :meth:`coefficient_ring` coerce into
            the :meth:`growth_group` and the :meth:`coefficient_ring` of this
            asymptotic ring, respectively.

        TESTS::

            sage: AR_ZZ = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ); AR_ZZ
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: x_ZZ = AR_ZZ.gen()
            sage: AR_QQ = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); AR_QQ
            Asymptotic Ring <x^QQ> over Rational Field
            sage: x_QQ = AR_QQ.gen()
            sage: AR_QQ.has_coerce_map_from(AR_ZZ)  # indirect doctest
            True
            sage: x_ZZ * x_QQ
            x^2

        ::

            sage: AR_QQ.has_coerce_map_from(QQ)
            True
            sage: AR_QQ.has_coerce_map_from(ZZ)
            True
        """
        from sage.data_structures.mutable_poset import MutablePoset
        if R == MutablePoset:
            return
        if self.coefficient_ring.has_coerce_map_from(R):
            return True
        if self.growth_group.has_coerce_map_from(R):
            return True
        elif isinstance(R, AsymptoticRing):
            if self.growth_group.has_coerce_map_from(R.growth_group) and \
                    self.coefficient_ring.has_coerce_map_from(R.coefficient_ring):
                return True


    def _repr_(self):
        r"""
        A representation string of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ',
            ....:                     coefficient_ring=ZZ)
            sage: repr(AR)  # indirect doctest
            'Asymptotic Ring <x^ZZ> over Integer Ring'
        """
        try:
            G = '<' + self.growth_group._repr_(condense=True) + '>'
        except TypeError:
            G = repr(self.growth_group)
        return 'Asymptotic Ring %s over %s' % (G, self.coefficient_ring)


    def _an_element_(self):
        r"""
        Return an element of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        An :class:`AsymptoticExpansion`.

        EXAMPLES::

            sage: AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ).an_element()
            z^(3/2) + O(z^(1/2))
            sage: AsymptoticRing(growth_group='z^ZZ', coefficient_ring=QQ).an_element()
            1/8*z^3 + O(z)
            sage: AsymptoticRing(growth_group='z^QQ', coefficient_ring=QQ).an_element()
            1/8*z^(3/2) + O(z^(1/2))
        """
        E = self.term_monoid('exact')
        O = self.term_monoid('O')
        return self(E.an_element(), simplify=False, convert=False)**3 + \
            self(O.an_element(), simplify=False, convert=False)


    def some_elements(self):
        r"""
        Return some elements of this term monoid.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from itertools import islice
            sage: A = AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ)
            sage: tuple(islice(A.some_elements(), int(10)))
            (z^(3/2) + O(z^(1/2)),
             O(z^(1/2)),
             z^(3/2) + O(z^(-1/2)),
             -z^(3/2) + O(z^(1/2)),
             O(z^(-1/2)),
             O(z^2),
             z^6 + O(z^(1/2)),
             -z^(3/2) + O(z^(-1/2)),
             O(z^2),
             z^(3/2) + O(z^(-2)))
        """
        from sage.misc.mrange import cantor_product
        E = self.term_monoid('exact')
        O = self.term_monoid('O')
        return iter(self(e, simplify=False, convert=False)**3 +
                    self(o, simplify=False, convert=False)
                    for e, o in cantor_product(
                            E.some_elements(), O.some_elements()))


    def gens(self):
        r"""
        Return a tuple with generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple of asymptotic expansions.

        .. NOTE::

            Generators do not necessarily exist. This depends on the
            underlying growth group. For example,
            :class:`monomial growth groups <sage.rings.asymptotic.growth_group.MonomialGrowthGroup>`
            have a generator, and exponential growth groups
            do not.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.gens()
            (x,)
            sage: B.<y,z> = AsymptoticRing(growth_group='y^ZZ * z^ZZ', coefficient_ring=QQ)
            sage: B.gens()
            (y, z)
        """
        return tuple(self.create_summand('exact',
                                         growth=g,
                                         coefficient=self.coefficient_ring(1))
                     for g in self.growth_group.gens_monomial())


    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of this asymptotic ring.

        INPUT:

        - ``n`` -- (default: `0`) a non-negative integer.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: R.gen()
            x
        """
        return self.gens()[n]


    def ngens(self):
        r"""
        Return the number of generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.ngens()
            1
        """
        return len(self.growth_group.gens_monomial())


    def coefficients_of_generating_function(self, function, singularities, precision=None,
                                            return_singular_expansions=False,
                                            error_term=None):
        r"""
        Return the asymptotic growth of the coefficients of some
        generating function by means of Singularity Analysis.

        INPUT:

        - ``function`` -- a callable function in one variable.

        - ``singularities`` -- list of dominant singularities of the function.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``return_singular_expansions`` -- (default: ``False``) a boolean.
          If set, the singular expansions are also returned.

        - ``error_term`` -- (default: ``None``) an asymptotic expansion.
          If ``None``, then this is interpreted as zero.
          The contributions of the coefficients are added to ``error_term``
          during Singularity Analysis.

        OUTPUT:

        - If ``return_singular_expansions=False``: An asymptotic expansion from
          this ring.

        - If ``return_singular_expansions=True``: A named tuple with
          components ``asymptotic_expansion`` and
          ``singular_expansions``. The former contains an asymptotic
          expansion from this ring, the latter is a dictionary which
          contains the singular expansions around the singularities.

        .. TODO::

            Make this method more usable by implementing the
            processing of symbolic expressions.

        EXAMPLES:

        Catalan numbers::

            sage: def catalan(z):
            ....:     return (1-(1-4*z)^(1/2))/(2*z)
            sage: B.<n> = AsymptoticRing('QQ^n * n^QQ', QQ)
            sage: B.coefficients_of_generating_function(catalan, (1/4,), precision=3)
            1/sqrt(pi)*4^n*n^(-3/2) - 9/8/sqrt(pi)*4^n*n^(-5/2)
            + 145/128/sqrt(pi)*4^n*n^(-7/2) + O(4^n*n^(-4))
            sage: B.coefficients_of_generating_function(catalan, (1/4,), precision=2,
            ....:                        return_singular_expansions=True)
            SingularityAnalysisResult(asymptotic_expansion=1/sqrt(pi)*4^n*n^(-3/2)
            - 9/8/sqrt(pi)*4^n*n^(-5/2) + O(4^n*n^(-3)),
            singular_expansions={1/4: 2 - 2*T^(-1/2)
            + 2*T^(-1) - 2*T^(-3/2) + O(T^(-2))})

        Unit fractions::

            sage: def logarithmic(z):
            ....:     return -log(1-z)
            sage: B.coefficients_of_generating_function(logarithmic, (1,), precision=5)
            n^(-1) + O(n^(-3))

        Harmonic numbers::

            sage: def harmonic(z):
            ....:     return -log(1-z)/(1-z)
            sage: B.<n> = AsymptoticRing('QQ^n * n^QQ * log(n)^QQ', QQ)
            sage: ex = B.coefficients_of_generating_function(harmonic, (1,), precision=13); ex
            log(n) + euler_gamma + 1/2*n^(-1) - 1/12*n^(-2) + 1/120*n^(-4)
            + O(n^(-6))
            sage: ex.has_same_summands(asymptotic_expansions.HarmonicNumber(
            ....:    'n', precision=5))
            True

        .. WARNING::

            Once singular expansions around points other than infinity
            are implemented (:trac:`20050`), the output in the case
            ``return_singular_expansions`` will change to return singular
            expansions around the singularities.

        In the following example, the result is an exact asymptotic expression
        for sufficiently large `n` (i.e., there might be finitely many
        exceptional values). This is encoded by an `O(0)` error term::

            sage: def f(z):
            ....:     return z/(1-z)
            sage: B.coefficients_of_generating_function(f, (1,), precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got 1 + O(0)
            The error term O(0) means 0 for sufficiently large n.

        In this case, we can manually intervene by adding an error term
        that suits us::

            sage: B.coefficients_of_generating_function(f, (1,), precision=3,
            ....:                                       error_term=O(n^-100))
            1 + O(n^(-100))
        """
        from sage.symbolic.ring import SR
        from .misc import NotImplementedOZero

        singular_expansions = {}

        OZeroEncountered = False

        A = AsymptoticRing('T^QQ * log(T)^QQ', coefficient_ring=SR,
                           default_prec=precision)
        T = A.gen()

        if error_term is None:
            result = self.zero()
        else:
            result = error_term
        for singularity in singularities:
            singular_expansion = A(function((1-1/T)*singularity))
            singular_expansions[singularity] = singular_expansion

            try:
                contribution = singular_expansion._singularity_analysis_(
                        var='Z', zeta=singularity,
                        precision=precision).subs(Z=self.gen())
            except NotImplementedOZero as ozero:
                OZeroEncountered = True
                result += ozero.exact_part.subs(Z=self.gen())
            else:
                result += contribution

        if OZeroEncountered and result.is_exact():
            raise NotImplementedOZero(self, exact_part=result)

        if return_singular_expansions:
            from collections import namedtuple
            SingularityAnalysisResult = namedtuple(
                'SingularityAnalysisResult',
                ['asymptotic_expansion', 'singular_expansions'])
            return SingularityAnalysisResult(
                asymptotic_expansion=result,
                singular_expansions=singular_expansions)
        else:
            return result


    def create_summand(self, type, data=None, **kwds):
        r"""
        Create a simple asymptotic expansion consisting of a single
        summand.

        INPUT:

        - ``type`` -- 'O' or 'exact'.

        - ``data`` -- the element out of which a summand has to be created.

        - ``growth`` -- an element of the :meth:`growth_group`.

        - ``coefficient`` -- an element of the :meth:`coefficient_ring`.

        .. NOTE::

            Either ``growth`` and ``coefficient`` or ``data`` have to
            be specified.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            This method calls the factory :class:`TermMonoid
            <sage.rings.asymptotic.term_monoid.TermMonoidFactory>`
            with the appropriate arguments.

        EXAMPLES::

            sage: R = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: R.create_summand('O', x^2)
            O(x^2)
            sage: R.create_summand('exact', growth=x^456, coefficient=123)
            123*x^456
            sage: R.create_summand('exact', data=12*x^13)
            12*x^13

        TESTS::

            sage: R.create_summand('exact', data='12*x^13')
            12*x^13
            sage: R.create_summand('exact', data='x^13 * 12')
            12*x^13
            sage: R.create_summand('exact', data='x^13')
            x^13
            sage: R.create_summand('exact', data='12')
            12
            sage: R.create_summand('exact', data=12)
            12

        ::

            sage: Z = R.change_parameter(coefficient_ring=Zmod(3))
            sage: Z.create_summand('exact', data=42)
            0

        ::

            sage: R.create_summand('O', growth=42*x^2, coefficient=1)
            Traceback (most recent call last):
            ...
            ValueError: Growth 42*x^2 is not valid in O-Term Monoid x^ZZ with implicit coefficients in Integer Ring.
            > *previous* ValueError: 42*x^2 is not in Growth Group x^ZZ.

        ::

            sage: AR.<z> = AsymptoticRing('z^QQ', QQ)
            sage: AR.create_summand('exact', growth='z^2')
            Traceback (most recent call last):
            ...
            TypeError: Cannot create exact term: only 'growth' but
            no 'coefficient' specified.
        """
        from .term_monoid import ZeroCoefficientError
        TM = self.term_monoid(type)

        if data is None:
            try:
                data = kwds.pop('growth')
            except KeyError:
                raise TypeError("Neither 'data' nor 'growth' are specified.")
            if type == 'exact' and kwds.get('coefficient') is None:
                raise TypeError("Cannot create exact term: only 'growth' "
                                "but no 'coefficient' specified.")

        try:
            return self(TM(data, **kwds), simplify=False, convert=False)
        except ZeroCoefficientError:
            return self.zero()


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ * QQ^y', coefficient_ring=QQ)
            sage: A.variable_names()
            ('x', 'y')
        """
        return self.growth_group.variable_names()


    def construction(self):
        r"""
        Return the construction of this asymptotic ring.

        OUTPUT:

        A pair whose first entry is an
        :class:`asymptotic ring construction functor <AsymptoticRingFunctor>`
        and its second entry the coefficient ring.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ * QQ^y', coefficient_ring=QQ)
            sage: A.construction()
            (AsymptoticRing<x^ZZ * QQ^y * Signs^y>, Rational Field)

        .. SEEALSO::

            :doc:`asymptotic_ring`,
            :class:`AsymptoticRing`,
            :class:`AsymptoticRingFunctor`.

        TESTS:

        :trac:`22392`::

            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRing
            sage: class MyAsymptoticRing(AsymptoticRing):
            ....:     pass
            sage: A = MyAsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: A.construction()[0].cls
            <class '__main__.MyAsymptoticRing'>
        """
        return (AsymptoticRingFunctor(self.growth_group,
                                      default_prec=self.default_prec,
                                      category=self.category(),
                                      cls=self._underlying_class()),
                self.coefficient_ring)

    @staticmethod
    def B(self, valid_from=0):
        r""""
        Create a B-term.

        INPUT:

        - ``valid_from`` -- dictionary mapping variable names to lower bounds
          for the corresponding variable. The bound implied by this term is valid when
          all variables are at least their corresponding lower bound. If a number
          is passed to ``valid_from``, then the lower bounds for all variables of
          the asymptotic expansion are set to this number

        OUTPUT:

        A B-term

        EXAMPLES::

            sage: A.<x> = AsymptoticRing(growth_group='x^ZZ * QQ^y', coefficient_ring=QQ)
            sage: A.B(2*x^3, {x: 5})
            B(2*x^3, x >= 5)
        """
        return self.B(valid_from)


from sage.categories.pushout import ConstructionFunctor


class AsymptoticRingFunctor(ConstructionFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`asymptotic rings <AsymptoticRing>`.

    INPUT:

    - ``growth_group`` -- a partially ordered group (see
      :class:`AsymptoticRing` or
      :doc:`growth_group` for details).

    - ``default_prec`` -- ``None`` (default) or an integer.

    - ``category`` -- ``None`` (default) or a category.

    - ``cls`` -- :class:`AsymptoticRing` (default) or a derived class.

    EXAMPLES::

        sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ).construction()  # indirect doctest
        (AsymptoticRing<x^ZZ>, Rational Field)

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`AsymptoticRing`,
        :class:`sage.rings.asymptotic.growth_group.AbstractGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.growth_group.ExponentialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.growth_group.MonomialGrowthGroupFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.

    TESTS::

        sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
        sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
        sage: cm = sage.structure.element.get_coercion_model()
        sage: cm.record_exceptions()
        sage: cm.common_parent(X, Y)
        Asymptotic Ring <x^ZZ * y^ZZ> over Rational Field
        sage: sage.structure.element.coercion_traceback()  # not tested

    ::

        sage: from sage.categories.pushout import pushout
        sage: pushout(AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ), QQ)
        Asymptotic Ring <x^ZZ> over Rational Field
    """

    rank = 13


    def __init__(self, growth_group,
                 default_prec=None, category=None,
                 term_monoid_factory=None, locals=None,
                 cls=None):
        r"""
        See :class:`AsymptoticRingFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRingFunctor
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: AsymptoticRingFunctor(GrowthGroup('x^ZZ'))
            AsymptoticRing<x^ZZ>
        """
        self.growth_group = growth_group
        if cls is None:
            self.cls = AsymptoticRing
        else:
            self.cls = cls
        self._default_prec_ = default_prec
        self._category_ = category
        self._term_monoid_factory_ = term_monoid_factory
        self._locals_ = locals

        from sage.categories.rings import Rings
        super(ConstructionFunctor, self).__init__(
            Rings(), Rings())


    def _repr_(self):
        r"""
        Return a representation string of this functor.

        OUTPUT:

        A string.

        TESTS::

            sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ).construction()[0]  # indirect doctest
            AsymptoticRing<x^ZZ>

        :trac:`22392`::

            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRing
            sage: class MyAsymptoticRing(AsymptoticRing):
            ....:     pass
            sage: A = MyAsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: A.construction()
            (MyAsymptoticRing<x^ZZ>, Rational Field)
        """
        return '{}<{}>'.format(self.cls.__name__,
                               self.growth_group._repr_(condense=True))


    def _apply_functor(self, coefficient_ring):
        r"""
        Apply this functor to the given ``coefficient_ring``.

        INPUT:

        - ``base`` - anything :class:`~sage.rings.asymptotic.growth_group.MonomialGrowthGroup` accepts.

        OUTPUT:

        An :class:`AsymptoticRing`.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: F, C = A.construction()
            sage: F(C)  # indirect doctest
            Asymptotic Ring <x^ZZ> over Rational Field

        TESTS:

        :trac:`22392`::

            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRing
            sage: class MyAsymptoticRing(AsymptoticRing):
            ....:     pass
            sage: A = MyAsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: type(A.construction()[0](ZZ))
            <class '__main__.MyAsymptoticRing_with_category'>


            sage: C = CyclotomicField(3)
            sage: P = C['z']
            sage: type(P(2) * A.gen())
            <class '...MyAsymptoticRing_with_category.element_class'>

        :trac:`22396`::

            sage: A.<n> = AsymptoticRing('ZZ^n * n^ZZ', ZZ, default_prec=3)
            sage: 1/(QQ(1)+n)
            n^(-1) - n^(-2) + n^(-3) + O(n^(-4))
        """
        kwds = {'growth_group': self.growth_group,
                'coefficient_ring': coefficient_ring}
        parameters = ('category', 'default_prec',
                      'term_monoid_factory', 'locals')
        for parameter in parameters:
            value = getattr(self, '_{}_'.format(parameter))
            if value is not None:
                kwds[parameter] = value
        return self.cls(**kwds)


    def merge(self, other):
        r"""
        Merge this functor with ``other`` if possible.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A functor or ``None``.

        EXAMPLES::

            sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_X = X.construction()[0]
            sage: F_Y = Y.construction()[0]
            sage: F_X.merge(F_X)
            AsymptoticRing<x^ZZ>
            sage: F_X.merge(F_Y)
            AsymptoticRing<x^ZZ * y^ZZ>

        TESTS:

        :trac:`22396`::

            sage: AN = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_AN = AN.construction()[0]; F_AN._default_prec_ = None
            sage: A3 = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ, default_prec=3)
            sage: F_A3 = A3.construction()[0]
            sage: A5 = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ, default_prec=5)
            sage: F_A5 = A5.construction()[0]

            sage: F_AN.merge(F_AN)(ZZ).default_prec
            20
            sage: F_AN.merge(F_A3)(ZZ).default_prec
            3
            sage: F_AN.merge(F_A5)(ZZ).default_prec
            5
            sage: F_A3.merge(F_AN)(ZZ).default_prec
            3
            sage: F_A3.merge(F_A3)(ZZ).default_prec
            3
            sage: F_A3.merge(F_A5)(ZZ).default_prec
            3
            sage: F_A5.merge(F_AN)(ZZ).default_prec
            5
            sage: F_A5.merge(F_A3)(ZZ).default_prec
            3
            sage: F_A5.merge(F_A5)(ZZ).default_prec
            5

            sage: A = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F1 = A.construction()[0]
            sage: F2 = A.construction()[0]; F2._category_ = Rings()
            sage: F1.merge(F2)._category_
            Category of rings
        """
        if self == other:
            return self

        if isinstance(other, AsymptoticRingFunctor) and self.cls == other.cls:
            from sage.structure.element import get_coercion_model
            cm = get_coercion_model()
            try:
                G = cm.common_parent(self.growth_group, other.growth_group)
            except TypeError:
                pass
            else:
                if (self._default_prec_ is None
                    and other._default_prec_ is None):
                    default_prec = None
                elif self._default_prec_ is None:
                    default_prec = other._default_prec_
                elif other._default_prec_ is None:
                    default_prec = self._default_prec_
                else:
                    default_prec = min(self._default_prec_,
                                       other._default_prec_)
                if (self._category_ is None
                    and other._category_ is None):
                    category = None
                elif self._category_ is None:
                    category = other._category_
                elif other._category_ is None:
                    category = self._category_
                else:
                    category = self._category_ | other._category_

                return AsymptoticRingFunctor(
                    G,
                    default_prec=default_prec,
                    category=category,
                    cls=self.cls)


    def __eq__(self, other):
        r"""
        Return whether this functor is equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_X = X.construction()[0]
            sage: F_Y = Y.construction()[0]
            sage: F_X == F_X
            True
            sage: F_X == F_Y
            False
        """
        return (type(self) == type(other)
                and self.growth_group == other.growth_group
                and self._default_prec_ == other._default_prec_
                and self._category_ == other._category_
                and self.cls == other.cls)


    def __ne__(self, other):
        r"""
        Return whether this functor is not equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_X = X.construction()[0]
            sage: F_Y = Y.construction()[0]
            sage: F_X != F_X
            False
            sage: F_X != F_Y
            True
        """
        return not self == other
