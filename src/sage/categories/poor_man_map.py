# sage_setup: distribution = sagemath-categories
r"""
Poor Man's map
"""
# ****************************************************************************
#       Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2016 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sage.structure.sage_object


class PoorManMap(sage.structure.sage_object.SageObject):
    """
    A class for maps between sets which are not (yet) modeled by parents.

    Could possibly disappear when all combinatorial classes / enumerated sets will be parents.

    INPUT:

    - ``function`` -- a callable or an iterable of callables. This represents
      the underlying function used to implement this map. If it is an iterable,
      then the callables will be composed to implement this map.

    - ``domain`` -- the domain of this map or ``None`` if the domain is not
      known or should remain unspecified

    - ``codomain`` -- the codomain of this map or ``None`` if the codomain is
      not known or should remain unspecified

    - ``name`` -- a name for this map or ``None`` if this map has no particular
      name

    EXAMPLES::

        sage: from sage.categories.poor_man_map import PoorManMap
        sage: f = PoorManMap(factorial, domain=(1, 2, 3), codomain=(1, 2, 6))
        sage: f
        A map from (1, 2, 3) to (1, 2, 6)
        sage: f(3)
        6

    The composition of several functions can be created by passing in a tuple
    of functions::

        sage: i = PoorManMap((factorial, sqrt), domain=(1, 4, 9), codomain=(1, 2, 6))

    However, the same effect can also be achieved by just composing maps::

        sage: g = PoorManMap(factorial, domain=(1, 2, 3), codomain=(1, 2, 6))
        sage: h = PoorManMap(sqrt, domain=(1, 4, 9), codomain=(1, 2, 3))
        sage: i == g*h
        True
    """
    def __init__(self, function, domain=None, codomain=None, name=None):
        """
        TESTS::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain=(1, 2, 3), codomain=(1, 2, 6))
            sage: g = PoorManMap(sqrt, domain=(1, 4, 9), codomain=(1, 2, 6))

            sage: TestSuite(f).run()
            sage: TestSuite(f*g).run()
        """
        from collections.abc import Iterable
        if not isinstance(function, Iterable):
            function = (function,)
        self._functions = tuple(function)
        self._domain = domain
        self._codomain = codomain
        self._name = name

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: PoorManMap(lambda x: x+2)  # indirect doctest
            A map
            sage: PoorManMap(lambda x: x+2, domain=(1,2,3))
            A map from (1, 2, 3)
            sage: PoorManMap(lambda x: x+2, domain=(1,2,3))
            A map from (1, 2, 3)
            sage: PoorManMap(lambda x: x+2, codomain=(3,4,5))
            A map to (3, 4, 5)
        """
        return ((self._name if self._name is not None else "A map") +
                (" from %s" % (self._domain,) if self._domain is not None else "") +
                (" to %s" % (self._codomain,) if self._codomain is not None else ""))

    def domain(self):
        """
        Return the domain of ``self``.

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: PoorManMap(lambda x: x+1, domain=(1,2,3), codomain=(2,3,4)).domain()
            (1, 2, 3)
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of ``self``.

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: PoorManMap(lambda x: x+1, domain=(1,2,3), codomain=(2,3,4)).codomain()
            (2, 3, 4)
        """
        return self._codomain

    def __eq__(self, other):
        r"""
        Return whether this map is equal to ``other``.

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6))
            sage: g = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6))
            sage: h1 = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6,8))
            sage: h2 = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6,8))
            sage: h3 = PoorManMap(factorial, domain=(1,2,3,4), codomain=(1,2,6))
            sage: h4 = PoorManMap(lambda x: x, domain=(1,2,3), codomain=(1,2,6))
            sage: f == g, f == h1, f == h2, f == h3, f == h4, f == 1, 1 == f
            (True, False, False, False, False, False, False)
        """
        if isinstance(other, PoorManMap):
            return (self._functions == other._functions
                    and self._domain == other._domain
                    and self._codomain == other._codomain
                    and self._name == other._name)
        else:
            return False

    def __ne__(self, other):
        r"""
        Return whether this map is not equal to ``other``.

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6))
            sage: g = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6))
            sage: h1 = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6,8))
            sage: h2 = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6,8))
            sage: h3 = PoorManMap(factorial, domain=(1,2,3,4), codomain=(1,2,6))
            sage: h4 = PoorManMap(lambda x: x, domain=(1,2,3), codomain=(1,2,6))
            sage: f != g, f != h1, f != h2, f != h3, f != h4, f != 1, 1 != f
            (False, True, True, True, True, True, True)
        """
        return not (self == other)

    def __hash__(self):
        r"""
        Return a hash value for this map.

        TESTS::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6))
            sage: g = PoorManMap(factorial, domain=(1,2,3), codomain=(1,2,6))
            sage: hash(f) == hash(g)
            True
        """
        return hash((self._functions, self._domain, self._codomain, self._name))

    def __mul__(self, other):
        r"""
        Composition.

        INPUT:

        - ``self`` -- a map `f`
        - ``other`` -- a map `g`

        Returns the composition map `f\circ g` of `f` and `g`

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(lambda x: x+1, domain=(1,2,3), codomain=(2,3,4))
            sage: g = PoorManMap(lambda x: -x,  domain=(2,3,4), codomain=(-2,-3,-4))
            sage: g*f
            A map from (1, 2, 3) to (-2, -3, -4)

        Note that the compatibility of the domains and codomains is for performance
        reasons only checked for proper parents. For example, the incompatibility
        is not detected here::

            sage: f*g
            A map from (2, 3, 4) to (2, 3, 4)

        But it is detected here::

            sage: g = PoorManMap(factorial, domain=ZZ, codomain=ZZ)
            sage: h = PoorManMap(sqrt, domain=RR, codomain=CC)                          # needs sage.rings.real_mpfr
            sage: g*h                                                                   # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            ValueError: the codomain Complex Field with 53 bits of precision
            does not coerce into the domain Integer Ring
            sage: h*g                                                                   # needs sage.rings.real_mpfr
            A map from Integer Ring to Complex Field with 53 bits of precision
        """
        self_domain = self.domain()

        try:
            other_codomain = other.codomain()
        except AttributeError:
            other_codomain = None

        if self_domain is not None and other_codomain is not None:
            from sage.structure.parent import Parent
            if isinstance(self_domain, Parent) and isinstance(other_codomain, Parent):
                if not self_domain.has_coerce_map_from(other_codomain):
                    raise ValueError("the codomain %r does not coerce into the domain %r" % (other_codomain, self_domain))

        codomain = self.codomain()
        try:
            domain = other.domain()
        except AttributeError:
            domain = None

        if isinstance(other, PoorManMap):
            other = other._functions
        else:
            other = (other,)

        return PoorManMap(self._functions + other, domain=domain, codomain=codomain)

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(lambda x: x+1, domain=(1,2,3), codomain=(2,3,4))
            sage: f(2)
            3

            sage: g = PoorManMap(lambda x: -x,  domain=(2,3,4), codomain=(-2,-3,-4))
            sage: (g*f)(2)
            -3
        """
        for function in reversed(self._functions):
            args = [function(*args)]
        return args[0]

    def _sympy_(self):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: h = PoorManMap(sin, domain=RR, codomain=RR)
            sage: h._sympy_()                                                           # needs sympy sage.symbolic
            sin
        """
        from sympy import Lambda, sympify
        if len(self._functions) == 1:
            return sympify(self._functions[0])
        raise NotImplementedError
