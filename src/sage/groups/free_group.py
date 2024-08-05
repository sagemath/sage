"""
Free Groups

Free groups and finitely presented groups are implemented as a wrapper
over the corresponding GAP objects.

A free group can be created by giving the number of generators, or their names.
It is also possible to create indexed generators::

    sage: G.<x,y,z> = FreeGroup();  G
    Free Group on generators {x, y, z}
    sage: FreeGroup(3)
    Free Group on generators {x0, x1, x2}
    sage: FreeGroup('a,b,c')
    Free Group on generators {a, b, c}
    sage: FreeGroup(3,'t')
    Free Group on generators {t0, t1, t2}

The elements can be created by operating with the generators, or by passing a list
with the indices of the letters to the group:

EXAMPLES::

    sage: G.<a,b,c> = FreeGroup()
    sage: a*b*c*a
    a*b*c*a
    sage: G([1,2,3,1])
    a*b*c*a
    sage: a * b / c * b^2
    a*b*c^-1*b^2
    sage: G([1,1,2,-1,-3,2])
    a^2*b*a^-1*c^-1*b

You can use call syntax to replace the generators with a set of
arbitrary ring elements::

    sage: g =  a * b / c * b^2
    sage: g(1,2,3)
    8/3
    sage: M1 = identity_matrix(2)
    sage: M2 = matrix([[1,1],[0,1]])
    sage: M3 = matrix([[0,1],[1,0]])
    sage: g([M1, M2, M3])
    [1 3]
    [1 2]

AUTHORS:

- Miguel Angel Marco Buzunariz
- Volker Braun
"""

# ****************************************************************************
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.groups import Groups
from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.structure.unique_representation import UniqueRepresentation
from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.structure.sequence import Sequence
from sage.structure.element import coercion_model, parent


def is_FreeGroup(x):
    """
    Test whether ``x`` is a :class:`FreeGroup_class`.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.groups.free_group import is_FreeGroup
        sage: is_FreeGroup('a string')                                                  # needs sage.combinat
        False
        sage: is_FreeGroup(FreeGroup(0))
        True
        sage: is_FreeGroup(FreeGroup(index_set=ZZ))                                     # needs sage.combinat
        True
    """
    if isinstance(x, FreeGroup_class):
        return True
    from sage.groups.indexed_free_group import IndexedFreeGroup
    return isinstance(x, IndexedFreeGroup)


def _lexi_gen(zeroes=False):
    """
    Return a generator object that produces variable names suitable for the
    generators of a free group.

    INPUT:

    - ``zeroes`` -- boolean (default: ``False``); if ``True``, the
      integers appended to the output string begin at zero at the
      first iteration through the alphabet

    OUTPUT:

    Python generator object which outputs a character from the alphabet on each
    ``next()`` call in lexicographical order. The integer `i` is appended
    to the output string on the `i`-th iteration through the alphabet.

    EXAMPLES::

        sage: from sage.groups.free_group import _lexi_gen
        sage: itr = _lexi_gen()
        sage: F = FreeGroup([next(itr) for i in [1..10]]); F
        Free Group on generators {a, b, c, d, e, f, g, h, i, j}
        sage: it = _lexi_gen()
        sage: [next(it) for i in range(10)]
        ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
        sage: itt = _lexi_gen(True)
        sage: [next(itt) for i in range(10)]
        ['a0', 'b0', 'c0', 'd0', 'e0', 'f0', 'g0', 'h0', 'i0', 'j0']
        sage: test = _lexi_gen()
        sage: ls = [next(test) for i in range(3*26)]
        sage: ls[2*26:3*26]
        ['a2', 'b2', 'c2', 'd2', 'e2', 'f2', 'g2', 'h2', 'i2', 'j2', 'k2', 'l2', 'm2',
         'n2', 'o2', 'p2', 'q2', 'r2', 's2', 't2', 'u2', 'v2', 'w2', 'x2', 'y2', 'z2']

    TESTS::

        sage: from sage.groups.free_group import _lexi_gen
        sage: test = _lexi_gen()
        sage: ls = [next(test) for i in range(500)]
        sage: ls[234], ls[260]
        ('a9', 'a10')
    """
    count = Integer(0)
    while True:
        mwrap, ind = count.quo_rem(26)
        if mwrap == 0 and not zeroes:
            name = ''
        else:
            name = str(mwrap)
        name = chr(ord('a') + ind) + name
        yield name
        count += 1


class FreeGroupElement(ElementLibGAP):
    """
    A wrapper of GAP's Free Group elements.

    INPUT:

    - ``x`` -- something that determines the group element; either a
      :class:`~sage.libs.gap.element.GapElement` or the Tietze list
      (see :meth:`Tietze`) of the group element

    - ``parent`` -- the parent :class:`FreeGroup`

    EXAMPLES::

        sage: G = FreeGroup('a, b')
        sage: x = G([1, 2, -1, -2])
        sage: x
        a*b*a^-1*b^-1
        sage: y = G([2, 2, 2, 1, -2, -2, -2])
        sage: y
        b^3*a*b^-3
        sage: x*y
        a*b*a^-1*b^2*a*b^-3
        sage: y*x
        b^3*a*b^-3*a*b*a^-1*b^-1
        sage: x^(-1)
        b*a*b^-1*a^-1
        sage: x == x*y*y^(-1)
        True
    """

    def __init__(self, parent, x):
        """
        The Python constructor.

        See :class:`FreeGroupElement` for details.

        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: x = G([1, 2, -1, -1])
            sage: x # indirect doctest
            a*b*a^-2
            sage: y = G([2, 2, 2, 1, -2, -2, -1])
            sage: y # indirect doctest
            b^3*a*b^-2*a^-1

            sage: TestSuite(G).run()
            sage: TestSuite(x).run()
        """
        if not isinstance(x, GapElement):
            try:
                l = x.Tietze()
            except AttributeError:
                l = list(x)
            if len(l) > 0:
                if min(l) < -parent.ngens() or parent.ngens() < max(l):
                    raise ValueError('generators not in the group')
            if 0 in l:
                raise ValueError('zero does not denote a generator')
            i = 0
            while i < len(l)-1:
                if l[i] == -l[i+1]:
                    l.pop(i)
                    l.pop(i)
                    if i > 0:
                        i = i-1
                else:
                    i = i+1
            AbstractWordTietzeWord = libgap.eval('AbstractWordTietzeWord')
            x = AbstractWordTietzeWord(l, parent.gap().GeneratorsOfGroup())
        ElementLibGAP.__init__(self, parent, x)

    def __hash__(self):
        r"""
        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: hash(a*b*b*~a) == hash((1, 2, 2, -1))
            True
        """
        return hash(self.Tietze())

    def _latex_(self):
        r"""
        Return a LaTeX representation.

        OUTPUT: string; a valid LaTeX math command sequence

        EXAMPLES::

            sage: F.<a,b,c> = FreeGroup()
            sage: f = F([1, 2, 2, -3, -1]) * c^15 * a^(-23)
            sage: f._latex_()
            'a\\cdot b^{2}\\cdot c^{-1}\\cdot a^{-1}\\cdot c^{15}\\cdot a^{-23}'

            sage: F = FreeGroup(3)
            sage: f = F([1, 2, 2, -3, -1]) * F.gen(2)^11 * F.gen(0)^(-12)
            sage: f._latex_()
            'x_{0}\\cdot x_{1}^{2}\\cdot x_{2}^{-1}\\cdot x_{0}^{-1}\\cdot x_{2}^{11}\\cdot x_{0}^{-12}'

            sage: F.<a,b,c> = FreeGroup()
            sage: G = F /  (F([1, 2, 1, -3, 2, -1]), F([2, -1]))
            sage: f = G([1, 2, 2, -3, -1]) * G.gen(2)^15 * G.gen(0)^(-23)
            sage: f._latex_()
            'a\\cdot b^{2}\\cdot c^{-1}\\cdot a^{-1}\\cdot c^{15}\\cdot a^{-23}'

            sage: F = FreeGroup(4)
            sage: G = F.quotient((F([1, 2, 4, -3, 2, -1]), F([2, -1])))
            sage: f = G([1, 2, 2, -3, -1]) * G.gen(3)^11 * G.gen(0)^(-12)
            sage: f._latex_()
            'x_{0}\\cdot x_{1}^{2}\\cdot x_{2}^{-1}\\cdot x_{0}^{-1}\\cdot x_{3}^{11}\\cdot x_{0}^{-12}'
        """
        import re
        s = self._repr_()
        s = re.sub('([a-z]|[A-Z])([0-9]+)', r'\g<1>_{\g<2>}', s)
        s = re.sub(r'(\^)(-)([0-9]+)', r'\g<1>{\g<2>\g<3>}', s)
        s = re.sub(r'(\^)([0-9]+)', r'\g<1>{\g<2>}', s)
        s = s.replace('*', r'\cdot ')
        return s

    def __reduce__(self):
        """
        Implement pickling.

        TESTS::

            sage: F.<a,b> = FreeGroup()
            sage: a.__reduce__()
            (Free Group on generators {a, b}, ((1,),))
            sage: (a*b*a^-1).__reduce__()
            (Free Group on generators {a, b}, ((1, 2, -1),))
        """
        return (self.parent(), (self.Tietze(),))

    @cached_method
    def Tietze(self):
        """
        Return the Tietze list of the element.

        The Tietze list of a word is a list of integers that represent
        the letters in the word.  A positive integer `i` represents
        the letter corresponding to the `i`-th generator of the group.
        Negative integers represent the inverses of generators.

        OUTPUT: tuple of integers

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: a.Tietze()
            (1,)
            sage: x = a^2 * b^(-3) * a^(-2)
            sage: x.Tietze()
            (1, 1, -2, -2, -2, -1, -1)

        TESTS::

            sage: type(a.Tietze())
            <... 'tuple'>
            sage: type(a.Tietze()[0])
            <class 'sage.rings.integer.Integer'>
        """
        tl = self.gap().TietzeWordAbstractWord()
        return tuple(tl.sage())

    def fox_derivative(self, gen, im_gens=None, ring=None):
        r"""
        Return the Fox derivative of ``self`` with respect to a given
        generator ``gen`` of the free group.

        Let `F` be a free group with free generators
        `x_1, x_2, \ldots, x_n`. Let `j \in \left\{ 1, 2, \ldots ,
        n \right\}`. Let `a_1, a_2, \ldots, a_n` be `n`
        invertible elements of a ring `A`. Let `a : F \to A^\times`
        be the (unique) homomorphism from `F` to the multiplicative
        group of invertible elements of `A` which sends each `x_i`
        to `a_i`. Then, we can define a map `\partial_j : F \to A`
        by the requirements that

        .. MATH::

            \partial_j (x_i) = \delta_{i, j}
            \qquad \qquad \text{ for all indices } i \text{ and } j

        and

        .. MATH::

            \partial_j (uv) = \partial_j(u) + a(u) \partial_j(v)
            \qquad \qquad \text{ for all } u, v \in F .

        This map `\partial_j` is called the `j`-th Fox derivative
        on `F` induced by `(a_1, a_2, \ldots, a_n)`.

        The most well-known case is when `A` is the group ring
        `\ZZ [F]` of `F` over `\ZZ`, and when `a_i = x_i \in A`.
        In this case, `\partial_j` is simply called the `j`-th
        Fox derivative on `F`.

        INPUT:

        - ``gen`` -- the generator with respect to which the
          derivative will be computed. If this is `x_j`, then the
          method will return `\partial_j`.

        - ``im_gens`` -- (optional) the images of the generators
          (given as a list or iterable). This is the list
          `(a_1, a_2, \ldots, a_n)`.
          If not provided, it defaults to
          `(x_1, x_2, \ldots, x_n)` in the group ring
          `\ZZ [F]`.

        - ``ring`` -- (optional) the ring in which the elements
          of the list  `(a_1, a_2, \ldots, a_n)` lie. If not
          provided, this ring is inferred from these elements.

        OUTPUT:

        The fox derivative of ``self`` with respect to ``gen``
        (induced by ``im_gens``).
        By default, it is an element of the group algebra with
        integer coefficients.
        If ``im_gens`` are provided, the result lives in the
        algebra where ``im_gens`` live.

        EXAMPLES::

            sage: G = FreeGroup(5)
            sage: G.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x0)
            -x0^-1 + x0^-1*x1 - x0^-1*x1*x0*x2*x0^-1
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x1)
            x0^-1
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x2)
            x0^-1*x1*x0
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x3)
            0

        If ``im_gens`` is given, the images of the generators are
        mapped to them::

            sage: F = FreeGroup(3)
            sage: a = F([2,1,3,-1,2])
            sage: a.fox_derivative(F([1]))
            x1 - x1*x0*x2*x0^-1
            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: a.fox_derivative(F([1]),[t,t,t])
            t - t^2
            sage: S.<t1,t2,t3> = LaurentPolynomialRing(ZZ)
            sage: a.fox_derivative(F([1]),[t1,t2,t3])
            -t2*t3 + t2
            sage: R.<x,y,z> = QQ[]
            sage: a.fox_derivative(F([1]),[x,y,z])
            -y*z + y
            sage: a.inverse().fox_derivative(F([1]),[x,y,z])
            (z - 1)/(y*z)

        The optional parameter ``ring`` determines the ring `A`::

            sage: u = a.fox_derivative(F([1]), [1,2,3], ring=QQ)
            sage: u
            -4
            sage: parent(u)
            Rational Field
            sage: u = a.fox_derivative(F([1]), [1,2,3], ring=R)
            sage: u
            -4
            sage: parent(u)
            Multivariate Polynomial Ring in x, y, z over Rational Field

        TESTS::

            sage: F = FreeGroup(3)
            sage: a = F([])
            sage: a.fox_derivative(F([1]))
            0
            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: a.fox_derivative(F([1]),[t,t,t])
            0
        """
        if gen not in self.parent().generators():
            raise ValueError("Fox derivative can only be computed with respect to generators of the group")
        l = list(self.Tietze())
        if im_gens is None:
            F = self.parent()
            R = F.algebra(IntegerRing())
            R_basis = R.basis()
            symb = [R_basis[a] for a in F.gens()]
            symb += reversed([R_basis[a.inverse()] for a in F.gens()])
            if ring is not None:
                R = ring
                symb = [R(i) for i in symb]
        else:
            if ring is None:
                R = Sequence(im_gens).universe()
            else:
                R = ring
            symb = list(im_gens)
            symb += reversed([a**(-1) for a in im_gens])
        i = gen.Tietze()[0]  # So ``gen`` is the `i`-th
                             # generator of the free group.
        a = R.zero()
        coef = R.one()
        while l:
            b = l.pop(0)
            if b == i:
                a += coef * R.one()
                coef *= symb[b-1]
            elif b == -i:
                a -= coef * symb[b]
                coef *= symb[b]
            elif b > 0:
                coef *= symb[b-1]
            else:
                coef *= symb[b]
        return a

    @cached_method
    def syllables(self):
        r"""
        Return the syllables of the word.

        Consider a free group element `g = x_1^{n_1} x_2^{n_2} \cdots
        x_k^{n_k}`. The uniquely-determined subwords `x_i^{e_i}`
        consisting only of powers of a single generator are called the
        syllables of `g`.

        OUTPUT:

        The tuple of syllables. Each syllable is given as a pair
        `(x_i, e_i)` consisting of a generator and a nonzero integer.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: w = a^2 * b^-1 * a^3
            sage: w.syllables()
            ((a, 2), (b, -1), (a, 3))
        """
        g = self.gap().UnderlyingElement()
        k = g.NumberSyllables().sage()
        exponent_syllable = libgap.eval('ExponentSyllable')
        generator_syllable = libgap.eval('GeneratorSyllable')
        result = []
        gen = self.parent().gen
        for i in range(k):
            exponent = exponent_syllable(g, i+1).sage()
            generator = gen(generator_syllable(g, i+1).sage() - 1)
            result.append( (generator, exponent) )
        return tuple(result)

    def __call__(self, *values):
        """
        Replace the generators of the free group by corresponding
        elements of the iterable ``values`` in the group element
        ``self``.

        INPUT:

        - ``*values`` -- a sequence of values, or a
          list/tuple/iterable of the same length as the number of
          generators of the free group

        OUTPUT: the product of ``values`` in the order and with exponents
        specified by ``self``

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: w = a^2 * b^-1 * a^3
            sage: w(1, 2)
            1/2
            sage: w(2, 1)
            32
            sage: w.subs(b=1, a=2)   # indirect doctest
            32

        TESTS::

            sage: w([1, 2])
            1/2
            sage: w((1, 2))
            1/2
            sage: w(i+1 for i in range(2))
            1/2

        Check that :issue:`25017` is fixed::

            sage: F = FreeGroup(2)
            sage: x0, x1 = F.gens()
            sage: u = F(1)
            sage: parent(u.subs({x1:x0})) is F
            True

            sage: F = FreeGroup(2)
            sage: x0, x1 = F.gens()
            sage: u = x0*x1
            sage: u.subs({x0:3, x1:2})
            6
            sage: u.subs({x0:1r, x1:2r})
            2
            sage: M0 = matrix(ZZ,2,[1,1,0,1])
            sage: M1 = matrix(ZZ,2,[1,0,1,1])
            sage: u.subs({x0: M0, x1: M1})
            [2 1]
            [1 1]

        TESTS::

            sage: F.<x,y> = FreeGroup()
            sage: F.one().subs(x=x, y=1)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Free Group on generators {x, y}' and 'Integer Ring'
        """
        if len(values) == 1:
            try:
                values = list(values[0])
            except TypeError:
                pass
        G = self.parent()
        if len(values) != G.ngens():
            raise ValueError('number of values has to match the number of generators')
        replace = dict(zip(G.gens(), values))
        new_parent = coercion_model.common_parent(*[parent(v) for v in values])
        try:
            return new_parent.prod(replace[gen] ** power
                                   for gen, power in self.syllables())
        except AttributeError:
            return prod(new_parent(replace[gen]) ** power
                        for gen, power in self.syllables())


def FreeGroup(n=None, names='x', index_set=None, abelian=False, **kwds):
    """
    Construct a Free Group.

    INPUT:

    - ``n`` -- integer (default: ``None``); the number of
      generators (if not specified the ``names`` are counted)

    - ``names`` -- string or list/tuple/iterable of strings (default:
      ``'x'``); the generator names or name prefix

    - ``index_set`` -- (optional) an index set for the generators; if
      specified then the optional keyword ``abelian`` can be used

    - ``abelian`` -- boolean (default: ``False``); whether to construct a free
      abelian group or a free group

    .. NOTE::

        If you want to create a free group, it is currently preferential to
        use ``Groups().free(...)`` as that does not load GAP.

    EXAMPLES::

        sage: G.<a,b> = FreeGroup();  G
        Free Group on generators {a, b}
        sage: H = FreeGroup('a, b')
        sage: G is H
        True
        sage: FreeGroup(0)
        Free Group on generators {}

    The entry can be either a string with the names of the generators,
    or the number of generators and the prefix of the names to be
    given. The default prefix is ``'x'`` ::

        sage: FreeGroup(3)
        Free Group on generators {x0, x1, x2}
        sage: FreeGroup(3, 'g')
        Free Group on generators {g0, g1, g2}
        sage: FreeGroup()
        Free Group on generators {x}

    We give two examples using the ``index_set`` option::

        sage: FreeGroup(index_set=ZZ)                                                   # needs sage.combinat
        Free group indexed by Integer Ring
        sage: FreeGroup(index_set=ZZ, abelian=True)                                     # needs sage.combinat
        Free abelian group indexed by Integer Ring

    TESTS::

        sage: G1 = FreeGroup(2, 'a,b')
        sage: G2 = FreeGroup('a,b')
        sage: G3.<a,b> = FreeGroup()
        sage: G1 is G2, G2 is G3
        (True, True)
    """
    # Support Freegroup('a,b') syntax
    if n is not None:
        try:
            n = Integer(n)
        except TypeError:
            names = n
            n = None
    # derive n from counting names
    if n is None:
        if isinstance(names, str):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)
    from sage.structure.category_object import normalize_names
    names = normalize_names(n, names)
    if index_set is not None or abelian:
        if abelian:
            from sage.groups.indexed_free_group import IndexedFreeAbelianGroup
            return IndexedFreeAbelianGroup(index_set, names=names, **kwds)

        from sage.groups.indexed_free_group import IndexedFreeGroup
        return IndexedFreeGroup(index_set, names=names, **kwds)
    return FreeGroup_class(names)


def wrap_FreeGroup(libgap_free_group):
    """
    Wrap a LibGAP free group.

    This function changes the comparison method of
    ``libgap_free_group`` to comparison by Python ``id``. If you want
    to put the LibGAP free group into a container (set, dict) then you
    should understand the implications of
    :meth:`~sage.libs.gap.element.GapElement._set_compare_by_id`. To
    be safe, it is recommended that you just work with the resulting
    Sage :class:`FreeGroup_class`.

    INPUT:

    - ``libgap_free_group`` -- a LibGAP free group

    OUTPUT: a Sage :class:`FreeGroup_class`

    EXAMPLES:

    First construct a LibGAP free group::

        sage: F = libgap.FreeGroup(['a', 'b'])
        sage: type(F)
        <class 'sage.libs.gap.element.GapElement'>

    Now wrap it::

        sage: from sage.groups.free_group import wrap_FreeGroup
        sage: wrap_FreeGroup(F)
        Free Group on generators {a, b}

    TESTS:

    Check that we can do it twice (see :issue:`12339`) ::

        sage: G = libgap.FreeGroup(['a', 'b'])
        sage: wrap_FreeGroup(G)
        Free Group on generators {a, b}
    """
    assert libgap_free_group.IsFreeGroup()
    libgap_free_group._set_compare_by_id()
    names = tuple( str(g) for g in libgap_free_group.GeneratorsOfGroup() )
    return FreeGroup_class(names, libgap_free_group)


class FreeGroup_class(UniqueRepresentation, Group, ParentLibGAP):
    """
    A class that wraps GAP's FreeGroup.

    See :func:`FreeGroup` for details.

    TESTS::

        sage: G = FreeGroup('a, b')
        sage: TestSuite(G).run()
        sage: G.category()
        Category of infinite groups
    """
    Element = FreeGroupElement

    def __init__(self, generator_names, libgap_free_group=None):
        """
        Python constructor.

        INPUT:

        - ``generator_names`` -- tuple of strings; the names of the
          generators

        - ``libgap_free_group`` -- a LibGAP free group (default: ``None``);
          the LibGAP free group to wrap (if ``None``, a suitable group will be
          constructed)

        TESTS::

            sage: G.<a,b> = FreeGroup() # indirect doctest
            sage: G
            Free Group on generators {a, b}
            sage: G.variable_names()
            ('a', 'b')
        """
        self._assign_names(generator_names)
        if libgap_free_group is None:
            libgap_free_group = libgap.FreeGroup(generator_names)
        ParentLibGAP.__init__(self, libgap_free_group)
        if not generator_names:
            cat = Groups().Finite()
        else:
            cat = Groups().Infinite()
        Group.__init__(self, category=cat)

    def _repr_(self):
        """
        TESTS::

            sage: G = FreeGroup('a, b')
            sage: G # indirect doctest
            Free Group on generators {a, b}
            sage: G._repr_()
            'Free Group on generators {a, b}'
        """
        return 'Free Group on generators {' + ', '.join(self.variable_names()) + '}'

    def rank(self):
        """
        Return the number of generators of ``self``.

        Alias for :meth:`ngens`.

        OUTPUT: integer

        EXAMPLES::

            sage: G = FreeGroup('a, b');  G
            Free Group on generators {a, b}
            sage: G.rank()
            2
            sage: H = FreeGroup(3, 'x')
            sage: H
            Free Group on generators {x0, x1, x2}
            sage: H.rank()
            3
        """
        return self.ngens()

    def _gap_init_(self):
        """
        Return the string used to construct the object in gap.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G._gap_init_()
            'FreeGroup(["x0", "x1", "x2"])'
        """
        gap_names = [ '"' + s + '"' for s in self.variable_names() ]
        gen_str = ', '.join(gap_names)
        return 'FreeGroup(['+gen_str+'])'

    def _element_constructor_(self, *args, **kwds):
        """
        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: G([1, 2, 1]) # indirect doctest
            a*b*a
            sage: G([1, 2, -2, 1, 1, -2]) # indirect doctest
            a^3*b^-1

            sage: G( a.gap() )
            a
            sage: type(_)
            <class 'sage.groups.free_group.FreeGroup_class_with_category.element_class'>

        Check that conversion between free groups follow the convention that
        names are preserved::

            sage: F = FreeGroup('a,b')
            sage: G = FreeGroup('b,a')
            sage: G(F.gen(0))
            a
            sage: F(G.gen(0))
            b
            sage: a,b = F.gens()
            sage: G(a^2*b^-3*a^-1)
            a^2*b^-3*a^-1

        Check that :issue:`17246` is fixed::

            sage: F = FreeGroup(0)
            sage: F([])
            1

        Check that 0 isn't considered the identity::

            sage: F = FreeGroup('x')
            sage: F(0)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        if len(args) != 1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        if x == 1 or x == [] or x == ():
            return self.one()
        try:
            P = x.parent()
        except AttributeError:
            return self.element_class(self, x, **kwds)
        if isinstance(P, FreeGroup_class):
            names = {P._names[abs(i)-1] for i in x.Tietze()}
            if names.issubset(self._names):
                return self([i.sign()*(self._names.index(P._names[abs(i)-1])+1)
                             for i in x.Tietze()])
            else:
                raise ValueError('generators of %s not in the group' % x)
        return self.element_class(self, x, **kwds)

    def abelian_invariants(self):
        r"""
        Return the Abelian invariants of ``self``.

        The Abelian invariants are given by a list of integers
        `i_1 \dots i_j`, such that the abelianization of the
        group is isomorphic to

        .. MATH::

            \ZZ / (i_1) \times \dots \times \ZZ / (i_j)

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: F.abelian_invariants()
            (0, 0)
        """
        return (0,) * self.ngens()

    def quotient(self, relations, **kwds):
        """
        Return the quotient of ``self`` by the normal subgroup generated
        by the given elements.

        This quotient is a finitely presented groups with the same
        generators as ``self``, and relations given by the elements of
        ``relations``.

        INPUT:

        - ``relations`` -- list/tuple/iterable with the elements of
          the free group
        - further named arguments, that are passed to the constructor
          of a finitely presented group

        OUTPUT:

        A finitely presented group, with generators corresponding to
        the generators of the free group, and relations corresponding
        to the elements in ``relations``.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: F.quotient([a*b^2*a, b^3])
            Finitely presented group < a, b | a*b^2*a, b^3 >

        Division is shorthand for :meth:`quotient` ::

            sage: F /  [a*b^2*a, b^3]
            Finitely presented group < a, b | a*b^2*a, b^3 >

        Relations are converted to the free group, even if they are not
        elements of it (if possible) ::

            sage: F1.<a,b,c,d> = FreeGroup()
            sage: F2.<a,b> = FreeGroup()
            sage: r = a*b/a
            sage: r.parent()
            Free Group on generators {a, b}
            sage: F1/[r]
            Finitely presented group < a, b, c, d | a*b*a^-1 >
        """
        from sage.groups.finitely_presented import FinitelyPresentedGroup
        return FinitelyPresentedGroup(self, tuple(map(self, relations) ), **kwds)

    __truediv__ = quotient
