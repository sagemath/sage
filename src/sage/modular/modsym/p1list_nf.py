# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.number_field
r"""
Lists of Manin symbols over number fields, elements of `\mathbb{P}^1(R/N)`

Lists of elements of `\mathbb{P}^1(R/N)` where `R` is the ring of integers of a number
field `K` and `N` is an integral ideal.

AUTHORS:

- Maite Aranes (2009): Initial version

EXAMPLES:

We define a P1NFList:

::

    sage: x = polygen(QQ, 'x')
    sage: k.<a> = NumberField(x^3 + 11)
    sage: N = k.ideal(5, a^2 - a + 1)
    sage: P = P1NFList(N); P
    The projective line over
     the ring of integers modulo the Fractional ideal (5, a^2 - a + 1)

List operations with the P1NFList:

::

    sage: len(P)
    26
    sage: [p for p in P]
    [M-symbol (0: 1) of level Fractional ideal (5, a^2 - a + 1),
    ...
    M-symbol (1: 2*a^2 + 2*a) of level Fractional ideal (5, a^2 - a + 1)]

The elements of the P1NFList are M-symbols:

::

    sage: type(P[2])
    <class 'sage.modular.modsym.p1list_nf.MSymbol'>

Definition of MSymbols:

::

    sage: alpha = MSymbol(N, 3, a^2); alpha
    M-symbol (3: a^2) of level Fractional ideal (5, a^2 - a + 1)

Find the index of the class of an M-Symbol `(c: d)` in the list:

::

    sage: i = P.index(alpha)
    sage: P[i].c*alpha.d - P[i].d*alpha.c in N
    True

Lift an MSymbol to a matrix in `SL(2, R)`:

::

    sage: alpha = MSymbol(N, a + 2, 3*a^2)
    sage: alpha.lift_to_sl2_Ok()
    [-1, 4*a^2 - 13*a + 23, a + 2, 5*a^2 + 3*a - 3]
    sage: Ok = k.ring_of_integers()
    sage: M = Matrix(Ok, 2, alpha.lift_to_sl2_Ok())
    sage: det(M)
    1
    sage: M[1][1] - alpha.d in N
    True

Lift an MSymbol from P1NFList to a matrix in `SL(2, R)`

::

    sage: P[3]
    M-symbol (1: -2*a) of level Fractional ideal (5, a^2 - a + 1)
    sage: P.lift_to_sl2_Ok(3)
    [0, -1, 1, -2*a]
"""
# ****************************************************************************
#       Copyright (C) 2009, Maite Aranes <M.T.Aranes@warwick.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.richcmp import richcmp_method, richcmp
from sage.structure.sage_object import SageObject

from sage.misc.search import search

_level_cache = {}  # The info stored here is used in the normalization of MSymbols.


def P1NFList_clear_level_cache():
    """
    Clear the global cache of data for the level ideals.

    EXAMPLES::

        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^3 + 11)
        sage: N = k.ideal(a+1)
        sage: alpha = MSymbol(N, 2*a^2, 5)
        sage: alpha.normalize()
        M-symbol (-4*a^2: 5*a^2) of level Fractional ideal (a + 1)
        sage: sage.modular.modsym.p1list_nf._level_cache
        {Fractional ideal (a + 1): (...)}
        sage: sage.modular.modsym.p1list_nf.P1NFList_clear_level_cache()
        sage: sage.modular.modsym.p1list_nf._level_cache
        {}
    """
    global _level_cache
    _level_cache = {}


@richcmp_method
class MSymbol(SageObject):
    r"""
    The constructor for an M-symbol over a number field.

    INPUT:

    - ``N`` -- integral ideal (the modulus or level)

    - ``c`` -- integral element of the underlying number field or an MSymbol of
      level N

    - ``d`` -- (optional) when present, it must be an integral element such
      that `\langle c\rangle + \langle d\rangle + N = R`, where `R` is the
      corresponding ring of integers

    - ``check`` -- boolean (default: ``True``); if ``check=False`` the
      constructor does not check the condition
      `\langle c\rangle + \langle d\rangle + N = R`

    OUTPUT:

    An M-symbol modulo the given ideal `N`, i.e. an element of the
    projective line `\mathbb{P}^1(R/N)`, where `R` is the ring of integers of
    the underlying number field.

    EXAMPLES::

        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^3 + 11)
        sage: N = k.ideal(a + 1, 2)
        sage: MSymbol(N, 3, a^2 + 1)
        M-symbol (3: a^2 + 1) of level Fractional ideal (2, a + 1)

    We can give a tuple as input:

    ::

        sage: MSymbol(N, (1, 0))
        M-symbol (1: 0) of level Fractional ideal (2, a + 1)

    We get an error if `\langle c\rangle`, `\langle d\rangle` and `N` are not coprime:

    ::

        sage: MSymbol(N, 2*a, a - 1)
        Traceback (most recent call last):
        ...
        ValueError: (2*a, a - 1) is not an element of P1(R/N).
        sage: MSymbol(N, (0, 0))
        Traceback (most recent call last):
        ...
        ValueError: (0, 0) is not an element of P1(R/N).

    Saving and loading works:

    ::

        sage: alpha = MSymbol(N, 3, a^2 + 1)
        sage: loads(dumps(alpha))==alpha
        True
    """
    def __init__(self, N, c, d=None, check=True):
        """
        See ``MSymbol`` for full documentation.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^4 + 13*x - 7)
            sage: N = k.ideal(5)
            sage: MSymbol(N, 0, 6*a)
            M-symbol (0: 6*a) of level Fractional ideal (5)
            sage: MSymbol(N, a^2 + 3, 7)
            M-symbol (a^2 + 3: 7) of level Fractional ideal (5)
        """
        k = N.number_field()
        R = k.ring_of_integers()
        self.__N = N
        if d is None:  # if we give a list (c, d) or an MSymbol as input
            if isinstance(c, MSymbol):
                if c.N() is N:
                    c1 = R(c[0])
                    d1 = R(c[1])
                else:
                    raise ValueError("Cannot change level of an MSymbol")
            else:
                try:
                    c1 = R(c[0])
                    d1 = R(c[1])
                except (ValueError, TypeError):
                    raise TypeError("Unable to create a Manin symbol from %s" % c)
        else:
            try:
                c1 = R(c)
                d1 = R(d)
            except (ValueError, TypeError):
                raise TypeError("Unable to create a Manin symbol from (%s, %s)" % (c, d))
        if check:
            if (c1.is_zero() and d1.is_zero()) or not N.is_coprime(k.ideal(c1, d1)):
                raise ValueError("(%s, %s) is not an element of P1(R/N)." % (c1, d1))
        self.__c, self.__d = (c1, d1)

    def __repr__(self):
        r"""
        Return the string representation of this MSymbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: MSymbol(N, 3, a)
            M-symbol (3: a) of level Fractional ideal (3, 1/2*a - 1/2)
        """
        return "M-symbol (%s: %s) of level %s" % (self.__c, self.__d, self.__N)

    def _latex_(self):
        r"""
        Return latex representation of ``self``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^4 + 13*x - 7)
            sage: N = k.ideal(a^3 - 1)
            sage: alpha = MSymbol(N, 3, 5*a^2 - 1)
            sage: latex(alpha) # indirect doctest
            \(3: 5 a^{2} - 1\)
        """
        return r"\(%s: %s\)" % (self.c._latex_(), self.d._latex_())

    def __richcmp__(self, other, op):
        r"""
        Comparison function for objects of the class MSymbol.

        The order is the same as for the underlying lists of lists.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: alpha = MSymbol(N, 3, a)
            sage: beta = MSymbol(N, 1, 0)
            sage: alpha < beta
            False
            sage: beta = MSymbol(N, 3, a + 1)
            sage: alpha < beta
            True
        """
        if not isinstance(other, MSymbol):
            raise ValueError("You can only compare with another M-symbol")
        return richcmp([self.__c.list(), self.__d.list()],
                       [other.__c.list(), other.__d.list()], op)

    def N(self):
        r"""
        Return the level or modulus of this MSymbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: alpha = MSymbol(N, 3, a)
            sage: alpha.N()
            Fractional ideal (3, 1/2*a - 1/2)
        """
        return self.__N

    def tuple(self):
        r"""
        Return the :class:`MSymbol` as a list `(c, d)`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: alpha = MSymbol(N, 3, a); alpha
            M-symbol (3: a) of level Fractional ideal (3, 1/2*a - 1/2)
            sage: alpha.tuple()
            (3, a)
        """
        return self.__c, self.__d

    def __getitem__(self, n):
        r"""
        Indexing function for the list defined by an M-symbol.

        INPUT:

        - ``n`` -- integer (0 or 1, since the list defined by an M-symbol has
          length 2)

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: alpha = MSymbol(N, 3, a); alpha
            M-symbol (3: a) of level Fractional ideal (3, 1/2*a - 1/2)
            sage: alpha[0]
            3
            sage: alpha[1]
            a
        """
        return self.tuple()[n]

    def __get_c(self):
        r"""
        Return the first coefficient of the M-symbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(a + 1, 2)
            sage: alpha = MSymbol(N, 3, a^2 + 1)
            sage: alpha.c # indirect doctest
            3
        """
        return self.__c
    c = property(__get_c)

    def __get_d(self):
        r"""
        Return the second coefficient of the M-symbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(a + 1, 2)
            sage: alpha = MSymbol(N, 3, a^2 + 1)
            sage: alpha.d # indirect doctest
            a^2 + 1
        """
        return self.__d
    d = property(__get_d)

    def lift_to_sl2_Ok(self):
        r"""
        Lift the :class:`MSymbol` to an element of `SL(2, O_k)`, where `O_k` is the ring
        of integers of the corresponding number field.

        OUTPUT:

        A list of integral elements `[a, b, c', d']` that are the entries of
        a `2\times 2` matrix with determinant 1. The lower two entries are congruent
        (modulo the level) to the coefficients `c`, `d` of the :class:`MSymbol` ``self``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: alpha = MSymbol(N, 3*a + 1, a)
            sage: alpha.lift_to_sl2_Ok()
            [0, -1, 1, a]
        """
        return lift_to_sl2_Ok(self.__N, self.__c, self.__d)

    def normalize(self, with_scalar=False):
        r"""
        Return a normalized :class:`MSymbol` (a canonical representative of an element
        of `\mathbb{P}^1(R/N)`) equivalent to ``self``.

        INPUT:

        - ``with_scalar`` -- boolean (default: ``False``)

        OUTPUT:

        - (only if ``with_scalar=True``) a transforming scalar `u`, such that
          `(u*c', u*d')` is congruent to `(c: d)` (mod `N`), where `(c: d)`
          are the coefficients of ``self`` and `N` is the level.

        -  a normalized :class:`MSymbol` `(c': d')` equivalent to ``self``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3, a - 1)
            sage: alpha1 = MSymbol(N, 3, a); alpha1
            M-symbol (3: a) of level Fractional ideal (3, 1/2*a - 1/2)
            sage: alpha1.normalize()
            M-symbol (0: 1) of level Fractional ideal (3, 1/2*a - 1/2)
            sage: alpha2 = MSymbol(N, 4, a + 1)
            sage: alpha2.normalize()
            M-symbol (1: -a) of level Fractional ideal (3, 1/2*a - 1/2)

        We get the scaling factor by setting ``with_scalar=True``:

        ::

            sage: alpha1.normalize(with_scalar=True)
            (a, M-symbol (0: 1) of level Fractional ideal (3, 1/2*a - 1/2))
            sage: r, beta1 = alpha1.normalize(with_scalar=True)
            sage: r*beta1.c - alpha1.c in N
            True
            sage: r*beta1.d - alpha1.d in N
            True
            sage: r, beta2 = alpha2.normalize(with_scalar=True)
            sage: r*beta2.c - alpha2.c in N
            True
            sage: r*beta2.d - alpha2.d in N
            True
        """
        N = self.__N
        k = N.number_field()
        R = k.ring_of_integers()

        if self.__c in N:
            if with_scalar:
                return N.reduce(self.d), MSymbol(N, 0, 1)
            else:
                return MSymbol(N, 0, 1)
        if self.d in N:
            if with_scalar:
                return N.reduce(self.c), MSymbol(N, 1, 0)
            else:
                return MSymbol(N, 1, 0)
        if N.is_coprime(self.c):
            cinv = R(self.c).inverse_mod(N)
            if with_scalar:
                return N.reduce(self.c), MSymbol(N, 1, N.reduce(self.d*cinv))
            else:
                return MSymbol(N, 1, N.reduce(self.d*cinv))

        if N in _level_cache:
            Lfacs, Lxs = _level_cache[N]
        else:
            Lfacs = [p**e for p, e in N.factor()]
            Lxs = [(N/p).element_1_mod(p) for p in Lfacs]
            # Lfacs, Lxs only depend of the ideal: same lists every time we
            # call normalize for a given level, so we store the lists.
            _level_cache[N] = (Lfacs, Lxs)
        u = 0  # normalizer factor
        p_i = 0
        for p in Lfacs:
            if p.is_coprime(self.c):
                inv = self.c.inverse_mod(p)
            else:
                inv = self.d.inverse_mod(p)
            u = u + inv*Lxs[p_i]
            p_i = p_i + 1
        c, d = (N.reduce(u*self.c), N.reduce(u*self.d))
        if (c - 1) in N:
            c = R(1)
        if with_scalar:
            return u.inverse_mod(N), MSymbol(N, c, d)
        else:
            return MSymbol(N, c, d)


# ************************************************************************
#        P1NFList class                                                  *
# ************************************************************************


@richcmp_method
class P1NFList(SageObject):
    r"""
    The class for `\mathbb{P}^1(R/N)`, the projective line modulo `N`, where
    `R` is the ring of integers of a number field `K` and `N` is an integral ideal.

    INPUT:

    - ``N`` -- integral ideal (the modulus or level)

    OUTPUT:

    A :class:`P1NFList` object representing `\mathbb{P}^1(R/N)`.

    EXAMPLES::

        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^3 + 11)
        sage: N = k.ideal(5, a + 1)
        sage: P = P1NFList(N); P
        The projective line over the ring of integers modulo the Fractional ideal (5, a + 1)

    Saving and loading works.

    ::

        sage: loads(dumps(P)) == P
        True
    """
    def __init__(self, N):
        r"""
        The constructor for the class P1NFList. See ``P1NFList`` for full
        documentation.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 5)
            sage: N = k.ideal(3, a - 1)
            sage: P = P1NFList(N); P
            The projective line over the ring of integers modulo the Fractional ideal (3, a + 2)
        """
        self.__N = N
        self.__list = p1NFlist(N)
        self.__list.sort()

    def __richcmp__(self, other, op):
        r"""
        Comparison function for objects of the class P1NFList.

        The order is the same as for the underlying modulus.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N1 = k.ideal(3, a + 1)
            sage: P1 = P1NFList(N1)
            sage: N2 = k.ideal(a + 2)
            sage: P2 = P1NFList(N2)
            sage: P1 < P2
            True
            sage: P1 > P2
            False
            sage: P1 == P1NFList(N1)
            True
        """
        if not isinstance(other, P1NFList):
            raise ValueError("You can only compare with another P1NFList")
        return richcmp(self.__N, other.__N, op)

    def __getitem__(self, n):
        r"""
        Standard indexing function for the class P1NFList.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(a)
            sage: P = P1NFList(N)
            sage: list(P) == P._P1NFList__list
            True
            sage: j = randint(0,len(P)-1)
            sage: P[j] == P._P1NFList__list[j]
            True
        """
        return self.__list[n]

    def __len__(self):
        r"""
        Return the length of this P1NFList.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a^2 - a + 1)
            sage: P = P1NFList(N)
            sage: len(P)
            26
        """
        return len(self.__list)

    def __repr__(self):
        r"""
        Return the string representation of this P1NFList.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a + 1)
            sage: P = P1NFList(N); P
            The projective line over the ring of integers modulo the Fractional ideal (5, a + 1)
        """
        return "The projective line over the ring of integers modulo the %s" % self.__N

    def list(self):
        r"""
        Return the underlying list of this :class:`P1NFList` object.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a+1)
            sage: P = P1NFList(N)
            sage: type(P)
            <class 'sage.modular.modsym.p1list_nf.P1NFList'>
            sage: type(P.list())
            <... 'list'>
        """
        return self.__list

    def normalize(self, c, d=None, with_scalar=False):
        r"""
        Return a normalised element of `\mathbb{P}^1(R/N)`.

        INPUT:

        - ``c`` -- integral element of the underlying number field, or an
          MSymbol

        - ``d`` -- (optional) when present, it must be an integral element of
          the number field such that `(c, d)` defines an M-symbol of level `N`

        - ``with_scalar`` -- boolean (default: ``False``)

        OUTPUT:

        - (only if ``with_scalar=True``) a transforming scalar `u`, such that
          `(u*c', u*d')` is congruent to `(c: d)` (mod `N`).

        - a normalized :class:`MSymbol` `(c': d')` equivalent to `(c: d)`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 31)
            sage: N = k.ideal(5, a + 3)
            sage: P = P1NFList(N)
            sage: P.normalize(3, a)
            M-symbol (1: 2*a) of level Fractional ideal (5, 1/2*a + 3/2)

        We can use an :class:`MSymbol` as input:

        ::

            sage: alpha = MSymbol(N, 3, a)
            sage: P.normalize(alpha)
            M-symbol (1: 2*a) of level Fractional ideal (5, 1/2*a + 3/2)

        If we are interested in the normalizing scalar:

        ::

            sage: P.normalize(alpha, with_scalar=True)
            (-a, M-symbol (1: 2*a) of level Fractional ideal (5, 1/2*a + 3/2))
            sage: r, beta = P.normalize(alpha, with_scalar=True)
            sage: (r*beta.c - alpha.c in N) and (r*beta.d - alpha.d in N)
            True
        """
        if d is None:
            try:
                c = MSymbol(self.__N, c)  # check that c is an MSymbol
            except ValueError:  # catch special case of wrong level
                raise ValueError("The MSymbol is of a different level")
            return c.normalize(with_scalar)
        return MSymbol(self.N(), c, d).normalize(with_scalar)

    def N(self):
        r"""
        Return the level or modulus of this :class:`P1NFList`.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 31)
            sage: N = k.ideal(5, a + 3)
            sage: P = P1NFList(N)
            sage: P.N()
            Fractional ideal (5, 1/2*a + 3/2)
        """
        return self.__N

    def index(self, c, d=None, with_scalar=False):
        r"""
        Return the index of the class of the pair `(c, d)` in the fixed list
        of representatives of `\mathbb{P}^1(R/N)`.

        INPUT:

        - ``c`` -- integral element of the corresponding number field, or an
          :class:`MSymbol`

        - ``d`` -- (optional) when present, it must be an integral element of
          the number field such that `(c, d)` defines an M-symbol of level `N`

        - ``with_scalar`` -- boolean (default: ``False``)

        OUTPUT:

        - ``u`` -- the normalizing scalar (only if ``with_scalar=True``)

        - ``i`` -- the index of `(c, d)` in the list

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 31)
            sage: N = k.ideal(5, a + 3)
            sage: P = P1NFList(N)
            sage: P.index(3,a)
            5
            sage: P[5]==MSymbol(N, 3, a).normalize()
            True

        We can give an :class:`MSymbol` as input:

        ::

            sage: alpha = MSymbol(N, 3, a)
            sage: P.index(alpha)
            5

        We cannot look for the class of an :class:`MSymbol` of a different level:

        ::

            sage: M = k.ideal(a + 1)
            sage: beta = MSymbol(M, 0, 1)
            sage: P.index(beta)
            Traceback (most recent call last):
            ...
            ValueError: The MSymbol is of a different level

        If we are interested in the transforming scalar:

        ::

            sage: alpha = MSymbol(N, 3, a)
            sage: P.index(alpha, with_scalar=True)
            (-a, 5)
            sage: u, i = P.index(alpha, with_scalar=True)
            sage: (u*P[i].c - alpha.c in N) and (u*P[i].d - alpha.d in N)
            True
        """
        if d is None:
            try:
                c = MSymbol(self.__N, c)  # check that c is an MSymbol
            except ValueError:  # catch special case of wrong level
                raise ValueError("The MSymbol is of a different level")
            if with_scalar:
                u, norm_c = c.normalize(with_scalar=True)
            else:
                norm_c = c.normalize()
        else:
            if with_scalar:
                u, norm_c = MSymbol(self.__N, c, d).normalize(with_scalar=True)
            else:
                norm_c = MSymbol(self.__N, c, d).normalize()
        t, i = search(self.__list, norm_c)
        if t:
            if with_scalar:
                return u, i
            else:
                return i
        return False

    def index_of_normalized_pair(self, c, d=None):
        r"""
        Return the index of the class `(c, d)` in the fixed list of
        representatives of `\mathbb(P)^1(R/N)`.

        INPUT:

        - ``c`` -- integral element of the corresponding number field, or a
          normalized :class:`MSymbol`

        - ``d`` -- (optional) when present, it must be an integral element of
          the number field such that `(c, d)` defines a normalized M-symbol of
          level `N`

        OUTPUT: ``i`` -- the index of `(c, d)` in the list

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 31)
            sage: N = k.ideal(5, a + 3)
            sage: P = P1NFList(N)
            sage: P.index_of_normalized_pair(1, 0)
            3
            sage: j = randint(0,len(P)-1)
            sage: P.index_of_normalized_pair(P[j])==j
            True
        """
        if d is None:
            try:
                c = MSymbol(self.__N, c)  # check that c is an MSymbol
            except ValueError:  # catch special case of wrong level
                raise ValueError("The MSymbol is of a different level")
            t, i = search(self.__list, c)
        else:
            t, i = search(self.__list, MSymbol(self.__N, c, d))
        if t:
            return i
        return False

    def lift_to_sl2_Ok(self, i):
        r"""
        Lift the `i`-th element of this :class:`P1NFList` to an element of `SL(2, R)`,
        where `R` is the ring of integers of the corresponding number field.

        INPUT:

        - ``i`` -- integer (index of the element to lift)

        OUTPUT:

        If the `i`-th element is `(c : d)`, the function returns a list of
        integral elements `[a, b, c', d']` that defines a `2\times 2` matrix with
        determinant 1 and such that `c=c'` (mod `N`) and `d=d'` (mod `N`).

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3)
            sage: P = P1NFList(N)
            sage: len(P)
            16
            sage: P[5]
            M-symbol (1/2*a + 1/2: -a) of level Fractional ideal (3)
            sage: P.lift_to_sl2_Ok(5)
            [-a, 2*a - 2, 1/2*a + 1/2, -a]

        ::

            sage: Ok = k.ring_of_integers()
            sage: L = [Matrix(Ok, 2, P.lift_to_sl2_Ok(i)) for i in range(len(P))]
            sage: all(det(L[i]) == 1 for i in range(len(L)))
            True
        """
        return self[i].lift_to_sl2_Ok()

    def apply_S(self, i):
        r"""
        Applies the matrix `S` = [0, -1, 1, 0] to the `i`-th M-Symbol of the list.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        integer -- the index of the M-Symbol obtained by the right action of
        the matrix `S` = [0, -1, 1, 0] on the `i`-th M-Symbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a + 1)
            sage: P = P1NFList(N)
            sage: j = P.apply_S(P.index_of_normalized_pair(1, 0))
            sage: P[j]
            M-symbol (0: 1) of level Fractional ideal (5, a + 1)

        We test that S has order 2:

        ::

            sage: j = randint(0,len(P)-1)
            sage: P.apply_S(P.apply_S(j))==j
            True
        """
        c, d = self.__list[i].tuple()
        t, j = search(self.__list, self.normalize(d, -c))
        return j

    def apply_TS(self, i):
        r"""
        Applies the matrix `TS` = [1, -1, 0, 1] to the `i`-th M-Symbol of the list.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        integer -- the index of the M-Symbol obtained by the right action of
        the matrix `TS` = [1, -1, 0, 1] on the `i`-th M-Symbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a + 1)
            sage: P = P1NFList(N)
            sage: P.apply_TS(3)
            2

        We test that `TS` has order 3:

        ::

            sage: j = randint(0,len(P)-1)
            sage: P.apply_TS(P.apply_TS(P.apply_TS(j)))==j
            True
        """
        c, d = self.__list[i].tuple()
        t, j = search(self.__list, self.normalize(c + d, -c))
        return j

    def apply_T_alpha(self, i, alpha=1):
        r"""
        Applies the matrix `T_{alpha}` = [1, `alpha`, 0, 1] to the `i`-th
        M-Symbol of the list.

        INPUT:

        - ``i`` -- integer

        - ``alpha`` -- (default: 1) element of the corresponding ring of integers

        OUTPUT:

        integer -- the index of the M-Symbol obtained by the right action of
        the matrix `T_{alpha}` = [1, `alpha`, 0, 1] on the `i`-th M-Symbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a + 1)
            sage: P = P1NFList(N)
            sage: P.apply_T_alpha(4, a^ 2 - 2)
            3

        We test that `T_a*T_b = T_{(a+b)}`:

        ::

            sage: P.apply_T_alpha(3, a^2 - 2)==P.apply_T_alpha(P.apply_T_alpha(3,a^2),-2)
            True
        """
        c, d = self.__list[i].tuple()
        t, j = search(self.__list, self.normalize(c, alpha*c + d))
        return j

    def apply_J_epsilon(self, i, e1, e2=1):
        r"""
        Apply the matrix `J_{\epsilon}` = [e1, 0, 0, e2] to the `i`-th
        M-Symbol of the list.

        e1, e2 are units of the underlying number field.

        INPUT:

        - ``i`` -- integer

        - ``e1`` -- unit

        - ``e2`` -- unit (default: 1)

        OUTPUT:

        integer -- the index of the M-Symbol obtained by the right action of
        the matrix `J_{\epsilon}` = [e1, 0, 0, e2] on the `i`-th M-Symbol.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: k.<a> = NumberField(x^3 + 11)
            sage: N = k.ideal(5, a + 1)
            sage: P = P1NFList(N)
            sage: u = k.unit_group().gens_values(); u
            [-1, -2*a^2 - 4*a + 1]
            sage: P.apply_J_epsilon(4, -1)
            2
            sage: P.apply_J_epsilon(4, u[0], u[1])
            5

        ::

            sage: k.<a> = NumberField(x^4 + 13*x - 7)
            sage: N = k.ideal(a + 1)
            sage: P = P1NFList(N)
            sage: u = k.unit_group().gens_values(); u
            [-1, -a^3 - a^2 - a - 12, -a^3 - 3*a^2 + 1]
            sage: P.apply_J_epsilon(3, u[2]^2)==P.apply_J_epsilon(P.apply_J_epsilon(3, u[2]),u[2])
            True
        """
        c, d = self.__list[i].tuple()
        t, j = search(self.__list, self.normalize(c*e1, d*e2))
        return j


# *************************************************************************
#  Global functions:
#    - p1NFList --compute list of M-symbols
#    - lift_to_sl2_Ok
#    - make_coprime -- need it for ``lift_to_sl2_Ok``
#    - psi -- useful to check cardinality of the M-symbols list
# *************************************************************************

def p1NFlist(N):
    r"""
    Return a list of the normalized elements of `\mathbb{P}^1(R/N)`, where
    `N` is an integral ideal.

    INPUT:

    - ``N`` -- integral ideal (the level or modulus)

    EXAMPLES::

        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^2 + 23)
        sage: N = k.ideal(3)
        sage: from sage.modular.modsym.p1list_nf import p1NFlist, psi
        sage: len(p1NFlist(N))==psi(N)
        True
    """
    k = N.number_field()

    L = [MSymbol(N, k(0), k(1), check=False)]
    # N.residues() = iterator through the residues mod N
    L = L + [MSymbol(N, k(1), r, check=False) for r in N.residues()]

    from sage.arith.misc import divisors
    for D in divisors(N):
        if not D.is_trivial() and D != N:
            # we find Dp ideal coprime to N, in inverse class to D
            if D.is_principal():
                Dp = k.ideal(1)
                c = D.gens_reduced()[0]
            else:
                it = k.primes_of_degree_one_iter()
                Dp = next(it)
                while not Dp.is_coprime(N) or not (Dp*D).is_principal():
                    Dp = next(it)
                c = (D*Dp).gens_reduced()[0]
            # now we find all the (c,d)'s which have associated divisor D
            I = D + N/D
            for d in (N/D).residues():
                if I.is_coprime(d):
                    M = D.prime_to_idealM_part(N/D)
                    u = (Dp*M).element_1_mod(N/D)
                    d1 = u*d + (1-u)
                    L.append(MSymbol(N, c, d1, check=False).normalize())
    return L


def lift_to_sl2_Ok(N, c, d):
    r"""
    Lift a pair (c, d) to an element of `SL(2, O_k)`, where `O_k` is the ring
    of integers of the corresponding number field.

    INPUT:

    - ``N`` -- number field ideal

    - ``c`` -- integral element of the number field

    - ``d`` -- integral element of the number field

    OUTPUT:

    A list `[a, b, c', d']` of integral elements that are the entries of
    a `2\times 2` matrix with determinant 1. The lower two entries are congruent to
    `c`, `d` modulo the ideal `N`.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list_nf import lift_to_sl2_Ok
        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^2 + 23)
        sage: Ok = k.ring_of_integers()
        sage: N = k.ideal(3)
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 1, a))
        sage: det(M)
        1
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 0, a))
        sage: det(M)
        1
        sage: (M[1][0] in N) and (M[1][1] - a in N)
        True
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 0, 0))
        Traceback (most recent call last):
        ...
        ValueError: Cannot lift (0, 0) to an element of Sl2(Ok).

    ::

        sage: k.<a> = NumberField(x^3 + 11)
        sage: Ok = k.ring_of_integers()
        sage: N = k.ideal(3, a - 1)
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 2*a, 0))
        sage: det(M)
        1
        sage: (M[1][0] - 2*a in N) and (M[1][1] in N)
        True
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 4*a^2, a + 1))
        sage: det(M)
        1
        sage: (M[1][0] - 4*a^2 in N) and (M[1][1] - (a+1) in N)
        True

    ::

        sage: k.<a> = NumberField(x^4 - x^3 -21*x^2 + 17*x + 133)
        sage: Ok = k.ring_of_integers()
        sage: N = k.ideal(7, a)
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 0, a^2 - 1))
        sage: det(M)
        1
        sage: (M[1][0] in N) and (M[1][1] - (a^2-1) in N)
        True
        sage: M = Matrix(Ok, 2, lift_to_sl2_Ok(N, 0, 7))
        Traceback (most recent call last):
        ...
        ValueError: <0> + <7> and the Fractional ideal (7, a) are not coprime.
    """
    k = N.number_field()
    # check the input
    if c.is_zero() and d.is_zero():
        raise ValueError("Cannot lift (%s, %s) to an element of Sl2(Ok)." % (c, d))
    if not N.is_coprime(k.ideal(c, d)):
        raise ValueError("<%s> + <%s> and the %s are not coprime." % (c, d, N))
    # a few special cases
    if c - 1 in N:
        return [k(0), k(-1), 1, d]
    if d - 1 in N:
        return [k(1), k(0), c, 1]
    if c.is_zero():  # and d!=1, so won't happen for normalized M-symbols (c: d)
        it = k.primes_of_degree_one_iter()
        q = k.ideal(1)
        while not (q.is_coprime(d) and (q*N).is_principal()):
            q = next(it)
        m = (q*N).gens_reduced()[0]
        B = k.ideal(m).element_1_mod(k.ideal(d))
        return [(1-B)/d, -B/m, m, d]
    if d.is_zero():  # and c!=1, so won't happen for normalized M-symbols (c: d)
        it = k.primes_of_degree_one_iter()
        q = k.ideal(1)
        while not (q.is_coprime(c) and (q*N).is_principal()):
            q = next(it)
        m = (q*N).gens_reduced()[0]
        B = k.ideal(c).element_1_mod(k.ideal(m))
        return [(1-B)/m, -B/c, c, m]

    c, d = make_coprime(N, c, d)

    B = k.ideal(c).element_1_mod(k.ideal(d))
    b = -B/c
    a = (1-B)/d
    return [a, b, c, d]


def make_coprime(N, c, d):
    r"""
    Return (c, d') so d' is congruent to d modulo N, and such that c and d' are
    coprime (`\langle c\rangle + \langle d'\rangle = R`).

    INPUT:

    - ``N`` -- number field ideal

    - ``c`` -- integral element of the number field

    - ``d`` -- integral element of the number field

    OUTPUT:

    A pair `(c, d')` where `c`, `d'` are integral elements of the corresponding
    number field, with `d'` congruent to `d` mod `N`, and such that `\langle c\rangle + \langle d'\rangle = R`
    (`R` being the corresponding ring of integers).

    EXAMPLES::

        sage: from sage.modular.modsym.p1list_nf import make_coprime
        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^2 + 23)
        sage: N = k.ideal(3, a - 1)
        sage: c = 2*a; d = a + 1
        sage: N.is_coprime(k.ideal(c, d))
        True
        sage: k.ideal(c).is_coprime(d)
        False
        sage: c, dp = make_coprime(N, c, d)
        sage: k.ideal(c).is_coprime(dp)
        True
    """
    k = N.number_field()
    if k.ideal(c).is_coprime(d):
        return c, d
    else:
        q = k.ideal(c).prime_to_idealM_part(d)
        it = k.primes_of_degree_one_iter()
        r = k.ideal(1)
        qN = q*N
        while not (r.is_coprime(c) and (r*qN).is_principal()):
            r = next(it)
        m = (r*qN).gens_reduced()[0]
        d1 = d + m
        return c, d1


def psi(N):
    r"""
    The index `[\Gamma : \Gamma_0(N)]`, where `\Gamma = GL(2, R)` for `R` the
    corresponding ring of integers, and `\Gamma_0(N)` standard congruence
    subgroup.

    EXAMPLES::

        sage: from sage.modular.modsym.p1list_nf import psi
        sage: x = polygen(QQ, 'x')
        sage: k.<a> = NumberField(x^2 + 23)
        sage: N = k.ideal(3, a - 1)
        sage: psi(N)
        4

    ::

        sage: k.<a> = NumberField(x^2 + 23)
        sage: N = k.ideal(5)
        sage: psi(N)
        26
    """
    if not N.is_integral():
        raise ValueError("psi only defined for integral ideals")

    from sage.misc.misc_c import prod
    return prod([(np + 1) * np**(e - 1)
                 for np, e in [(p.absolute_norm(), e)
                               for p, e in N.factor()]])
