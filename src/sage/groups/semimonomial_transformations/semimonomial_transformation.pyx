# sage_setup: distribution = sagemath-groups
# sage.doctest: needs sage.rings.finite_rings
r"""
Elements of a semimonomial transformation group

The semimonomial transformation group of degree `n` over a ring `R` is
the semidirect product of the monomial transformation group of degree `n`
(also known as the complete monomial group over the group of units
`R^{\times}` of `R`) and the group of ring automorphisms.

The multiplication of two elements `(\phi, \pi, \alpha)(\psi, \sigma, \beta)`
with

- `\phi, \psi \in  {R^{\times}}^n`

- `\pi, \sigma \in S_n` (with the multiplication `\pi\sigma`
  done from left to right (like in GAP) --
  that is, `(\pi\sigma)(i) = \sigma(\pi(i))` for all `i`.)

- `\alpha, \beta \in Aut(R)`

is defined by

.. MATH::

    (\phi, \pi, \alpha)(\psi, \sigma, \beta) =
    (\phi \cdot \psi^{\pi, \alpha}, \pi\sigma, \alpha \circ \beta)

with
`\psi^{\pi, \alpha} = (\alpha(\psi_{\pi(1)-1}), \ldots, \alpha(\psi_{\pi(n)-1}))`
and an elementwisely defined multiplication of vectors. (The indexing
of vectors is `0`-based here, so `\psi = (\psi_0, \psi_1, \ldots, \psi_{n-1})`.)

The parent is
:class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialTransformationGroup`.

AUTHORS:

- Thomas Feulner (2012-11-15): initial version
- Thomas Feulner (2013-12-27): :issue:`15576` dissolve dependency on
    Permutations.options.mul

EXAMPLES::

    sage: S = SemimonomialTransformationGroup(GF(4, 'a'), 4)
    sage: G = S.gens()
    sage: G[0]*G[1]
    ((a, 1, 1, 1); (1,2,3,4), Ring endomorphism of Finite Field in a of size 2^2
      Defn: a |--> a)

TESTS::

    sage: TestSuite(G[0]).run()
"""
from cpython.object cimport PyObject_RichCompare


def _is_id(f, R):
    """
    Test some automorphism `f` of a ring `R` if it is the identity
    automorphism.

    EXAMPLES::

        sage: from sage.groups.semimonomial_transformations.semimonomial_transformation import _is_id
        sage: F.<a> = GF(8)
        sage: f = F.hom([a**2])
        sage: _is_id(f, F)
        False
    """
    for r in R.gens():
        if r != f(r):
            return False
    return True


def _inverse(f, R):
    """
    Return the inverse to the automorphism `f` of a ring `R`.

    EXAMPLES::

        sage: from sage.groups.semimonomial_transformations.semimonomial_transformation import _inverse
        sage: F.<a> = GF(8)
        sage: f = F.hom([a**2])
        sage: _inverse(f, F)*f == F.hom([a])
        True
    """
    g = f
    while not _is_id(g*f, R):
        g *= f
    return g


cdef class SemimonomialTransformation(MultiplicativeGroupElement):
    r"""
    An element in the semimonomial group over a ring `R`. See
    :class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialTransformationGroup`
    for the details on the multiplication of two elements.

    The init method should never be called directly. Use the call via the
    parent
    :class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialTransformationGroup`.
    instead.

    EXAMPLES::

        sage: F.<a> = GF(9)
        sage: S = SemimonomialTransformationGroup(F, 4)
        sage: g = S(v = [2, a, 1, 2])
        sage: h = S(perm = Permutation('(1,2,3,4)'), autom=F.hom([a**3]))
        sage: g*h
        ((2, a, 1, 2); (1,2,3,4), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1)
        sage: h*g
        ((2*a + 1, 1, 2, 2); (1,2,3,4), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1)
        sage: S(g)
        ((2, a, 1, 2); (), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
        sage: S(1)  # the one element in the group
        ((1, 1, 1, 1); (), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
    """
    def __init__(self, parent, v, perm, alpha):
        r"""
        The init method should never be called directly. Use the call via the
        parent instead. See
        :meth:`sage.groups.semimonomial_transformations.semimonomial_transformation.SemimonomialTransformation.__call__`.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: S = SemimonomialTransformationGroup(F, 4)
            sage: g = S(v = [2, a, 1, 2])  #indirect doctest
        """
        MultiplicativeGroupElement.__init__(self, parent)
        self.v = tuple(v)
        self.perm = perm
        self.alpha = alpha

    cdef _new_c(self):
        # Create a copy of self.
        cdef SemimonomialTransformation x
        x = SemimonomialTransformation.__new__(SemimonomialTransformation)
        x._parent = self._parent
        x.v = self.v
        x.perm = self.perm
        x.alpha = self.alpha
        return x

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: s = SemimonomialTransformationGroup(F, 4).an_element()
            sage: t = copy(s)  # indirect doctest
            sage: t is s
            False
            sage: t == s
            True
        """
        return self._new_c()

    def __hash__(self):
        """
        Return hash of this element.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: hash( SemimonomialTransformationGroup(F, 4).an_element() )  #random #indirect doctest
            6279637968393375107
        """
        return hash(self.v) + hash(self.perm) + hash(self.get_autom())

    cpdef _mul_(left, _right):
        r"""
        Multiplication of elements.

        The multiplication of two elements `(\phi, \pi, \alpha)` and
        `(\psi, \sigma, \beta)` with

            - `\phi, \psi \in  {R^{\times}}^n`

            - `\pi, \sigma \in S_n`

            - `\alpha, \beta \in Aut(R)`

        is defined by:

        .. MATH::

            (\phi, \pi, \alpha)(\psi, \sigma, \beta) =
            (\phi \cdot \psi^{\pi, \alpha}, \pi\sigma, \alpha \circ \beta)

        with
        `\psi^{\pi, \alpha} = (\alpha(\psi_{\pi(1)-1}), \ldots, \alpha(\psi_{\pi(n)-1}))`
        and an elementwisely defined multiplication of vectors. (The indexing
        of vectors is `0`-based here, so `\psi = (\psi_0, \psi_1, \ldots, \psi_{n-1})`.)
        Furthermore, the multiplication `\pi\sigma` is done from left to right
        (like in GAP) -- that is, `(\pi\sigma)(i) = \sigma(\pi(i))` for all `i`.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: s = SemimonomialTransformationGroup(F, 4).an_element()
            sage: s*s  # indirect doctest
            ((a, 2*a + 1, 1, 1); (1,3)(2,4), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
        """
        cdef SemimonomialTransformation right = <SemimonomialTransformation> _right
        cdef i
        v = left.perm.action(right.v)
        alpha = left.get_autom()
        v = [left.v[i]*alpha(v[i]) for i in range(left.parent().degree())]
        return left.parent()(v=v, perm=left.perm.right_action_product(right.perm),
                             autom=alpha*right.get_autom(), check=False)

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: S = SemimonomialTransformationGroup(F, 4)
            sage: s = S.an_element()
            sage: s*s**(-1) == S(1)  # indirect doctest
            True
        """
        cdef i
        alpha = _inverse(self.get_autom(), self.get_autom().domain())
        inv_perm = self.perm.inverse()
        v = [alpha(self.v[i]**(-1)) for i in range(len(self.v))]
        return self.parent()(v=inv_perm.action(v), perm=inv_perm, autom=alpha,
                             check=False)

    def __repr__(self):
        """
        String representation of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: SemimonomialTransformationGroup(F, 4).an_element()  # indirect doctest
            ((a, 1, 1, 1); (1,4,3,2), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1)
        """
        return "(%s; %s, %s)"%(self.v, self.perm.cycle_string(),
                               self.get_autom())

    cpdef _richcmp_(left, _right, int op):
        """
        Compare group elements ``self`` and ``right``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: g = SemimonomialTransformationGroup(F, 4).gens()
            sage: g[0] > g[1]  # indirect doctest
            True
            sage: g[1] != g[2]  # indirect doctest
            True
        """
        cdef SemimonomialTransformation right = <SemimonomialTransformation> _right
        return PyObject_RichCompare([left.v, left.perm, left.get_autom()],
                                    [right.v, right.perm, right.get_autom()],
                                    op)

    def __reduce__(self):
        """
        Return a function and its arguments needed to create this
        semimonomial group element. This is used in pickling.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: SemimonomialTransformationGroup(F, 4).an_element().__reduce__()
            (Semimonomial transformation group over Finite Field in a of size 3^2 of degree 4, (0, (a, 1, 1, 1), [4, 1, 2, 3], Ring endomorphism of Finite Field in a of size 3^2
              Defn: a |--> 2*a + 1))
        """
        return (self.parent(), (0, self.v, self.perm, self.get_autom()))

    def get_v(self):
        """
        Return the component corresponding to `{R^{\times}}^n` of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: SemimonomialTransformationGroup(F, 4).an_element().get_v()
            (a, 1, 1, 1)
        """
        return self.v

    def get_v_inverse(self):
        """
        Return the (elementwise) inverse of the component corresponding to
        `{R^{\times}}^n` of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: SemimonomialTransformationGroup(F, 4).an_element().get_v_inverse()
            (a + 2, 1, 1, 1)
        """
        return tuple(x**(-1) for x in self.v)

    def get_perm(self):
        """
        Return the component corresponding to `S_n` of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: SemimonomialTransformationGroup(F, 4).an_element().get_perm()
            [4, 1, 2, 3]
        """
        return self.perm

    def get_autom(self):
        """
        Return the component corresponding to `Aut(R)` of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: SemimonomialTransformationGroup(F, 4).an_element().get_autom()
            Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1
        """
        return self.alpha

    def invert_v(self):
        """
        Elementwisely invert all entries of ``self`` which
        correspond to the component `{R^{\times}}^n`.

        The other components of ``self`` keep unchanged.

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: x = copy(SemimonomialTransformationGroup(F, 4).an_element())
            sage: x.invert_v()
            sage: x.get_v() == SemimonomialTransformationGroup(F, 4).an_element().get_v_inverse()
            True
        """
        self.v = tuple([x**(-1) for x in self.v])
