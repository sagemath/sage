# sage.doctest: needs sage.rings.finite_rings
r"""
Semimonomial transformation group

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

where
`\psi^{\pi, \alpha} = (\alpha(\psi_{\pi(1)-1}), \ldots, \alpha(\psi_{\pi(n)-1}))`
and the multiplication of vectors is defined elementwisely. (The indexing
of vectors is `0`-based here, so `\psi = (\psi_0, \psi_1, \ldots, \psi_{n-1})`.)

.. TODO::

    Up to now, this group is only implemented for finite fields because of
    the limited support of automorphisms for arbitrary rings.

AUTHORS:

- Thomas Feulner (2012-11-15): initial version

EXAMPLES::

    sage: S = SemimonomialTransformationGroup(GF(4, 'a'), 4)
    sage: G = S.gens()
    sage: G[0]*G[1]
    ((a, 1, 1, 1); (1,2,3,4), Ring endomorphism of Finite Field in a of size 2^2
      Defn: a |--> a)

TESTS::

    sage: TestSuite(S).run()
    sage: TestSuite(S.an_element()).run()
"""
from __future__ import annotations
from sage.rings.integer import Integer

from sage.groups.group import FiniteGroup
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.action import Action
from sage.combinat.permutation import Permutation
from sage.groups.semimonomial_transformations.semimonomial_transformation import SemimonomialTransformation


class SemimonomialTransformationGroup(FiniteGroup, UniqueRepresentation):
    r"""
    A semimonomial transformation group over a ring.

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

    where
    `\psi^{\pi, \alpha} = (\alpha(\psi_{\pi(1)-1}), \ldots, \alpha(\psi_{\pi(n)-1}))`
    and the multiplication of vectors is defined elementwisely. (The indexing
    of vectors is `0`-based here, so `\psi = (\psi_0, \psi_1, \ldots, \psi_{n-1})`.)

    .. TODO::

        Up to now, this group is only implemented for finite fields because of
        the limited support of automorphisms for arbitrary rings.

    EXAMPLES::

        sage: F.<a> = GF(9)
        sage: S = SemimonomialTransformationGroup(F, 4)
        sage: g = S(v = [2, a, 1, 2])
        sage: h = S(perm = Permutation('(1,2,3,4)'), autom=F.hom([a**3]))
        sage: g*h
        ((2, a, 1, 2); (1,2,3,4),
         Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1)
        sage: h*g
        ((2*a + 1, 1, 2, 2); (1,2,3,4),
         Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1)
        sage: S(g)
        ((2, a, 1, 2); (),
         Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
        sage: S(1)
        ((1, 1, 1, 1); (),
         Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
    """
    Element = SemimonomialTransformation

    def __init__(self, R, len):
        r"""
        Initialization.

        INPUT:

        - ``R`` -- a ring

        - ``len`` -- the degree of the monomial group

        OUTPUT: the complete semimonomial group

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: S = SemimonomialTransformationGroup(F, 4)
        """
        if not R.is_field():
            raise NotImplementedError('the ring must be a field')
        self._R = R
        self._len = len

        from sage.categories.finite_groups import FiniteGroups
        super().__init__(category=FiniteGroups())

    def _element_constructor_(self, arg1, v=None, perm=None, autom=None, check=True):
        r"""
        Coerce ``arg1`` into this permutation group, if ``arg1`` is 0,
        then we will try to coerce ``(v, perm, autom)``.

        INPUT:

        - ``arg1`` -- (optional) either the integers 0, 1 or an element of ``self``

        - ``v`` -- (optional) a vector of length ``self.degree()``

        - ``perm`` -- (optional) a permutation of degree ``self.degree()``

        - ``autom`` -- (optional) an automorphism of the ring

        EXAMPLES::

            sage: F.<a> = GF(9)
            sage: S = SemimonomialTransformationGroup(F, 4)
            sage: S(1)
            ((1, 1, 1, 1); (), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
            sage: g = S(v=[1,1,1,a])
            sage: S(g)
            ((1, 1, 1, a); (), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
            sage: S(perm=Permutation('(1,2)(3,4)'))
            ((1, 1, 1, 1); (1,2)(3,4), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> a)
            sage: S(autom=F.hom([a**3]))
            ((1, 1, 1, 1); (), Ring endomorphism of Finite Field in a of size 3^2 Defn: a |--> 2*a + 1)
        """
        from sage.categories.homset import End
        R = self.base_ring()
        if arg1 == 0:
            if v is None:
                v = [R.one()] * self.degree()
            if perm is None:
                perm = Permutation(range(1, self.degree() + 1))
            if autom is None:
                autom = R.hom(R.gens())

            if check:
                try:
                    v = [R(x) for x in v]
                except TypeError:
                    raise TypeError('the vector attribute %s ' % v +
                                    'should be iterable')
                if len(v) != self.degree():
                    raise ValueError('the length of the vector is %s,' % len(v) +
                                     ' should be %s' % self.degree())
                if not all(x.parent() is R and x.is_unit() for x in v):
                    raise ValueError('there is at least one element in the ' +
                                     'list %s not lying in %s ' % (v, R) +
                                     'or which is not invertible')
                try:
                    perm = Permutation(perm)
                except TypeError:
                    raise TypeError('the permutation attribute %s ' % perm +
                                    'could not be converted to a permutation')
                if len(perm) != self.degree():
                    txt = 'the permutation length is {}, should be {}'
                    raise ValueError(txt.format(len(perm), self.degree()))

                try:
                    if autom.parent() != End(R):
                        autom = End(R)(autom)
                except TypeError:
                    raise TypeError('%s of type %s' % (autom, type(autom)) +
                                    ' is not coerceable to an automorphism')
            return self.Element(self, v, perm, autom)
        else:
            try:
                if arg1.parent() is self:
                    return arg1
            except AttributeError:
                pass
            try:
                from sage.rings.integer import Integer
                if Integer(arg1) == 1:
                    return self()
            except TypeError:
                pass
            raise TypeError('the first argument must be an integer' +
                            ' or an element of this group')

    def base_ring(self):
        r"""
        Return the underlying ring of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: SemimonomialTransformationGroup(F, 3).base_ring() is F
            True
        """
        return self._R

    def degree(self) -> Integer:
        r"""
        Return the degree of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: SemimonomialTransformationGroup(F, 3).degree()
            3
        """
        return self._len

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: SemimonomialTransformationGroup(F, 3).an_element()  # indirect doctest
            ((a, 1, 1); (1,3,2), Ring endomorphism of Finite Field in a of size 2^2 Defn: a |--> a + 1)
        """
        R = self.base_ring()
        v = [R.primitive_element()] + [R.one()] * (self.degree() - 1)
        p = Permutation([self.degree()] + list(range(1, self.degree())))

        if not R.is_prime_field():
            f = R.hom([R.gen()**R.characteristic()])
        else:
            f = R.Hom(R).identity()
        return self(0, v, p, f)

    def __contains__(self, item) -> bool:
        r"""
        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: S = SemimonomialTransformationGroup(F, 3)
            sage: 1 in S  # indirect doctest
            True
            sage: a in S  # indirect doctest
            False
        """
        try:
            self(item, check=True)
        except TypeError:
            return False
        return True

    def gens(self) -> tuple:
        r"""
        Return a tuple of generators of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: SemimonomialTransformationGroup(F, 3).gens()
            (((a, 1, 1); (),
              Ring endomorphism of Finite Field in a of size 2^2 Defn: a |--> a),
             ((1, 1, 1); (1,2,3),
              Ring endomorphism of Finite Field in a of size 2^2 Defn: a |--> a),
             ((1, 1, 1); (1,2),
              Ring endomorphism of Finite Field in a of size 2^2 Defn: a |--> a),
             ((1, 1, 1); (),
              Ring endomorphism of Finite Field in a of size 2^2 Defn: a |--> a + 1))
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        R = self.base_ring()
        l = [self(v=([R.primitive_element()] + [R.one()] * (self.degree() - 1)))]
        l.extend(self(perm=Permutation(g))
                 for g in SymmetricGroup(self.degree()).gens())
        if R.is_field() and not R.is_prime_field():
            l.append(self(autom=R.hom([R.primitive_element()**R.characteristic()])))
        return tuple(l)

    def order(self) -> Integer:
        r"""
        Return the number of elements of ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: SemimonomialTransformationGroup(F, 5).order() == (4-1)**5 * factorial(5) * 2
            True
        """
        from sage.arith.misc import factorial
        from sage.categories.homset import End
        n = self.degree()
        R = self.base_ring()
        if R.is_field():
            multgroup_size = len(R) - 1
            autgroup_size = R.degree()
        else:
            multgroup_size = R.unit_group_order()
            autgroup_size = len([x for x in End(R) if x.is_injective()])
        return multgroup_size**n * factorial(n) * autgroup_size

    def _get_action_(self, X, op, self_on_left):
        r"""
        If ``self`` is the semimonomial group of degree `n` over `R`, then
        there is the natural action on `R^n` and on matrices `R^{m \times n}`
        for arbitrary integers `m` from the left.

        See also:
        :class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialActionVec` and
        :class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialActionMat`

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: s = SemimonomialTransformationGroup(F, 3).an_element()
            sage: v = (F**3).0
            sage: s*v   # indirect doctest
            (0, 1, 0)
            sage: M = MatrixSpace(F, 3).one()
            sage: s*M   # indirect doctest
            [    0     1     0]
            [    0     0     1]
            [a + 1     0     0]
        """
        if self_on_left:
            try:
                A = SemimonomialActionVec(self, X)
                return A
            except ValueError:
                pass

            try:
                A = SemimonomialActionMat(self, X)
                return A
            except ValueError:
                pass

        return None

    def _repr_(self) -> str:
        r"""
        Return a string describing ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: SemimonomialTransformationGroup(F, 3)  # indirect doctest
            Semimonomial transformation group over Finite Field in a of size 2^2 of degree 3
        """
        return ('Semimonomial transformation group over %s' % self.base_ring() +
                ' of degree %s' % self.degree())

    def _latex_(self) -> str:
        r"""
        Method for describing ``self`` in LaTeX.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: latex(SemimonomialTransformationGroup(F, 3))  # indirect doctest
            \left(\Bold{F}_{2^{2}}^3\wr\langle (1,2,3), (1,2) \rangle \right) \rtimes \operatorname{Aut}(\Bold{F}_{2^{2}})
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        ring_latex = self.base_ring()._latex_()
        return ('\\left(' + ring_latex + '^' + str(self.degree()) + '\\wr' +
                SymmetricGroup(self.degree())._latex_() +
                ' \\right) \\rtimes \\operatorname{Aut}(' + ring_latex + ')')


class SemimonomialActionVec(Action):
    r"""
    The natural left action of the semimonomial group on vectors.

    The action is defined by:
    `(\phi, \pi, \alpha)*(v_0, \ldots, v_{n-1}) :=
    (\alpha(v_{\pi(1)-1}) \cdot \phi_0^{-1}, \ldots, \alpha(v_{\pi(n)-1}) \cdot \phi_{n-1}^{-1})`.
    (The indexing of vectors is `0`-based here, so
    `\psi = (\psi_0, \psi_1, \ldots, \psi_{n-1})`.)
    """
    def __init__(self, G, V, check=True):
        r"""
        Initialization.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: s = SemimonomialTransformationGroup(F, 3).an_element()
            sage: v = (F**3).1
            sage: s*v   # indirect doctest
            (0, 0, 1)
        """
        if check:
            from sage.modules.free_module import FreeModule_generic
            if not isinstance(G, SemimonomialTransformationGroup):
                raise ValueError('%s is not a semimonomial group' % G)
            if not isinstance(V, FreeModule_generic):
                raise ValueError('%s is not a free module' % V)
            if V.ambient_module() != V:
                raise ValueError('%s is not equal to its ambient module' % V)
            if V.dimension() != G.degree():
                raise ValueError('%s has a dimension different to the degree of %s' % (V, G))
            if V.base_ring() != G.base_ring():
                raise ValueError('%s and %s have different base rings' % (V, G))

        Action.__init__(self, G, V.dense_module())

    def _act_(self, a, b):
        r"""
        Apply the semimonomial group element `a` to the vector `b`.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: s = SemimonomialTransformationGroup(F, 3).an_element()
            sage: v = (F**3).1
            sage: s*v   # indirect doctest
            (0, 0, 1)
        """
        b = b.apply_map(a.get_autom())
        b = self.codomain()(a.get_perm().action(b))
        return b.pairwise_product(self.codomain()(a.get_v_inverse()))


class SemimonomialActionMat(Action):
    r"""
    The left action of
    :class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialTransformationGroup`
    on matrices over the same ring whose number of columns is equal to the degree.
    See :class:`~sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialActionVec`
    for the definition of the action on the row vectors of such a matrix.
    """
    def __init__(self, G, M, check=True):
        r"""
        Initialization.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: s = SemimonomialTransformationGroup(F, 3).an_element()
            sage: M = MatrixSpace(F, 3).one()
            sage: s*M  # indirect doctest
            [    0     1     0]
            [    0     0     1]
            [a + 1     0     0]
        """
        if check:
            from sage.matrix.matrix_space import MatrixSpace
            if not isinstance(G, SemimonomialTransformationGroup):
                raise ValueError('%s is not a semimonomial group' % G)
            if not isinstance(M, MatrixSpace):
                raise ValueError('%s is not a matrix space' % M)
            if M.ncols() != G.degree():
                raise ValueError('the number of columns of %s' % M +
                                 ' and the degree of %s are different' % G)
            if M.base_ring() != G.base_ring():
                raise ValueError('%s and %s have different base rings' % (M, G))
        Action.__init__(self, G, M)

    def _act_(self, a, b):
        r"""
        Apply the semimonomial group element `a` to the matrix `b`.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: s = SemimonomialTransformationGroup(F, 3).an_element()
            sage: M = MatrixSpace(F, 3).one()
            sage: s*M  # indirect doctest
            [    0     1     0]
            [    0     0     1]
            [a + 1     0     0]
        """
        return self.codomain()([a * x for x in b.rows()])
