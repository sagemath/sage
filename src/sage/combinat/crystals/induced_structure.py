# sage.doctest: needs sage.combinat sage.modules
r"""
Induced Crystals

We construct a crystal structure on a set induced by a bijection `\Phi`.

AUTHORS:

- Travis Scrimshaw (2014-05-15): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper


class InducedCrystal(UniqueRepresentation, Parent):
    r"""
    A crystal induced from an injection.

    Let `X` be a set and let `C` be crystal and consider any injection
    `\Phi : X \to C`. We induce a crystal structure on `X` by considering
    `\Phi` to be a crystal morphism.

    Alternatively we can induce a crystal structure on some (sub)set of `X`
    by considering an injection `\Phi : C \to X` considered as a crystal
    morphism. This form is also useful when the set `X` is not explicitly
    known.

    INPUT:

    - ``X`` -- the base set
    - ``phi`` -- the map `\Phi`
    - ``inverse`` -- (optional) the inverse map `\Phi^{-1}`
    - ``from_crystal`` -- boolean (default: ``False``); if the induced
      structure is of the second type `\Phi : C \to X`

    EXAMPLES:

    We construct a crystal structure of Gelfand-Tsetlin patterns by going
    through their bijection with semistandard tableaux::

        sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
        sage: G = GelfandTsetlinPatterns(4, 3)
        sage: phi = lambda x: D(x.to_tableau())
        sage: phi_inv = lambda x: G(x.to_tableau())
        sage: I = crystals.Induced(G, phi, phi_inv)
        sage: I.digraph().is_isomorphic(D.digraph(), edge_labels=True)
        True

    Now we construct the above example but inducing the structure going the
    other way (from tableaux to Gelfand-Tsetlin patterns). This can also
    give us more information coming from the crystal. ::

        sage: D2 = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
        sage: G2 = GelfandTsetlinPatterns(4, 1)
        sage: phi2 = lambda x: D2(x.to_tableau())
        sage: phi2_inv = lambda x: G2(x.to_tableau())
        sage: I2 = crystals.Induced(D2, phi2_inv, phi2, from_crystal=True)
        sage: I2.module_generators
        ([[0, 0, 0, 0], [0, 0, 0], [0, 0], [0]],
         [[1, 0, 0, 0], [1, 0, 0], [1, 0], [1]],
         [[1, 1, 0, 0], [1, 1, 0], [1, 1], [1]],
         [[1, 1, 1, 0], [1, 1, 1], [1, 1], [1]],
         [[1, 1, 1, 1], [1, 1, 1], [1, 1], [1]])

    We check an example when the codomain is larger than the domain
    (although here the crystal structure is trivial)::

        sage: P = Permutations(4)
        sage: D = crystals.Tableaux(['A',3], shapes=Partitions(4))
        sage: T = crystals.TensorProduct(D, D)
        sage: phi = lambda p: T(D(RSK(p)[0]), D(RSK(p)[1]))
        sage: phi_inv = lambda d: RSK_inverse(d[0].to_tableau(), d[1].to_tableau(), output='permutation')
        sage: all(phi_inv(phi(p)) == p for p in P) # Check it really is the inverse
        True
        sage: I = crystals.Induced(P, phi, phi_inv)
        sage: I.digraph()
        Digraph on 24 vertices

    We construct an example without a specified inverse map::

        sage: X = Words(2,4)
        sage: L = crystals.Letters(['A',1])
        sage: T = crystals.TensorProduct(*[L]*4)
        sage: Phi = lambda x : T(*[L(i) for i in x])
        sage: I = crystals.Induced(X, Phi)
        sage: I.digraph()
        Digraph on 16 vertices
    """
    @staticmethod
    def __classcall_private__(cls, X, phi, inverse=None, from_crystal=False):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
            sage: G = GelfandTsetlinPatterns(4, 3)
            sage: phi = lambda x: D(x.to_tableau())
            sage: phi_inv = lambda x: G(x.to_tableau())
            sage: I1 = crystals.Induced(G, phi, phi_inv)
            sage: I2 = crystals.Induced(G, phi, phi_inv)
            sage: I1 is I2
            True
        """
        if from_crystal:
            return InducedFromCrystal(X, phi, inverse)

        return super().__classcall__(cls, X, phi, inverse)

    def __init__(self, X, phi, inverse):
        """
        Initialize ``self``.

        TESTS:

        Note that pickling only works when the input functions
        can be pickled::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: import __main__
            sage: __main__.phi = phi
            sage: __main__.phi_inv = phi_inv
            sage: I = crystals.Induced(G, phi, phi_inv)
            sage: TestSuite(I).run()
        """
        try:
            codomain = phi.codomain()
        except AttributeError:
            codomain = phi(X.an_element()).parent()

        self._set = X
        self._phi = phi

        if inverse is None:
            try:
                inverse = ~self._phi
            except (TypeError, ValueError):
                try:
                    inverse = self._phi.section()
                except AttributeError:
                    if X.cardinality() == float('inf'):
                        raise ValueError("the inverse map must be defined for infinite sets")
                    self._preimage = {}
                    for x in X:
                        y = phi(x)
                        if y in self._preimage:
                            raise ValueError("the map is not injective")
                        self._preimage[y] = x
                    inverse = self._preimage.__getitem__
        self._inverse = inverse

        self._cartan_type = codomain.cartan_type()
        Parent.__init__(self, category=codomain.category())

        self.module_generators = self

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: crystals.Induced(G, phi, phi_inv)
            Crystal of Gelfand-Tsetlin patterns of width 4 and max value 1
             induced by <function phi at 0x...>
        """
        return "Crystal of {} induced by {}".format(self._set, self._phi)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: I = crystals.Induced(G, phi, phi_inv)
            sage: I([[1,1,0,0],[1,0,0],[1,0],[1]])
            [[1, 1, 0, 0], [1, 0, 0], [1, 0], [1]]
            sage: I(D(3,2,1))
            [[1, 1, 1, 0], [1, 1, 1], [1, 1], [1]]
            sage: I([[1,1,0,0],[1,0,0],[0,1],[1]])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert [[1, 1, 0, 0], [1, 0, 0], [0, 1], [1]] to Crystal of Gelfand-Tsetlin patterns of width 4 and max value 1 induced by <function phi at 0x...>
        """
        if x in self._set:
            return self.element_class(self, self._set(x))

        try:
            return self.element_class(self, self._inverse(x))
        except (TypeError, ValueError, AttributeError):
            raise TypeError("unable to convert {!r} to {}".format(x, self))

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: I = crystals.Induced(G, phi, phi_inv)
            sage: all(g in I for g in G)
            True
            sage: [[1,1,0,0],[1,0,0],[1,0],[1]] in I
            True
            sage: [[1,1,0,0],[1,0,0],[0,1],[1]] in I
            False
        """
        if isinstance(x, InducedCrystal.Element):
            return x.parent() == self

        return x in self._set

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: I = crystals.Induced(G, phi, phi_inv)
            sage: sorted([x for x in I])
            [[[0, 0, 0, 0], [0, 0, 0], [0, 0], [0]],
             [[1, 0, 0, 0], [0, 0, 0], [0, 0], [0]],
             [[1, 0, 0, 0], [1, 0, 0], [0, 0], [0]],
             [[1, 0, 0, 0], [1, 0, 0], [1, 0], [0]],
             [[1, 0, 0, 0], [1, 0, 0], [1, 0], [1]],
             [[1, 1, 0, 0], [1, 0, 0], [0, 0], [0]],
             [[1, 1, 0, 0], [1, 0, 0], [1, 0], [0]],
             [[1, 1, 0, 0], [1, 0, 0], [1, 0], [1]],
             [[1, 1, 0, 0], [1, 1, 0], [1, 0], [0]],
             [[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]],
             [[1, 1, 0, 0], [1, 1, 0], [1, 1], [1]],
             [[1, 1, 1, 0], [1, 1, 0], [1, 0], [0]],
             [[1, 1, 1, 0], [1, 1, 0], [1, 0], [1]],
             [[1, 1, 1, 0], [1, 1, 0], [1, 1], [1]],
             [[1, 1, 1, 0], [1, 1, 1], [1, 1], [1]],
             [[1, 1, 1, 1], [1, 1, 1], [1, 1], [1]]]
        """
        for x in self._set:
            yield self.element_class(self, x)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: P = Permutations(4)
            sage: D = crystals.Tableaux(['A',3], shapes=Partitions(4))
            sage: T = crystals.TensorProduct(D, D)
            sage: phi = lambda p: T(D(RSK(p)[0]), D(RSK(p)[1]))
            sage: phi_inv = lambda d: RSK_inverse(d[0].to_tableau(), d[1].to_tableau(), output='permutation')
            sage: I = crystals.Induced(P, phi, phi_inv)
            sage: I.cardinality() == factorial(4)
            True
        """
        return self._set.cardinality()

    class Element(ElementWrapper):
        """
        An element of an induced crystal.
        """

        def e(self, i):
            """
            Return `e_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
                sage: G = GelfandTsetlinPatterns(4, 3)
                sage: phi = lambda x: D(x.to_tableau())
                sage: phi_inv = lambda x: G(x.to_tableau())
                sage: I = crystals.Induced(G, phi, phi_inv)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: elt.e(1)
                sage: elt.e(2)
                [[1, 1, 0, 0], [1, 1, 0], [1, 1], [1]]
                sage: elt.e(3)
            """
            P = self.parent()
            ret = P._phi(self.value).e(i)
            if ret is None:
                return None
            try:
                return self.__class__(P, P._inverse(ret))
            except (ValueError, TypeError, AttributeError):
                return None

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
                sage: G = GelfandTsetlinPatterns(4, 3)
                sage: phi = lambda x: D(x.to_tableau())
                sage: phi_inv = lambda x: G(x.to_tableau())
                sage: I = crystals.Induced(G, phi, phi_inv)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: elt.f(1)
                [[1, 1, 0, 0], [1, 1, 0], [1, 0], [0]]
                sage: elt.f(2)
                sage: elt.f(3)
                [[1, 1, 0, 0], [1, 0, 0], [1, 0], [1]]
            """
            P = self.parent()
            ret = P._phi(self.value).f(i)
            if ret is None:
                return None
            try:
                return self.__class__(P, P._inverse(ret))
            except (ValueError, TypeError, AttributeError):
                return None

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
                sage: G = GelfandTsetlinPatterns(4, 3)
                sage: phi = lambda x: D(x.to_tableau())
                sage: phi_inv = lambda x: G(x.to_tableau())
                sage: I = crystals.Induced(G, phi, phi_inv)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: [elt.epsilon(i) for i in I.index_set()]
                [0, 1, 0]
            """
            return self.parent()._phi(self.value).epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
                sage: G = GelfandTsetlinPatterns(4, 3)
                sage: phi = lambda x: D(x.to_tableau())
                sage: phi_inv = lambda x: G(x.to_tableau())
                sage: I = crystals.Induced(G, phi, phi_inv)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: [elt.phi(i) for i in I.index_set()]
                [1, 0, 1]
            """
            return self.parent()._phi(self.value).phi(i)

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,3))
                sage: G = GelfandTsetlinPatterns(4, 3)
                sage: phi = lambda x: D(x.to_tableau())
                sage: phi_inv = lambda x: G(x.to_tableau())
                sage: I = crystals.Induced(G, phi, phi_inv)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: elt.weight()
                (1, 0, 1, 0)
            """
            return self.parent()._phi(self.value).weight()


class InducedFromCrystal(UniqueRepresentation, Parent):
    r"""
    A crystal induced from an injection.

    Alternatively we can induce a crystal structure on some (sub)set of `X`
    by considering an injection `\Phi : C \to X` considered as a crystal
    morphism.

    .. SEEALSO::

        :class:`InducedCrystal`

    INPUT:

    - ``X`` -- the base set
    - ``phi`` -- the map `\Phi`
    - ``inverse`` -- (optional) the inverse map `\Phi^{-1}`

    EXAMPLES:

    We construct a crystal structure on generalized permutations with a
    fixed first row by using RSK::

        sage: C = crystals.Tableaux(['A',3], shape=[2,1])
        sage: def psi(x):
        ....:     ret = RSK_inverse(x.to_tableau(), Tableau([[1,1],[2]]))
        ....:     return (tuple(ret[0]), tuple(ret[1]))
        sage: psi_inv = lambda x: C(RSK(*x)[0])
        sage: I = crystals.Induced(C, psi, psi_inv, from_crystal=True)
    """

    def __init__(self, X, phi, inverse):
        """
        Initialize ``self``.

        TESTS:

        Note that pickling only works when the input functions
        can be pickled::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: import __main__
            sage: __main__.phi = phi
            sage: __main__.phi_inv = phi_inv
            sage: I = crystals.Induced(D, phi_inv, phi, from_crystal=True)
            sage: TestSuite(I).run()
        """
        self._crystal = X
        self._phi = phi

        if inverse is None:
            try:
                inverse = ~self._phi
            except (TypeError, ValueError):
                try:
                    inverse = self._phi.section()
                except AttributeError:
                    if X.cardinality() == float('inf'):
                        raise ValueError("the inverse map must be defined for infinite sets")
                    self._preimage = {}
                    for x in X:
                        y = phi(x)
                        if y in self._preimage:
                            raise ValueError("the map is not injective")
                        self._preimage[y] = x
                    inverse = self._preimage.__getitem__
        self._inverse = inverse

        self._cartan_type = X.cartan_type()
        Parent.__init__(self, category=X.category())
        self.module_generators = tuple(self.element_class(self, phi(mg))
                                       for mg in X.module_generators)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
            sage: G = GelfandTsetlinPatterns(4, 1)
            sage: def phi(x): return D(x.to_tableau())
            sage: def phi_inv(x): return G(x.to_tableau())
            sage: crystals.Induced(D, phi_inv, phi, from_crystal=True)
            Crystal induced by <function phi_inv at 0x...> from
             The crystal of tableaux of type ['A', 3] and shape(s)
             [[], [1], [1, 1], [1, 1, 1], [1, 1, 1, 1]]
        """
        return "Crystal induced by {} from {}".format(self._phi, self._crystal)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[2,1])
            sage: def psi(x):
            ....:     ret = RSK_inverse(x.to_tableau(), Tableau([[1,1],[2]]))
            ....:     return (tuple(ret[0]), tuple(ret[1]))
            sage: psi_inv = lambda x: C(RSK(*x)[0])
            sage: I = crystals.Induced(C, psi, psi_inv, from_crystal=True)
            sage: I([[1, 1, 2], [2, 2, 1]])
            ((1, 1, 2), (2, 2, 1))
            sage: I(C(2,1,3))
            ((1, 1, 2), (2, 3, 1))
        """
        if x in self._crystal:
            return self.element_class(self, self._phi(self._crystal(x)))

        try:
            return self.element_class(self, self._phi(self._inverse(x)))
        except (TypeError, ValueError, AttributeError):
            raise TypeError("unable to convert {!r} to {}".format(x, self))

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[2,1])
            sage: def psi(x):
            ....:     ret = RSK_inverse(x.to_tableau(), Tableau([[1,1],[2]]))
            ....:     return (tuple(ret[0]), tuple(ret[1]))
            sage: psi_inv = lambda x: C(RSK(*x)[0])
            sage: I = crystals.Induced(C, psi, psi_inv, from_crystal=True)
            sage: ((1, 1, 2), (2, 2, 1)) in I
            True
            sage: ((1, 2, 2), (1, 1, 2)) in I
            False
            sage: ((1, 2, 3), (1, 2, 3)) in I
            False
            sage: ((1, 2, 2), (1, 3, 2)) in I
            False
        """
        if isinstance(x, InducedFromCrystal.Element):
            return x.parent() == self

        try:
            y = self._inverse(x)
            return y in self._crystal and self._phi(y) == x
        except (ValueError, TypeError):
            return False

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',2], shape=[2,1])
            sage: def psi(x):
            ....:     ret = RSK_inverse(x.to_tableau(), Tableau([[1,1],[2]]))
            ....:     return (tuple(ret[0]), tuple(ret[1]))
            sage: psi_inv = lambda x: C(RSK(*x)[0])
            sage: I = crystals.Induced(C, psi, psi_inv, from_crystal=True)
            sage: sorted(x for x in I)
            [((1, 1, 2), (1, 2, 1)),
             ((1, 1, 2), (2, 2, 1)),
             ((1, 1, 2), (2, 3, 1)),
             ((1, 1, 2), (3, 3, 1)),
             ((1, 1, 2), (3, 3, 2)),
             ((1, 1, 2), (1, 3, 1)),
             ((1, 1, 2), (1, 3, 2)),
             ((1, 1, 2), (2, 3, 2))]
        """
        for x in self._crystal:
            yield self.element_class(self, self._phi(x))

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[2,1])
            sage: def psi(x):
            ....:     ret = RSK_inverse(x.to_tableau(), Tableau([[1,1],[2]]))
            ....:     return (tuple(ret[0]), tuple(ret[1]))
            sage: psi_inv = lambda x: C(RSK(*x)[0])
            sage: I = crystals.Induced(C, psi, psi_inv, from_crystal=True)
            sage: I.cardinality() == C.cardinality()
            True
        """
        return self._crystal.cardinality()

    class Element(ElementWrapper):
        """
        An element of an induced crystal.
        """

        def e(self, i):
            """
            Return `e_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
                sage: G = GelfandTsetlinPatterns(4, 1)
                sage: def phi(x): return G(x.to_tableau())
                sage: def phi_inv(x): return D(G(x).to_tableau())
                sage: I = crystals.Induced(D, phi, phi_inv, from_crystal=True)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: elt.e(1)
                sage: elt.e(2)
                [[1, 1, 0, 0], [1, 1, 0], [1, 1], [1]]
                sage: elt.e(3)
            """
            P = self.parent()
            ret = P._inverse(self.value).e(i)
            if ret is None:
                return None
            return self.__class__(P, P._phi(ret))

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
                sage: G = GelfandTsetlinPatterns(4, 1)
                sage: def phi(x): return G(x.to_tableau())
                sage: def phi_inv(x): return D(G(x).to_tableau())
                sage: I = crystals.Induced(D, phi, phi_inv, from_crystal=True)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: elt.f(1)
                [[1, 1, 0, 0], [1, 1, 0], [1, 0], [0]]
                sage: elt.f(2)
                sage: elt.f(3)
                [[1, 1, 0, 0], [1, 0, 0], [1, 0], [1]]
            """
            P = self.parent()
            ret = P._inverse(self.value).f(i)
            if ret is None:
                return None
            return self.__class__(P, P._phi(ret))

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
                sage: G = GelfandTsetlinPatterns(4, 1)
                sage: def phi(x): return G(x.to_tableau())
                sage: def phi_inv(x): return D(G(x).to_tableau())
                sage: I = crystals.Induced(D, phi, phi_inv, from_crystal=True)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: [elt.epsilon(i) for i in I.index_set()]
                [0, 1, 0]
            """
            return self.parent()._inverse(self.value).epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
                sage: G = GelfandTsetlinPatterns(4, 1)
                sage: def phi(x): return G(x.to_tableau())
                sage: def phi_inv(x): return D(G(x).to_tableau())
                sage: I = crystals.Induced(D, phi, phi_inv, from_crystal=True)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: [elt.epsilon(i) for i in I.index_set()]
                [0, 1, 0]
            """
            return self.parent()._inverse(self.value).phi(i)

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: D = crystals.Tableaux(['A',3], shapes=PartitionsInBox(4,1))
                sage: G = GelfandTsetlinPatterns(4, 1)
                sage: def phi(x): return G(x.to_tableau())
                sage: def phi_inv(x): return D(G(x).to_tableau())
                sage: I = crystals.Induced(D, phi, phi_inv, from_crystal=True)
                sage: elt = I([[1, 1, 0, 0], [1, 1, 0], [1, 0], [1]])
                sage: elt.weight()
                (1, 0, 1, 0)
            """
            return self.parent()._inverse(self.value).weight()
