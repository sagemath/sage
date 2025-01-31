r"""
Valuations which are defined as limits of valuations.

The discrete valuation of a complete field extends uniquely to a finite field
extension. This is not the case anymore for fields which are not complete with
respect to their discrete valuation. In this case, the extensions essentially
correspond to the factors of the defining polynomial of the extension over the
completion. However, these factors only exist over the completion and this
makes it difficult to write down such valuations with a representation of
finite length.

More specifically, let `v` be a discrete valuation on `K` and let `L=K[x]/(G)`
a finite extension thereof. An extension of `v` to `L` can be represented as a
discrete pseudo-valuation `w'` on `K[x]` which sends `G` to infinity.
However, such `w'` might not be described by an :mod:`augmented valuation <sage.rings.valuation.augmented_valuation>`
over a :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>` anymore. Instead, we may need to write is as a
limit of augmented valuations.

The classes in this module provide the means of writing down such limits and
resulting valuations on quotients.

AUTHORS:

- Julian Rüth (2016-10-19): initial version

EXAMPLES:

In this function field, the unique place of ``K`` which corresponds to the zero
point has two extensions to ``L``. The valuations corresponding to these
extensions can only be approximated::

    sage: # needs sage.rings.function_field
    sage: K.<x> = FunctionField(QQ)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)
    sage: v = K.valuation(1)
    sage: w = v.extensions(L); w
    [[ (x - 1)-adic valuation, v(y + 1) = 1 ]-adic valuation,
     [ (x - 1)-adic valuation, v(y - 1) = 1 ]-adic valuation]

The same phenomenon can be observed for valuations on number fields::

    sage: # needs sage.rings.number_field
    sage: K = QQ
    sage: R.<t> = K[]
    sage: L.<t> = K.extension(t^2 + 1)
    sage: v = QQ.valuation(5)
    sage: w = v.extensions(L); w
    [[ 5-adic valuation, v(t + 2) = 1 ]-adic valuation,
     [ 5-adic valuation, v(t + 3) = 1 ]-adic valuation]

.. NOTE::

    We often rely on approximations of valuations even if we could represent the
    valuation without using a limit. This is done to improve performance as many
    computations already can be done correctly with an approximation::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: v = K.valuation(1/x)
        sage: w = v.extension(L); w
        Valuation at the infinite place
        sage: w._base_valuation._base_valuation._improve_approximation()
        sage: w._base_valuation._base_valuation._approximation
        [ Gauss valuation induced by Valuation at the infinite place,
            v(y) = 1/2, v(y^2 - 1/x) = +Infinity ]

REFERENCES:

Limits of inductive valuations are discussed in [Mac1936I]_ and [Mac1936II]_. An
overview can also be found in Section 4.6 of [Rüt2014]_.
"""
# ****************************************************************************
#       Copyright (C) 2016-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.abstract_method import abstract_method
from .valuation import DiscretePseudoValuation, InfiniteDiscretePseudoValuation
from sage.structure.factory import UniqueFactory


class LimitValuationFactory(UniqueFactory):
    r"""
    Return a limit valuation which sends the polynomial ``G`` to infinity and
    is greater than or equal than ``base_valuation``.

    INPUT:

    - ``base_valuation`` -- a discrete (pseudo-)valuation on a polynomial ring
      which is a discrete valuation on the coefficient ring which can be
      uniquely augmented (possibly only in the limit) to a pseudo-valuation
      that sends ``G`` to infinity.

    - ``G`` -- a squarefree polynomial in the domain of ``base_valuation``

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, QQ.valuation(2))
        sage: w = valuations.LimitValuation(v, x)
        sage: w(x)
        +Infinity
    """
    def create_key(self, base_valuation, G):
        r"""
        Create a key from the parameters of this valuation.

        EXAMPLES:

        Note that this does not normalize ``base_valuation`` in any way. It is
        easily possible to create the same limit in two different ways::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: w = valuations.LimitValuation(v, x)  # indirect doctest
            sage: v = v.augmentation(x, infinity)
            sage: u = valuations.LimitValuation(v, x)
            sage: u == w
            False

        The point here is that this is not meant to be invoked from user code.
        But mostly from other factories which have made sure that the
        parameters are normalized already.
        """
        if not base_valuation.restriction(G.parent().base_ring()).is_discrete_valuation():
            raise ValueError("base_valuation must be discrete on the coefficient ring.")
        return base_valuation, G

    def create_object(self, version, key):
        r"""
        Create an object from ``key``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: w = valuations.LimitValuation(v, x^2 + 1)  # indirect doctest
        """
        base_valuation, G = key
        from .valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(base_valuation.domain())
        return parent.__make_element_class__(MacLaneLimitValuation)(parent, base_valuation, G)


LimitValuation = LimitValuationFactory("sage.rings.valuation.limit_valuation.LimitValuation")


class LimitValuation_generic(DiscretePseudoValuation):
    r"""
    Base class for limit valuations.

    A limit valuation is realized as an approximation of a valuation and means
    to improve that approximation when necessary.

    EXAMPLES::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: v = K.valuation(0)
        sage: w = v.extension(L)
        sage: w._base_valuation
        [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 , … ]

    The currently used approximation can be found in the ``_approximation``
    field::

        sage: w._base_valuation._approximation                                          # needs sage.rings.function_field
        [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]

    TESTS::

        sage: from sage.rings.valuation.limit_valuation import LimitValuation_generic
        sage: isinstance(w._base_valuation, LimitValuation_generic)                     # needs sage.rings.function_field
        True
        sage: TestSuite(w._base_valuation).run()        # long time                     # needs sage.rings.function_field
    """
    def __init__(self, parent, approximation):
        r"""
        TESTS::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<i> = QQ.extension(x^2 + 1)
            sage: v = K.valuation(2)
            sage: from sage.rings.valuation.limit_valuation import LimitValuation_generic
            sage: isinstance(v._base_valuation, LimitValuation_generic)
            True
        """
        DiscretePseudoValuation.__init__(self, parent)

        self._initial_approximation = approximation
        self._approximation = approximation

    def reduce(self, f, check=True):
        r"""
        Return the reduction of ``f`` as an element of the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring`.

        INPUT:

        - ``f`` -- an element in the domain of this valuation of nonnegative
          valuation

        - ``check`` -- whether or not to check that ``f`` has nonnegative
          valuation (default: ``True``)

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 1))
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: w.reduce(y)  # indirect doctest
            u1
        """
        f = self.domain().coerce(f)
        self._improve_approximation_for_reduce(f)
        F = self._approximation.reduce(f, check=check)
        return self.residue_ring()(F)

    def _call_(self, f):
        r"""
        Return the valuation of ``f``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: w(y)  # indirect doctest
            1/2
        """
        self._improve_approximation_for_call(f)
        return self._approximation(f)

    @abstract_method
    def _improve_approximation_for_reduce(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the reduction of ``f``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 1337))

        For the unique extension over the place at 1337, the initial
        approximation is sufficient to compute the reduction of ``y``::

            sage: # needs sage.rings.function_field
            sage: v = K.valuation(1337)
            sage: w = v.extension(L)
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by (x - 1337)-adic valuation, v(y) = 1/2 ]
            sage: w.reduce(y)
            0
            sage: u._approximation
            [ Gauss valuation induced by (x - 1337)-adic valuation, v(y) = 1/2 ]

        However, at a place over 1341, the initial approximation is not sufficient
        for some values (note that 1341-1337 is a square)::

            sage: # needs sage.rings.function_field
            sage: v = K.valuation(1341)
            sage: w = v.extensions(L)[1]
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by (x - 1341)-adic valuation, v(y - 2) = 1 ]
            sage: w.reduce((y - 2) / (x - 1341))  # indirect doctest
            1/4
            sage: u._approximation
            [ Gauss valuation induced by (x - 1341)-adic valuation, v(y - 1/4*x + 1333/4) = 2 ]
            sage: w.reduce((y - 1/4*x + 1333/4) / (x - 1341)^2)  # indirect doctest
            -1/64
            sage: u._approximation
            [ Gauss valuation induced by (x - 1341)-adic valuation,
                v(y + 1/64*x^2 - 1349/32*x + 1819609/64) = 3 ]
        """

    @abstract_method
    def _improve_approximation_for_call(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the valuation of ``f``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 23))

        For the unique extension over the place at 23, the initial
        approximation is sufficient to compute all valuations::

            sage: # needs sage.rings.function_field
            sage: v = K.valuation(23)
            sage: w = v.extension(L)
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by (x - 23)-adic valuation, v(y) = 1/2 ]
            sage: w(x - 23)
            1
            sage: u._approximation
            [ Gauss valuation induced by (x - 23)-adic valuation, v(y) = 1/2 ]

        However, due to performance reasons, sometimes we improve the
        approximation though it would not have been necessary (performing the
        improvement step is faster in this case than checking whether the
        approximation is sufficient)::

            sage: # needs sage.rings.function_field
            sage: w(y)  # indirect doctest
            1/2
            sage: u._approximation
            [ Gauss valuation induced by (x - 23)-adic valuation, v(y) = 1/2, v(y^2 - x + 23) = +Infinity ]
        """

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: w = v.extension(L)
            sage: w._base_valuation # indirect doctest
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2 , … ]
        """
        from sage.rings.infinity import infinity
        from .augmented_valuation import AugmentedValuation_base
        if self._initial_approximation(self._G) is not infinity:
            if isinstance(self._initial_approximation, AugmentedValuation_base):
                return repr(self._initial_approximation)[:-1] + ", … ]"
        return repr(self._initial_approximation)


class MacLaneLimitValuation(LimitValuation_generic, InfiniteDiscretePseudoValuation):
    r"""
    A limit valuation that is a pseudo-valuation on polynomial ring `K[x]`
    which sends a square-free polynomial `G` to infinity.

    This uses the MacLane algorithm to compute the next element in the limit.

    It starts from a first valuation ``approximation`` which has a unique
    augmentation that sends `G` to infinity and whose uniformizer must be a
    uniformizer of the limit and whose residue field must contain the residue
    field of the limit.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<i> = QQ.extension(x^2 + 1)
        sage: v = K.valuation(2)
        sage: u = v._base_valuation; u
        [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 , … ]
    """
    def __init__(self, parent, approximation, G):
        r"""
        TESTS::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<i> = QQ.extension(x^2 + 1)
            sage: v = K.valuation(2)
            sage: u = v._base_valuation
            sage: from sage.rings.valuation.limit_valuation import MacLaneLimitValuation
            sage: isinstance(u, MacLaneLimitValuation)
            True
        """
        LimitValuation_generic.__init__(self, parent, approximation)
        InfiniteDiscretePseudoValuation.__init__(self, parent)

        self._G = G
        self._next_coefficients = None
        self._next_valuations = None

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: v = GaussianIntegers().valuation(2)
            sage: u = v._base_valuation
            sage: u.extensions(QQ['x'])
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 , … ]]
        """
        if self.domain() is ring:
            return [self]
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
        if isinstance(ring, PolynomialRing_generic) and self.domain().base_ring().is_subring(ring.base_ring()):
            if self.domain().base_ring().fraction_field() is ring.base_ring():
                return [LimitValuation(self._initial_approximation.change_domain(ring),
                        self._G.change_ring(ring.base_ring()))]
            else:
                # we need to recompute the mac lane approximants over this base
                # ring because it could split differently
                pass
        return super().extensions(ring)

    def lift(self, F):
        r"""
        Return a lift of ``F`` from the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring` to the domain of
        this valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^4 - x^2 - 2*x - 1)
            sage: v = K.valuation(1)
            sage: w = v.extensions(L)[1]; w
            [ (x - 1)-adic valuation, v(y^2 - 2) = 1 ]-adic valuation
            sage: s = w.reduce(y); s
            u1
            sage: w.lift(s)  # indirect doctest
            y
        """
        F = self.residue_ring().coerce(F)
        return self._initial_approximation.lift(F)

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: w.uniformizer()  # indirect doctest
            y
        """
        return self._initial_approximation.uniformizer()

    def _call_(self, f):
        r"""
        Return the valuation of ``f``.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: K = QQ
            sage: R.<x> = K[]
            sage: vK = K.valuation(2)
            sage: f = (x^2 + 7) * (x^2 + 9)
            sage: V = vK.mac_lane_approximants(f, require_incomparability=True)
            sage: w = valuations.LimitValuation(V[0], f)
            sage: w((x^2 + 7) * (x + 3))
            3/2
            sage: w = valuations.LimitValuation(V[1], f)
            sage: w((x^2 + 7) * (x + 3))
            +Infinity
            sage: w = valuations.LimitValuation(V[2], f)
            sage: w((x^2 + 7) * (x + 3))
            +Infinity
        """
        self._improve_approximation_for_call(f)
        if self._G.divides(f):
            from sage.rings.infinity import infinity
            return infinity
        return self._approximation(f)

    def _improve_approximation(self):
        r"""
        Perform one step of the Mac Lane algorithm to improve our approximation.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: w = v.extension(L)
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2 ]
            sage: u._improve_approximation()
            sage: u._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2, v(t^2 + 1) = +Infinity ]

        This method has no effect, if the approximation is already an infinite
        valuation::

            sage: u._improve_approximation()                                            # needs sage.rings.number_field
            sage: u._approximation                                                      # needs sage.rings.number_field
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2, v(t^2 + 1) = +Infinity ]
        """
        from sage.rings.infinity import infinity
        if self._approximation(self._G) is infinity:
            # an infinite valuation can not be improved further
            return

        approximations = self._approximation.mac_lane_step(self._G,
                          assume_squarefree=True,
                          assume_equivalence_irreducible=True,
                          check=False,
                          principal_part_bound=1 if self._approximation.E() * self._approximation.F() == self._approximation.phi().degree() else None,
                          report_degree_bounds_and_caches=True)
        assert (len(approximations) == 1)
        self._approximation, _, _, self._next_coefficients, self._next_valuations = approximations[0]

    def _improve_approximation_for_call(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the valuation of ``f``.

        EXAMPLES:

        In this examples, the approximation is increased unnecessarily. The
        first approximation would have been precise enough to compute the
        valuation of ``t + 2``. However, it is faster to improve the
        approximation (perform one step of the Mac Lane algorithm) than to
        check for this::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(5)
            sage: w = v.extensions(L)[0]
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 2) = 1 ]
            sage: w(t + 2) # indirect doctest
            1
            sage: u._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 7) = 2 ]

        ALGORITHM:

            Write `L=K[x]/(G)` and consider `g` a representative of the class
            of ``f`` in `K[x]` (of minimal degree.) Write `v` for
            ``self._approximation`` and `\phi` for the last key polynomial of
            `v`. With repeated quotient and remainder `g` has a unique
            expansion as `g=\sum a_i\phi^i`.  Suppose that `g` is an
            equivalence-unit with respect to ``self._approximation``, i.e.,
            `v(a_0) < v(a_i\phi^i)` for all `i\ne 0`. If we denote the limit
            valuation as `w`, then `v(a_i\phi^i)=w(a_i\phi^i)` since the
            valuation of key polynomials does not change during augmentations
            (Theorem 6.4 in [Mac1936II]_.) By the strict triangle inequality,
            `w(g)=v(g)`.
            Note that any `g` which is coprime to `G` is an equivalence-unit
            after finitely many steps of the Mac Lane algorithm. Indeed,
            otherwise the valuation of `g` would be infinite (follows from
            Theorem 5.1 in [Mac1936II]_) since the valuation of the key
            polynomials increases.
            When `f` is not coprime to `G`, consider `s=gcd(f,G)` and write
            `G=st`. Since `G` is squarefree, either `s` or `t` have finite
            valuation. With the above algorithm, this can be decided in
            finitely many steps. From this we can deduce the valuation of `s`
            (and in fact replace `G` with the factor with infinite valuation
            for all future computations.)
        """
        from sage.rings.infinity import infinity
        if self._approximation(self._approximation.phi()) is infinity:
            # an infinite valuation can not be improved further
            return

        if f == 0:
            # zero always has infinite valuation (actually, this might
            # not be desirable for inexact zero elements with leading
            # zero coefficients.)
            return

        while not self._approximation.is_equivalence_unit(f):
            # TODO: I doubt that this really works over inexact fields
            s = self._G.gcd(f)
            if s.is_constant():
                self._improve_approximation()
            else:
                t = self._G // s

                while True:
                    if self._approximation.is_equivalence_unit(s):
                        # t has infinite valuation
                        self._G = t
                        return self._improve_approximation_for_call(f // s)
                    if self._approximation.is_equivalence_unit(t):
                        # s has infinite valuation
                        self._G = s
                        return

                    self._improve_approximation()

    def _improve_approximation_for_reduce(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the reduction of ``f``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(13)
            sage: w = v.extensions(L)[0]
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by 13-adic valuation, v(t + 5) = 1 ]
            sage: w.reduce((t + 5) / 13) # indirect doctest
            8
            sage: u._approximation
            [ Gauss valuation induced by 13-adic valuation, v(t - 29/2) = 2 ]

        ALGORITHM:

            The reduction produced by the approximation is correct for an
            equivalence-unit, see :meth:`_improve_approximation_for_call`.
        """
        if self._approximation(f) > 0:
            return
        self._improve_approximation_for_call(f)

    def residue_ring(self):
        r"""
        Return the residue ring of this valuation, which is always a field.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: w = v.extension(L)
            sage: w.residue_ring()
            Finite Field of size 2
        """
        R = self._initial_approximation.residue_ring()
        from sage.categories.fields import Fields
        if R in Fields():
            # the approximation ends in v(phi)=infty
            return R
        else:
            from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
            assert (isinstance(R, PolynomialRing_generic))
            return R.base_ring()

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: R.<x> = QQ[]
            sage: F = (x^2 + 7) * (x^2 + 9)
            sage: G = (x^2 + 7)
            sage: V = QQ.valuation(2).mac_lane_approximants(F, require_incomparability=True)
            sage: valuations.LimitValuation(V[0], F) >= valuations.LimitValuation(V[1], F)
            False
            sage: valuations.LimitValuation(V[1], F) >= valuations.LimitValuation(V[1], G)
            True
            sage: valuations.LimitValuation(V[2], F) >= valuations.LimitValuation(V[2], G)
            True
        """
        if other.is_trivial():
            return other.is_discrete_valuation()
        if isinstance(other, MacLaneLimitValuation):
            if self._approximation.restriction(self._approximation.domain().base_ring()) == other._approximation.restriction(other._approximation.domain().base_ring()):
                # Two MacLane limit valuations v,w over the same constant
                # valuation are either equal or incomparable; neither v>w nor
                # v<w can hold everywhere.
                # They are equal iff they approximate the same factor of their
                # defining G. Note that they can be equal even if the defining
                # G is different, so we need to make sure that this can not be
                # the case.
                self._improve_approximation_for_call(other._G)
                other._improve_approximation_for_call(self._G)
                if self._G != other._G:
                    assert self._G.gcd(other._G).is_one()
                    return False

                # If the valuations are comparable, they must approximate the
                # same factor of G (see the documentation of LimitValuation:
                # the approximation must *uniquely* single out a valuation.)
                return (self._initial_approximation >= other._initial_approximation
                        or self._initial_approximation <= other._initial_approximation)

        return super()._ge_(other)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: w = v.extension(L)
            sage: w._base_valuation.restriction(K)
            2-adic valuation
        """
        if ring.is_subring(self.domain().base()):
            return self._initial_approximation.restriction(ring)
        return super().restriction(ring)

    def _weakly_separating_element(self, other):
        r"""
        Return an element in the domain of this valuation which has
        positive valuation with respect to this valuation and higher
        valuation with respect to this valuation than with respect to
        ``other``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: w = v.extension(L)
            sage: v = QQ.valuation(5)
            sage: u,uu = v.extensions(L)
            sage: w._base_valuation._weakly_separating_element(u._base_valuation)   # long time
            2
            sage: u._base_valuation._weakly_separating_element(uu._base_valuation)  # long time
            t + 2

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: v = K.valuation(1/x)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - 1/(x^2 + 1))
            sage: u,uu = v.extensions(L)
            sage: v = K.valuation(x)
            sage: w,ww = v.extensions(L)
            sage: v = K.valuation(1)
            sage: v = v.extension(L)
            sage: u.separating_element([uu,ww,w,v])  # random output                # long time
            ((8*x^4 + 12*x^2 + 4)/(x^2 - x))*y + (8*x^4 + 8*x^2 + 1)/(x^3 - x^2)

        The underlying algorithm is quite naive and might not terminate in
        reasonable time. In particular, the order of the arguments sometimes
        has a huge impact on the runtime::

            sage: u.separating_element([ww,w,v,uu])  # not tested, takes forever
        """
        from .scaled_valuation import ScaledValuation_generic
        v = self.restriction(self.domain().base())
        if isinstance(v, ScaledValuation_generic):
            v = v._base_valuation
        u = other.restriction(self.domain().base())
        if isinstance(u, ScaledValuation_generic):
            u = u._base_valuation

        if u == v:
            # phi of the initial approximant must be good enough to separate it
            # from any other approximant of an extension
            ret = self._initial_approximation.phi()
            assert (self(ret) > other(ret))  # I could not come up with an example where this fails
            return ret
        else:
            # if the valuations are sane, it should be possible to separate
            # them with constants
            return self.domain()(v._weakly_separating_element(u))

    def value_semigroup(self):
        r"""
        Return the value semigroup of this valuation.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(5)
            sage: u,uu = v.extensions(L)
            sage: u.value_semigroup()
            Additive Abelian Semigroup generated by -1, 1
        """
        return self._initial_approximation.value_semigroup()

    def element_with_valuation(self, s):
        r"""
        Return an element with valuation ``s``.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: u = v.extension(L)
            sage: u.element_with_valuation(1/2)
            t + 1
        """
        return self._initial_approximation.element_with_valuation(s)

    def _relative_size(self, f):
        r"""
        Return an estimate on the coefficient size of ``f``.

        The number returned is an estimate on the factor between the number of
        bits used by ``f`` and the minimal number of bits used by an element
        congruent to ``f``.

        This is used by :meth:`simplify` to decide whether simplification of
        coefficients is going to lead to a significant shrinking of the
        coefficients of ``f``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: u = v.extension(L)
            sage: u._relative_size(1024*t + 1024)
            6
        """
        return self._initial_approximation._relative_size(f)

    def simplify(self, f, error=None, force=False):
        r"""
        Return a simplified version of ``f``.

        Produce an element which differs from ``f`` by an element of valuation
        strictly greater than the valuation of ``f`` (or strictly greater than
        ``error`` if set.)

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: u = v.extension(L)
            sage: u.simplify(t + 1024, force=True)
            t
        """
        f = self.domain().coerce(f)

        self._improve_approximation_for_call(f)
        # now _approximation is sufficiently precise to compute a valid
        # simplification of f

        if error is None:
            error = self(f) if force else self.upper_bound(f)

        return self._approximation.simplify(f, error=error, force=force)

    def lower_bound(self, f):
        r"""
        Return a lower bound of this valuation at ``x``.

        Use this method to get an approximation of the valuation of ``x``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: u = v.extension(L)
            sage: u.lower_bound(1024*t + 1024)
            10
            sage: u(1024*t + 1024)
            21/2
        """
        f = self.domain().coerce(f)
        return self._approximation.lower_bound(f)

    def upper_bound(self, f):
        r"""
        Return an upper bound of this valuation at ``x``.

        Use this method to get an approximation of the valuation of ``x``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: u = v.extension(L)
            sage: u.upper_bound(1024*t + 1024)
            21/2
            sage: u(1024*t + 1024)
            21/2
        """
        f = self.domain().coerce(f)
        self._improve_approximation_for_call(f)
        return self._approximation.upper_bound(f)

    def is_negative_pseudo_valuation(self):
        r"""
        Return whether this valuation attains `-\infty`.

        EXAMPLES:

        For a Mac Lane limit valuation, this is never the case, so this
        method always returns ``False``::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = QQ.valuation(2)
            sage: u = v.extension(L)
            sage: u.is_negative_pseudo_valuation()
            False
        """
        return False
