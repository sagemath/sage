from sage.misc.latex import latex
from sage.misc.prandom import shuffle
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.hyperelliptic_curves_smooth_model.weighted_projective_homset import (
    SchemeHomset_points_weighted_projective_ring,
)
from sage.schemes.projective.projective_space import ProjectiveSpace, _CommRings
from sage.structure.all import UniqueRepresentation
from sage.structure.category_object import normalize_names


def WeightedProjectiveSpace(weights, R=None, names=None):
    r"""
    Return a weighted projective space with the given ``weights`` over the ring ``R``.

    EXAMPLES::

        sage: # TODO: add example of point on this space (it doesn't work right now)
        sage: WP = WeightedProjectiveSpace([1, 3, 1]); WP
        Weighted Projective Space of dimension 2 with weights (1, 3, 1) over Integer Ring
    """
    if (
        isinstance(weights, (MPolynomialRing_base, PolynomialRing_general))
        and R is None
    ):
        if names is not None:
            # Check for the case that the user provided a variable name
            # That does not match what we wanted to use from R
            names = normalize_names(weights.ngens(), names)
            if weights.variable_names() != names:
                # The provided name doesn't match the name of R's variables
                raise NameError(
                    "variable names passed to ProjectiveSpace conflict with names in ring"
                )
        A = WeightedProjectiveSpace(
            weights.ngens() - 1, weights.base_ring(), names=weights.variable_names()
        )
        A._coordinate_ring = weights
        return A

    if isinstance(R, (int, Integer, list, tuple)):
        weights, R = R, weights
    elif R is None:
        R = ZZ

    # WeightedProjectiveSpace(5) -> just return unweighted version
    if isinstance(weights, (int, Integer)):
        return ProjectiveSpace(weights, R=R, names=names)
    elif isinstance(weights, (list, tuple)):
        # Make it hashable
        weights = tuple(map(Integer, weights))
        if any(w <= 0 for w in weights):
            raise TypeError(
                f"weights(={weights}) should only consist of positive integers"
            )
    else:
        raise TypeError(f"weights={weights} must be an integer, a list or a tuple")

    if names is None:
        names = "x"

    # TODO: Specialise implementation to projective spaces over non-rings. But
    # since we don't really implement extra functionalities, I don't think we
    # care.
    if R in _CommRings:
        return WeightedProjectiveSpace_ring(weights, R=R, names=names)

    raise TypeError(f"R (={R}) must be a commutative ring")


class WeightedProjectiveSpace_ring(UniqueRepresentation, AmbientSpace):
    @staticmethod
    def __classcall__(cls, weights: tuple[Integer], R=ZZ, names=None):
        # __classcall_ is the "preprocessing" step for UniqueRepresentation
        # see docs of CachedRepresentation
        # weights should be a tuple, also because it should be hashable
        if not isinstance(weights, tuple):
            raise TypeError(
                f"weights(={weights}) is not a tuple. Please use the `WeightedProjectiveSpace`"
                " constructor"
            )

        # TODO: Do we normalise the weights to make it coprime?
        normalized_names = normalize_names(len(weights), names)
        return super().__classcall__(cls, weights, R, normalized_names)

    def __init__(self, weights: tuple[Integer], R=ZZ, names=None):
        """
        Initialization function.

        EXAMPLES::

            sage: WeightedProjectiveSpace(Zp(5), [1, 3, 1], 'y')                        # needs sage.rings.padics
            Weighted Projective Space of dimension 2 with weights (1, 3, 1) over 5-adic Ring with
            capped relative precision 20
            sage: WeightedProjectiveSpace(QQ, 5, 'y')
            Projective Space of dimension 5 over Rational Field
            sage: _ is ProjectiveSpace(QQ, 5, 'y')
            True
        """
        AmbientSpace.__init__(self, len(weights) - 1, R)
        self._weights = weights
        self._assign_names(names)

    def weights(self) -> tuple[Integer]:
        """
        Return the tuple of weights of this weighted projective space.

        EXAMPLES::

            sage: WeightedProjectiveSpace(QQ, [1, 3, 1]).weights()
            (1, 3, 1)
        """
        return self._weights

    def ngens(self) -> Integer:
        """
        Return the number of generators of this weighted projective space.

        This is the number of variables in the coordinate ring of ``self``.

        EXAMPLES::

            sage: WeightedProjectiveSpace(QQ, [1, 3, 1]).ngens()
            3
            sage: WeightedProjectiveSpace(ZZ, 5).ngens()
            6
        """
        return self.dimension_relative() + 1

    def _check_satisfies_equations(self, v: list[Integer] | tuple[Integer]) -> bool:
        """
        Return ``True`` if ``v`` defines a point on the weighted projective
        plane; raise a :class:`TypeError` otherwise.

        EXAMPLES::

            sage: P = WeightedProjectiveSpace(ZZ, [1, 3, 1])
            sage: P._check_satisfies_equations([1, 1, 0])
            True

            sage: P._check_satisfies_equations([1, 0])
            Traceback (most recent call last):
            ...
            TypeError: the list v=[1, 0] must have 3 components

            sage: P._check_satisfies_equations([1/2, 0, 1])
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
        """
        if not isinstance(v, (list, tuple)):
            raise TypeError(f"the argument v={v} must be a list or tuple")

        n = self.ngens()
        if len(v) != n:
            raise TypeError(f"the list v={v} must have {n} components")

        v = list(map(self.base_ring(), v))
        if all(vi.is_zero() for vi in v):
            raise TypeError("the zero vector is not a point in projective space")

        return True

    def coordinate_ring(self) -> PolynomialRing_general:
        """
        Return the coordinate ring of this weighted projective space.

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace(GF(19^2, 'α'), [1, 3, 4, 1], 'abcd')
            sage: # needs sage.rings.finite_rings
            sage: R = WP.coordinate_ring(); R
            Multivariate Polynomial Ring in a, b, c, d over Finite Field in α of size 19^2
            sage: R.term_order()
            Weighted degree reverse lexicographic term order with weights (1, 3, 4, 1)

        ::

            sage: WP = WeightedProjectiveSpace(QQ, [1, 1, 1], ['alpha', 'beta', 'gamma'])
            sage: R = WP.coordinate_ring(); R
            Multivariate Polynomial Ring in alpha, beta, gamma over Rational Field
            sage: R.term_order()
            Weighted degree reverse lexicographic term order with weights (1, 1, 1)
        """
        if not hasattr(self, "_coordinate_ring"):
            term_order = TermOrder("wdegrevlex", self.weights())
            self._coordinate_ring = PolynomialRing(
                self.base_ring(),
                self.dimension_relative() + 1,
                names=self.variable_names(),
                order=term_order,
            )
        return self._coordinate_ring

    def _validate(self, polynomials):
        """
        If ``polynomials`` is a tuple of valid polynomial functions on ``self``,
        return ``polynomials``, otherwise raise ``TypeError``.

        Since this is a weighted projective space, polynomials must be
        homogeneous with respect to the grading of this space.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of
            this space.

        OUTPUT:

        - tuple of polynomials in the coordinate ring of this space.

        EXAMPLES::

            sage: P.<x, y, z> = WeightedProjectiveSpace(QQ, [1, 3, 1])
            sage: P._validate([x*y - z^4, x])
            [x*y - z^4, x]
            sage: P._validate([x*y - z^2, x])
            Traceback (most recent call last):
            ...
            TypeError: x*y - z^2 is not homogeneous with weights (1, 3, 1)
            sage: P._validate(x*y - z)
            Traceback (most recent call last):
            ...
            TypeError: the argument polynomials=x*y - z must be a list or tuple
        """
        if not isinstance(polynomials, (list, tuple)):
            raise TypeError(
                f"the argument polynomials={polynomials} must be a list or tuple"
            )

        R = self.coordinate_ring()
        for f in map(R, polynomials):
            if not f.is_homogeneous():
                raise TypeError(f"{f} is not homogeneous with weights {self.weights()}")

        return polynomials

    def _latex_(self):
        r"""
        Return a LaTeX representation of this projective space.

        EXAMPLES::

            sage: print(latex(WeightedProjectiveSpace(ZZ, [1, 3, 1], 'x')))
            {\mathbf P}_{\Bold{Z}}^{[1, 3, 1]}

        TESTS::

            sage: WeightedProjectiveSpace(Zp(5), [2, 1, 3], 'y')._latex_()              # needs sage.rings.padics
            '{\\mathbf P}_{\\Bold{Z}_{5}}^{[2, 1, 3]}'
        """
        return (
            f"{{\\mathbf P}}_{{{latex(self.base_ring())}}}^{{{list(self.weights())}}}"
        )

    def _morphism(self, *_, **__):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.
        """
        raise NotImplementedError(
            "_morphism not implemented for weighted projective space"
        )

    def _homset(self, *_, **__):
        """
        Construct the Hom-set.
        """
        raise NotImplementedError(
            "_homset not implemented for weighted projective space"
        )

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.
        """
        # raise NotImplementedError("_point_homset not implemented for weighted projective space")
        return SchemeHomset_points_weighted_projective_ring(*args, **kwds)

    def point(self, v, check=True):
        """
        Create a point on this weighted projective space.

        INPUT:

        INPUT:

        - ``v`` -- anything that defines a point

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: A point of this weighted projective space.

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace(QQ, [1, 3, 1])
            sage: WP.point([2, 3, 1])
            (2 : 3 : 1)
        """
        from sage.rings.infinity import infinity

        if v is infinity or (
            isinstance(v, (list, tuple)) and len(v) == 1 and v[0] is infinity
        ):
            if self.dimension_relative() > 1:
                raise ValueError("%s not well defined in dimension > 1" % v)
            v = [1, 0]

        return self.point_homset()(v, check=check)

    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.
        """
        from sage.schemes.hyperelliptic_curves_smooth_model.weighted_projective_point import (
            SchemeMorphism_point_weighted_projective_ring,
        )
        return SchemeMorphism_point_weighted_projective_ring(*args, **kwds)

    def _repr_(self) -> str:
        """
        Return a string representation of this projective space.

        EXAMPLES::

            sage: WeightedProjectiveSpace(Qp(5), [1, 3, 1], 'x')                        # needs sage.rings.padics
            Weighted Projective Space of dimension 2 with weights (1, 3, 1)
            over 5-adic Field with capped relative precision 20
        """
        return (
            f"Weighted Projective Space of dimension {self.dimension_relative()} with weights"
            f" {self.weights()} over {self.base_ring()}"
        )

    def _an_element_(self):
        r"""
        Return a (preferably typical) element of this space.

        This is used both for illustration and testing purposes.

        OUTPUT: a point in this projective space.

        EXAMPLES::

            sage: # TODO: Enable this
            sage: WeightedProjectiveSpace(ZZ, [1, 3, 1], 'x').an_element()  # random
            (7 : 6 : 5 : 1)
            sage: WeightedProjectiveSpace(ZZ["y"], [2, 3, 1], 'x').an_element()  # random
            (7*y : 6*y : 5*y : 1)
        """
        n = self.dimension_relative()
        R = self.base_ring()
        coords = [(n + 1 - i) * R.an_element() for i in range(n)] + [R.one()]
        shuffle(coords)
        return self(coords)

    def subscheme(self, *_, **__):
        raise NotImplementedError("subscheme of weighted projective space has not been implemented")
