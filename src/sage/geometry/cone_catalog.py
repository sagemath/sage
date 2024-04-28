r"""
Catalog of common polyhedral convex cones

This module provides shortcut functions, grouped under the
globally-available ``cones`` prefix, to create some common cones:

- The nonnegative orthant,
- The rearrangement cone of order ``p``,
- The Schur cone,
- The trivial cone.

At the moment, only convex rational polyhedral cones are
supported---specifically, those cones that can be built using the
:func:`Cone` constructor. As a result, each shortcut method can be
passed either an ambient dimension ``ambient_dim``, or a toric
``lattice`` (from which the dimension can be inferred) to determine
the ambient space.

Here are some typical usage examples::

    sage: cones.nonnegative_orthant(2).rays()
    N(1, 0),
    N(0, 1)
    in 2-d lattice N

::

    sage: cones.rearrangement(2,2).rays()
    N( 1,  0),
    N( 1, -1),
    N(-1,  1)
    in 2-d lattice N

::

    sage: cones.schur(3).rays()
    N(1, -1,  0),
    N(0,  1, -1)
    in 3-d lattice N

::

    sage: cones.trivial(3).rays()
    Empty collection
    in 3-d lattice N

To specify some other lattice, pass it as an argument to the function::

    sage: K = cones.nonnegative_orthant(3)
    sage: cones.schur(lattice=K.dual().lattice())
    2-d cone in 3-d lattice M

For more information about these cones, see the documentation for the
individual functions and the references therein.
"""

#
# WARNING:
#
# Don't import anything into global scope here. If you do, it will
# show up in the cones.<tab> list, and that list should contain only
# the top-level non-underscore functions defined in this module.
#

def _preprocess_args(ambient_dim, lattice):
    r"""
    Preprocess arguments for cone-constructing functions.

    Each cone-constructing function in this module performs some
    preprocessing on its ``ambient_dim`` and ``lattice`` arguments to
    ensure that they both wind up defined and compatible with one
    another. In particular, the ``ambient_dim`` can be inferred from
    a ``lattice`` and a default lattice can be chosen based on an
    ``ambient_dim``. Both can be specified so long as the rank of
    the ``lattice`` matches the ``ambient_dim``; it is an error
    if the two numbers are different.

    Since each function performs exactly the same steps,

    - Check that either ``ambient_dim`` or ``lattice`` was given,
    - Infer the ``ambient_dim`` from the ``lattice`` if ``lattice``
      was given and ``ambient_dim`` was not,
    - Use the default ``lattice`` if ``ambient_dim`` was given and
      ``lattice`` was not,
    - Ensure that the rank of the ``lattice`` equals ``ambient_dim``,

    we collect them here.

    INPUT:

    - ``ambient_dim`` -- a nonnegative integer; the dimension of the
      ambient space in which the cone will live

    - ``lattice`` -- a toric lattice; the lattice in which the cone
      will live

    OUTPUT:

    A pair ``(ambient_dim, lattice)`` containing the processed
    arguments. Both are guaranteed to be defined (that is, not
    ``None``) and the rank of the ``lattice`` will equal
    ``ambient_dim``.

    EXAMPLES:

    If the ``lattice`` argument is ``None``, then the default lattice
    with rank equal to ``ambient_dim`` will be used::

        sage: from sage.geometry.cone_catalog import _preprocess_args
        sage: _preprocess_args(3, None)
        (3, 3-d lattice N)

    If the ``ambient_dim`` argument is ``None``, then the rank of the
    ``lattice`` will be used as the ambient dimension::

        sage: from sage.geometry.cone_catalog import _preprocess_args
        sage: _preprocess_args(None, ToricLattice(3))
        (3, 3-d lattice N)

    So long as the rank of ``lattice`` equals ``ambient_dim``, it is
    not an error for both to be defined::

        sage: from sage.geometry.cone_catalog import _preprocess_args
        sage: _preprocess_args(3, ToricLattice(3))
        (3, 3-d lattice N)

    .. NOTE::

        The two error conditions are tested in the cone-constructing
        functions themselves, which introduces some duplication but
        serves to ensure that ``_preprocess_args()`` is actually called
        in each such function.
    """
    from sage.geometry.toric_lattice import ToricLattice

    if ambient_dim is None and lattice is None:
        raise ValueError("either the ambient dimension or the lattice "
                         "must be specified")

    if ambient_dim is None:
        ambient_dim = lattice.rank()

    if lattice is None:
        lattice = ToricLattice(ambient_dim)

    if lattice.rank() != ambient_dim:
        raise ValueError("lattice rank=%d and ambient_dim=%d "
                         "are incompatible" % (lattice.rank(), ambient_dim))

    return (ambient_dim, lattice)


def nonnegative_orthant(ambient_dim=None, lattice=None):
    r"""
    The nonnegative orthant in ``ambient_dim`` dimensions, or living
    in ``lattice``.

    The nonnegative orthant consists of all componentwise-nonnegative
    vectors. It is the convex-conic hull of the standard basis.

    INPUT:

    - ``ambient_dim`` -- a nonnegative integer (default: ``None``); the
      dimension of the ambient space

    - ``lattice`` -- a toric lattice (default: ``None``); the lattice in
      which the cone will live

    If ``ambient_dim`` is omitted, then it will be inferred from the
    rank of ``lattice``. If the ``lattice`` is omitted, then the
    default lattice of rank ``ambient_dim`` will be used.

    A ``ValueError`` is raised if neither ``ambient_dim`` nor
    ``lattice`` are specified. It is also a ``ValueError`` to specify
    both ``ambient_dim`` and ``lattice`` unless the rank of
    ``lattice`` is equal to ``ambient_dim``.

    OUTPUT:

    A :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone` living in ``lattice``
    and having ``ambient_dim`` standard basis vectors as its
    generators. Each generating ray has the integer ring as its
    base ring.

    A ``ValueError`` can be raised if the inputs are incompatible or
    insufficient. See the INPUT documentation for details.

    REFERENCES:

    - Chapter 2 in [BV2009]_ (Examples 2.4, 2.14, and 2.23 in particular)

    EXAMPLES::

        sage: cones.nonnegative_orthant(3).rays()
        N(1, 0, 0),
        N(0, 1, 0),
        N(0, 0, 1)
        in 3-d lattice N

    TESTS:

    We can construct the trivial cone as the nonnegative orthant in a
    trivial vector space::

        sage: cones.nonnegative_orthant(0)
        0-d cone in 0-d lattice N

    The nonnegative orthant is a proper cone::

        sage: ambient_dim = ZZ.random_element(10)
        sage: K = cones.nonnegative_orthant(ambient_dim)
        sage: K.is_proper()
        True

    If a ``lattice`` was given, it is actually used::

        sage: L = ToricLattice(3, 'M')
        sage: cones.nonnegative_orthant(lattice=L)
        3-d cone in 3-d lattice M

    Unless the rank of the lattice disagrees with ``ambient_dim``::

        sage: L = ToricLattice(1, 'M')
        sage: cones.nonnegative_orthant(3, lattice=L)
        Traceback (most recent call last):
        ...
        ValueError: lattice rank=1 and ambient_dim=3 are incompatible

    We also get an error if no arguments are given::

        sage: cones.nonnegative_orthant()
        Traceback (most recent call last):
        ...
        ValueError: either the ambient dimension or the lattice must
        be specified
    """
    from sage.geometry.cone import Cone
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ

    (ambient_dim, lattice) = _preprocess_args(ambient_dim, lattice)

    I = matrix.identity(ZZ, ambient_dim)
    return Cone(I.rows(), lattice)


def rearrangement(p, ambient_dim=None, lattice=None):
    r"""
    The rearrangement cone of order ``p`` in ``ambient_dim``
    dimensions, or living in ``lattice``.

    The rearrangement cone of order ``p`` in ``ambient_dim``
    dimensions consists of all vectors of length ``ambient_dim``
    whose smallest ``p`` components sum to a nonnegative number.

    For example, the rearrangement cone of order one has its single
    smallest component nonnegative. This implies that all components
    are nonnegative, and that therefore the rearrangement cone of
    order one is the nonnegative orthant in its ambient space.

    When ``p`` and ``ambient_dim`` are equal, all components of the
    cone's elements must sum to a nonnegative number. In other
    words, the rearrangement cone of order ``ambient_dim`` is a
    half-space.

    INPUT:

    - ``p`` -- a nonnegative integer; the number of components to
      "rearrange", between ``1`` and ``ambient_dim`` inclusive

    - ``ambient_dim`` -- a nonnegative integer (default: ``None``); the
      dimension of the ambient space

    - ``lattice`` -- a toric lattice (default: ``None``); the lattice in
      which the cone will live

    If ``ambient_dim`` is omitted, then it will be inferred from the
    rank of ``lattice``. If the ``lattice`` is omitted, then the
    default lattice of rank ``ambient_dim`` will be used.

    A ``ValueError`` is raised if neither ``ambient_dim`` nor
    ``lattice`` are specified. It is also a ``ValueError`` to specify
    both ``ambient_dim`` and ``lattice`` unless the rank of
    ``lattice`` is equal to ``ambient_dim``.

    It is also a ``ValueError`` to specify a non-integer ``p``.

    OUTPUT:

    A :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone` representing the
    rearrangement cone of order ``p`` living in ``lattice``, with
    ambient dimension ``ambient_dim``. Each generating ray has the
    integer ring as its base ring.

    A ``ValueError`` can be raised if the inputs are incompatible or
    insufficient. See the INPUT documentation for details.

    ALGORITHM:

    Suppose that the ambient space is of dimension `n`. The extreme
    directions of the rearrangement cone for `1 \le p \le n-1` are
    given by [Jeong2017]_ Theorem 5.2.3. When `2 \le p \le n-2` (that
    is, if we ignore `p = 1` and `p = n-1`), they consist of

    - the standard basis `\left\{e_{1},e_{2},\ldots,e_{n}\right\}` for
      the ambient space, and

    - the `n` vectors `\left(1,1,\ldots,1\right)^{T} - pe_{i}` for
      `i = 1,2,\ldots,n`.

    Special cases are then given for `p = 1` and `p = n-1` in the
    theorem. However in SageMath we don't need conically-independent
    extreme directions. We only need a generating set, because the
    :func:`Cone` function will eliminate any redundant generators. And
    one can easily verify that the special-case extreme directions for
    `p = 1` and `p = n-1` are contained in the conic hull of the `2n`
    generators just described. The half space resulting from `p = n`
    is also covered by this set of generators, so for all valid `p` we
    simply take the conic hull of those `2n` vectors.

    REFERENCES:

    - [GJ2016]_, Section 4

    - [HS2010]_, Example 2.21

    - [Jeong2017]_, Section 5.2

    EXAMPLES:

    The rearrangement cones of order one are nonnegative orthants::

        sage: orthant = cones.nonnegative_orthant(6)
        sage: cones.rearrangement(1,6).is_equivalent(orthant)
        True

    When ``p`` and ``ambient_dim`` are equal, the rearrangement cone
    is a half-space, so we expect its lineality to be one less than
    ``ambient_dim`` because it will contain a hyperplane but is not
    the entire space::

        sage: cones.rearrangement(5,5).lineality()
        4

    Jeong's Proposition 5.2.1 [Jeong2017]_ states that all rearrangement
    cones are proper when ``p`` is less than ``ambient_dim``::

        sage: all( cones.rearrangement(p, ambient_dim).is_proper()
        ....:              for ambient_dim in range(10)
        ....:              for p in range(1, ambient_dim) )
        True

    Jeong's Corollary 5.2.4 [Jeong2017]_ states that if `p = n-1` in
    an `n`-dimensional ambient space, then the Lyapunov rank of the
    rearrangement cone is `n`, and that for all other `p > 1` its
    Lyapunov rank is one::

        sage: all( cones.rearrangement(p, ambient_dim).lyapunov_rank()
        ....:      ==
        ....:      ambient_dim
        ....:              for ambient_dim in range(2, 10)
        ....:              for p in [ ambient_dim-1 ] )
        True
        sage: all( cones.rearrangement(p, ambient_dim).lyapunov_rank() == 1
        ....:              for ambient_dim in range(3, 10)
        ....:              for p in range(2, ambient_dim-1) )
        True

    TESTS:

    Jeong's Proposition 5.2.1 [Jeong2017]_ states that rearrangement
    cones are permutation-invariant::

        sage: ambient_dim = ZZ.random_element(2,10).abs()
        sage: p = ZZ.random_element(1, ambient_dim)
        sage: K = cones.rearrangement(p, ambient_dim)
        sage: P = SymmetricGroup(ambient_dim).random_element().matrix()                 # needs sage.groups
        sage: all(K.contains(P*r) for r in K)                                           # needs sage.groups
        True

    The smallest ``p`` components of every element of the rearrangement
    cone should sum to a nonnegative number. In other words, the
    generators really are what we think they are::

        sage: def _has_rearrangement_property(v,p):
        ....:     return sum( sorted(v)[0:p] ) >= 0
        sage: all(
        ....:   _has_rearrangement_property(
        ....:     cones.rearrangement(p, ambient_dim).random_element(),
        ....:     p
        ....:   )
        ....:   for ambient_dim in range(2, 10)
        ....:   for p in range(1, ambient_dim+1)
        ....: )
        True

    The rearrangement cone of order ``p`` is, almost by definition,
    contained in the rearrangement cone of order ``p + 1``::

        sage: ambient_dim = ZZ.random_element(2,10)
        sage: p = ZZ.random_element(1, ambient_dim)
        sage: K1 = cones.rearrangement(p, ambient_dim)
        sage: K2 = cones.rearrangement(p+1, ambient_dim)
        sage: all( x in K2 for x in K1 )
        True

    Jeong's Proposition 5.2.1 [Jeong2017]_ states that the rearrangement
    cone of order ``p`` is linearly isomorphic to the rearrangement
    cone of order ``ambient_dim - p`` when ``p`` is less than
    ``ambient_dim``::

        sage: ambient_dim = ZZ.random_element(2,10)
        sage: p = ZZ.random_element(1, ambient_dim)
        sage: K1 = cones.rearrangement(p, ambient_dim)
        sage: K2 = cones.rearrangement(ambient_dim-p, ambient_dim)
        sage: Mp = ((1/p)*matrix.ones(QQ, ambient_dim)
        ....:    - matrix.identity(QQ, ambient_dim))
        sage: Cone( (Mp*K2.rays()).columns() ).is_equivalent(K1)
        True

    The order ``p`` should be an integer between ``1`` and
    ``ambient_dim``, inclusive::

        sage: cones.rearrangement(0,3)
        Traceback (most recent call last):
        ...
        ValueError: order p=0 should be an integer between 1 and
        ambient_dim=3, inclusive
        sage: cones.rearrangement(5,3)
        Traceback (most recent call last):
        ...
        ValueError: order p=5 should be an integer between 1 and
        ambient_dim=3, inclusive
        sage: cones.rearrangement(3/2, 3)
        Traceback (most recent call last):
        ...
        ValueError: order p=3/2 should be an integer between 1 and
        ambient_dim=3, inclusive

    If a ``lattice`` was given, it is actually used::

        sage: L = ToricLattice(3, 'M')
        sage: cones.rearrangement(2, 3, lattice=L)
        3-d cone in 3-d lattice M

    Unless the rank of the lattice disagrees with ``ambient_dim``::

        sage: L = ToricLattice(1, 'M')
        sage: cones.rearrangement(2, 3, lattice=L)
        Traceback (most recent call last):
        ...
        ValueError: lattice rank=1 and ambient_dim=3 are incompatible

    We also get an error if neither the ambient dimension nor lattice
    are specified::

        sage: cones.rearrangement(3)
        Traceback (most recent call last):
        ...
        ValueError: either the ambient dimension or the lattice must
        be specified
    """
    from sage.geometry.cone import Cone
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ

    (ambient_dim, lattice) = _preprocess_args(ambient_dim, lattice)

    if p < 1 or p > ambient_dim or p not in ZZ:
        raise ValueError("order p=%s should be an integer between 1 "
                         "and ambient_dim=%d, inclusive" % (p, ambient_dim))

    I = matrix.identity(ZZ, ambient_dim)
    M = matrix.ones(ZZ, ambient_dim) - p * I
    G = matrix.identity(ZZ, ambient_dim).rows() + M.rows()
    return Cone(G, lattice=lattice)


def schur(ambient_dim=None, lattice=None):
    r"""
    The Schur cone in ``ambient_dim`` dimensions, or living
    in ``lattice``.

    The Schur cone in `n` dimensions induces the majorization ordering
    on the ambient space. If `\left\{e_{1}, e_{2}, \ldots,
    e_{n}\right\}` is the standard basis for the space, then its
    generators are `\left\{e_{i} - e_{i+1}\ |\ 1 \le i \le
    n-1\right\}`. Its dual is the downward monotonic cone.

    INPUT:

    - ``ambient_dim`` -- a nonnegative integer (default: ``None``); the
      dimension of the ambient space

    - ``lattice`` -- a toric lattice (default: ``None``); the lattice in
      which the cone will live

    If ``ambient_dim`` is omitted, then it will be inferred from the
    rank of ``lattice``. If the ``lattice`` is omitted, then the
    default lattice of rank ``ambient_dim`` will be used.

    A ``ValueError`` is raised if neither ``ambient_dim`` nor
    ``lattice`` are specified. It is also a ``ValueError`` to specify
    both ``ambient_dim`` and ``lattice`` unless the rank of
    ``lattice`` is equal to ``ambient_dim``.

    OUTPUT:

    A :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone` representing the Schur
    cone living in ``lattice``, with ambient dimension ``ambient_dim``.
    Each generating ray has the integer ring as its base ring.

    A ``ValueError`` can be raised if the inputs are incompatible or
    insufficient. See the INPUT documentation for details.

    REFERENCES:

    - [GS2010]_, Section 3.1

    - [IS2005]_, Example 7.3

    - [SS2016]_, Example 7.4

    EXAMPLES:

    Verify the claim [SS2016]_ that the maximal angle between any two
    generators of the Schur cone and the nonnegative orthant in
    dimension five is `\left(3/4\right)\pi`::

        sage: # needs sage.rings.number_fields
        sage: P = cones.schur(5)
        sage: Q = cones.nonnegative_orthant(5)
        sage: G = ( g.change_ring(QQbar).normalized() for g in P )
        sage: H = ( h.change_ring(QQbar).normalized() for h in Q )
        sage: actual = max(arccos(u.inner_product(v)) for u in G for v in H)
        sage: expected = 3*pi/4
        sage: abs(actual - expected).n() < 1e-12
        True

    The dual of the Schur cone is the "downward monotonic cone"
    [GS2010]_, whose elements' entries are in non-increasing order::

        sage: ambient_dim = ZZ.random_element(10)
        sage: K = cones.schur(ambient_dim).dual()
        sage: x = K.random_element()
        sage: all( x[i] >= x[i+1] for i in range(ambient_dim-1) )
        True

    TESTS:

    We get the trivial cone when ``ambient_dim`` is zero::

        sage: cones.schur(0).is_trivial()
        True

    The Schur cone induces the majorization ordering, as in Iusem
    and Seeger's [IS2005]_ Example 7.3::

        sage: def majorized_by(x,y):
        ....:     return (all(sum(x[0:i]) <= sum(y[0:i])
        ....:                 for i in range(x.degree()-1))
        ....:             and sum(x) == sum(y))
        sage: ambient_dim = ZZ.random_element(10)
        sage: V = VectorSpace(QQ, ambient_dim)
        sage: S = cones.schur(ambient_dim)
        sage: majorized_by(V.zero(), S.random_element())
        True
        sage: x = V.random_element()
        sage: y = V.random_element()
        sage: majorized_by(x,y) == ( (y-x) in S )
        True

    If a ``lattice`` was given, it is actually used::

        sage: L = ToricLattice(3, 'M')
        sage: cones.schur(3, lattice=L)
        2-d cone in 3-d lattice M

    Unless the rank of the lattice disagrees with ``ambient_dim``::

        sage: L = ToricLattice(1, 'M')
        sage: cones.schur(3, lattice=L)
        Traceback (most recent call last):
        ...
        ValueError: lattice rank=1 and ambient_dim=3 are incompatible

    We also get an error if no arguments are given::

        sage: cones.schur()
        Traceback (most recent call last):
        ...
        ValueError: either the ambient dimension or the lattice must
        be specified
    """
    from sage.geometry.cone import Cone
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ

    (ambient_dim, lattice) = _preprocess_args(ambient_dim, lattice)

    def _f(i,j):
        if i == j:
            return 1
        elif j - i == 1:
            return -1
        else:
            return 0

    # The "max" below catches the trivial case where ambient_dim == 0.
    S = matrix(ZZ, max(0, ambient_dim-1), ambient_dim, _f)

    return Cone(S.rows(), lattice)


def trivial(ambient_dim=None, lattice=None):
    r"""
    The trivial cone with no nonzero generators in ``ambient_dim``
    dimensions, or living in ``lattice``.

    INPUT:

    - ``ambient_dim`` -- a nonnegative integer (default: ``None``); the
      dimension of the ambient space

    - ``lattice`` -- a toric lattice (default: ``None``); the lattice in
      which the cone will live

    If ``ambient_dim`` is omitted, then it will be inferred from the
    rank of ``lattice``. If the ``lattice`` is omitted, then the
    default lattice of rank ``ambient_dim`` will be used.

    A ``ValueError`` is raised if neither ``ambient_dim`` nor
    ``lattice`` are specified. It is also a ``ValueError`` to specify
    both ``ambient_dim`` and ``lattice`` unless the rank of
    ``lattice`` is equal to ``ambient_dim``.

    OUTPUT:

    A :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone` representing the
    trivial cone with no nonzero generators living in ``lattice``,
    with ambient dimension ``ambient_dim``.

    A ``ValueError`` can be raised if the inputs are incompatible or
    insufficient. See the INPUT documentation for details.

    EXAMPLES:

    Construct the trivial cone, containing only the origin, in three
    dimensions::

        sage: cones.trivial(3)
        0-d cone in 3-d lattice N

    If a ``lattice`` is given, the trivial cone will live in that
    lattice::

        sage: L = ToricLattice(3, 'M')
        sage: cones.trivial(3, lattice=L)
        0-d cone in 3-d lattice M

    TESTS:

    We can construct the trivial cone in a trivial ambient space::

        sage: cones.trivial(0)
        0-d cone in 0-d lattice N

    An error is raised if the rank of the lattice disagrees with
    ``ambient_dim``::

        sage: L = ToricLattice(1, 'M')
        sage: cones.trivial(3, lattice=L)
        Traceback (most recent call last):
        ...
        ValueError: lattice rank=1 and ambient_dim=3 are incompatible

    We also get an error if no arguments are given::

        sage: cones.trivial()
        Traceback (most recent call last):
        ...
        ValueError: either the ambient dimension or the lattice must
        be specified
    """
    from sage.geometry.cone import Cone

    (ambient_dim, lattice) = _preprocess_args(ambient_dim, lattice)

    return Cone([], lattice)
