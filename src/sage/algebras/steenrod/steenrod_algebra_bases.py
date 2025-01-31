"""
Steenrod algebra bases

AUTHORS:

- John H. Palmieri (2008-07-30): version 0.9
- John H. Palmieri (2010-06-30): version 1.0
- Simon King (2011-10-25): Fix the use of cached functions

This package defines functions for computing various bases of the
Steenrod algebra, and for converting between the Milnor basis and
any other basis.

This packages implements a number of different bases, at least at
the prime 2. The Milnor and Serre-Cartan bases are the most
familiar and most standard ones, and all of the others are defined
in terms of one of these. The bases are described in the
documentation for the function
:func:`steenrod_algebra_basis`; also see the papers by
Monks [Mon1998]_ and Wood [Woo1998]_ for more information about them. For
commutator bases, see the preprint by Palmieri and Zhang [PZ2008]_.

- ``'milnor'`` -- Milnor basis

- ``'serre-cartan'`` or ``'adem'`` or ``'admissible'`` -- Serre-Cartan basis

Most of the rest of the bases are only defined when `p=2`.  The only
exceptions are the `P^s_t`-bases and the commutator bases, which are
defined at all primes.

- ``'wood_y'`` -- Wood's Y basis

- ``'wood_z'`` -- Wood's Z basis

- ``'wall'``, ``'wall_long'`` -- Wall's basis

- ``'arnon_a'``, ``'arnon_a_long'`` -- Arnon's A basis

- ``'arnon_c'`` -- Arnon's C basis

- ``'pst'``, ``'pst_rlex'``, ``'pst_llex'``, ``'pst_deg'``, ``'pst_revz'`` --
  various `P^s_t`-bases

- ``'comm'``, ``'comm_rlex'``, ``'comm_llex'``, ``'comm_deg'``, ``'comm_revz'``,
   or these with ``'_long'`` appended -- various commutator bases

The main functions provided here are

- :func:`steenrod_algebra_basis`.  This computes a tuple representing
  basis elements for the Steenrod algebra in a given degree, at a
  given prime, with respect to a given basis. It is a cached function.

- :func:`convert_to_milnor_matrix`.  This returns the change-of-basis
  matrix, in a given degree, from any basis to the Milnor basis. It is
  a cached function.

- :func:`convert_from_milnor_matrix`.  This returns the inverse of the
  previous matrix.

INTERNAL DOCUMENTATION:

If you want to implement a new basis for the Steenrod algebra:

In the file :file:`steenrod_algebra.py`:

For the class :class:`SteenrodAlgebra_generic
<sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic>`, add functionality to the
methods:

- :meth:`_repr_term <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic._repr_term>`

- :meth:`degree_on_basis <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic.degree_on_basis>`

- :meth:`_milnor_on_basis <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic._milnor_on_basis>`

- :meth:`an_element <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic.an_element>`

In the file :file:`steenrod_algebra_misc.py`:

- add functionality to :func:`get_basis_name
  <sage.algebras.steenrod.steenrod_algebra_misc.get_basis_name>`: this
  should accept as input various synonyms for the basis, and its
  output should be a canonical name for the basis.

- add a function ``BASIS_mono_to_string`` like
  :func:`milnor_mono_to_string
  <sage.algebras.steenrod.steenrod_algebra_misc.milnor_mono_to_string>`
  or one of the other similar functions.

In this file :file:`steenrod_algebra_bases.py`:

- add appropriate lines to :func:`steenrod_algebra_basis`.

- add a function to compute the basis in a given dimension (to be
  called by :func:`steenrod_algebra_basis`).

- modify :func:`steenrod_basis_error_check` so it checks the new
  basis.

If the basis has an intrinsic way of defining a product, implement it
in the file :file:`steenrod_algebra_mult.py` and also in the
:meth:`product_on_basis
<sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic.product_on_basis>`
method for :class:`SteenrodAlgebra_generic
<sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic>` in
:file:`steenrod_algebra.py`.
"""

# ****************************************************************************
#  Copyright (C) 2008-2010 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
# ****************************************************************************

from sage.misc.cachefunc import cached_function


@cached_function
def convert_to_milnor_matrix(n, basis, p=2, generic='auto'):
    r"""
    Change-of-basis matrix, 'basis' to Milnor, in dimension
    `n`, at the prime `p`.

    INPUT:

    - ``n`` -- nonnegative integer, the dimension
    - ``basis`` -- string, the basis from which to convert
    - ``p`` -- positive prime number (default: 2)

    OUTPUT:

    ``matrix`` -- change-of-basis matrix, a square matrix over `\GF{p}`

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.algebras.steenrod.steenrod_algebra_bases import convert_to_milnor_matrix
        sage: convert_to_milnor_matrix(5, 'adem')  # indirect doctest
        [0 1]
        [1 1]
        sage: convert_to_milnor_matrix(45, 'milnor')
        111 x 111 dense matrix over Finite Field of size 2 (use the '.str()' method to see the entries)
        sage: convert_to_milnor_matrix(12, 'wall')
        [1 0 0 1 0 0 0]
        [1 1 0 0 0 1 0]
        [0 1 0 1 0 0 0]
        [0 0 0 1 0 0 0]
        [1 1 0 0 1 0 0]
        [0 0 1 1 1 0 1]
        [0 0 0 0 1 0 1]

    The function takes an optional argument, the prime `p` over
    which to work::

        sage: # needs sage.modules
        sage: convert_to_milnor_matrix(17, 'adem', 3)
        [0 0 1 1]
        [0 0 0 1]
        [1 1 1 1]
        [0 1 0 1]
        sage: convert_to_milnor_matrix(48, 'adem', 5)
        [0 1]
        [1 1]
        sage: convert_to_milnor_matrix(36, 'adem', 3)
        [0 0 1]
        [0 1 0]
        [1 2 0]
    """
    from sage.matrix.constructor import matrix
    from sage.rings.finite_rings.finite_field_constructor import GF
    from .steenrod_algebra import SteenrodAlgebra
    if generic == 'auto':
        generic = p != 2
    if n == 0:
        return matrix(GF(p), 1, 1, [[1]])
    milnor_base = steenrod_algebra_basis(n, 'milnor', p, generic=generic)
    rows = []
    A = SteenrodAlgebra(basis=basis, p=p, generic=generic)
    for poly in A.basis(n):
        d = poly.milnor().monomial_coefficients()
        for v in milnor_base:
            entry = d.get(v, 0)
            rows = rows + [entry]
    d = len(milnor_base)
    return matrix(GF(p), d, d, rows)


def convert_from_milnor_matrix(n, basis, p=2, generic='auto'):
    r"""
    Change-of-basis matrix, Milnor to ``basis``, in dimension `n`.

    INPUT:

    - ``n`` -- nonnegative integer, the dimension

    - ``basis`` -- string, the basis to which to convert

    - ``p`` -- positive prime number (default: 2)

    OUTPUT:

    ``matrix`` -- change-of-basis matrix, a square matrix over `\GF{p}`

    .. NOTE::

        This is called internally.  It is not intended for casual
        users, so no error checking is made on the integer `n`, the
        basis name, or the prime.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.algebras.steenrod.steenrod_algebra_bases import convert_from_milnor_matrix, convert_to_milnor_matrix
        sage: convert_from_milnor_matrix(12, 'wall')
        [1 0 0 1 0 0 0]
        [0 0 1 1 0 0 0]
        [0 0 0 1 0 1 1]
        [0 0 0 1 0 0 0]
        [1 0 1 0 1 0 0]
        [1 1 1 0 0 0 0]
        [1 0 1 0 1 0 1]
        sage: convert_from_milnor_matrix(38, 'serre_cartan')
        72 x 72 dense matrix over Finite Field of size 2 (use the '.str()' method to see the entries)
        sage: x = convert_to_milnor_matrix(20, 'wood_y')
        sage: y = convert_from_milnor_matrix(20, 'wood_y')
        sage: x*y
        [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]

    The function takes an optional argument, the prime `p` over
    which to work::

        sage: convert_from_milnor_matrix(17, 'adem', 3)                                 # needs sage.modules
        [2 1 1 2]
        [0 2 0 1]
        [1 2 0 0]
        [0 1 0 0]
    """
    mat = convert_to_milnor_matrix(n, basis, p, generic)
    if mat.nrows():
        return convert_to_milnor_matrix(n, basis, p, generic).inverse()
    return mat


@cached_function
def steenrod_algebra_basis(n, basis='milnor', p=2, **kwds):
    r"""
    Basis for the Steenrod algebra in degree `n`.

    INPUT:

    - ``n`` -- nonnegative integer
    - ``basis`` -- string, which basis to use (default: ``'milnor'``)
    - ``p`` -- positive prime number (default: 2)
    - ``profile`` -- profile function (default: ``None``); this
      is just passed on to the functions :func:`milnor_basis` and
      :func:`pst_basis`
    - ``truncation_type`` -- truncation type, either 0 or ``Infinity``
      (default: ``Infinity`` if no profile function is specified,
      0 otherwise).  This is just passed on to the function
      :func:`milnor_basis`.
    - ``generic`` -- boolean (default: ``None``)

    OUTPUT:

    Tuple of objects representing basis elements for the Steenrod algebra
    in dimension n.

    The choices for the string ``basis`` are as follows; see the
    documentation for :mod:`sage.algebras.steenrod.steenrod_algebra`
    for details on each basis:

    - ``'milnor'`` -- Milnor basis
    - ``'serre-cartan'`` or ``'adem'`` or ``'admissible'`` -- Serre-Cartan basis
    - ``'pst'``, ``'pst_rlex'``, ``'pst_llex'``, ``'pst_deg'``, ``'pst_revz'`` --
      various `P^s_t`-bases
    - ``'comm'``, ``'comm_rlex'``, ``'comm_llex'``, ``'comm_deg'``, ``'comm_revz'``, or
      any of these with ``'_long'`` appended -- various commutator bases

    The rest of these bases are only defined when `p=2`.

    - ``'wood_y'`` -- Wood's Y basis
    - ``'wood_z'`` -- Wood's Z basis
    - ``'wall'`` or ``'wall_long'`` -- Wall's basis
    - ``'arnon_a'`` or ``'arnon_a_long'`` -- Arnon's A basis
    - ``'arnon_c'`` -- Arnon's C basis

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import steenrod_algebra_basis
        sage: steenrod_algebra_basis(7, 'milnor') # indirect doctest
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: steenrod_algebra_basis(5)   # milnor basis is the default
        ((2, 1), (5,))

    Bases in negative dimensions are empty::

        sage: steenrod_algebra_basis(-2, 'wall')
        ()

    The third (optional) argument to 'steenrod_algebra_basis' is the
    prime p::

        sage: steenrod_algebra_basis(9, 'milnor', p=3)
        (((1,), (1,)), ((0,), (2,)))
        sage: steenrod_algebra_basis(9, 'milnor', 3)
        (((1,), (1,)), ((0,), (2,)))
        sage: steenrod_algebra_basis(17, 'milnor', 3)
        (((2,), ()), ((1,), (3,)), ((0,), (0, 1)), ((0,), (4,)))

    Other bases::

        sage: steenrod_algebra_basis(7, 'admissible')
        ((7,), (6, 1), (4, 2, 1), (5, 2))
        sage: steenrod_algebra_basis(13, 'admissible', p=3)
        ((1, 3, 0), (0, 3, 1))
        sage: steenrod_algebra_basis(5, 'wall')
        (((2, 2), (0, 0)), ((1, 1), (1, 0)))
        sage: steenrod_algebra_basis(5, 'wall_long')
        (((2, 2), (0, 0)), ((1, 1), (1, 0)))
        sage: steenrod_algebra_basis(5, 'pst-rlex')
        (((0, 1), (2, 1)), ((1, 1), (0, 2)))
    """
    from .steenrod_algebra_misc import get_basis_name
    try:
        if n < 0 or int(n) != n:
            return ()
    except TypeError:
        return ()

    generic = kwds.get("generic", p != 2)

    basis_name = get_basis_name(basis, p, generic=generic)
    if basis_name.find('long') >= 0:
        basis_name = basis_name.rsplit('_', 1)[0]

    profile = kwds.get("profile", None)
    if (profile is not None and profile != () and profile != ((), ())
            and basis != 'milnor' and basis.find('pst') == -1):
        raise ValueError("Profile functions may only be used with the Milnor or pst bases")

    # Milnor basis
    if basis_name == 'milnor':
        return milnor_basis(n, p, **kwds)
    # Serre-Cartan basis
    elif basis_name == 'serre-cartan':
        return serre_cartan_basis(n, p, **kwds)
    # Atomic bases, p odd:
    elif generic and (basis_name.find('pst') >= 0
                    or basis_name.find('comm') >= 0):
        return atomic_basis_odd(n, basis_name, p, **kwds)
    # Atomic bases, p=2
    elif not generic and (basis_name == 'woody' or basis_name == 'woodz'
                     or basis_name == 'wall' or basis_name == 'arnona'
                     or basis_name.find('pst') >= 0
                     or basis_name.find('comm') >= 0):
        return atomic_basis(n, basis_name, **kwds)
    # Arnon 'C' basis
    elif not generic and basis == 'arnonc':
        return arnonC_basis(n)
    else:
        raise ValueError("Unknown basis: %s at the prime %s" % (basis, p))


# helper functions for producing bases

def restricted_partitions(n, l, no_repeats=False):
    """
    Iterator over 'restricted' partitions of `n`: partitions with parts taken
    from list `l`.

    INPUT:

    - ``n`` -- nonnegative integer
    - ``l`` -- list of positive integers
    - ``no_repeats`` -- boolean (default: ``False``); if ``True``,
      only return partitions with no repeated parts

    OUTPUT: iterator of lists

    One could also use ``Partitions(n, parts_in=l)``, but this
    function may be faster.  Also, while ``Partitions(n, parts_in=l,
    max_slope=-1)`` should in theory return the partitions of `n` with
    parts in ``l`` with no repetitions, the ``max_slope=-1`` argument
    is ignored, so it does not work.  (At the moment, the
    ``no_repeats=True`` case is the only one used in the code.)

    .. TODO::

        This should be re-implemented in a non-recursive way.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import restricted_partitions
        sage: list(restricted_partitions(10, [7,5,1]))
        [[7, 1, 1, 1], [5, 5], [5, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        sage: list(restricted_partitions(10, [6,5,4,3,2,1], no_repeats=True))
        [[6, 4], [6, 3, 1], [5, 4, 1], [5, 3, 2], [4, 3, 2, 1]]
        sage: list(restricted_partitions(10, [6,4,2]))
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
        sage: list(restricted_partitions(10, [6,4,2], no_repeats=True))
        [[6, 4]]

    ``l`` may have repeated elements. If ``no_repeats`` is ``False``, this
    has no effect. If ``no_repeats`` is ``True``, and if the repeated
    elements appear consecutively in ``l``, then each element may be
    used only as many times as it appears in ``l``::

        sage: list(restricted_partitions(10, [6,4,2,2], no_repeats=True))
        [[6, 4], [6, 2, 2]]
        sage: list(restricted_partitions(10, [6,4,2,2,2], no_repeats=True))
        [[6, 4], [6, 2, 2], [4, 2, 2, 2]]

    (If the repeated elements don't appear consecutively, the results
    are likely meaningless, containing several partitions more than
    once, for example.)

    In the following examples, ``no_repeats`` is ``False``::

        sage: list(restricted_partitions(10, [6,4,2]))
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
        sage: list(restricted_partitions(10, [6,4,2,2,2]))
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
        sage: list(restricted_partitions(10, [6,4,4,4,2,2,2,2,2,2]))
        [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    """
    if n < 0:
        return

    if n == 0:
        yield []
        return

    index = 1 if no_repeats else 0
    old_i = 0
    for i in l:
        if old_i != i:
            for sigma in restricted_partitions(n - i, l[index:], no_repeats):
                yield [i] + sigma
        index += 1
        old_i = i


def xi_degrees(n, p=2, reverse=True):
    r"""
    Decreasing list of degrees of the `\xi_i`'s, starting in degree `n`.

    INPUT:

    - ``n`` -- integer
    - ``p`` -- prime number (default: 2)
    - ``reverse`` -- boolean (default: ``True``)

    OUTPUT: list of integers

    When `p=2`: decreasing list of the degrees of the `\xi_i`'s with
    degree at most `n`.

    At odd primes: decreasing list of these degrees, each divided by
    `2(p-1)`.

    If ``reverse`` is ``False``, then return an increasing list rather
    than a decreasing one.

    EXAMPLES::

        sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17)
        [15, 7, 3, 1]
        sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17, reverse=False)
        [1, 3, 7, 15]
        sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17, p=3)
        [13, 4, 1]
        sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(400, p=17)
        [307, 18, 1]
    """
    from sage.rings.integer import Integer
    if n <= 0:
        return []
    N = Integer(n*(p-1) + 1)
    l = [(p**d-1)//(p-1) for d in range(1, N.exact_log(p)+1)]
    if reverse:
        l.reverse()
    return l


########################################################
# Functions for defining bases.

# These should each return a tuple of tuples of the appropriate form
# for the basis.  For example, at the prime 2, the Milnor basis
# element Sq(a,b,c,...) corresponds to the tuple (a, b, c, ...), while
# at odd primes, the element Q_i Q_j ... P(a, b, ...) corresponds to
# the pair ((i, j, ...), (a, b, ...)).  See each function for more
# information.

def milnor_basis(n, p=2, **kwds):
    r"""
    Milnor basis in dimension `n` with profile function ``profile``.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``p`` -- positive prime number (default: 2)

    - ``profile`` -- profile function (default: ``None``).
      Together with ``truncation_type``, specify the profile function
      to be used; ``None`` means the profile function for the entire
      Steenrod algebra.  See
      :mod:`sage.algebras.steenrod.steenrod_algebra` and
      :func:`SteenrodAlgebra <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
      for information on profile functions.

    - ``truncation_type`` -- truncation type, either 0 or ``Infinity``
      (default: ``Infinity`` if no profile function is specified,
      0 otherwise)

    OUTPUT: tuple of mod `p` Milnor basis elements in dimension `n`

    At the prime 2, the Milnor basis consists of symbols of the form
    `\text{Sq}(m_1, m_2, ..., m_t)`, where each
    `m_i` is a nonnegative integer and if `t>1`, then
    `m_t \neq 0`. At odd primes, it consists of symbols of the
    form `Q_{e_1} Q_{e_2} ... P(m_1, m_2, ..., m_t)`,
    where `0 \leq e_1 < e_2 < ...`, each `m_i` is a
    nonnegative integer, and if `t>1`, then
    `m_t \neq 0`.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import milnor_basis
        sage: milnor_basis(7)
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: milnor_basis(7, 2)
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: milnor_basis(4, 2)
        ((1, 1), (4,))
        sage: milnor_basis(4, 2, profile=[2,1])
        ((1, 1),)
        sage: milnor_basis(4, 2, profile=(), truncation_type=0)
        ()
        sage: milnor_basis(4, 2, profile=(), truncation_type=Infinity)
        ((1, 1), (4,))
        sage: milnor_basis(9, 3)
        (((1,), (1,)), ((0,), (2,)))
        sage: milnor_basis(17, 3)
        (((2,), ()), ((1,), (3,)), ((0,), (0, 1)), ((0,), (4,)))
        sage: milnor_basis(48, p=5)
        (((), (0, 1)), ((), (6,)))
        sage: len(milnor_basis(100,3))
        13
        sage: len(milnor_basis(200,7))
        0
        sage: len(milnor_basis(240,7))
        3
        sage: len(milnor_basis(240,7, profile=((),()), truncation_type=Infinity))
        3
        sage: len(milnor_basis(240,7, profile=((),()), truncation_type=0))
        0
    """
    generic = kwds.get('generic', p != 2)

    if n == 0:
        if not generic:
            return ((),)
        else:
            return (((), ()),)

    from sage.rings.infinity import Infinity
    from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
    profile = kwds.get("profile", None)
    trunc = kwds.get("truncation_type", None)
    if trunc is None:
        if profile is not None:
            trunc = 0
        else:
            trunc = Infinity

    result = []
    if not generic:
        for mono in WeightedIntegerVectors(n, xi_degrees(n, reverse=False)):
            exponents = list(mono)
            while exponents and exponents[-1] == 0:
                exponents.pop(-1)
            # check profile:
            okay = True
            if profile is not None and len(profile) > 0:
                for i in range(len(exponents)):
                    if ((len(profile) > i and exponents[i] >= 2**profile[i])
                        or (len(profile) <= i and trunc < Infinity
                            and exponents[i] >= 2**trunc)):
                        okay = False
                        break
            else:
                # profile is empty
                okay = (trunc == Infinity)
            if okay:
                result.append(tuple(exponents))
    else:  # p odd
        # first find the P part of each basis element.
        # in this part of the code (the P part), all dimensions are
        # divided by 2(p-1).
        for dim in range(n//(2*(p-1)) + 1):
            if dim == 0:
                P_result = [[0]]
            else:
                P_result = []
            for mono in WeightedIntegerVectors(dim, xi_degrees(dim, p=p, reverse=False)):
                p_mono = list(mono)
                while p_mono and p_mono[-1] == 0:
                    p_mono.pop(-1)
                if p_mono:
                    P_result.append(p_mono)
            # now find the Q part of the basis element.
            # dimensions here are back to normal.
            for p_mono in P_result:
                deg = n - 2*dim*(p-1)
                q_degrees = [1+2*(p-1)*d for d in
                             xi_degrees(int((deg - 1)//(2*(p-1))), p)] + [1]
                q_degrees_decrease = q_degrees
                q_degrees.reverse()
                if deg % (2*(p-1)) <= len(q_degrees):
                    # if this inequality fails, no way to have a partition
                    # with distinct parts.
                    for sigma in restricted_partitions(deg,
                                                       q_degrees_decrease,
                                                       no_repeats=True):
                        index = 0
                        q_mono = []
                        for q in q_degrees:
                            if q in sigma:
                                q_mono.append(index)
                            index += 1
                        # check profile:
                        okay = True
                        if profile is not None and (len(profile[0]) > 0
                                                    or len(profile[1]) > 0):
                            # check profile function for q_mono
                            for i in q_mono:
                                if ((len(profile[1]) > i and profile[1][i] == 1)
                                    or (len(profile[1]) <= i and trunc == 0)):
                                    okay = False
                                    break
                            # check profile function for p_mono
                            for i in range(len(p_mono)):
                                if okay and ((len(profile[0]) > i and p_mono[i] >= p**profile[0][i])
                                             or (len(profile[0]) <= i and trunc < Infinity
                                                 and p_mono[i] >= p**trunc)):
                                    okay = False
                                    break
                        else:
                            # profile is empty
                            okay = (trunc == Infinity)
                        if okay:
                            if list(p_mono) == [0]:
                                p_mono = []
                            result.append((tuple(q_mono), tuple(p_mono)))
    return tuple(result)


def serre_cartan_basis(n, p=2, bound=1, **kwds):
    r"""
    Serre-Cartan basis in dimension `n`.

    INPUT:

    - ``n`` -- nonnegative integer
    - ``bound`` -- positive integer (optional)
    - ``prime`` -- positive prime number (default: 2)

    OUTPUT: tuple of mod `p` Serre-Cartan basis elements in dimension `n`

    The Serre-Cartan basis consists of 'admissible monomials in the
    Steenrod squares'. Thus at the prime 2, it consists of monomials
    `\text{Sq}^{m_1} \text{Sq}^{m_2} ... \text{Sq}^{m_t}` with `m_i
    \geq 2m_{i+1}` for each `i`. At odd primes, it consists of
    monomials `\beta^{e_0} P^{s_1} \beta^{e_1} P^{s_2} ...  P^{s_k}
    \beta^{e_k}` with each `e_i` either 0 or 1, `s_i \geq p s_{i+1} +
    e_i` for all `i`, and `s_k \geq 1`.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import serre_cartan_basis
        sage: serre_cartan_basis(7)
        ((7,), (6, 1), (4, 2, 1), (5, 2))
        sage: serre_cartan_basis(13,3)
        ((1, 3, 0), (0, 3, 1))
        sage: serre_cartan_basis(50,5)
        ((1, 5, 0, 1, 1), (1, 6, 1))

    If optional argument ``bound`` is present, include only those monomials
    whose last term is at least ``bound`` (when `p=2`), or those for which
    `s_k - e_k \geq bound` (when `p` is odd). ::

        sage: serre_cartan_basis(7, bound=2)
        ((7,), (5, 2))
        sage: serre_cartan_basis(13, 3, bound=3)
        ((1, 3, 0),)
    """
    generic = kwds.get('generic', p != 2)

    if n == 0:
        return ((),)
    else:
        if not generic:
            # Build basis recursively.  last = last term.
            # last is >= bound, and we will append (last,) to the end of
            # elements from serre_cartan_basis (n - last, bound=2 * last).
            # This means that 2 last <= n - last, or 3 last <= n.
            result = [(n,)]
            for last in range(bound, 1+n//3):
                for vec in serre_cartan_basis(n - last, bound=2 * last):
                    new = vec + (last,)
                    result.append(new)
        else: # p odd
            if n % (2 * (p-1)) == 0 and n//(2 * (p-1)) >= bound:
                result = [(0, int(n//(2 * (p-1))), 0)]
            elif n == 1:
                result = [(1,)]
            else:
                result = []
            # 2 cases: append P^{last}, or append P^{last} beta
            # case 1: append P^{last}
            for last in range(bound, 1+n//(2*(p - 1))):
                if n - 2*(p-1)*last > 0:
                    for vec in serre_cartan_basis(n - 2*(p-1)*last,
                                                  p, p*last, generic=generic):
                        result.append(vec + (last,0))
            # case 2: append P^{last} beta
            if bound == 1:
                bound = 0
            for last in range(bound+1, 1+n//(2*(p - 1))):
                basis = serre_cartan_basis(n - 2*(p-1)*last - 1,
                                           p, p*last, generic=generic)
                for vec in basis:
                    if vec == ():
                        vec = (0,)
                    new = vec + (last, 1)
                    result.append(new)
    return tuple(result)


def atomic_basis(n, basis, **kwds):
    r"""
    Basis for dimension `n` made of elements in 'atomic' degrees:
    degrees of the form `2^i (2^j - 1)`.

    This works at the prime 2 only.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``basis`` -- string, the name of the basis

    - ``profile`` -- profile function (default: ``None``).
      Together with ``truncation_type``, specify the profile function
      to be used; ``None`` means the profile function for the entire
      Steenrod algebra.  See
      :mod:`sage.algebras.steenrod.steenrod_algebra` and
      :func:`SteenrodAlgebra` for information on profile functions.

    - ``truncation_type`` -- truncation type, either 0 or ``Infinity``
      (default: ``Infinity`` if no profile function is specified,
      0 otherwise).

    OUTPUT: tuple of basis elements in dimension `n`

    The atomic bases include Wood's Y and Z bases, Wall's basis,
    Arnon's A basis, the `P^s_t`-bases, and the commutator
    bases. (All of these bases are constructed similarly, hence their
    constructions have been consolidated into a single function. Also,
    see the documentation for :func:`steenrod_algebra_basis` for
    descriptions of them.)  For `P^s_t`-bases, you may also specify a
    profile function and truncation type; profile functions are ignored
    for the other bases.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import atomic_basis
        sage: atomic_basis(6,'woody')
        (((1, 0), (0, 1), (0, 0)), ((2, 0), (1, 0)), ((1, 1),))
        sage: atomic_basis(8,'woodz')
        (((2, 0), (0, 1), (0, 0)), ((0, 2), (0, 0)), ((1, 1), (1, 0)), ((3, 0),))
        sage: atomic_basis(6,'woodz') == atomic_basis(6, 'woody')
        True
        sage: atomic_basis(9,'woodz') == atomic_basis(9, 'woody')
        False

    Wall's basis::

        sage: atomic_basis(8,'wall')
        (((2, 2), (1, 0), (0, 0)), ((2, 0), (0, 0)), ((2, 1), (1, 1)), ((3, 3),))

    Arnon's A basis::

        sage: atomic_basis(7,'arnona')
        (((0, 0), (1, 1), (2, 2)), ((0, 0), (2, 1)), ((1, 0), (2, 2)), ((2, 0),))

    `P^s_t`-bases::

        sage: atomic_basis(7,'pst_rlex')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((2, 1), (0, 2)), ((0, 3),))
        sage: atomic_basis(7,'pst_llex')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((0, 2), (2, 1)), ((0, 3),))
        sage: atomic_basis(7,'pst_deg')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((0, 2), (2, 1)), ((0, 3),))
        sage: atomic_basis(7,'pst_revz')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((0, 2), (2, 1)), ((0, 3),))

    Commutator bases::

        sage: atomic_basis(7,'comm_rlex')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((2, 1), (0, 2)), ((0, 3),))
        sage: atomic_basis(7,'comm_llex')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((0, 2), (2, 1)), ((0, 3),))
        sage: atomic_basis(7,'comm_deg')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((0, 2), (2, 1)), ((0, 3),))
        sage: atomic_basis(7,'comm_revz')
        (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 2)), ((0, 2), (2, 1)), ((0, 3),))
    """
    def degree_dictionary(n, basis):
        """
        Dictionary of atomic degrees for basis up to degree `n`.

        The keys for the dictionary are the atomic degrees -- the numbers of
        the form `2^i (2^j - 1)` -- which are less than or equal to `n`. The value
        associated to such a degree depends on basis; it has the form
        `(s,t)`, where `(s,t)` is a pair of integers which indexes the
        corresponding element.
        """
        dict = {}
        if basis.find('wood') >= 0:
            k = 0
            m = 0
            deg = 2**m * (2**(k+1) - 1)
            while deg <= n:
                dict[deg] = (m,k)
                if m > 0:
                    m = m - 1
                    k = k + 1
                else:
                    m = k + 1
                    k = 0
                deg = 2**m * (2**(k+1) - 1)
        elif basis.find('wall') >= 0 or basis.find('arnon') >= 0:
            k = 0
            m = 0
            deg = 2**k * (2**(m-k+1) - 1)
            while deg <= n:
                dict[deg] = (m,k)
                if k == 0:
                    m = m + 1
                    k = m
                else:
                    k = k - 1
                deg = 2**k * (2**(m-k+1) - 1)
        elif basis.find('pst') >= 0 or basis.find('comm') >= 0:
            s = 0
            t = 1
            deg = 2**s * (2**t - 1)
            while deg <= n:
                if basis.find('pst') >= 0:
                    dict[deg] = (s,t)
                else:  # comm
                    dict[deg] = (s,t)
                if s == 0:
                    s = t
                    t = 1
                else:
                    s = s - 1
                    t = t + 1
                deg = 2**s * (2**t - 1)
        return dict

    def sorting_pair(s, t, basis):   # pair used for sorting the basis
        if basis.find('wood') >= 0 and basis.find('z') >= 0:
            return (-s-t,-s)
        elif basis.find('wood') >= 0 or basis.find('wall') >= 0 or \
                basis.find('arnon') >= 0:
            return (-s,-t)
        elif basis.find('rlex') >= 0:
            return (t,s)
        elif basis.find('llex') >= 0:
            return (s,t)
        elif basis.find('deg') >= 0:
            return (s+t,t)
        elif basis.find('revz') >= 0:
            return (s+t,s)

    from sage.rings.infinity import Infinity
    profile = kwds.get("profile", None)
    trunc = kwds.get("truncation_type", None)
    if profile is not None and trunc is None:
        trunc = 0

    if n == 0:
        return ((),)
    else:
        result = []
        degrees_etc = degree_dictionary(n, basis)
        degrees = list(degrees_etc)
        for sigma in restricted_partitions(n, degrees, no_repeats=True):
            big_list = [degrees_etc[part] for part in sigma]
            big_list.sort(key=lambda x: sorting_pair(x[0], x[1], basis))
            # reverse = True)
            # arnon: sort like wall, then reverse end result
            if basis.find('arnon') >= 0:
                big_list.reverse()

            # check profile:
            okay = True
            if basis.find('pst') >= 0:
                if profile is not None and len(profile) > 0:
                    for (s,t) in big_list:
                        if ((len(profile) > t-1 and profile[t-1] <= s)
                            or (len(profile) <= t-1 and trunc < Infinity)):
                            okay = False
                            break
            if okay:
                result.append(tuple(big_list))
        return tuple(result)


@cached_function
def arnonC_basis(n, bound=1):
    r"""
    Arnon's C basis in dimension `n`.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``bound`` -- positive integer (optional)

    OUTPUT: tuple of basis elements in dimension `n`

    The elements of Arnon's C basis are monomials of the form
    `\text{Sq}^{t_1} ... \text{Sq}^{t_m}` where for each
    `i`, we have `t_i \leq 2t_{i+1}` and
    `2^i | t_{m-i}`.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import arnonC_basis
        sage: arnonC_basis(7)
        ((7,), (2, 5), (4, 3), (4, 2, 1))

    If optional argument ``bound`` is present, include only those monomials
    whose first term is at least as large as ``bound``::

        sage: arnonC_basis(7,3)
        ((7,), (4, 3), (4, 2, 1))
    """
    if n == 0:
        return ((),)

    # Build basis recursively.  first = first term.
    # first is >= bound, and we will prepend (first,) to the
    # elements from arnonC_basis (n - first, first / 2).
    # first also must be divisible by 2**(len(old-basis-elt))
    # This means that 3 first <= 2 n.
    result = [(n,)]
    for first in range(bound, 1 + 2 * n // 3):
        tup = (first,)
        result.extend(tup + vec
                      for vec in arnonC_basis(n - first, max(first // 2, 1))
                      if not first % 2**len(vec))
    return tuple(result)


def atomic_basis_odd(n, basis, p, **kwds):
    r"""
    `P^s_t`-bases and commutator basis in dimension `n` at odd primes.

    This function is called ``atomic_basis_odd`` in analogy with
    :func:`atomic_basis`.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``basis`` -- string, the name of the basis

    - ``p`` -- positive prime number

    - ``profile`` -- profile function (default: ``None``).
      Together with ``truncation_type``, specify the profile function
      to be used; ``None`` means the profile function for the entire
      Steenrod algebra.  See
      :mod:`sage.algebras.steenrod.steenrod_algebra` and
      :func:`SteenrodAlgebra` for information on profile functions.

    - ``truncation_type`` -- truncation type, either 0 or ``Infinity``
      (default: ``Infinity`` if no profile function is specified,
      0 otherwise).

    OUTPUT: tuple of basis elements in dimension `n`

    The only possible difference in the implementations for `P^s_t`
    bases and commutator bases is that the former make sense, and
    require filtering, if there is a nontrivial profile function.
    This function is called by :func:`steenrod_algebra_basis`, and it
    will not be called for commutator bases if there is a profile
    function, so we treat the two bases exactly the same.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import atomic_basis_odd
        sage: atomic_basis_odd(8, 'pst_rlex', 3)
        (((), (((0, 1), 2),)),)

        sage: atomic_basis_odd(18, 'pst_rlex', 3)
        (((0, 2), ()), ((0, 1), (((1, 1), 1),)))
        sage: atomic_basis_odd(18, 'pst_rlex', 3, profile=((), (2,2,2)))
        (((0, 2), ()),)
    """
    def sorting_pair(s, t, basis):   # pair used for sorting the basis
        if basis.find('rlex') >= 0:
            return (t, s)
        elif basis.find('llex') >= 0:
            return (s, t)
        elif basis.find('deg') >= 0:
            return (s + t, t)
        elif basis.find('revz') >= 0:
            return (s + t, s)

    generic = kwds.get('generic', p != 2)
    if n == 0:
        if not generic:
            return ((),)
        else:
            return (((), ()),)

    from sage.rings.integer import Integer
    from sage.rings.infinity import Infinity
    from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
    profile = kwds.get("profile", None)
    trunc = kwds.get("truncation_type", 0)

    result = []
    for dim in range(n//(2*p-2) + 1):
        P_result = []
        for v in WeightedIntegerVectors(dim, xi_degrees(dim, p=p, reverse=False)):
            mono = []
            for t, a in enumerate(v):
                for s, pow in enumerate(Integer(a).digits(p)):
                    if pow > 0:
                        mono.append(((s, t+1), pow))
            P_result.append(mono)
        for p_mono in P_result:
            p_mono.sort(key=lambda x: sorting_pair(x[0][0], x[0][1], basis))
            deg = n - 2*dim*(p-1)
            q_degrees = [1+2*(p-1)*d for d in
                         xi_degrees((deg - 1)//(2*(p-1)), p)] + [1]
            q_degrees_decrease = q_degrees
            q_degrees.reverse()
            if deg % (2*(p-1)) <= len(q_degrees):
                # if this inequality fails, no way to have a partition
                # with distinct parts.
                for sigma in restricted_partitions(deg,
                                                   q_degrees_decrease,
                                                   no_repeats=True):
                    q_mono = [index for index, q in enumerate(q_degrees)
                              if q in sigma]
                    # check profile:
                    okay = True
                    if profile is not None and profile != ((), ()):
                        # check profile function for q_mono
                        for i in q_mono:
                            if ((len(profile[1]) > i and profile[1][i] == 1)
                                or (len(profile[1]) <= i and trunc == 0)):
                                okay = False
                                break

                        for ((s, t), _) in p_mono:
                            if ((len(profile[0]) > t-1 and profile[0][t-1] <= s)
                                or (len(profile[0]) <= t-1 and trunc < Infinity)):
                                okay = False
                                break

                    if okay:
                        if list(p_mono) == [0]:
                            p_mono = []
                        result.append((tuple(q_mono), tuple(p_mono)))
    return tuple(result)


#############################################################################
def steenrod_basis_error_check(dim, p, **kwds):
    """
    This performs crude error checking.

    INPUT:

    - ``dim`` -- nonnegative integer
    - ``p`` -- positive prime number

    OUTPUT: none

    This checks to see if the different bases have the same length, and
    if the change-of-basis matrices are invertible. If something goes
    wrong, an error message is printed.

    This function checks at the prime `p` as the dimension goes up
    from 0 to ``dim``.

    If you set the Sage verbosity level to a positive integer (using
    ``set_verbose(n)``), then some extra messages will be printed.

    EXAMPLES::

        sage: # long time
        sage: from sage.algebras.steenrod.steenrod_algebra_bases import steenrod_basis_error_check
        sage: steenrod_basis_error_check(15, 2)
        sage: steenrod_basis_error_check(15, 2, generic=True)
        sage: steenrod_basis_error_check(40, 3)
        sage: steenrod_basis_error_check(80, 5)
    """
    from sage.misc.verbose import verbose
    generic = kwds.get('generic', p != 2)

    if not generic:
        bases = ('adem','woody', 'woodz', 'wall', 'arnona', 'arnonc',
                 'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz',
                 'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz')
    else:
        bases = ('adem',
                 'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz',
                 'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz')

    for i in range(dim):
        if i % 5 == 0:
            verbose("up to dimension %s" % i)
        milnor_dim = len(steenrod_algebra_basis.f(i,'milnor',p=p,generic=generic))
        for B in bases:
            if milnor_dim != len(steenrod_algebra_basis.f(i,B,p,generic=generic)):
                print("problem with milnor/{} in dimension {}".format(B, i))
            mat = convert_to_milnor_matrix.f(i,B,p,generic=generic)
            if mat.nrows() != 0 and not mat.is_invertible():
                print("%s invertibility problem in dim %s at p=%s" % (B, i, p))

    verbose("done checking, no profiles")

    bases = ('pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz')
    if not generic:
        profiles = [(4,3,2,1), (2,2,3,1,1), (0,0,0,2)]
    else:
        profiles = [((3,2,1), ()), ((), (2,1,2)), ((3,2,1), (2,2,2,2))]

    for i in range(dim):
        if i % 5 == 0:
            verbose("up to dimension %s" % i)
        for pro in profiles:
            milnor_dim = len(steenrod_algebra_basis.f(i,'milnor',p=p,profile=pro,generic=generic))
            for B in bases:
                if milnor_dim != len(steenrod_algebra_basis.f(i,B,p,profile=pro,generic=generic)):
                    print("problem with milnor/%s in dimension %s with profile %s" % (B, i, pro))

    verbose("done checking with profiles")
