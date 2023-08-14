# -*- coding: utf-8 -*-
r"""
Chain complexes

This module implements bounded chain complexes of free `R`-modules,
for any commutative ring `R` (although the interesting things, like
homology, only work if `R` is the integers or a field).

Fix a ring `R`.  A chain complex over `R` is a collection of
`R`-modules `\{C_n\}` indexed by the integers, with `R`-module maps
`d_n : C_n \rightarrow C_{n+1}` such that `d_{n+1} \circ d_n = 0` for
all `n`.  The maps `d_n` are called *differentials*.

One can vary this somewhat: the differentials may decrease degree by
one instead of increasing it: sometimes a chain complex is defined
with `d_n : C_n \rightarrow C_{n-1}` for each `n`. Indeed, the
differentials may change dimension by any fixed integer.

Also, the modules may be indexed over an abelian group other than the
integers, e.g., `\ZZ^{m}` for some integer `m \geq 1`, in which case
the differentials may change the grading by any element of that
grading group. The elements of the grading group are generally called
degrees, so `C_n` is the module in degree `n` and so on.

In this implementation, the ring `R` must be commutative and the
modules `C_n` must be free `R`-modules.  As noted above, homology
calculations will only work if the ring `R` is either `\ZZ` or a
field.  The modules may be indexed by any free abelian group.  The
differentials may increase degree by 1 or decrease it, or indeed
change it by any fixed amount: this is controlled by the
``degree_of_differential`` parameter used in defining the chain
complex.

AUTHORS:

- John H. Palmieri (2009-04): initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2013 John H. Palmieri <palmieri@math.washington.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from functools import reduce

from sage.structure.parent import Parent
from sage.structure.element import ModuleElement, is_Vector, coercion_model
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.misc.latex import latex
from sage.misc.superseded import deprecation
from sage.rings.fast_arith import prime_range
from sage.homology.homology_group import HomologyGroup
from sage.misc.persist import register_unpickle_override


def _latex_module(R, m):
    """
    LaTeX string representing a free module over ``R`` of rank ``m``.

    INPUT:

    - ``R`` -- a commutative ring
    - ``m`` -- non-negative integer

    This is used by the ``_latex_`` method for chain complexes.

    EXAMPLES::

        sage: from sage.homology.chain_complex import _latex_module
        sage: _latex_module(ZZ, 3)
        '\\Bold{Z}^{3}'
        sage: _latex_module(ZZ, 0)
        '0'
        sage: _latex_module(GF(3), 1)
        '\\Bold{F}_{3}^{1}'
    """
    if m == 0:
        return str(latex(0))
    return str(latex(FreeModule(R, m)))


def ChainComplex(data=None, base_ring=None, grading_group=None,
                 degree_of_differential=1, degree=1,
                 check=True):
    r"""
    Define a chain complex.

    INPUT:

    - ``data`` -- the data defining the chain complex; see below for
      more details.

    The following keyword arguments are supported:

    - ``base_ring`` -- a commutative ring (optional), the ring over
      which the chain complex is defined. If this is not specified,
      it is determined by the data defining the chain complex.

    - ``grading_group`` -- a additive free abelian group (optional,
      default ``ZZ``), the group over which the chain complex is
      indexed.

    - ``degree_of_differential`` -- element of grading_group
      (optional, default ``1``). The degree of the differential.

    - ``degree`` -- alias for ``degree_of_differential``.

    - ``check`` -- boolean (optional, default ``True``). If ``True``,
      check that each consecutive pair of differentials are
      composable and have composite equal to zero.

    OUTPUT:

    A chain complex.

    .. WARNING::

       Right now, homology calculations will only work if the base
       ring is either `\ZZ` or a field, so please take this into account
       when defining a chain complex.

    Use data to define the chain complex.  This may be in any of the
    following forms.

    1. a dictionary with integers (or more generally, elements of
       grading_group) for keys, and with ``data[n]`` a matrix representing
       (via left multiplication) the differential coming from degree
       `n`.  (Note that the shape of the matrix then determines the
       rank of the free modules `C_n` and `C_{n+d}`.)

    2. a list/tuple/iterable of the form `[C_0, d_0, C_1, d_1, C_2,
       d_2, ...]`, where each `C_i` is a free module and each `d_i` is
       a matrix, as above.  This only makes sense if ``grading_group``
       is `\ZZ` and ``degree`` is 1.

    3. a list/tuple/iterable of the form `[r_0, d_0, r_1, d_1, r_2,
       d_2, \ldots]`, where `r_i` is the rank of the free module `C_i`
       and each `d_i` is a matrix, as above.  This only makes sense if
       ``grading_group`` is `\ZZ` and ``degree`` is 1.

    4. a list/tuple/iterable of the form `[d_0, d_1, d_2, \ldots]` where
       each `d_i` is a matrix, as above.  This only makes sense if
       ``grading_group`` is `\ZZ` and ``degree`` is 1.

    .. NOTE::

       In fact, the free modules `C_i` in case 2 and the ranks `r_i`
       in case 3 are ignored: only the matrices are kept, and from
       their shapes, the ranks of the modules are determined.
       (Indeed, if ``data`` is a list or tuple, then any element which
       is not a matrix is discarded; thus the list may have any number
       of different things in it, and all of the non-matrices will be
       ignored.)  No error checking is done to make sure, for
       instance, that the given modules have the appropriate ranks for
       the given matrices.  However, as long as ``check`` is True, the
       code checks to see if the matrices are composable and that each
       appropriate composite is zero.

    If the base ring is not specified, then the matrices are examined
    to determine a ring over which they are all naturally defined, and
    this becomes the base ring for the complex.  If no such ring can
    be found, an error is raised.  If the base ring is specified, then
    the matrices are converted automatically to this ring when
    defining the chain complex.  If some matrix cannot be converted,
    then an error is raised.

    EXAMPLES::

        sage: ChainComplex()
        Trivial chain complex over Integer Ring

        sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
        sage: C
        Chain complex with at most 2 nonzero terms over Integer Ring

        sage: m = matrix(ZZ, 2, 2, [0, 1, 0, 0])
        sage: D = ChainComplex([m, m], base_ring=GF(2)); D
        Chain complex with at most 3 nonzero terms over Finite Field of size 2
        sage: D == loads(dumps(D))
        True
        sage: D.differential(0)==m, m.is_immutable(), D.differential(0).is_immutable()
        (True, False, True)

    Note that when a chain complex is defined in Sage, new
    differentials may be created: every nonzero module in the chain
    complex must have a differential coming from it, even if that
    differential is zero::

        sage: IZ = ChainComplex({0: identity_matrix(ZZ, 1)})
        sage: diff = IZ.differential()  # the differentials in the chain complex
        sage: diff[-1], diff[0], diff[1]
        ([], [1], [])
        sage: IZ.differential(1).parent()
        Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        sage: mat = ChainComplex({0: matrix(ZZ, 3, 4)}).differential(1)
        sage: mat.nrows(), mat.ncols()
        (0, 3)

    Defining the base ring implicitly::

        sage: ChainComplex([matrix(QQ, 3, 1), matrix(ZZ, 4, 3)])
        Chain complex with at most 3 nonzero terms over Rational Field
        sage: ChainComplex([matrix(GF(125, 'a'), 3, 1), matrix(ZZ, 4, 3)])              # needs sage.rings.finite_rings
        Chain complex with at most 3 nonzero terms over Finite Field in a of size 5^3

    If the matrices are defined over incompatible rings, an error results::

        sage: ChainComplex([matrix(GF(125, 'a'), 3, 1), matrix(QQ, 4, 3)])              # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        TypeError: no common canonical parent for objects with parents:
        'Finite Field in a of size 5^3' and 'Rational Field'

    If the base ring is given explicitly but is not compatible with
    the matrices, an error results::

        sage: ChainComplex([matrix(GF(125, 'a'), 3, 1)], base_ring=QQ)                  # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 0 to a rational
    """
    if grading_group is None:
        grading_group = ZZ
    if degree_of_differential != 1 and degree != 1:
        raise ValueError('specify only one of degree_of_differential or degree, not both')
    if degree_of_differential != 1:
        degree = degree_of_differential
    try:
        degree = grading_group(degree)
    except Exception:
        raise ValueError('degree is not an element of the grading group')

    # transform data into data_dict
    if data is None or (isinstance(data, (list, tuple)) and len(data) == 0):
        data_dict = {}
    elif isinstance(data, dict):  # data is dictionary
        data_dict = data
    else:  # data is list/tuple/iterable
        data_matrices = [x for x in data if isinstance(x, Matrix)]
        if degree != 1:
            raise ValueError('degree must be +1 if the data argument is a list or tuple')
        if grading_group != ZZ:
            raise ValueError('grading_group must be ZZ if the data argument is a list or tuple')
        data_dict = {grading_group(i): m for i, m in enumerate(data_matrices)}

    if base_ring is None:
        if not data_dict:
            base_ring = ZZ
        else:
            bases = tuple(x.base_ring() for x in data_dict.values())
            base_ring = coercion_model.common_parent(*bases)

    # make sure values in data_dict are appropriate matrices
    for n in list(data_dict):
        if n not in grading_group:
            raise ValueError('one of the dictionary keys is not an element of the grading group')
        mat = data_dict[n]
        if not isinstance(mat, Matrix):
            raise TypeError('one of the differentials in the data is not a matrix')
        if mat.base_ring() is base_ring:
            if not mat.is_immutable():
                mat = copy(mat)  # do not make any arguments passed immutable
                mat.set_immutable()
        else:
            mat = mat.change_ring(base_ring)
            mat.set_immutable()
        data_dict[n] = mat

    # include any "obvious" zero matrices that are not 0x0
    for n in list(data_dict):  # note: data_dict will be mutated in this loop
        mat1 = data_dict[n]
        if (mat1.nrows(), mat1.ncols()) == (0, 0):
            del data_dict[n]
        if (mat1.nrows() != 0) and (n+degree not in data_dict):
            if n+2*degree in data_dict:
                mat2 = matrix(base_ring, data_dict[n+2*degree].ncols(), mat1.nrows())
            else:
                mat2 = matrix(base_ring, 0, mat1.nrows())
            mat2.set_immutable()
            data_dict[n+degree] = mat2
        if (mat1.ncols() != 0) and (n-degree not in data_dict):
            if n-2*degree in data_dict:
                mat0 = matrix(base_ring, mat1.ncols(), data_dict[n-2*degree].nrows())
            else:
                mat0 = matrix(base_ring, mat1.ncols(), 0)
            mat0.set_immutable()
            data_dict[n-degree] = mat0

    # check that this is a complex: going twice is zero
    if check:
        for n in data_dict:
            mat0 = data_dict[n]
            try:
                mat1 = data_dict[n+degree]
            except KeyError:
                continue
            try:
                prod = mat1 * mat0
            except TypeError:
                raise TypeError('the differentials d_{{{}}} and d_{{{}}} are not compatible: '
                                'their product is not defined'.format(n, n+degree))
            if not prod.is_zero():
                raise ValueError('the differentials d_{{{}}} and d_{{{}}} are not compatible: '
                                 'their composition is not zero.'.format(n, n+degree))

    return ChainComplex_class(grading_group, degree, base_ring, data_dict)


class Chain_class(ModuleElement):

    def __init__(self, parent, vectors, check=True):
        r"""
        A Chain in a Chain Complex

        A chain is collection of module elements for each module `C_n`
        of the chain complex `(C_n, d_n)`. There is no restriction on
        how the differentials `d_n` act on the elements of the chain.

        .. NOTE::

            You must use the chain complex to construct chains.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])},
            ....:                  base_ring=GF(7))
            sage: C.category()
            Category of chain complexes over Finite Field of size 7

        TESTS::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: TestSuite(c).run()
        """
        # only nonzero vectors shall be stored, ensuring this is the
        # job of the _element constructor_
        assert all(v.is_immutable() and not v.is_zero()
                   and v.base_ring() is parent.base_ring()
                   for v in vectors.values())
        self._vec = vectors
        super(Chain_class, self).__init__(parent)

    def vector(self, degree):
        """
        Return the free module element in ``degree``.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([1, 2, 3]), 1: vector([4, 5])})
            sage: c.vector(0)
            (1, 2, 3)
            sage: c.vector(1)
            (4, 5)
            sage: c.vector(2)
            ()
        """
        try:
            return self._vec[degree]
        except KeyError:
            return self.parent().free_module(degree).zero()

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C()
            Trivial chain
            sage: C({0: vector([1, 2, 3])})
            Chain(0:(1, 2, 3))
            sage: c = C({0: vector([1, 2, 3]), 1: vector([4, 5])});  c
            Chain with 2 nonzero terms over Integer Ring
            sage: c._repr_()
            'Chain with 2 nonzero terms over Integer Ring'
        """
        n = len(self._vec)
        if n == 0:
            return 'Trivial chain'

        if n == 1:
            deg, vec = next(iter(self._vec.items()))
            return 'Chain({0}:{1})'.format(deg, vec)

        return 'Chain with {0} nonzero terms over {1}'.format(
            n, self.parent().base_ring())

    def _ascii_art_(self):
        """
        Return an ascii art representation.

        Note that arrows go to the left so that composition of
        differentials is the usual matrix multiplication.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]),
            ....:                   1: zero_matrix(1,2)})
            sage: c = C({0: vector([1, 2, 3]), 1: vector([4, 5])})
            sage: ascii_art(c)
               d_2       d_1       d_0  [1]  d_-1
            0 <---- [0] <---- [4] <---- [2] <----- 0
                              [5]       [3]
        """
        from sage.typeset.ascii_art import AsciiArt

        def arrow_art(d):
            d_str = ['  d_{0}  '.format(d)]
            arrow = ' <' + '-'*(len(d_str[0])-3) + ' '
            d_str.append(arrow)
            return AsciiArt(d_str, baseline=0)

        def vector_art(d):
            v = self.vector(d)
            if v.degree() == 0:
                return AsciiArt(['0'])
            v = str(v.column()).splitlines()
            return AsciiArt(v, baseline=len(v)//2)

        result = []
        chain_complex = self.parent()
        for ordered in chain_complex.ordered_degrees():
            ordered = list(reversed(ordered))
            if len(ordered) == 0:
                return AsciiArt(['0'])
            result_ordered = vector_art(ordered[0] + chain_complex.degree_of_differential())
            for n in ordered:
                result_ordered += arrow_art(n) + vector_art(n)
            result = [result_ordered] + result
        concatenated = result[0]
        for r in result[1:]:
            concatenated += AsciiArt([' ... ']) + r
        return concatenated

    def _unicode_art_(self):
        """
        Return a unicode art representation.

        Note that arrows go to the left so that composition of
        differentials is the usual matrix multiplication.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]),
            ....:                   1: zero_matrix(1,2)})
            sage: c = C({0: vector([1, 2, 3]), 1: vector([4, 5])})
            sage: unicode_art(c)
                                        ⎛1⎞
               d_2       d_1  ⎛4⎞  d_0  ⎜2⎟  d_-1
            0 <──── (0) <──── ⎝5⎠ <──── ⎝3⎠ <───── 0
        """
        from sage.typeset.unicode_art import UnicodeArt

        def arrow_art(d):
            d_str = [u'  d_{0}  '.format(d)]
            arrow = u' <' + u'─' * (len(d_str[0]) - 3) + u' '
            d_str.append(arrow)
            return UnicodeArt(d_str, baseline=0)

        def vector_art(d):
            v = self.vector(d)
            if not v.degree():
                return UnicodeArt([u'0'])
            w = matrix(v).transpose()
            return w._unicode_art_()

        result = []
        chain_complex = self.parent()
        for ordered in chain_complex.ordered_degrees():
            ordered = list(reversed(ordered))
            if not ordered:
                return UnicodeArt([u'0'])
            result_ordered = vector_art(ordered[0] +
                                        chain_complex.degree_of_differential())
            for n in ordered:
                result_ordered += arrow_art(n) + vector_art(n)
            result = [result_ordered] + result
        concatenated = result[0]
        for r in result[1:]:
            concatenated += UnicodeArt([u' ... ']) + r
        return concatenated

    def is_cycle(self):
        """
        Return whether the chain is a cycle.

        OUTPUT:

        Boolean. Whether the elements of the chain are in the kernel
        of the differentials.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: c.is_cycle()
            True
        """
        chain_complex = self.parent()
        for d, v in self._vec.items():
            dv = chain_complex.differential(d) * v
            if not dv.is_zero():
                return False
        return True

    def is_boundary(self):
        """
        Return whether the chain is a boundary.

        OUTPUT:

        Boolean. Whether the elements of the chain are in the image of
        the differentials.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: c.is_boundary()
            False
            sage: z3 = C({1:(1, 0)})
            sage: z3.is_cycle()
            True
            sage: (2*z3).is_boundary()
            False
            sage: (3*z3).is_boundary()
            True
        """
        chain_complex = self.parent()
        for d, v in self._vec.items():
            d = chain_complex.differential(d - chain_complex.degree_of_differential()).transpose()
            if v not in d.image():
                return False
        return True

    def _add_(self, other):
        """
        Module addition

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: c + c
            Chain with 2 nonzero terms over Integer Ring
            sage: ascii_art(c + c)
               d_1       d_0  [0]  d_-1
            0 <---- [6] <---- [2] <----- 0
                    [8]       [4]
        """
        vectors = {}
        for d in set(list(self._vec) + list(other._vec)):
            v = self.vector(d) + other.vector(d)
            if not v.is_zero():
                v.set_immutable()
                vectors[d] = v
        parent = self.parent()
        return parent.element_class(parent, vectors)

    def _lmul_(self, scalar):
        """
        Scalar multiplication

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: 2 * c
            Chain with 2 nonzero terms over Integer Ring
            sage: 2 * c == c + c == c * 2
            True
        """
        vectors = dict()
        for d, v in self._vec.items():
            v = scalar * v
            if not v.is_zero():
                v.set_immutable()
                vectors[d] = v
        parent = self.parent()
        return parent.element_class(parent, vectors)

    def __eq__(self, other):
        """
        Return ``True`` if this chain is equal to ``other``.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: c == c
            True
            sage: c == C(0)
            False
        """
        if type(self) is not type(other) or self.parent() != other.parent():
            return False
        return self._vec == other._vec

    def __ne__(self, other):
        """
        Return ``True`` if this chain is not equal to ``other``.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: c = C({0: vector([0, 1, 2]), 1: vector([3, 4])})
            sage: c != c
            False
            sage: c != C(0)
            True
        """
        return not self == other


class ChainComplex_class(Parent):
    r"""
    See :func:`ChainComplex` for full documentation.

    The differentials are required to be in the following canonical form:

    * All differentials that are not `0 \times 0` must be specified
      (even if they have zero rows or zero columns), and

    * Differentials that are `0 \times 0` must not be specified.

    * Immutable matrices over the ``base_ring``

    This and more is ensured by the assertions in the
    constructor. The :func:`ChainComplex` factory function must
    ensure that only valid input is passed.

    EXAMPLES::

        sage: C = ChainComplex(); C
        Trivial chain complex over Integer Ring

        sage: D = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
        sage: D
        Chain complex with at most 2 nonzero terms over Integer Ring
    """
    def __init__(self, grading_group, degree_of_differential, base_ring, differentials):
        """
        Initialize ``self``.

        TESTS::

            sage: ChainComplex().base_ring()
            Integer Ring

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: TestSuite(C).run()
        """
        if any(d.base_ring() != base_ring or not d.is_immutable() or
               (d.ncols(), d.nrows()) == (0, 0)
               for d in differentials.values()):
            raise ValueError('invalid differentials')
        if degree_of_differential.parent() is not grading_group:
            raise ValueError('the degree_of_differential.parent() must be grading_group')
        if grading_group is not ZZ and grading_group.is_multiplicative():
            raise ValueError('grading_group must be either ZZ or multiplicative')
        # all differentials (excluding the 0x0 ones) must be specified to the constructor
        if any(dim+degree_of_differential not in differentials and d.nrows() != 0
               for dim, d in differentials.items()):
            raise ValueError('invalid differentials')
        if any(dim-degree_of_differential not in differentials and d.ncols() != 0
               for dim, d in differentials.items()):
            raise ValueError('invalid differentials')
        self._grading_group = grading_group
        self._degree_of_differential = degree_of_differential
        self._diff = differentials

        from sage.categories.chain_complexes import ChainComplexes
        category = ChainComplexes(base_ring)
        super(ChainComplex_class, self).__init__(base=base_ring, category=category)

    Element = Chain_class

    def _element_constructor_(self, vectors, check=True):
        """
        The element constructor.

        This is part of the Parent/Element framework. Calling the
        parent uses this method to construct elements.

        TESTS::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D._element_constructor_(0)
            Trivial chain
            sage: D({0:[2, 3]})
            Chain(0:(2, 3))
        """
        if not vectors:  # special case: the zero chain
            return self.element_class(self, {})
        if isinstance(vectors, Chain_class):
            vectors = vectors._vec
        data = dict()
        for degree, vec in vectors.items():
            if not is_Vector(vec):
                vec = vector(self.base_ring(), vec)
                vec.set_immutable()
            if check and vec.degree() != self.free_module_rank(degree):
                raise ValueError('vector dimension does not match module dimension')
            if vec.is_zero():
                continue
            if vec.base_ring() != self.base_ring():
                vec = vec.change_ring(self.base_ring())
            if not vec.is_immutable():
                vec = copy(vec)
                vec.set_immutable()
            data[degree] = vec
        return self.element_class(self, data)

    def random_element(self):
        """
        Return a random element.

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D.random_element()    # random output
            Chain with 1 nonzero terms over Integer Ring
       """
        vec = dict()
        for d in self.nonzero_degrees():
            vec[d] = self.free_module(d).random_element()
        return self(vec)

    _an_element_ = random_element

    @cached_method
    def rank(self, degree, ring=None):
        r"""
        Return the rank of a differential

        INPUT:

        - ``degree`` -- an element `\delta` of the grading
          group. Which differential `d_{\delta}` we want to know the
          rank of

        - ``ring`` -- (optional) a commutative ring `S`;
          if specified, the rank is computed after changing to this ring

        OUTPUT:

        The rank of the differential `d_{\delta} \otimes_R S`, where
        `R` is the base ring of the chain complex.

        EXAMPLES::

            sage: C = ChainComplex({0:matrix(ZZ, [[2]])})
            sage: C.differential(0)
            [2]
            sage: C.rank(0)
            1
            sage: C.rank(0, ring=GF(2))
            0
        """
        degree = self.grading_group()(degree)
        try:
            d = self._diff[degree]
        except IndexError:
            return ZZ.zero()
        if d.nrows() == 0 or d.ncols() == 0:
            return ZZ.zero()
        if ring is None:
            return d.rank()
        return d.change_ring(ring).rank()

    def grading_group(self):
        r"""
        Return the grading group.

        OUTPUT:

        The discrete abelian group that indexes the individual modules
        of the complex. Usually `\ZZ`.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([0, 3])
            sage: C = ChainComplex(grading_group=G, degree=G(vector([1,2])))
            sage: C.grading_group()
            Additive abelian group isomorphic to Z + Z/3
            sage: C.degree_of_differential()
            (1, 2)
        """
        return self._grading_group

    @cached_method
    def nonzero_degrees(self):
        r"""
        Return the degrees in which the module is non-trivial.

        See also :meth:`ordered_degrees`.

        OUTPUT:

        The tuple containing all degrees `n` (grading group elements)
        such that the module `C_n` of the chain is non-trivial.

        EXAMPLES::

            sage: one = matrix(ZZ, [[1]])
            sage: D = ChainComplex({0: one, 2: one, 6:one})
            sage: ascii_art(D)
                       [1]                             [1]       [0]       [1]
            0 <-- C_7 <---- C_6 <-- 0  ...  0 <-- C_3 <---- C_2 <---- C_1 <---- C_0 <-- 0
            sage: D.nonzero_degrees()
            (0, 1, 2, 3, 6, 7)
        """
        return tuple(sorted(n for n, d in self._diff.items()
                            if d.ncols()))

    @cached_method
    def ordered_degrees(self, start=None, exclude_first=False):
        r"""
        Sort the degrees in the order determined by the differential

        INPUT:

        - ``start`` -- (default: ``None``) a degree (element of the grading
          group) or ``None``

        - ``exclude_first`` -- boolean (optional; default:
          ``False``); whether to exclude the lowest degree -- this is a
          handy way to just get the degrees of the non-zero modules,
          as the domain of the first differential is zero.

        OUTPUT:

        If ``start`` has been specified, the longest tuple of degrees

        * containing ``start`` (unless ``start`` would be the first
          and ``exclude_first=True``),

        * in ascending order relative to :meth:`degree_of_differential`, and

        * such that none of the corresponding differentials are `0\times 0`.

        If ``start`` has not been specified, a tuple of such tuples of
        degrees. One for each sequence of non-zero differentials. They
        are returned in sort order.

        EXAMPLES::

            sage: one = matrix(ZZ, [[1]])
            sage: D = ChainComplex({0: one, 2: one, 6:one})
            sage: ascii_art(D)
                       [1]                             [1]       [0]       [1]
            0 <-- C_7 <---- C_6 <-- 0  ...  0 <-- C_3 <---- C_2 <---- C_1 <---- C_0 <-- 0
            sage: D.ordered_degrees()
            ((-1, 0, 1, 2, 3), (5, 6, 7))
            sage: D.ordered_degrees(exclude_first=True)
            ((0, 1, 2, 3), (6, 7))
            sage: D.ordered_degrees(6)
            (5, 6, 7)
            sage: D.ordered_degrees(5, exclude_first=True)
            (6, 7)
        """
        if start is None:
            result = []
            degrees = set(self._diff)
            while degrees:
                ordered = self.ordered_degrees(degrees.pop())
                degrees.difference_update(ordered)
                if exclude_first:
                    ordered = tuple(ordered[1:])
                result.append(ordered)
            result.sort()
            return tuple(result)

        from collections import deque
        result = deque()
        result.append(start)

        next_deg = start + self.degree_of_differential()
        while next_deg in self._diff:
            result.append(next_deg)
            next_deg += self.degree_of_differential()

        prev_deg = start - self.degree_of_differential()
        while prev_deg in self._diff:
            result.appendleft(prev_deg)
            prev_deg -= self.degree_of_differential()

        if exclude_first:
            result.popleft()
        return tuple(result)

    def degree_of_differential(self):
        """
        Return the degree of the differentials of the complex

        OUTPUT:

        An element of the grading group.

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D.degree_of_differential()
            1
        """
        return self._degree_of_differential

    def differential(self, dim=None):
        """
        The differentials which make up the chain complex.

        INPUT:

        - ``dim`` -- element of the grading group (optional, default
          ``None``); if this is ``None``, return a dictionary of all
          of the differentials, or if this is a single element, return
          the differential starting in that dimension

        OUTPUT:

        Either a dictionary of all of the differentials or a single
        differential (i.e., a matrix).

        EXAMPLES::

            sage: D = ChainComplex({0: matrix(ZZ, 2, 2, [1,0,0,2])})
            sage: D.differential(0)
            [1 0]
            [0 2]
            sage: D.differential(-1)
            []
            sage: C = ChainComplex({0: identity_matrix(ZZ, 40)})
            sage: diff = C.differential()
            sage: diff[-1]
            40 x 0 dense matrix over Integer Ring (use the '.str()' method to see the entries)
            sage: diff[0]
            40 x 40 dense matrix over Integer Ring (use the '.str()' method to see the entries)
            sage: diff[1]
            []
        """
        if dim is None:
            return copy(self._diff)
        dim = self.grading_group()(dim)
        try:
            return self._diff[dim]
        except KeyError:
            pass
        # all differentials that are not 0x0 are in self._diff
        # TODO: turn differentials into morphisms between free modules?
        return matrix(self.base_ring(), 0, 0)

    def dual(self):
        """
        The dual chain complex to ``self``.

        Since all modules in ``self`` are free of finite rank, the
        dual in dimension `n` is isomorphic to the original chain
        complex in dimension `n`, and the corresponding boundary
        matrix is the transpose of the matrix in the original complex.
        This converts a chain complex to a cochain complex and vice versa.

        EXAMPLES::

            sage: C = ChainComplex({2: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.degree_of_differential()
            1
            sage: C.differential(2)
            [3 0 0]
            [0 0 0]
            sage: C.dual().degree_of_differential()
            -1
            sage: C.dual().differential(3)
            [3 0]
            [0 0]
            [0 0]
        """
        data = {}
        deg = self.degree_of_differential()
        for d in self.differential():
            data[(d+deg)] = self.differential()[d].transpose()
        return ChainComplex(data, degree=-deg)

    def free_module_rank(self, degree):
        r"""
        Return the rank of the free module at the given ``degree``.

        INPUT:

        - ``degree`` -- an element of the grading group

        OUTPUT:

        Integer. The rank of the free module `C_n` at the given degree
        `n`.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]), 1: matrix(ZZ, [[0, 1]])})
            sage: [C.free_module_rank(i) for i in range(-2, 5)]
            [0, 0, 3, 2, 1, 0, 0]
        """
        try:
            return self._diff[degree].ncols()
        except KeyError:
            return ZZ.zero()

    def free_module(self, degree=None):
        r"""
        Return the free module at fixed ``degree``, or their sum.

        INPUT:

        - ``degree`` -- an element of the grading group or ``None`` (default).

        OUTPUT:

        The free module `C_n` at the given degree `n`. If the degree
        is not specified, the sum `\bigoplus C_n` is returned.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]), 1: matrix(ZZ, [[0, 1]])})
            sage: C.free_module()
            Ambient free module of rank 6 over the principal ideal domain Integer Ring
            sage: C.free_module(0)
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: C.free_module(1)
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: C.free_module(2)
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
        """
        if degree is None:
            rank = sum([mat.ncols() for mat in self.differential().values()])
        else:
            rank = self.free_module_rank(degree)
        return FreeModule(self.base_ring(), rank)

    def __hash__(self):
        """
        The hash is formed by combining the hashes of

        - the base ring
        - the differentials -- the matrices and their degrees
        - the degree of the differential of the chain complex

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: D = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: hash(C) == hash(D)
            True
        """
        return (hash(self.base_ring())
                ^ hash(tuple(self.differential().items()))
                ^ hash(self.degree_of_differential()))

    def __eq__(self, other):
        """
        Return ``True`` iff this chain complex is the same as other: that
        is, if the base rings and the matrices of the two are the
        same.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])},
            ....:                  base_ring=GF(2))
            sage: D = ChainComplex({0: matrix(GF(2), 2, 3, [1, 0, 0, 0, 0, 0]),
            ....:                   1: matrix(ZZ, 0, 2),
            ....:                   3: matrix(ZZ, 0, 0)})  # base_ring determined from the matrices
            sage: C == D
            True
        """
        if not isinstance(other, ChainComplex_class) or self.base_ring() != other.base_ring():
            return False
        R = self.base_ring()
        equal = True
        for d, mat in self.differential().items():
            if d not in other.differential():
                equal = equal and mat.ncols() == 0 and mat.nrows() == 0
            else:
                equal = (equal and
                         other.differential()[d].change_ring(R) == mat.change_ring(R))
        for d, mat in other.differential().items():
            if d not in self.differential():
                equal = equal and mat.ncols() == 0 and mat.nrows() == 0
        return equal

    def __ne__(self, other):
        """
        Return ``True`` iff this chain complex is not the same as other.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])},
            ....:                  base_ring=GF(2))
            sage: D = ChainComplex({0: matrix(GF(2), 2, 3, [1, 0, 0, 0, 0, 0]),
            ....:                   1: matrix(ZZ, 0, 2),
            ....:                   3: matrix(ZZ, 0, 0)})  # base_ring determined from the matrices
            sage: C != D
            False
            sage: E = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])},
            ....:                  base_ring=ZZ)
            sage: C != E
            True
        """
        return not self == other

    def _homology_chomp(self, deg, base_ring, verbose, generators):
        """
        Helper function for :meth:`homology`.

        This function is deprecated.

        INPUT:

        - ``deg`` -- integer (one specific homology group) or ``None``
          (all of those that can be non-zero)

        - ``base_ring`` -- the base ring (must be the integers
          or a prime field)

        - ``verbose`` -- boolean, whether to print some messages

        - ``generators`` --  boolean, whether to also return generators
          for homology

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])}, base_ring=GF(2))
            sage: C._homology_chomp(None, GF(2), False, False)  # optional - chomp, needs sage.rings.finite_rings
            doctest:...: DeprecationWarning: the CHomP interface is deprecated; hence so is this function
            See https://github.com/sagemath/sage/issues/33777 for details.
            {0: Vector space of dimension 2 over Finite Field of size 2, 1: Vector space of dimension 1 over Finite Field of size 2}

            sage: D = ChainComplex({0: matrix(ZZ,1,0,[]), 1: matrix(ZZ,1,1,[0]),
            ....:   2: matrix(ZZ,0,1,[])})
            sage: D._homology_chomp(None, GF(2), False, False)  # optional - chomp, needs sage.rings.finite_rings
            {1: Vector space of dimension 1 over Finite Field of size 2,
            2: Vector space of dimension 1 over Finite Field of size 2}
        """
        deprecation(33777, "the CHomP interface is deprecated; hence so is this function")
        from sage.interfaces.chomp import homchain
        H = homchain(self, base_ring=base_ring, verbose=verbose,
                     generators=generators)
        if H is None:
            raise RuntimeError('ran CHomP, but no output')
        if deg is None:
            # all the homology groups that could be non-zero
            # one has to complete the answer of chomp
            result = H
            for idx in self.nonzero_degrees():
                if idx not in H:
                    result[idx] = HomologyGroup(0, base_ring)
            return result
        if deg in H:
            return H[deg]
        else:
            return HomologyGroup(0, base_ring)

    def homology(self, deg=None, base_ring=None, generators=False,
                 verbose=False, algorithm='pari'):
        r"""
        The homology of the chain complex.

        INPUT:

        - ``deg`` -- an element of the grading group for the chain
          complex (default: ``None``); the degree in which
          to compute homology -- if this is ``None``, return the
          homology in every degree in which the chain complex is
          possibly nonzero.

        - ``base_ring`` -- a commutative ring (optional, default is the
          base ring for the chain complex); must be either the
          integers `\ZZ` or a field

        - ``generators`` -- boolean (optional, default ``False``); if
          ``True``, return generators for the homology groups along with
          the groups. See :trac:`6100`

        - ``verbose`` - boolean (optional, default ``False``); if
          ``True``, print some messages as the homology is computed

        - ``algorithm`` - string (optional, default ``'pari'``); the
          options are:

          * ``'auto'``
          * ``'dhsw'``
          * ``'pari'``
          * ``'chomp'`` (this option is deprecated)

          See below for descriptions.

        OUTPUT:

        If the degree is specified, the homology in degree ``deg``.
        Otherwise, the homology in every dimension as a dictionary
        indexed by dimension.

        ALGORITHM:

        Over a
        field, just compute ranks and nullities, thus obtaining
        dimensions of the homology groups as vector spaces.  Over the
        integers, compute Smith normal form of the boundary matrices
        defining the chain complex according to the value of
        ``algorithm``.  If ``algorithm`` is ``'auto'``,
        then for each relatively small matrix, use the standard Sage
        method, which calls the Pari package.  For any large matrix,
        reduce it using the Dumas, Heckenbach, Saunders, and Welker
        elimination algorithm [DHSW2003]_: see
        :func:`~sage.homology.matrix_utils.dhsw_snf` for details.

        ``'no_chomp'`` is a synonym for ``'auto'``, maintained for
        backward-compatibility.

        ``algorithm`` may also be ``'pari'`` or ``'dhsw'``, which
        forces the named algorithm to be used regardless of the size
        of the matrices.

        Finally, if ``algorithm`` is set to ``'chomp'``, then use
        CHomP. CHomP is available at the web page
        http://chomp.rutgers.edu/, although the software has not been
        tested recently in Sage. The use of this option is deprecated;
        see :trac:`33777`.

        As of this writing, ``'pari'`` is the fastest standard option.

        .. WARNING::

           This only works if the base ring is the integers or a
           field.  Other values will return an error.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.homology()
            {0: Z x Z, 1: Z x C3}
            sage: C.homology(deg=1, base_ring=GF(3))
            Vector space of dimension 2 over Finite Field of size 3
            sage: D = ChainComplex({0: identity_matrix(ZZ, 4), 4: identity_matrix(ZZ, 30)})
            sage: D.homology()
            {0: 0, 1: 0, 4: 0, 5: 0}

        Generators: generators are given as a list of cycles, each of
        which is an element in the appropriate free module, and hence
        is represented as a vector. Each summand of the homology is
        listed separately, with a corresponding generator::

            sage: C.homology(1, generators=True)
            [(C3, Chain(1:(1, 0))), (Z, Chain(1:(0, 1)))]

        Tests for :trac:`6100`, the Klein bottle with generators::

            sage: d0 = matrix(ZZ, 0,1)
            sage: d1 = matrix(ZZ, 1,3, [[0,0,0]])
            sage: d2 = matrix(ZZ, 3,2, [[1,1], [1,-1], [-1,1]])
            sage: C_k = ChainComplex({0:d0, 1:d1, 2:d2}, degree=-1)
            sage: C_k.homology(generators=true)
            {0: [(Z, Chain(0:(1)))],
             1: [(C2, Chain(1:(0, 1, -1))), (Z, Chain(1:(0, 1, 0)))],
             2: []}

        From a torus using a field::

            sage: T = simplicial_complexes.Torus()                                      # needs sage.graphs
            sage: C_t = T.chain_complex()                                               # needs sage.graphs
            sage: C_t.homology(base_ring=QQ, generators=True)                           # needs sage.graphs
            {0: [(Vector space of dimension 1 over Rational Field,
               Chain(0:(0, 0, 0, 0, 0, 0, 1)))],
             1: [(Vector space of dimension 1 over Rational Field,
               Chain(1:(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1))),
              (Vector space of dimension 1 over Rational Field,
               Chain(1:(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, -1, 0)))],
             2: [(Vector space of dimension 1 over Rational Field,
               Chain(2:(1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1)))]}
        """
        if deg is not None and deg not in self.grading_group():
            raise ValueError('degree is not an element of the grading group')

        if base_ring is None:
            base_ring = self.base_ring()
        if not (base_ring.is_field() or base_ring is ZZ):
            raise NotImplementedError('can only compute homology if the base ring is the integers or a field')

        if algorithm not in ['dhsw', 'pari', 'auto', 'no_chomp', 'chomp']:
            raise NotImplementedError('algorithm not recognized')
        if algorithm == 'chomp':
            return self._homology_chomp(deg, base_ring, verbose, generators)

        if deg is None:
            deg = self.nonzero_degrees()
        if isinstance(deg, (list, tuple)):
            answer = {}
            for deg in self.nonzero_degrees():
                answer[deg] = self._homology_in_degree(deg, base_ring, verbose, generators, algorithm)
            return answer
        else:
            return self._homology_in_degree(deg, base_ring, verbose, generators, algorithm)

    def _homology_in_degree(self, deg, base_ring, verbose, generators, algorithm):
        """
        Helper method for :meth:`homology`.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.homology(1) == C._homology_in_degree(1, ZZ, False, False, 'auto')
            True
        """
        if deg not in self.nonzero_degrees():
            zero_homology = HomologyGroup(0, base_ring)
            if generators:
                return (zero_homology, vector(base_ring, []))
            else:
                return zero_homology
        if verbose:
            print('Computing homology of the chain complex in dimension %s...' % deg)

        fraction_field = base_ring.fraction_field()

        def change_ring(X):
            if X.base_ring() is base_ring:
                return X
            return X.change_ring(base_ring)

        # d_out is the differential going out of degree deg,
        # d_in is the differential entering degree deg
        differential = self.degree_of_differential()
        d_in = change_ring(self.differential(deg - differential))
        d_out = change_ring(self.differential(deg))
        d_out_rank = self.rank(deg, ring=fraction_field)
        d_out_nullity = d_out.ncols() - d_out_rank

        if d_in.is_zero():
            if generators:  # Include the generators of the nullspace
                return [(HomologyGroup(1, base_ring), self({deg: gen}))
                        for gen in d_out.right_kernel().basis()]
            else:
                return HomologyGroup(d_out_nullity, base_ring)

        if generators:
            orders, gens = self._homology_generators_snf(d_in, d_out, d_out_rank)
            answer = [(HomologyGroup(1, base_ring, [order]), self({deg: gen}))
                      for order, gen in zip(orders, gens)]
        else:
            if base_ring.is_field():
                d_in_rank = self.rank(deg-differential, ring=base_ring)
                answer = HomologyGroup(d_out_nullity - d_in_rank, base_ring)
            elif base_ring == ZZ:
                if d_in.ncols() == 0:
                    all_divs = [0] * d_out_nullity
                else:
                    if algorithm in ['auto', 'no_chomp']:
                        if ((d_in.ncols() > 300 and d_in.nrows() > 300)
                            or (min(d_in.ncols(), d_in.nrows()) > 100 and
                                d_in.ncols() + d_in.nrows() > 600)):
                            algorithm = 'dhsw'
                        else:
                            algorithm = 'pari'
                    if algorithm == 'dhsw':
                        from sage.homology.matrix_utils import dhsw_snf
                        all_divs = dhsw_snf(d_in, verbose=verbose)
                    elif algorithm == 'pari':
                        all_divs = d_in.elementary_divisors(algorithm)
                    else:
                        raise ValueError('unsupported algorithm')
                all_divs = all_divs[:d_out_nullity]
                # divisors equal to 1 produce trivial
                # summands, so filter them out
                divisors = [x for x in all_divs if x != 1]
                answer = HomologyGroup(len(divisors), base_ring, divisors)
            else:
                raise NotImplementedError('only base rings ZZ and fields are supported')
        return answer

    def _homology_generators_snf(self, d_in, d_out, d_out_rank):
        """
        Compute the homology generators using the Smith normal form.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.homology(1)
            Z x C3
            sage: C._homology_generators_snf(C.differential(0), C.differential(1), 0)
            ([3, 0], [(1, 0), (0, 1)])
        """
        # Find the kernel of the out-going differential.
        K = d_out.right_kernel().matrix().transpose().change_ring(d_out.base_ring())

        # Compute the induced map to the kernel
        S = K.augment(d_in).hermite_form()
        d_in_induced = S.submatrix(row=0, nrows=d_in.nrows()-d_out_rank,
                                   col=d_in.nrows()-d_out_rank, ncols=d_in.ncols())

        # Find the SNF of the induced matrix and appropriate generators
        (N, P, Q) = d_in_induced.smith_form()
        all_divs = [0]*N.nrows()
        non_triv = 0
        for i in range(0, N.nrows()):
            if i >= N.ncols():
                break
            all_divs[i] = N[i][i]
            if N[i][i] == 1:
                non_triv = non_triv + 1
        divisors = [x for x in all_divs if x != 1]
        gens = (K * P.inverse().submatrix(col=non_triv)).columns()
        return divisors, gens

    def betti(self, deg=None, base_ring=None):
        """
        The Betti number of the chain complex.

        That is, write the homology in this degree as a direct sum
        of a free module and a torsion module; the Betti number is the
        rank of the free summand.

        INPUT:

        - ``deg`` -- an element of the grading group for the chain
          complex or None (default ``None``); if ``None``,
          then return every Betti number, as a dictionary indexed by
          degree, or if an element of the grading group, then return
          the Betti number in that degree

        - ``base_ring`` -- a commutative ring (optional, default is the
          base ring for the chain complex); compute homology with
          these coefficients -- must be either the integers or a
          field

        OUTPUT:

        The Betti number in degree ``deg`` -- the rank of the free
        part of the homology module in this degree.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.betti(0)
            2
            sage: [C.betti(n) for n in range(5)]
            [2, 1, 0, 0, 0]
            sage: C.betti()
            {0: 2, 1: 1}

            sage: D = ChainComplex({0: matrix(GF(5), [[3, 1],[1, 2]])})
            sage: D.betti()
            {0: 1, 1: 1}
        """
        if base_ring is None:
            base_ring = self.base_ring()
        try:
            base_ring = base_ring.fraction_field()
        except AttributeError:
            raise NotImplementedError('only implemented if the base ring is ZZ or a field')
        H = self.homology(deg, base_ring=base_ring)
        if isinstance(H, dict):
            return {deg: homology_group.dimension()
                    for deg, homology_group in H.items()}
        else:
            return H.dimension()

    def torsion_list(self, max_prime, min_prime=2):
        r"""
        Look for torsion in this chain complex by computing its mod `p`
        homology for a range of primes `p`.

        INPUT:

        -  ``max_prime`` -- prime number; search for torsion mod `p` for
           all `p` strictly less than this number

        -  ``min_prime`` -- prime (optional, default 2); search for
           torsion mod `p` for primes at least as big as this

        Return a list of pairs `(p, d)` where `p` is a prime at which
        there is torsion and `d` is a list of dimensions in which this
        torsion occurs.

        The base ring for the chain complex must be the integers; if
        not, an error is raised.

        ALGORITHM:

        Let `C` denote the chain complex.  Let `P` equal
        ``max_prime``.  Compute the mod `P` homology of `C`, and use
        this as the base-line computation: the assumption is that this
        is isomorphic to the integral homology tensored with
        `\GF{P}`.  Then compute the mod `p` homology for a range of
        primes `p`, and record whenever the answer differs from the
        base-line answer.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.homology()
            {0: Z x Z, 1: Z x C3}
            sage: C.torsion_list(11)                                                    # needs sage.rings.finite_rings
            [(3, [1])]
            sage: C = ChainComplex([matrix(ZZ, 1, 1, [2]), matrix(ZZ, 1, 1), matrix(1, 1, [3])])
            sage: C.homology(1)
            C2
            sage: C.homology(3)
            C3
            sage: C.torsion_list(5)                                                     # needs sage.rings.finite_rings
            [(2, [1]), (3, [3])]
        """
        if self.base_ring() != ZZ:
            raise NotImplementedError('only implemented for base ring the integers')

        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

        answer = []
        torsion_free = self.betti(base_ring=GF(max_prime))
        for p in prime_range(min_prime, max_prime):
            mod_p_betti = self.betti(base_ring=GF(p))
            if mod_p_betti != torsion_free:
                diff_dict = {}
                temp_diff = {}
                D = self.degree_of_differential()
                for i in torsion_free:
                    temp_diff[i] = mod_p_betti.get(i, 0) - torsion_free[i]
                for i in temp_diff:
                    if temp_diff[i] > 0:
                        if i+D in diff_dict:
                            lower = diff_dict[i+D]
                        else:
                            lower = 0
                        current = temp_diff[i]
                        if current > lower:
                            diff_dict[i] = current - lower
                            if i-D in diff_dict:
                                diff_dict[i-D] -= current - lower
                differences = []
                for i in diff_dict:
                    if diff_dict[i] != 0:
                        differences.append(i)
                answer.append((p, differences))
        return answer

    def _Hom_(self, other, category=None):
        """
        Return the set of chain maps between chain complexes ``self``
        and ``other``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: S = simplicial_complexes.Sphere(2)
            sage: T = simplicial_complexes.Torus()
            sage: C = S.chain_complex(augmented=True, cochain=True)
            sage: D = T.chain_complex(augmented=True, cochain=True)
            sage: Hom(C, D)  # indirect doctest
            Set of Morphisms from Chain complex with at most 4 nonzero terms over
            Integer Ring to Chain complex with at most 4 nonzero terms over Integer
            Ring in Category of chain complexes over Integer Ring
        """
        from sage.homology.chain_complex_homspace import ChainComplexHomspace
        return ChainComplexHomspace(self, other)

    def _flip_(self):
        """
        Flip chain complex upside down (degree `n` gets changed to
        degree `-n`), thus turning a chain complex into a cochain complex
        without changing the homology (except for flipping it, too).

        EXAMPLES::

            sage: C = ChainComplex({2: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C.degree_of_differential()
            1
            sage: C.differential(2)
            [3 0 0]
            [0 0 0]
            sage: C._flip_().degree_of_differential()
            -1
            sage: C._flip_().differential(-2)
            [3 0 0]
            [0 0 0]
        """
        data = {}
        deg = self.degree_of_differential()
        for d in self.differential():
            data[-d] = self.differential()[d]
        return ChainComplex(data, degree=-deg)

    def shift(self, n=1):
        """
        Shift this chain complex `n` times.

        INPUT:

        - ``n`` -- an integer (optional, default 1)

        The *shift* operation is also sometimes called *translation* or
        *suspension*.

        To shift a chain complex by `n`, shift its entries up by `n`
        (if it is a chain complex) or down by `n` (if it is a cochain
        complex); that is, shifting by 1 always shifts in the opposite
        direction of the differential. In symbols, if `C` is a chain
        complex and `C[n]` is its `n`-th shift, then `C[n]_j =
        C_{j-n}`. The differential in the shift `C[n]` is obtained by
        multiplying each differential in `C` by `(-1)^n`.

        Caveat: different sources use different conventions for
        shifting: what we call `C[n]` might be called `C[-n]` in some
        places. See for example.
        https://ncatlab.org/nlab/show/suspension+of+a+chain+complex
        (which uses `C[n]` as we do but acknowledges `C[-n]`) or 1.2.8
        in [Wei1994]_ (which uses `C[-n]`).

        EXAMPLES::

            sage: # needs sage.graphs
            sage: S1 = simplicial_complexes.Sphere(1).chain_complex()
            sage: S1.shift(1).differential(2) == -S1.differential(1)
            True
            sage: S1.shift(2).differential(3) == S1.differential(1)
            True
            sage: S1.shift(3).homology(4)
            Z

        For cochain complexes, shifting goes in the other
        direction. Topologically, this makes sense if we grade the
        cochain complex for a space negatively::

            sage: # needs sage.graphs
            sage: T = simplicial_complexes.Torus()
            sage: co_T = T.chain_complex()._flip_()
            sage: co_T.homology()
            {-2: Z, -1: Z x Z, 0: Z}
            sage: co_T.degree_of_differential()
            1
            sage: co_T.shift(2).homology()
            {-4: Z, -3: Z x Z, -2: Z}

        You can achieve the same result by tensoring (on the left, to
        get the signs right) with a rank one free module in degree
        ``-n * deg``, if ``deg`` is the degree of the differential::

            sage: C = ChainComplex({-2: matrix(ZZ, 0, 1)})
            sage: C.tensor(co_T).homology()                                             # needs sage.graphs
            {-4: Z, -3: Z x Z, -2: Z}
        """
        deg = self.degree_of_differential()
        shift = n * deg
        sgn = (-1)**n
        return ChainComplex({k-shift: sgn * self._diff[k] for k in self._diff},
                            degree_of_differential=deg)

    def _chomp_repr_(self):
        r"""
        String representation of ``self`` suitable for use by the CHomP
        program.

        This function is deprecated.

        Since CHomP can only handle chain complexes, not cochain
        complexes, and since it likes its complexes to start in degree
        0, flip the complex over if necessary, and shift it to start
        in degree 0.  Note also that CHomP only works over the
        integers or a finite prime field.

        EXAMPLES::

            sage: C = ChainComplex({-2: matrix(ZZ, 1, 3, [3, 0, 0])}, degree=-1)
            sage: C._chomp_repr_()
            doctest:...: DeprecationWarning: the CHomP interface is deprecated; hence so is this function
            See https://github.com/sagemath/sage/issues/33777 for details.
            'chain complex\n\nmax dimension = 1\n\ndimension 0\n   boundary a1 = 0\n\ndimension 1\n   boundary a1 = + 3 * a1 \n   boundary a2 = 0\n   boundary a3 = 0\n\n'
            sage: C = ChainComplex({-2: matrix(ZZ, 1, 3, [3, 0, 0])}, degree=1)
            sage: C._chomp_repr_()
            'chain complex\n\nmax dimension = 1\n\ndimension 0\n   boundary a1 = 0\n\ndimension 1\n   boundary a1 = + 3 * a1 \n   boundary a2 = 0\n   boundary a3 = 0\n\n'
        """
        deprecation(33777, "the CHomP interface is deprecated; hence so is this function")
        deg = self.degree_of_differential()
        if (self.grading_group() != ZZ or
                (deg != 1 and deg != -1)):
            raise ValueError('CHomP only works on Z-graded chain complexes with '
                             'differential of degree 1 or -1')
        base_ring = self.base_ring()
        if (base_ring == QQ) or (base_ring != ZZ and not (base_ring.is_prime_field())):
            raise ValueError('CHomP doesn\'t compute over the rationals, only over Z or F_p')
        if deg == -1:
            diffs = self.differential()
        else:
            diffs = self._flip_().differential()

        if len(diffs) == 0:
            diffs = {0: matrix(ZZ, 0, 0)}

        maxdim = max(diffs)
        mindim = min(diffs)
        # will shift chain complex by subtracting mindim from
        # dimensions, so its bottom dimension is zero.
        s = "chain complex\n\nmax dimension = %s\n\n" % (maxdim - mindim - 1,)

        for i in range(0, maxdim - mindim):
            s += "dimension %s\n" % i
            mat = diffs.get(i + mindim, matrix(base_ring, 0, 0))
            for idx in range(mat.ncols()):
                s += "   boundary a%s = " % (idx + 1)
                # construct list of bdries
                col = mat.column(idx)
                nonzero_pos = col.nonzero_positions()
                if nonzero_pos:
                    for j in nonzero_pos:
                        entry = col[j]
                        if entry > 0:
                            sgn = "+"
                        else:
                            sgn = "-"
                            entry = -entry
                        s += "%s %s * a%s " % (sgn, entry, j+1)
                else:
                    s += "0"
                s += "\n"
            s += "\n"
        return s

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C
            Chain complex with at most 2 nonzero terms over Integer Ring
        """
        diffs = [mat for mat in self._diff.values() if mat.nrows() + mat.ncols() > 0]
        if len(diffs) == 0:
            s = 'Trivial chain complex'
        else:
            s = 'Chain complex with at most {0} nonzero terms'.format(len(diffs)-1)
        s += ' over {0}'.format(self.base_ring())
        return s

    def _ascii_art_(self):
        """
        Return an ascii art representation.

        Note that arrows go to the left so that composition of
        differentials is the usual matrix multiplication.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]), 1:zero_matrix(1,2)})
            sage: ascii_art(C)
                                    [3 0 0]
                        [0 0]       [0 0 0]
             0 <-- C_2 <------ C_1 <-------- C_0 <-- 0

            sage: one = matrix(ZZ, [[1]])
            sage: D = ChainComplex({0: one, 2: one, 6:one})
            sage: ascii_art(D)
                        [1]                             [1]       [0]       [1]
             0 <-- C_7 <---- C_6 <-- 0  ...  0 <-- C_3 <---- C_2 <---- C_1 <---- C_0 <-- 0
        """
        from sage.typeset.ascii_art import AsciiArt

        def arrow_art(n):
            d_n = self.differential(n)
            if d_n.nrows() == 0 or d_n.ncols() == 0:
                return AsciiArt(['<--'])
            d_str = [' '+line+' ' for line in str(d_n).splitlines()]
            arrow = '<' + '-'*(len(d_str[0])-1)
            d_str.append(arrow)
            return AsciiArt(d_str)

        def module_art(n):
            C_n = self.free_module(n)
            if C_n.rank() == 0:
                return AsciiArt([' 0 '])
            else:
                return AsciiArt([' C_{0} '.format(n)])

        result = []
        for ordered in self.ordered_degrees():
            ordered = list(reversed(ordered))
            if len(ordered) == 0:
                return AsciiArt(['0'])
            result_ordered = module_art(ordered[0] + self.degree_of_differential())
            for n in ordered:
                result_ordered += arrow_art(n) + module_art(n)
            result = [result_ordered] + result
        concatenated = result[0]
        for r in result[1:]:
            concatenated += AsciiArt([' ... ']) + r
        return concatenated

    def _unicode_art_(self):
        """
        Return a unicode art representation.

        Note that arrows go to the left so that composition of
        differentials is the usual matrix multiplication.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0]), 1:zero_matrix(1,2)})
            sage: unicode_art(C)
                                ⎛3 0 0⎞
                      (0 0)     ⎝0 0 0⎠
            0 <── C_2 <──── C_1 <────── C_0 <── 0

            sage: one = matrix(ZZ, [[1]])
            sage: D = ChainComplex({0: one, 2: one, 6:one})
            sage: unicode_art(D)
                      (1)                           (1)     (0)     (1)
            0 <── C_7 <── C_6 <── 0  ...  0 <── C_3 <── C_2 <── C_1 <── C_0 <── 0
        """
        from sage.typeset.unicode_art import UnicodeArt

        def arrow_art(n):
            d_n = self.differential(n)
            if not d_n.nrows() or not d_n.ncols():
                return UnicodeArt([u'<──'])
            d_str = list(d_n._unicode_art_())
            arrow = u'<' + u'─' * (len(d_str[0]) - 1)
            d_str.append(arrow)
            return UnicodeArt(d_str)

        def module_art(n):
            C_n = self.free_module(n)
            if not C_n.rank():
                return UnicodeArt([u' 0 '])
            else:
                return UnicodeArt([u' C_{0} '.format(n)])

        result = []
        for ordered in self.ordered_degrees():
            ordered = list(reversed(ordered))
            if not ordered:
                return UnicodeArt([u'0'])
            result_ordered = module_art(ordered[0] + self.degree_of_differential())
            for n in ordered:
                result_ordered += arrow_art(n) + module_art(n)
            result = [result_ordered] + result
        concatenated = result[0]
        for r in result[1:]:
            concatenated += UnicodeArt([u' ... ']) + r
        return concatenated

    def _latex_(self):
        """
        LaTeX print representation.

        EXAMPLES::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: C._latex_()
            '\\Bold{Z}^{3} \\xrightarrow{d_{0}} \\Bold{Z}^{2}'

            sage: ChainComplex()._latex_()
            '0'

            sage: G = AdditiveAbelianGroup([0, 0])
            sage: m = matrix([0])
            sage: C = ChainComplex(grading_group=G, degree=G(vector([1,2])), data={G.zero(): m})
            sage: C._latex_()
            '\\Bold{Z}^{1} \\xrightarrow{d_{\\text{\\texttt{(0,{ }0)}}}} \\Bold{Z}^{1}'
        """
#         Warning: this is likely to screw up if, for example, the
#         degree of the differential is 2 and there are nonzero terms
#         in consecutive dimensions (e.g., in dimensions 0 and 1).  In
#         such cases, the representation might show a differential
#         connecting these terms, although the differential goes from
#         dimension 0 to dimension 2, and from dimension 1 to
#         dimension 3, etc.  I don't know how much effort should be
#         put into trying to fix this.
        string = ""
        diffs = self._diff
        if len(diffs) == 0:
            return "0"
        deg = self.degree_of_differential()
        ring = self.base_ring()
        backwards = bool(deg < 0)
        sorted_list = sorted(diffs.keys(), reverse=backwards)
        if len(diffs) <= 6:
            for n in sorted_list[1:-1]:
                mat = diffs[n]
                string += _latex_module(ring, mat.ncols())
                string += " \\xrightarrow{d_{%s}} " % latex(n)
            mat = diffs[sorted_list[-1]]
            string += _latex_module(ring, mat.ncols())
        else:
            for n in sorted_list[:2]:
                mat = diffs[n]
                string += _latex_module(ring, mat.ncols())
                string += " \\xrightarrow{d_{%s}} " % latex(n)
            string += "\\dots "
            n = sorted_list[-2]
            string += "\\xrightarrow{d_{%s}} " % latex(n)
            mat = diffs[sorted_list[-1]]
            string += _latex_module(ring, mat.ncols())
        return string

    def cartesian_product(self, *factors, **kwds):
        r"""
        Return the direct sum (Cartesian product) of ``self`` with ``D``.

        Let `C` and `D` be two chain complexes with differentials
        `\partial_C` and `\partial_D`, respectively, of the same degree (so
        they must also have the same grading group).
        The direct sum `S = C \oplus D` is a chain complex given by
        `S_i = C_i \oplus D_i` with differential
        `\partial = \partial_C \oplus \partial_D`.

        INPUT:

        - ``subdivide`` -- (default: ``False``) whether to subdivide the
          the differential matrices

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: C = ChainComplex([matrix([[-y],[x]]), matrix([[x, y]])])
            sage: D = ChainComplex([matrix([[x-y]]), matrix([[0], [0]])])
            sage: ascii_art(C.cartesian_product(D))
                        [x y 0]       [   -y     0]
                        [0 0 0]       [    x     0]
                        [0 0 0]       [    0 x - y]
             0 <-- C_2 <-------- C_1 <-------------- C_0 <-- 0

            sage: D = ChainComplex({1:matrix([[x-y]]), 4:matrix([[x], [y]])})
            sage: ascii_art(D)
                        [x]
                        [y]                     [x - y]
             0 <-- C_5 <---- C_4 <-- 0 <-- C_2 <-------- C_1 <-- 0
            sage: ascii_art(cartesian_product([C, D]))
                                                                          [-y]
                        [x]                     [    x     y     0]       [ x]
                        [y]                     [    0     0 x - y]       [ 0]
             0 <-- C_5 <---- C_4 <-- 0 <-- C_2 <-------------------- C_1 <----- C_0 <-- 0

        The degrees of the differentials must agree::

            sage: C = ChainComplex({1:matrix([[x]])}, degree_of_differential=-1)
            sage: D = ChainComplex({1:matrix([[x]])}, degree_of_differential=1)
            sage: C.cartesian_product(D)
            Traceback (most recent call last):
            ...
            ValueError: the degrees of the differentials must match

        TESTS::

            sage: C = ChainComplex({2:matrix([[-1],[2]]), 1:matrix([[2, 1]])},
            ....:                  degree_of_differential=-1)
            sage: ascii_art(C.cartesian_product(C, subdivide=True))
                                        [-1| 0]
                                        [ 2| 0]
                        [2 1|0 0]       [--+--]
                        [---+---]       [ 0|-1]
                        [0 0|2 1]       [ 0| 2]
             0 <-- C_0 <---------- C_1 <-------- C_2 <-- 0

        ::

            sage: R.<x,y,z> = QQ[]
            sage: C1 = ChainComplex({1:matrix([[x]])})
            sage: C2 = ChainComplex({1:matrix([[y]])})
            sage: C3 = ChainComplex({1:matrix([[z]])})
            sage: ascii_art(cartesian_product([C1, C2, C3]))
                        [x 0 0]
                        [0 y 0]
                        [0 0 z]
             0 <-- C_2 <-------- C_1 <-- 0
            sage: ascii_art(C1.cartesian_product([C2, C3], subdivide=True))
                        [x|0|0]
                        [-+-+-]
                        [0|y|0]
                        [-+-+-]
                        [0|0|z]
             0 <-- C_2 <-------- C_1 <-- 0

        ::

            sage: R.<x> = ZZ[]
            sage: G = AdditiveAbelianGroup([0,7])
            sage: d = {G(vector([1,1])):matrix([[x]])}
            sage: C = ChainComplex(d, grading_group=G, degree=G(vector([2,1])))
            sage: ascii_art(C.cartesian_product(C))
                             [x 0]
                             [0 x]
             0 <-- C_(3, 2) <------ C_(1, 1) <-- 0
        """
        if not factors:
            return self
        if isinstance(factors[0], (list, tuple)):
            factors = factors[0]
        deg_diff = self.degree_of_differential()
        if any(D.degree_of_differential() != deg_diff for D in factors):
            raise ValueError("the degrees of the differentials must match")
        if any(D.grading_group() != self._grading_group for D in factors):
            raise ValueError("the grading groups must match")

        factors = [self] + list(factors)
        R = self.base_ring()
        zero = matrix(R, [])
        subdivide = kwds.get('subdivide', False)
        diffs = [D.differential() for D in factors]
        keys = reduce(lambda X, d: X.union(d.keys()), diffs, set())
        ret = {k: matrix.block_diagonal([d.get(k, zero) for d in diffs],
                                        subdivide=subdivide)
               for k in keys}
        return ChainComplex(ret, degree_of_differential=deg_diff,
                            grading_group=self._grading_group)

    def tensor(self, *factors, **kwds):
        r"""
        Return the tensor product of ``self`` with ``D``.

        Let `C` and `D` be two chain complexes with differentials
        `\partial_C` and `\partial_D`, respectively, of the same degree (so
        they must also have the same grading group).
        The tensor product `S = C \otimes D` is a chain complex given by

        .. MATH::

            S_i = \bigoplus_{a+b=i} C_a \otimes D_b

        with differential

        .. MATH::

            \partial(x \otimes y) = \partial_C x \otimes y
            + (-1)^{|a| \cdot |\partial_D|} x \otimes \partial_D y

        for `x \in C_a` and `y \in D_b`, where `|a|` is the degree of `a` and
        `|\partial_D|` is the degree of `\partial_D`.

        .. WARNING::

            If the degree of the differential is even, then this may not
            result in a valid chain complex.

        INPUT:

        - ``subdivide`` -- (default: ``False``) whether to subdivide the
          the differential matrices

        .. TODO::

            Make subdivision work correctly on multiple factors.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: C1 = ChainComplex({1:matrix([[x]])}, degree_of_differential=-1)
            sage: C2 = ChainComplex({1:matrix([[y]])}, degree_of_differential=-1)
            sage: C3 = ChainComplex({1:matrix([[z]])}, degree_of_differential=-1)
            sage: ascii_art(C1.tensor(C2))
                                    [ x]
                        [y x]       [-y]
             0 <-- C_0 <------ C_1 <----- C_2 <-- 0
            sage: ascii_art(C1.tensor(C2).tensor(C3))
                                      [ y  x  0]       [ x]
                                      [-z  0  x]       [-y]
                        [z y x]       [ 0 -z -y]       [ z]
             0 <-- C_0 <-------- C_1 <----------- C_2 <----- C_3 <-- 0

        ::

            sage: C = ChainComplex({2:matrix([[-y],[x]]), 1:matrix([[x, y]])},
            ....:                  degree_of_differential=-1); ascii_art(C)
                                    [-y]
                        [x y]       [ x]
             0 <-- C_0 <------ C_1 <----- C_2 <-- 0
            sage: T = C.tensor(C)
            sage: T.differential(1)
            [x y x y]
            sage: T.differential(2)
            [-y  x  0  y  0  0]
            [ x  0  x  0  y  0]
            [ 0 -x -y  0  0 -y]
            [ 0  0  0 -x -y  x]
            sage: T.differential(3)
            [ x  y  0  0]
            [ y  0 -y  0]
            [-x  0  0 -y]
            [ 0  y  x  0]
            [ 0 -x  0  x]
            [ 0  0  x  y]
            sage: T.differential(4)
            [-y]
            [ x]
            [-y]
            [ x]

        The degrees of the differentials must agree::

            sage: C1p = ChainComplex({1:matrix([[x]])}, degree_of_differential=1)
            sage: C1.tensor(C1p)
            Traceback (most recent call last):
            ...
            ValueError: the degrees of the differentials must match

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: C1 = ChainComplex({1:matrix([[x]])})
            sage: C2 = ChainComplex({1:matrix([[y]])})
            sage: C3 = ChainComplex({1:matrix([[z]])})
            sage: ascii_art(tensor([C1, C2, C3]))
                                      [-y -z  0]       [ z]
                                      [ x  0 -z]       [-y]
                        [x y z]       [ 0  x  y]       [ x]
             0 <-- C_6 <-------- C_5 <----------- C_4 <----- C_3 <-- 0

        ::

            sage: R.<x,y> = ZZ[]
            sage: G = AdditiveAbelianGroup([0,7])
            sage: d1 = {G(vector([1,1])):matrix([[x]])}
            sage: C1 = ChainComplex(d1, grading_group=G, degree=G(vector([2,1])))
            sage: d2 = {G(vector([3,0])):matrix([[y]])}
            sage: C2 = ChainComplex(d2, grading_group=G, degree=G(vector([2,1])))
            sage: ascii_art(C1.tensor(C2))
                                                [y]
                             [ x -y]            [x]
             0 <-- C_(8, 3) <-------- C_(6, 2) <---- C_(4, 1) <-- 0

        Check that :trac:`21760` is fixed::

            sage: C = ChainComplex({0: matrix(ZZ, 0, 2)}, degree=-1)
            sage: ascii_art(C)
             0 <-- C_0 <-- 0
            sage: T = C.tensor(C)
            sage: ascii_art(T)
             0 <-- C_0 <-- 0
            sage: T.free_module_rank(0)
            4
        """
        if not factors:
            return self
        if isinstance(factors[0], (list, tuple)):
            factors = factors[0]
        deg_diff = self.degree_of_differential()
        if any(D.degree_of_differential() != deg_diff for D in factors):
            raise ValueError("the degrees of the differentials must match")
        if any(D.grading_group() != self._grading_group for D in factors):
            raise ValueError("the grading groups must match")

        R = self.base_ring()
        zero = R.zero()
        subdivide = kwds.get('subdivide', False)
        ret = self

        if self._grading_group is ZZ:
            def scalar(a):
                return (-1)**(a * deg_diff)
        else:
            def scalar(a):
                return (-1)**(sum(a) * sum(deg_diff))

        for D in factors:
            # Setup
            d = ret.differential()
            dD = D.differential()
            deg = sorted((k, ret.free_module_rank(k)) for k in d
                         if ret.free_module_rank(k) > 0)
            degD = sorted((k, D.free_module_rank(k)) for k in dD
                          if D.free_module_rank(k) > 0)
            diff = {}

            # Our choice for tensor products will be x # y = x1 * y + x2 * y + ...

            # Generate the data for the differential
            for a, r in deg:
                for b, s in degD:
                    rp = d[a].nrows()
                    sp = dD[b].nrows()
                    if a+b not in diff:
                        diff[a+b] = {}
                    mor = diff[a+b]
                    cur = {}
                    cur[(a+deg_diff, b)] = []
                    cur[(a, b+deg_diff)] = []

                    for i in range(r):
                        for j in range(s):
                            # \partial x_i \otimes y_j
                            vec = [zero]*(rp*s)
                            for k, val in enumerate(d[a].column(i)):
                                vec[s*k+j] += val
                            cur[(a+deg_diff, b)].append(vec)

                            # (-1)^a x_i \otimes \partial y_j
                            vec = [zero]*(r*sp)
                            for k, val in enumerate(dD[b].column(j)):
                                vec[sp*i+k] += scalar(a) * val
                            cur[(a, b+deg_diff)].append(vec)

                    mor[a, b] = cur

            # Parse the data into matrices
            to_delete = []
            for k in diff:
                # Get the data and interchange the indices
                mor = diff[k]
                row_keys = sorted(mor.keys())
                cols = {}
                col_widths = {}
                for dom in mor:
                    c = mor[dom]
                    for im in c:
                        if im not in cols:
                            cols[im] = {}
                            col_widths[im] = len(c[im])
                        cols[im][dom] = c[im]
                col_keys = sorted(cols.keys())
                # Now build the matrix
                M = []
                for ck in col_keys:
                    M.append([])
                    col = cols[ck]
                    for rk in row_keys:
                        if rk in col:
                            M[-1].append(matrix(R, col[rk]).transpose())
                        else:
                            M[-1].append(zero)
                diff[k] = matrix.block(M, subdivide=subdivide)

                # Flag for removal any 0x0 matrices
                if diff[k].nrows() == 0 and diff[k].ncols() == 0:
                    to_delete.append(k)

            # Delete the 0x0 matrices
            for k in to_delete:
                del diff[k]

            ret = ChainComplex(diff, degree_of_differential=deg_diff,
                               grading_group=self._grading_group)

        return ret


register_unpickle_override('sage.homology.chain_complex', 'ChainComplex', ChainComplex_class)
