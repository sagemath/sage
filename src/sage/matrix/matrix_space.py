r"""
Matrix Spaces

You can create any space `\text{Mat}_{n\times m}(R)` of
either dense or sparse matrices with given number of rows and
columns over any commutative or noncommutative ring.

EXAMPLES::

    sage: MS = MatrixSpace(QQ, 6,6, sparse=True); MS
    Full MatrixSpace of 6 by 6 sparse matrices over Rational Field
    sage: MS.base_ring()
    Rational Field
    sage: MS = MatrixSpace(ZZ, 3,5, sparse=False); MS
    Full MatrixSpace of 3 by 5 dense matrices over Integer Ring

TESTS::

    sage: matrix(RR,2,2,sparse=True)
    [0.000000000000000 0.000000000000000]
    [0.000000000000000 0.000000000000000]
    sage: matrix(GF(11), 2, 2, sparse=True)
    [0 0]
    [0 0]
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# System imports
import sys
import operator

# Sage matrix imports see :issue:`34283`

# Sage imports
import sage.structure.coerce_actions
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
import sage.rings.integer as integer
from sage.rings.finite_rings.finite_field_base import FiniteField
import sage.misc.latex as latex
import sage.modules.free_module

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.superseded import deprecated_function_alias
from sage.misc.persist import register_unpickle_override
from sage.categories.rings import Rings
from sage.categories.fields import Fields
from sage.categories.enumerated_sets import EnumeratedSets

from sage.misc.lazy_import import lazy_import
from sage.features.meataxe import Meataxe
lazy_import('sage.matrix.matrix_gfpn_dense', ['Matrix_gfpn_dense'],
            feature=Meataxe())
lazy_import('sage.groups.matrix_gps.matrix_group', ['MatrixGroup_base'])

_Rings = Rings()
_Fields = Fields()


def is_MatrixSpace(x):
    """
    Return whether ``self`` is an instance of ``MatrixSpace``.

    EXAMPLES::

        sage: from sage.matrix.matrix_space import is_MatrixSpace
        sage: MS = MatrixSpace(QQ,2)
        sage: A = MS.random_element()
        sage: is_MatrixSpace(MS)
        doctest:warning...
        DeprecationWarning: the function is_MatrixSpace is deprecated;
        use 'isinstance(..., MatrixSpace)' instead
        See https://github.com/sagemath/sage/issues/37924 for details.
        True
        sage: is_MatrixSpace(A)
        False
        sage: is_MatrixSpace(5)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37924, "the function is_MatrixSpace is deprecated; use 'isinstance(..., MatrixSpace)' instead")
    return isinstance(x, MatrixSpace)


def get_matrix_class(R, nrows, ncols, sparse, implementation):
    r"""
    Return a matrix class according to the input.

    .. NOTE::

        This returns the base class without the category.

    INPUT:

    - ``R`` -- a base ring

    - ``nrows`` -- number of rows

    - ``ncols`` -- number of columns

    - ``sparse`` -- boolean; whether the matrix class should be sparse

    - ``implementation`` -- ``None`` or string or a matrix class; a possible
      implementation. See the documentation of the constructor of :class:`MatrixSpace`.

    EXAMPLES::

        sage: from sage.matrix.matrix_space import get_matrix_class

        sage: get_matrix_class(ZZ, 4, 5, False, None)                                   # needs sage.libs.linbox
        <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        sage: get_matrix_class(ZZ, 4, 5, True, None)                                    # needs sage.libs.linbox
        <class 'sage.matrix.matrix_integer_sparse.Matrix_integer_sparse'>

        sage: get_matrix_class(ZZ, 3, 3, False, 'flint')                                # needs sage.libs.linbox
        <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        sage: get_matrix_class(ZZ, 3, 3, False, 'gap')                                  # needs sage.libs.gap
        <class 'sage.matrix.matrix_gap.Matrix_gap'>
        sage: get_matrix_class(ZZ, 3, 3, False, 'generic')
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>

        sage: get_matrix_class(GF(2^15), 3, 3, False, None)                             # needs sage.rings.finite_rings
        <class 'sage.matrix.matrix_gf2e_dense.Matrix_gf2e_dense'>
        sage: get_matrix_class(GF(2^17), 3, 3, False, None)                             # needs sage.rings.finite_rings
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>

        sage: get_matrix_class(GF(2), 2, 2, False, 'm4ri')                              # needs sage.libs.m4ri
        <class 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>
        sage: get_matrix_class(GF(4), 2, 2, False, 'm4ri')                              # needs sage.libs.m4ri sage.rings.finite_rings
        <class 'sage.matrix.matrix_gf2e_dense.Matrix_gf2e_dense'>
        sage: get_matrix_class(GF(7), 2, 2, False, 'linbox-float')                      # needs sage.libs.linbox
        <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
        sage: get_matrix_class(GF(7), 2, 2, False, 'linbox-double')                     # needs sage.libs.linbox
        <class 'sage.matrix.matrix_modn_dense_double.Matrix_modn_dense_double'>

        sage: get_matrix_class(RDF, 2, 2, False, 'numpy')                               # needs numpy
        <class 'sage.matrix.matrix_real_double_dense.Matrix_real_double_dense'>
        sage: get_matrix_class(CDF, 2, 3, False, 'numpy')                               # needs numpy sage.rings.complex_double
        <class 'sage.matrix.matrix_complex_double_dense.Matrix_complex_double_dense'>

        sage: get_matrix_class(GF(25,'x'), 4, 4, False, 'meataxe')          # optional - meataxe, needs sage.rings.finite_rings
        <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
        sage: get_matrix_class(IntegerModRing(3), 4, 4, False, 'meataxe')   # optional - meataxe
        <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
        sage: get_matrix_class(IntegerModRing(4), 4, 4, False, 'meataxe')
        Traceback (most recent call last):
        ...
        ValueError: 'meataxe' matrix can only deal with finite fields of order < 256
        sage: get_matrix_class(GF(next_prime(255)), 4, 4, False, 'meataxe')             # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: 'meataxe' matrix can only deal with finite fields of order < 256

        sage: get_matrix_class(ZZ, 3, 5, False, 'crazy_matrix')
        Traceback (most recent call last):
        ...
        ValueError: unknown matrix implementation 'crazy_matrix' over Integer Ring
        sage: get_matrix_class(GF(3), 2, 2, False, 'm4ri')
        Traceback (most recent call last):
        ...
        ValueError: 'm4ri' matrices are only available for fields of characteristic 2
        and order <= 65536
        sage: get_matrix_class(Zmod(2**30), 2, 2, False, 'linbox-float')                # needs sage.libs.linbox
        Traceback (most recent call last):
        ...
        ValueError: 'linbox-float' matrices can only deal with order < 256
        sage: get_matrix_class(Zmod(2**30), 2, 2, False, 'linbox-double')               # needs sage.libs.linbox
        Traceback (most recent call last):
        ...
        ValueError: 'linbox-double' matrices can only deal with order < 94906266

        sage: type(matrix(SR, 2, 2, 0))                                                 # needs sage.symbolic
        <class 'sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense'>
        sage: type(matrix(SR, 2, 2, 0, sparse=True))                                    # needs sage.symbolic
        <class 'sage.matrix.matrix_symbolic_sparse.Matrix_symbolic_sparse'>
        sage: type(matrix(GF(7), 2, range(4)))                                          # needs sage.libs.linbox
        <class 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
        sage: type(matrix(GF(16007), 2, range(4)))                                      # needs sage.libs.linbox
        <class 'sage.matrix.matrix_modn_dense_double.Matrix_modn_dense_double'>
        sage: type(matrix(CBF, 2, range(4)))                                            # needs sage.libs.flint
        <class 'sage.matrix.matrix_complex_ball_dense.Matrix_complex_ball_dense'>
        sage: type(matrix(GF(2), 2, range(4)))                                          # needs sage.libs.m4ri
        <class 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>
        sage: type(matrix(GF(64, 'z'), 2, range(4)))                                    # needs sage.libs.m4ri sage.rings.finite_rings
        <class 'sage.matrix.matrix_gf2e_dense.Matrix_gf2e_dense'>
        sage: type(matrix(GF(125, 'z'), 2, range(4)))                       # optional - meataxe, needs sage.rings.finite_rings
        <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
    """
    if isinstance(implementation, type):
        return implementation

    if not sparse:
        if implementation is None:
            # Choose default implementation:
            if R is sage.rings.integer_ring.ZZ:
                try:
                    from . import matrix_integer_dense
                except ImportError:
                    pass
                else:
                    return matrix_integer_dense.Matrix_integer_dense

            elif R is sage.rings.rational_field.QQ:
                try:
                    from . import matrix_rational_dense
                except ImportError:
                    pass
                else:
                    return matrix_rational_dense.Matrix_rational_dense

            elif isinstance(R, sage.rings.abc.RealDoubleField):
                try:
                    from . import matrix_real_double_dense
                except ImportError:
                    pass
                else:
                    return matrix_real_double_dense.Matrix_real_double_dense

            elif isinstance(R, sage.rings.abc.ComplexDoubleField):
                if implementation is None or implementation == 'numpy':
                    try:
                        from . import matrix_complex_double_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_complex_double_dense.Matrix_complex_double_dense

            elif isinstance(R, FiniteField):
                if R.order() == 2:
                    try:
                        from . import matrix_mod2_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_mod2_dense.Matrix_mod2_dense

                if R.characteristic() == 2 and R.order() <= 65536:  # 65536 == 2^16
                    try:
                        from . import matrix_gf2e_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_gf2e_dense.Matrix_gf2e_dense

                if (not R.is_prime_field()) and R.order() < 256:
                    try:
                        from . import matrix_gfpn_dense
                        return matrix_gfpn_dense.Matrix_gfpn_dense
                    except ImportError:
                        pass

            if isinstance(R, sage.rings.abc.IntegerModRing):
                try:
                    from . import matrix_modn_dense_double, matrix_modn_dense_float
                except ImportError:
                    pass
                else:
                    if R.order() < matrix_modn_dense_float.MAX_MODULUS:
                        return matrix_modn_dense_float.Matrix_modn_dense_float
                    if R.order() < matrix_modn_dense_double.MAX_MODULUS:
                        return matrix_modn_dense_double.Matrix_modn_dense_double

            if isinstance(R, sage.rings.abc.NumberField_cyclotomic):
                from . import matrix_cyclo_dense
                return matrix_cyclo_dense.Matrix_cyclo_dense

            try:
                from sage.symbolic.ring import SR
            except ImportError:
                pass
            else:
                if R is SR:
                    try:
                        from . import matrix_symbolic_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_symbolic_dense.Matrix_symbolic_dense

            if isinstance(R, sage.rings.abc.ComplexBallField):
                try:
                    from . import matrix_complex_ball_dense
                except ImportError:
                    pass
                else:
                    return matrix_complex_ball_dense.Matrix_complex_ball_dense

            try:
                from sage.rings.polynomial import polynomial_ring, multi_polynomial_ring_base
            except ImportError:
                pass
            else:
                if isinstance(R, polynomial_ring.PolynomialRing_generic) and R.base_ring() in _Fields:
                    try:
                        from . import matrix_polynomial_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_polynomial_dense.Matrix_polynomial_dense

                elif isinstance(R, multi_polynomial_ring_base.MPolynomialRing_base) and R.base_ring() in _Fields:
                    try:
                        from . import matrix_mpolynomial_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_mpolynomial_dense.Matrix_mpolynomial_dense
                elif isinstance(R, sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_mpair) and R.base_ring() in _Fields:
                    try:
                        from . import matrix_laurent_mpolynomial_dense
                    except ImportError:
                        pass
                    else:
                        return matrix_laurent_mpolynomial_dense.Matrix_laurent_mpolynomial_dense

            # The fallback
            from sage.matrix.matrix_generic_dense import Matrix_generic_dense
            return Matrix_generic_dense

        # Deal with request for a specific implementation
        if implementation == 'flint':
            if R is sage.rings.integer_ring.ZZ:
                from . import matrix_integer_dense
                return matrix_integer_dense.Matrix_integer_dense
            if R is sage.rings.rational_field.QQ:
                from . import matrix_rational_dense
                return matrix_rational_dense.Matrix_rational_dense
            raise ValueError("'flint' matrices are only available over the integers or the rationals")

        if implementation == 'm4ri':
            if R.is_field() and R.characteristic() == 2 and R.order() <= 65536:
                if R.order() == 2:
                    from . import matrix_mod2_dense
                    return matrix_mod2_dense.Matrix_mod2_dense
                from . import matrix_gf2e_dense
                return matrix_gf2e_dense.Matrix_gf2e_dense
            raise ValueError("'m4ri' matrices are only available for fields of characteristic 2 and order <= 65536")

        if implementation == 'meataxe':
            if R.is_field() and R.order() < 256:
                return Matrix_gfpn_dense
            raise ValueError("'meataxe' matrix can only deal with finite fields of order < 256")

        if implementation == 'numpy':
            if R is sage.rings.real_double.RDF:
                from . import matrix_real_double_dense
                return matrix_real_double_dense.Matrix_real_double_dense
            if R is sage.rings.complex_double.CDF:
                from . import matrix_complex_double_dense
                return matrix_complex_double_dense.Matrix_complex_double_dense
            if R is sage.rings.integer_ring.ZZ:
                from . import matrix_numpy_integer_dense
                return matrix_numpy_integer_dense.Matrix_numpy_integer_dense
            raise ValueError("'numpy' matrices are only available over RDF, CDF, and ZZ")

        if implementation == 'rational':
            if isinstance(R, sage.rings.abc.NumberField_cyclotomic):
                from . import matrix_cyclo_dense
                return matrix_cyclo_dense.Matrix_cyclo_dense
            raise ValueError("'rational' matrices are only available over a cyclotomic field")

        if implementation == 'linbox-float':
            from . import matrix_modn_dense_float
            if R.order() < matrix_modn_dense_float.MAX_MODULUS:
                return matrix_modn_dense_float.Matrix_modn_dense_float
            raise ValueError("'linbox-float' matrices can only deal with order < %s" % matrix_modn_dense_float.MAX_MODULUS)

        if implementation == 'linbox-double':
            from . import matrix_modn_dense_double
            if R.order() < matrix_modn_dense_double.MAX_MODULUS:
                return matrix_modn_dense_double.Matrix_modn_dense_double
            raise ValueError("'linbox-double' matrices can only deal with order < %s" % matrix_modn_dense_double.MAX_MODULUS)

        if implementation == 'generic':
            from sage.matrix.matrix_generic_dense import Matrix_generic_dense
            return Matrix_generic_dense

        if implementation == 'gap':
            from sage.matrix.matrix_gap import Matrix_gap
            return Matrix_gap

        raise ValueError("unknown matrix implementation %r over %r" % (implementation, R))

    # By now, we are dealing with sparse matrices
    if implementation is not None:
        raise ValueError("cannot choose an implementation for sparse matrices")

    if isinstance(R, sage.rings.abc.IntegerModRing):
        try:
            from . import matrix_modn_sparse
        except ImportError:
            pass
        else:
            if R.order() < matrix_modn_sparse.MAX_MODULUS:
                return matrix_modn_sparse.Matrix_modn_sparse

    if isinstance(R, sage.rings.rational_field.RationalField):
        try:
            from . import matrix_rational_sparse
        except ImportError:
            pass
        else:
            return matrix_rational_sparse.Matrix_rational_sparse

    if isinstance(R, sage.rings.integer_ring.IntegerRing_class):
        try:
            from . import matrix_integer_sparse
        except ImportError:
            pass
        else:
            return matrix_integer_sparse.Matrix_integer_sparse

    if isinstance(R, (sage.rings.abc.RealDoubleField, sage.rings.abc.ComplexDoubleField)):
        from . import matrix_double_sparse
        return matrix_double_sparse.Matrix_double_sparse
    try:
        from sage.symbolic.ring import SR
    except ImportError:
        pass
    else:
        if R is SR:
            try:
                from . import matrix_symbolic_sparse
            except ImportError:
                pass
            else:
                return matrix_symbolic_sparse.Matrix_symbolic_sparse

    # the fallback
    from sage.matrix.matrix_generic_sparse import Matrix_generic_sparse
    return Matrix_generic_sparse


class MatrixSpace(UniqueRepresentation, Parent):
    """
    The space of matrices of given size and base ring.

    INPUT:

    - ``base_ring`` -- a ring

    - ``nrows`` or ``row_keys`` -- nonnegative integer; the number of rows, or
      a finite family of arbitrary objects that index the rows of the matrix

    - ``ncols`` or ``column_keys`` -- nonnegative integer (default: ``nrows``);
      the number of columns, or a finite family of arbitrary objects that index
      the columns of the matrix

    - ``sparse`` -- boolean (default: ``False``); whether or not matrices
      are given a sparse representation

    - ``implementation`` -- (optional)  string or matrix class; a possible
      implementation. Depending on the base ring, the string can be

      - ``'generic'`` -- on any base rings

      - ``'flint'`` -- for integers and rationals

      - ``'meataxe'`` -- finite fields using the optional package :ref:`spkg_meataxe`

      - ``'m4ri'`` -- for characteristic 2 using the :ref:`spkg_m4ri` library

      - ``'linbox-float'`` -- for integer mod rings up to `2^8 = 256`

      - ``'linbox-double'`` -- for integer mod rings up to
        `floor(2^26*sqrt(2) + 1/2) = 94906266`

      - ``'numpy'`` -- for real and complex floating point numbers

    OUTPUT: a matrix space or, more generally, a homspace between free modules

    This factory function creates instances of various specialized classes
    depending on the input.  Not all combinations of options are
    implemented.

    - If the parameters ``row_keys`` or ``column_keys`` are provided, they
      must be finite families of objects. In this case, instances of
      :class:`CombinatorialFreeModule` are created via the factory function
      :func:`FreeModule`. Then the homspace between these modules is returned.

    EXAMPLES::

        sage: MatrixSpace(QQ, 2)
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: MatrixSpace(ZZ, 3, 2)
        Full MatrixSpace of 3 by 2 dense matrices over Integer Ring
        sage: MatrixSpace(ZZ, 3, sparse=False)
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

        sage: MatrixSpace(ZZ, 10, 5)
        Full MatrixSpace of 10 by 5 dense matrices over Integer Ring
        sage: MatrixSpace(ZZ, 10, 5).category()
        Category of infinite enumerated finite dimensional modules with basis over
         (Dedekind domains and euclidean domains
          and noetherian rings
          and infinite enumerated sets and metric spaces)
        sage: MatrixSpace(ZZ, 10, 10).category()
        Category of infinite enumerated finite dimensional algebras with basis over
         (Dedekind domains and euclidean domains
          and noetherian rings
          and infinite enumerated sets and metric spaces)
        sage: MatrixSpace(QQ, 10).category()
        Category of infinite finite dimensional algebras with basis over
         (number fields and quotient fields and metric spaces)

    Some examples of square 2 by 2 rational matrices::

        sage: MS = MatrixSpace(QQ, 2)
        sage: MS.dimension()
        4
        sage: MS.dims()
        (2, 2)
        sage: B = MS.basis()
        sage: list(B)
        [
        [1 0]  [0 1]  [0 0]  [0 0]
        [0 0], [0 0], [1 0], [0 1]
        ]
        sage: B[0,0]
        [1 0]
        [0 0]
        sage: B[0,1]
        [0 1]
        [0 0]
        sage: B[1,0]
        [0 0]
        [1 0]
        sage: B[1,1]
        [0 0]
        [0 1]
        sage: A = MS.matrix([1,2,3,4]); A
        [1 2]
        [3 4]

    The above matrix ``A`` can be multiplied by a 2 by 3 integer matrix::

        sage: MS2 = MatrixSpace(ZZ, 2, 3)
        sage: B = MS2.matrix([1,2,3,4,5,6])
        sage: A * B
        [ 9 12 15]
        [19 26 33]

    Using ``row_keys`` and ``column_keys``::

        sage: MS = MatrixSpace(ZZ, ['u', 'v'], ['a', 'b', 'c']); MS
        Set of Morphisms
         from Free module generated by {'a', 'b', 'c'} over Integer Ring
           to Free module generated by {'u', 'v'} over Integer Ring
           in Category of finite dimensional modules with basis over Integer Ring

    Check categories::

        sage: MatrixSpace(ZZ, 10, 5)
        Full MatrixSpace of 10 by 5 dense matrices over Integer Ring
        sage: MatrixSpace(ZZ, 10, 5).category()
        Category of infinite enumerated finite dimensional modules with basis over
         (Dedekind domains and euclidean domains
          and noetherian rings
          and infinite enumerated sets and metric spaces)
        sage: MatrixSpace(ZZ, 10, 10).category()
        Category of infinite enumerated finite dimensional algebras with basis over
         (Dedekind domains and euclidean domains
          and noetherian rings
          and infinite enumerated sets and metric spaces)
        sage: MatrixSpace(QQ, 10).category()
        Category of infinite finite dimensional algebras with basis over
         (number fields and quotient fields and metric spaces)

    TESTS::

        sage: MatrixSpace(ZZ, 1, 2^63)
        Traceback (most recent call last):
        ...
        OverflowError: number of rows and columns may be at most...
        sage: MatrixSpace(ZZ, 2^100, 10)
        Traceback (most recent call last):
        ...
        OverflowError: number of rows and columns may be at most...

    Check that different implementations play together as expected::

        sage: # needs sage.libs.linbox
        sage: M1 = MatrixSpace(ZZ, 2, implementation='flint')
        sage: M2 = MatrixSpace(ZZ, 2, implementation='generic')
        sage: type(M1(range(4)))
        <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        sage: type(M2(range(4)))
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        sage: M1(M2.an_element())
        [ 0  1]
        [-1  2]
        sage: M2(M1.an_element())
        [ 0  1]
        [-1  2]
        sage: all((A.get_action(B) is not None) == (A is B)
        ....:     for A in [M1, M2] for B in [M1, M2])
        True

    Check that libgap matrices over finite fields are working properly::

        sage: # needs sage.libs.gap
        sage: M2 = MatrixSpace(GF(2), 5, implementation='gap')
        sage: M2.one()
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        sage: m = M2.random_element()
        sage: M1 = MatrixSpace(GF(2), 5)
        sage: M1(m * m) == M1(m) * M1(m)
        True

    Check various combinations of dimensions and row/column keys::

        sage: M_ab_4 = MatrixSpace(QQ, ['a','b'], 4); M_ab_4
        Set of Morphisms (Linear Transformations)
         from Vector space of dimension 4 over Rational Field
           to Free module generated by {'a', 'b'} over Rational Field
        sage: TestSuite(M_ab_4).run()  # known bug
        sage: M_4_ab = MatrixSpace(QQ, 4, ['a','b']); M_4_ab
        Set of Morphisms
         from Free module generated by {'a', 'b'} over Rational Field
           to Vector space of dimension 4 over Rational Field
           in Category of finite dimensional vector spaces with basis
              over (number fields and quotient fields and metric spaces)
        sage: TestSuite(M_4_ab).run()  # known bug
        sage: M_ab_xy = MatrixSpace(QQ, ['a','b'], ['x','y'], nrows=2); M_ab_xy
        Set of Morphisms
         from Free module generated by {'x', 'y'} over Rational Field
           to Free module generated by {'a', 'b'} over Rational Field
           in Category of finite dimensional vector spaces with basis over Rational Field
        sage: TestSuite(M_ab_xy).run()  # known bug
        sage: MatrixSpace(QQ, ['a','b'], ['x','y'], nrows=4)
        Traceback (most recent call last):
        ...
        ValueError: inconsistent number of rows:
        should be cardinality of ['a', 'b'] but got 4
        sage: MatrixSpace(QQ, ['a','b'], ['x','y'], ncols=2)
        Set of Morphisms
         from Free module generated by {'x', 'y'} over Rational Field
           to Free module generated by {'a', 'b'} over Rational Field
           in Category of finite dimensional vector spaces with basis over Rational Field
        sage: MatrixSpace(QQ, ['a','b'], ['x','y'], ncols=4)
        Traceback (most recent call last):
        ...
        ValueError: inconsistent number of columns:
        should be cardinality of ['x', 'y'] but got 4
        sage: MatrixSpace(QQ, ['a','b'], ['x','y'], nrows=2, ncols=2)
        Set of Morphisms
         from Free module generated by {'x', 'y'} over Rational Field
           to Free module generated by {'a', 'b'} over Rational Field
           in Category of finite dimensional vector spaces with basis over Rational Field
        sage: MatrixSpace(QQ, ['a','b'], ['x','y'], nrows=2, ncols=4)
        Traceback (most recent call last):
        ...
        ValueError: inconsistent number of columns:
        should be cardinality of ['x', 'y'] but got 4
        sage: MatrixSpace(QQ, ['a','b'], ['x','y'], nrows=4, ncols=4)
        Traceback (most recent call last):
        ...
        ValueError: inconsistent number of columns:
        should be cardinality of ['x', 'y'] but got 4
        sage: MatrixSpace(QQ, 4, ['a','b'], nrows=4, ncols=2)
        Traceback (most recent call last):
        ...
        ValueError: duplicate values for nrows
    """

    @staticmethod
    def __classcall__(cls, base_ring,
                      nrows_or_row_keys=None, ncols_or_column_keys=None,
                      sparse=False, implementation=None, *,
                      nrows=None, ncols=None,
                      row_keys=None, column_keys=None,
                      **kwds):
        """
        Normalize the arguments to call the ``__init__`` constructor or delegate to another class.

        TESTS::

            sage: M1 = MatrixSpace(QQ, 2)
            sage: M2 = MatrixSpace(QQ, 2)
            sage: M1 is M2
            True
            sage: M3 = MatrixSpace(QQ, 2, implementation='flint')                       # needs sage.libs.flint
            sage: M1 is M3                                                              # needs sage.libs.flint
            True

        ::

            sage: M = MatrixSpace(ZZ, 10, implementation='flint')                       # needs sage.libs.flint
            sage: M                                                                     # needs sage.libs.flint
            Full MatrixSpace of 10 by 10 dense matrices over Integer Ring
            sage: loads(M.dumps()) is M                                                 # needs sage.libs.flint
            True

            sage: MatrixSpace(ZZ, 10, implementation='foobar')
            Traceback (most recent call last):
            ...
            ValueError: unknown matrix implementation 'foobar' over Integer Ring

        Check that :issue:`29466` is fixed::

            sage: class MyMatrixSpace(MatrixSpace):
            ....:     @staticmethod
            ....:     def __classcall__(cls, base_ring, nrows, ncols=None, my_option=True, sparse=False, implementation=None):
            ....:         return super().__classcall__(cls, base_ring, nrows, ncols=ncols, my_option=my_option, sparse=sparse, implementation=implementation)
            ....:
            ....:     def __init__(self, base_ring, nrows, ncols, sparse,  implementation, my_option=True):
            ....:         super().__init__(base_ring, nrows, ncols, sparse, implementation)
            ....:         self._my_option = my_option

            sage: MS1 = MyMatrixSpace(ZZ, 2)
            sage: MS1._my_option
            True
            sage: MS2 = MyMatrixSpace(ZZ, 2, my_option=False)
            sage: MS2._my_option
            False
        """
        if base_ring not in _Rings:
            raise TypeError("base_ring (=%s) must be a ring" % base_ring)

        if ncols_or_column_keys is not None:
            try:
                n = int(ncols_or_column_keys)
            except (TypeError, ValueError):
                if column_keys is not None:
                    raise ValueError("duplicate values for column_keys")
                column_keys = ncols_or_column_keys
            else:
                if ncols is not None:
                    raise ValueError("duplicate values for ncols")
                ncols = n
        if column_keys is not None and ncols is not None and ncols != len(column_keys):
            raise ValueError(f"inconsistent number of columns: should be cardinality of {column_keys} "
                             f"but got {ncols}")

        if nrows_or_row_keys is not None:
            try:
                n = int(nrows_or_row_keys)
            except (TypeError, ValueError):
                if row_keys is not None:
                    raise ValueError("duplicate values for row_keys")
                row_keys = nrows_or_row_keys
            else:
                if nrows is not None:
                    raise ValueError("duplicate values for nrows")
                nrows = n
        if row_keys is not None and nrows is not None and nrows != len(row_keys):
            raise ValueError(f"inconsistent number of rows: should be cardinality of {row_keys} "
                             f"but got {nrows}")

        if ncols is None and column_keys is None:
            ncols = nrows
            column_keys = row_keys

        sparse = bool(sparse)

        if row_keys is not None or column_keys is not None:
            from sage.categories.homset import Hom
            from sage.modules.free_module import FreeModule

            domain = FreeModule(base_ring, rank=ncols, basis_keys=column_keys,
                                sparse=sparse, **kwds)
            codomain = FreeModule(base_ring, rank=nrows, basis_keys=row_keys,
                                  sparse=sparse, **kwds)
            return Hom(domain, codomain)

        if nrows < 0:
            raise ArithmeticError("nrows must be nonnegative")
        if ncols < 0:
            raise ArithmeticError("ncols must be nonnegative")
        if nrows > sys.maxsize or ncols > sys.maxsize:
            raise OverflowError("number of rows and columns may be at most %s" % sys.maxsize)

        matrix_cls = get_matrix_class(base_ring, nrows, ncols, sparse, implementation)
        return super().__classcall__(cls, base_ring, nrows,
                                     ncols, sparse, matrix_cls, **kwds)

    def __init__(self, base_ring, nrows, ncols, sparse, implementation):
        r"""
        INPUT:

        - ``base_ring``

        - ``nrows`` -- positive integer; the number of rows

        - ``ncols`` -- positive integer (default: ``nrows``); the number of
          columns

        - ``sparse`` -- boolean (default: ``False``); whether or not matrices
          are given a sparse representation

        - ``implementation`` -- (optional) string or matrix class; a possible
          implementation. Depending on the base ring the string can be

           - ``'generic'`` -- on any base rings

           - ``'flint'`` -- for integers and rationals

           - ``'meataxe'`` -- finite fields, needs to install the optional package meataxe

           - ``m4ri`` -- for characteristic 2 using M4RI library

           - ``linbox-float`` -- for integer mod rings up to `2^8 = 256`

           - ``linbox-double`` -- for integer mod rings up to
             `floor(2^26*sqrt(2) + 1/2) = 94906266`

           - ``numpy`` -- for real and complex floating point numbers

        EXAMPLES::

            sage: MatrixSpace(QQ, 2)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: MatrixSpace(ZZ, 3, 2)
            Full MatrixSpace of 3 by 2 dense matrices over Integer Ring
            sage: MatrixSpace(ZZ, 3, sparse=False)
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

            sage: MatrixSpace(ZZ,10,5)
            Full MatrixSpace of 10 by 5 dense matrices over Integer Ring
            sage: MatrixSpace(ZZ,10,5).category()
            Category of infinite enumerated finite dimensional modules with basis over
             (Dedekind domains and euclidean domains
              and noetherian rings
              and infinite enumerated sets and metric spaces)
            sage: MatrixSpace(ZZ,10,10).category()
            Category of infinite enumerated finite dimensional algebras with basis over
             (Dedekind domains and euclidean domains
              and noetherian rings
              and infinite enumerated sets and metric spaces)
            sage: MatrixSpace(QQ,10).category()
            Category of infinite finite dimensional algebras with basis over
             (number fields and quotient fields and metric spaces)

        TESTS:

        We test that in the real or complex double dense case,
        conversion from the base ring is done by a call morphism.
        Note that by :issue:`9138`, other algebras usually
        get a conversion map by multiplication with the one element.
        ::

            sage: MS = MatrixSpace(RDF, 2, 2)
            sage: MS.convert_map_from(RDF)
            Coercion map:
              From: Real Double Field
              To:   Full MatrixSpace of 2 by 2 dense matrices over Real Double Field
            sage: MS = MatrixSpace(CDF, 2, 2)
            sage: MS.convert_map_from(CDF)
            Coercion map:
              From: Complex Double Field
              To:   Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field

        We check that :issue:`10095` is fixed::

            sage: M = Matrix(QQ, [[1 for dummy in range(125)]])
            sage: V = M.right_kernel()
            sage: V
            Vector space of degree 125 and dimension 124 over Rational Field
            Basis matrix:
            124 x 125 dense matrix over Rational Field
            sage: MatrixSpace(ZZ,20,20)(1).solve_right(MatrixSpace(ZZ,20,1).random_element())
            20 x 1 dense matrix over Rational Field (use the '.str()' method to see the entries)
            sage: MatrixSpace(ZZ,200,200)(1).solve_right(MatrixSpace(ZZ,200,1).random_element())
            200 x 1 dense matrix over Rational Field (use the '.str()' method to see the entries)
            sage: A = MatrixSpace(RDF,1000,1000).random_element()
            sage: B = MatrixSpace(RDF,1000,1000).random_element()

            sage: # needs numpy (otherwise timeout)
            sage: C = A * B

        We check that :issue:`18186` is fixed::

            sage: MatrixSpace(ZZ,0,3) in FiniteSets()
            True
            sage: MatrixSpace(Zmod(4),2) in FiniteSets()
            True
            sage: MatrixSpace(ZZ,2) in Sets().Infinite()
            True
        """
        # Checks of input data are supposed to be done in __classcall__
        assert isinstance(implementation, type)

        self.Element = implementation
        self.__nrows = nrows
        self.__ncols = ncols
        self.__is_sparse = sparse

        from sage.categories.modules import Modules
        from sage.categories.algebras import Algebras
        if nrows == ncols:
            category = Algebras(base_ring.category())
        else:
            category = Modules(base_ring.category())

        category = category.WithBasis().FiniteDimensional()

        if not self.__nrows or not self.__ncols:
            is_finite = True
        else:
            try:
                is_finite = base_ring.is_finite()
            except (AttributeError, NotImplementedError):
                is_finite = None

        if is_finite is True:
            category = category.Finite()
        elif is_finite is False:
            category = category.Infinite()

        if base_ring in EnumeratedSets():
            category = category.Enumerated()

        Parent.__init__(self, base_ring, category=category)

    def cardinality(self):
        r"""
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: MatrixSpace(GF(3), 2, 3).cardinality()
            729
            sage: MatrixSpace(ZZ, 2).cardinality()
            +Infinity
            sage: MatrixSpace(ZZ, 0, 3).cardinality()
            1
        """
        if not self.__nrows or not self.__ncols:
            from sage.rings.integer_ring import ZZ
            return ZZ.one()
        else:
            return self.base_ring().cardinality() ** (self.__nrows * self.__ncols)

    def characteristic(self):
        r"""
        Return the characteristic.

        EXAMPLES::

            sage: MatrixSpace(ZZ, 2).characteristic()
            0
            sage: MatrixSpace(GF(9), 0).characteristic()                                # needs sage.rings.finite_rings
            3
        """
        return self.base_ring().characteristic()

    def _has_default_implementation(self):
        r"""
        EXAMPLES::

            sage: MatrixSpace(ZZ, 2, implementation='generic')._has_default_implementation()
            False
            sage: MatrixSpace(ZZ, 2, implementation='flint')._has_default_implementation()          # needs sage.libs.linbox
            True
        """
        default = get_matrix_class(self.base_ring(), self.nrows(), self.ncols(), self.is_sparse(), None)
        return self.Element is default

    @lazy_attribute
    def transposed(self):
        """
        The transposed matrix space, having the same base ring and sparseness,
        but number of columns and rows is swapped.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(3), 7, 10)
            sage: MS.transposed
            Full MatrixSpace of 10 by 7 dense matrices over Finite Field of size 3
            sage: MS = MatrixSpace(GF(3), 7, 7)
            sage: MS.transposed is MS
            True

            sage: M = MatrixSpace(ZZ, 2, 3)
            sage: M.transposed
            Full MatrixSpace of 3 by 2 dense matrices over Integer Ring
        """
        return MatrixSpace(self._base, self.__ncols, self.__nrows,
                self.__is_sparse, self.Element)

    @lazy_attribute
    def _copy_zero(self):
        """
        Is it faster to copy a zero matrix or is it faster to create a
        new matrix from scratch?

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2), 20, 20)
            sage: MS._copy_zero
            False

            sage: MS = MatrixSpace(GF(3), 20, 20)
            sage: MS._copy_zero
            True
            sage: MS = MatrixSpace(GF(3), 200, 200)
            sage: MS._copy_zero
            False

            sage: MS = MatrixSpace(ZZ,200,200)
            sage: MS._copy_zero
            False
            sage: MS = MatrixSpace(ZZ,30,30)
            sage: MS._copy_zero
            True

            sage: MS = MatrixSpace(QQ,200,200)
            sage: MS._copy_zero
            False
            sage: MS = MatrixSpace(QQ,20,20)
            sage: MS._copy_zero
            False
        """
        if self.__is_sparse:
            return False
        elif self.Element is sage.matrix.matrix_mod2_dense.Matrix_mod2_dense:
            return False
        elif self.Element is sage.matrix.matrix_rational_dense.Matrix_rational_dense:
            return False
        elif self.__nrows > 40 and self.__ncols > 40:
            return False
        else:
            return True

    def _element_constructor_(self, entries, **kwds):
        """
        Construct an element of ``self`` from ``entries``.

        EXAMPLES::

            sage: k = GF(7)
            sage: G = MatrixGroup([matrix(k, 2, [1,1,0,1]), matrix(k, 2, [1,0,0,2])])
            sage: g = G.0
            sage: MatrixSpace(k, 2)(g)
            [1 1]
            [0 1]

        ::

            sage: MS = MatrixSpace(ZZ,2,4)
            sage: M2 = MS(range(8)); M2
            [0 1 2 3]
            [4 5 6 7]
            sage: M2 == MS(M2.rows())
            True

        ::

            sage: MS = MatrixSpace(ZZ,2,4, sparse=True)
            sage: M2 = MS(range(8)); M2
            [0 1 2 3]
            [4 5 6 7]
            sage: M2 == MS(M2.rows())
            True

        ::

            sage: MS = MatrixSpace(ZZ,2,2, sparse=True)
            sage: MS([1,2,3,4])
            [1 2]
            [3 4]

            sage: # needs sage.modular
            sage: MS = MatrixSpace(ZZ, 2)
            sage: g = Gamma0(5)([1,1,0,1])
            sage: MS(g)
            [1 1]
            [0 1]

        ::

            sage: MS = MatrixSpace(ZZ,2,2, sparse=True)
            sage: mat = MS(); mat
            [0 0]
            [0 0]
            sage: mat.is_mutable()
            True
            sage: mat2 = mat.change_ring(QQ); mat2.is_mutable()
            True

        TESTS:

        Ensure that :issue:`12020` is fixed::

            sage: rings = [ZZ, QQ, RDF]
            sage: rings.extend([RealField(100), ComplexField(100)])                     # needs sage.rings.real_mpfr
            sage: rings.append(CDF)                                                     # needs sage.rings.complex_double
            sage: rings.append(PolynomialRing(QQ, 'x'))
            sage: rings.append(PolynomialRing(CC, 2, 'x'))                              # needs sage.rings.real_mpfr
            sage: rings.append(SR)                                                      # needs sage.symbolic
            sage: rings.extend([GF(2), GF(11)])
            sage: rings.extend([GF(2^8,'a'), GF(3^19,'a')])                             # needs sage.rings.finite_rings
            sage: x = polygen(QQ)
            sage: rings.extend([NumberField(x^3 + 2, 'a'), CyclotomicField(4)])         # needs sage.rings.number_field
            sage: for R in rings:
            ....:     A = MatrixSpace(R, 60, 30, sparse=False)(0)
            ....:     B = A.augment(A)
            ....:     A = MatrixSpace(R, 60, 30, sparse=True)(0)
            ....:     B = A.augment(A)

        Check that :issue:`13012` is fixed::

            sage: m = zero_matrix(2, 3)
            sage: m
            [0 0 0]
            [0 0 0]
            sage: M = MatrixSpace(ZZ, 3, 5)
            sage: M.zero()
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            sage: M(m)
            Traceback (most recent call last):
            ...
            ValueError: inconsistent number of rows: should be 3 but got 2
            sage: M.matrix(m)
            Traceback (most recent call last):
            ...
            ValueError: inconsistent number of rows: should be 3 but got 2

        Check that :issue:`15110` is fixed::

            sage: S.<t> = LaurentSeriesRing(ZZ)
            sage: MS = MatrixSpace(S,1,1)
            sage: MS([[t]])   # given as a list of lists
            [t]
            sage: MS([t])     # given as a list of coefficients
            [t]
            sage: MS(t)       # given as a scalar matrix
            [t]

        Calling a matrix space `M` with a matrix in `M` as argument
        returns the original matrix unless ``copy=True`` is specified
        (:issue:`31078`)::

            sage: m = Matrix([[0, 1], [2, 3]])
            sage: M = m.parent()
            sage: M(m) is m
            True
            sage: M(m, copy=True) is m
            False
        """
        return self.element_class(self, entries, **kwds)

    def change_ring(self, R):
        """
        Return matrix space over R with otherwise same parameters as ``self``.

        INPUT:

        - ``R`` -- ring

        OUTPUT: a matrix space

        EXAMPLES::

            sage: Mat(QQ, 3, 5).change_ring(GF(7))
            Full MatrixSpace of 3 by 5 dense matrices
             over Finite Field of size 7
        """
        try:
            return self.__change_ring[R]
        except AttributeError:
            self.__change_ring = {}
        except KeyError:
            pass
        M = MatrixSpace(R, self.__nrows, self.__ncols, self.__is_sparse)
        self.__change_ring[R] = M
        return M

    def base_extend(self, R):
        """
        Return base extension of this matrix space to R.

        INPUT:

        - ``R`` -- ring

        OUTPUT: a matrix space

        EXAMPLES::

            sage: Mat(ZZ, 3, 5).base_extend(QQ)
            Full MatrixSpace of 3 by 5 dense matrices over Rational Field
            sage: Mat(QQ, 3, 5).base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        raise TypeError("no base extension defined")

    def construction(self):
        """
        EXAMPLES::

            sage: A = matrix(ZZ, 2, [1..4], sparse=True)
            sage: A.parent().construction()
            (MatrixFunctor, Integer Ring)
            sage: A.parent().construction()[0](QQ['x'])
            Full MatrixSpace of 2 by 2 sparse matrices over
             Univariate Polynomial Ring in x over Rational Field
            sage: parent(A/2)
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        """
        from sage.categories.pushout import MatrixFunctor
        return MatrixFunctor(self.__nrows, self.__ncols, is_sparse=self.is_sparse()), self.base_ring()

    def _get_action_(self, S, op, self_on_left):
        r"""
        Return the action of S on ``self``.

        INPUT:

        - ``S`` -- a parent

        - ``op`` -- an operator

        - ``self_on_left`` -- whether the operation is on left or on right

        EXAMPLES::

            sage: V = QQ^(2,3)
            sage: W1 = QQ^(3,4); W2 = QQ^(2,2)
            sage: V.get_action(W1, operator.mul)
            Left action by Full MatrixSpace of 2 by 3 dense matrices over Rational Field on Full MatrixSpace of 3 by 4 dense matrices over Rational Field
            sage: V.get_action(W2, operator.mul)
            sage: V.get_action(W1, operator.mul, self_on_left=False)
            sage: V.get_action(W2, operator.mul, self_on_left=False)
            Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field on Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        ::

            sage: V2 = QQ^2; V3 = QQ^3
            sage: V.get_action(V3, operator.mul)
            Left action by Full MatrixSpace of 2 by 3 dense matrices over Rational Field on Vector space of dimension 3 over Rational Field
            sage: V.get_action(V2, operator.mul)
            sage: V.get_action(V3, operator.mul, self_on_left=False)
            sage: V.get_action(V2, operator.mul, self_on_left=False)
            Right action by Full MatrixSpace of 2 by 3 dense matrices over Rational Field on Vector space of dimension 2 over Rational Field

        ::

            sage: V.get_action(ZZ, operator.mul)
            Right scalar multiplication by Integer Ring on Full MatrixSpace of 2 by 3 dense matrices over Rational Field
            sage: V.get_action(ZZ, operator.mul, self_on_left=False)
            Left scalar multiplication by Integer Ring on Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        """
        try:
            try:
                from sage.schemes.generic.homset import SchemeHomset_generic
                from sage.schemes.generic.homset import SchemeHomset_points
            except ImportError:
                SchemeHomset_generic = SchemeHomset_points = None
            if op is operator.mul:
                from . import action as matrix_action
                if self_on_left:
                    if isinstance(S, MatrixSpace):
                        # matrix multiplications
                        return matrix_action.MatrixMatrixAction(self, S)
                    elif isinstance(S, sage.modules.free_module.FreeModule_generic):
                        return matrix_action.MatrixVectorAction(self, S)
                    elif isinstance(S, SchemeHomset_points):
                        return matrix_action.MatrixSchemePointAction(self, S)
                    elif isinstance(S, SchemeHomset_generic):
                        return matrix_action.MatrixPolymapAction(self, S)
                    else:
                        # action of base ring
                        return sage.structure.coerce_actions.RightModuleAction(S, self)
                else:
                    if isinstance(S, MatrixSpace):
                        # matrix multiplications
                        return matrix_action.MatrixMatrixAction(S, self)
                    elif isinstance(S, sage.modules.free_module.FreeModule_generic):
                        return matrix_action.VectorMatrixAction(self, S)
                    elif isinstance(S, SchemeHomset_generic):
                        return matrix_action.PolymapMatrixAction(self, S)
                    else:
                        # action of base ring
                        return sage.structure.coerce_actions.LeftModuleAction(S, self)
        except TypeError:
            return None

    def _coerce_map_from_base_ring(self):
        """
        Return a coercion map from the base ring of ``self``.

        .. NOTE::

            This is only called for algebras of square matrices.

        EXAMPLES::

            sage: MS1 = MatrixSpace(QQ, 3)
            sage: MS1.coerce_map_from(QQ)
            Coercion map:
              From: Rational Field
              To:   Full MatrixSpace of 3 by 3 dense matrices over Rational Field
            sage: MS1.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Full MatrixSpace of 3 by 3 dense matrices over Rational Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Coercion map:
                      From: Rational Field
                      To:   Full MatrixSpace of 3 by 3 dense matrices over Rational Field

            sage: MS2 = MatrixSpace(ZZ, 3)
            sage: MS2.coerce_map_from(QQ)
            sage: MS2.coerce_map_from(ZZ)
            Coercion map:
              From: Integer Ring
              To:   Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

            sage: MatrixSpace(QQ, 1, 3).coerce_map_from(QQ)
        """
        return self._generic_coerce_map(self.base_ring())

    def _coerce_map_from_(self, S):
        r"""
        Canonical coercion from ``S`` to this matrix space.

        EXAMPLES::

            sage: MS1 = MatrixSpace(QQ, 3)
            sage: MS2 = MatrixSpace(ZZ, 3)
            sage: MS1.coerce_map_from(MS2)
            Coercion map:
              From: Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
              To:   Full MatrixSpace of 3 by 3 dense matrices over Rational Field
            sage: MS2.coerce_map_from(MS1)

        There are also coercions possible from matrix group and
        arithmetic subgroups::

            sage: MS = MatrixSpace(GF(3), 2, 2)
            sage: MS.coerce_map_from(GL(2, 3))
            Coercion map:
              From: General Linear Group of degree 2 over Finite Field of size 3
              To:   Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 3
            sage: MS.coerce_map_from(GL(2, 2))
            sage: MS.coerce_map_from(Gamma1(5))                                         # needs sage.modular
            Coercion map:
              From: Congruence Subgroup Gamma1(5)
              To:   Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 3

        TESTS:

        Check that :issue:`22091` is fixed::

            sage: A = Zmod(4)
            sage: R = MatrixSpace(A, 2)
            sage: G = GL(2, A)
            sage: R.coerce_map_from(G)
            Coercion map:
              From: General Linear Group of degree 2 over Ring of integers modulo 4
              To:   Full MatrixSpace of 2 by 2 dense matrices over Ring of integers modulo 4
            sage: R.coerce_map_from(GL(2, ZZ))
            Coercion map:
              From: General Linear Group of degree 2 over Integer Ring
              To:   Full MatrixSpace of 2 by 2 dense matrices over Ring of integers modulo 4

            sage: m = R([[1, 0], [0, 1]])
            sage: m in G
            True
            sage: m in list(G)                                                          # needs sage.libs.gap
            True
            sage: m == G(m)
            True

            sage: G = SL(3, QQ)
            sage: M = MatrixSpace(QQ, 3)
            sage: G.one() == M.identity_matrix()
            True
            sage: G.one() + M.identity_matrix()
            [2 0 0]
            [0 2 0]
            [0 0 2]

        Verify which coercion maps are allowed (this should form a
        poset)::

            sage: S = []
            sage: S += [MatrixSpace(ZZ, 3, implementation='flint')]                     # needs sage.libs.linbox
            sage: S += [MatrixSpace(ZZ, 3, implementation='generic')]
            sage: S += [MatrixSpace(ZZ, 3, implementation='gap')]                       # needs sage.libs.gap
            sage: S += [MatrixSpace(ZZ, 3, sparse=True)]
            sage: mult = ''
            sage: for A in S:
            ....:     for B in S:
            ....:         if A.has_coerce_map_from(B):
            ....:             mult += 'X'
            ....:         else:
            ....:             mult += ' '
            ....:     mult += '\n'
            sage: print(mult)                                                           # needs sage.libs.linbox sage.libs.gap
            XXXX
             X X
              XX
               X

        Thanks to the coercion model, arithmetic is allowed between all
        these parents::

            sage: for A in S:
            ....:     for B in S:
            ....:         a = A.an_element()
            ....:         b = B.an_element()
            ....:         dummy = (a * b) + (a - b)
        """
        B = self.base()

        if isinstance(S, MatrixSpace):
            # Disallow coercion if dimensions do not match
            if self.nrows() != S.nrows() or self.ncols() != S.ncols():
                return False
            T = S.base()
            if B is not T:
                # Matrix spaces over different base rings.
                # TODO: make this an actual map induced by the map
                # on the bases, see Issue #25540
                return B.has_coerce_map_from(T)

            # Base ring and dimensions are the same. So the only
            # difference can be the implementation and sparseness.
            if self.is_sparse() != S.is_sparse():
                # Allow coercion sparse -> dense
                return S.is_sparse()

            # Allow coercion to the default implementation.
            # As a consequence, the default implementation is considered
            # the "common parent" when mixing implementations.
            return self._has_default_implementation()

        # Check for other parents whose elements are some kind of matrices
        try:
            meth_matrix_space = S.matrix_space
        except AttributeError:
            pass
        else:
            MS = meth_matrix_space()
            if isinstance(S, MatrixGroup_base):
                return self.has_coerce_map_from(MS)

            try:
                from sage.modular.arithgroup.arithgroup_generic import ArithmeticSubgroup
            except ImportError:
                pass
            else:
                if isinstance(S, ArithmeticSubgroup):
                    return self.has_coerce_map_from(MS)

            return False

        # The parent is not matrix-like: coerce via base ring
        return (self.nrows() == self.ncols()) and self._coerce_map_via([B], S)

    def _repr_(self):
        """
        Return the string representation of a MatrixSpace.

        EXAMPLES::

            sage: MS = MatrixSpace(ZZ,2,4,true)
            sage: repr(MS)
            'Full MatrixSpace of 2 by 4 sparse matrices over Integer Ring'
            sage: MS
            Full MatrixSpace of 2 by 4 sparse matrices over Integer Ring

            sage: MatrixSpace(ZZ, 2, implementation='flint')                            # needs sage.libs.linbox
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: MatrixSpace(ZZ, 2, implementation='generic')
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring (using Matrix_generic_dense)
        """
        if self.is_sparse():
            s = "sparse"
        else:
            s = "dense"
        s = "Full MatrixSpace of %s by %s %s matrices over %s" % (
            self.__nrows, self.__ncols, s, self.base_ring())

        if not self._has_default_implementation():
            s += " (using {})".format(self.Element.__name__)

        return s

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: MS = MatrixSpace(ZZ,2,4,true)
            sage: MS._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return self.__nrows > 1
        return super()._repr_option(key)

    def _latex_(self):
        r"""
        Return the latex representation of a MatrixSpace.

        EXAMPLES::

            sage: MS3 = MatrixSpace(QQ,6,6,true)
            sage: latex(MS3)
            \mathrm{Mat}_{6\times 6}(\Bold{Q})
        """
        return "\\mathrm{Mat}_{%s\\times %s}(%s)" % (self.nrows(), self.ncols(),
                                                     latex.latex(self.base_ring()))

    def __len__(self):
        """
        Return number of elements of this matrix space if it fits in
        an int; raise a :exc:`TypeError` if there are infinitely many
        elements, and raise an :exc:`OverflowError` if there are finitely
        many but more than the size of an int.

        EXAMPLES::

            sage: len(MatrixSpace(GF(3), 3, 2))
            729
            sage: len(MatrixSpace(GF(3), 2, 3))
            729
            sage: 3^(2*3)
            729

            sage: len(MatrixSpace(GF(2003), 3, 2))                                      # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            OverflowError: cannot fit 'int' into an index-sized integer

            sage: len(MatrixSpace(QQ,3,2))
            Traceback (most recent call last):
            ...
            TypeError: len() of unsized object
        """
        return len(self.base_ring())**(self.nrows() * self.ncols())

    def __iter__(self):
        r"""
        Return a generator object which iterates through the elements of
        ``self``. The order in which the elements are generated is based on a
        'weight' of a matrix which is the number of iterations on the base
        ring that are required to reach that matrix.

        The ordering is similar to a degree negative lexicographic order in
        monomials in a multivariate polynomial ring.

        EXAMPLES: Consider the case of 2 x 2 matrices over GF(5).

        ::

            sage: list(GF(5))
            [0, 1, 2, 3, 4]
            sage: MS = MatrixSpace(GF(5), 2, 2)
            sage: l = list(MS)

        Then, consider the following matrices::

            sage: A = MS([2,1,0,1]); A
            [2 1]
            [0 1]
            sage: B = MS([1,2,1,0]); B
            [1 2]
            [1 0]
            sage: C = MS([1,2,0,0]); C
            [1 2]
            [0 0]

        A appears before B since the weight of one of A's entries exceeds
        the weight of the corresponding entry in B earliest in the list.

        ::

            sage: l.index(A)
            41
            sage: l.index(B)
            46

        However, A would come after the matrix C since C has a lower weight
        than A.

        ::

            sage: l.index(A)
            41
            sage: l.index(C)
            19

        The weights of matrices over other base rings are not as obvious.
        For example, the weight of

        ::

            sage: MS = MatrixSpace(ZZ, 2, 2)
            sage: MS([-1,0,0,0])
            [-1  0]
            [ 0  0]

        is 2 since

        ::

            sage: i = iter(ZZ)
            sage: next(i)
            0
            sage: next(i)
            1
            sage: next(i)
            -1

        Some more examples::

            sage: MS = MatrixSpace(GF(2), 2)
            sage: a = list(MS)
            sage: len(a)
            16
            sage: for m in a:
            ....:     print(m)
            ....:     print('-')
            [0 0]
            [0 0]
            -
            [1 0]
            [0 0]
            -
            [0 1]
            [0 0]
            -
            [0 0]
            [1 0]
            -
            [0 0]
            [0 1]
            -
            [1 1]
            [0 0]
            -
            [1 0]
            [1 0]
            -
            [1 0]
            [0 1]
            -
            [0 1]
            [1 0]
            -
            [0 1]
            [0 1]
            -
            [0 0]
            [1 1]
            -
            [1 1]
            [1 0]
            -
            [1 1]
            [0 1]
            -
            [1 0]
            [1 1]
            -
            [0 1]
            [1 1]
            -
            [1 1]
            [1 1]
            -

        ::

            sage: MS = MatrixSpace(GF(2), 2, 3)
            sage: a = list(MS)
            sage: len(a)
            64
            sage: a[0]
            [0 0 0]
            [0 0 0]

        ::

            sage: MS = MatrixSpace(ZZ, 2, 3)
            sage: i = iter(MS)
            sage: a = [ next(i) for _ in range(6) ]
            sage: a[0]
            [0 0 0]
            [0 0 0]
            sage: a[4]
            [0 0 0]
            [1 0 0]

        For degenerate cases, where either the number of rows or columns
        (or both) are zero, then the single element of the space is
        returned.

        ::

            sage: list(MatrixSpace(GF(2), 2, 0))
            [[]]
            sage: list(MatrixSpace(GF(2), 0, 2))
            [[]]
            sage: list(MatrixSpace(GF(2), 0, 0))
            [[]]

        If the base ring does not support iteration (for example, with the
        reals), then the matrix space over that ring does not support
        iteration either.

        ::

            sage: MS = MatrixSpace(RR, 2)
            sage: a = list(MS)
            Traceback (most recent call last):
            ...
            NotImplementedError: len() of an infinite set
        """
        # Make sure that we can iterate over the base ring
        base_ring = self.base_ring()
        base_iter = iter(base_ring)

        number_of_entries = (self.__nrows * self.__ncols)

        # If the number of entries is zero, then just
        # yield the empty matrix in that case and return
        if number_of_entries == 0:
            yield self(0)
            return

        import sage.combinat.integer_vector

        if not base_ring.is_finite():
            # When the base ring is not finite, then we should go
            # through and yield the matrices by "weight", which is
            # the total number of iterations that need to be done
            # on the base ring to reach the matrix.
            base_elements = [next(base_iter)]
            weight = 0
            while True:
                for iv in sage.combinat.integer_vector.IntegerVectors(weight, number_of_entries):
                    yield self([base_elements[i] for i in iv])
                weight += 1
                base_elements.append(next(base_iter))
        else:
            # In the finite case, we do a similar thing except that
            # the "weight" of each entry is bounded by the number
            # of elements in the base ring
            order = base_ring.order()
            base_elements = list(base_ring)
            for weight in range((order - 1) * number_of_entries + 1):
                for iv in sage.combinat.integer_vector.IntegerVectors(weight, number_of_entries, max_part=(order - 1)):
                    yield self([base_elements[i] for i in iv])

    def __getitem__(self, x):
        """
        Return a polynomial ring over this ring or the `n`-th element of this ring.

        This method implements the syntax ``R['x']`` to define polynomial rings
        over matrix rings, while still allowing to get the `n`-th element of a
        finite matrix ring with ``R[n]`` for backward compatibility.

        (If this behaviour proves desirable for all finite enumerated rings, it
        should eventually be implemented in the corresponding category rather
        than here.)

        .. SEEALSO::

            :meth:`sage.categories.rings.Rings.ParentMethod.__getitem__`,
            :meth:`sage.structure.parent.Parent.__getitem__`

        EXAMPLES::

            sage: MS = MatrixSpace(GF(3), 2, 2)
            sage: MS['x']
            Univariate Polynomial Ring in x
             over Full MatrixSpace of 2 by 2 dense matrices
              over Finite Field of size 3
            sage: MS[0]
            [0 0]
            [0 0]
            sage: MS[9]
            [0 2]
            [0 0]

            sage: MS = MatrixSpace(QQ, 7)
            sage: MS['x']
            Univariate Polynomial Ring in x over Full MatrixSpace of 7 by 7 dense matrices over Rational Field
            sage: MS[2]
            Traceback (most recent call last):
            ...
            AttributeError: 'MatrixSpace_with_category' object has no attribute 'list'...
        """
        if isinstance(x, (integer.Integer, int)):
            return self.list()[x]
        return super().__getitem__(x)

    def basis(self):
        """
        Return a basis for this matrix space.

        .. WARNING::

            This will of course compute every generator of this matrix
            space. So for large dimensions, this could take a long time,
            waste a massive amount of memory (for dense matrices), and
            is likely not very useful. Don't use this on large matrix
            spaces.

        EXAMPLES::

            sage: list(Mat(ZZ,2,2).basis())
            [
            [1 0]  [0 1]  [0 0]  [0 0]
            [0 0], [0 0], [1 0], [0 1]
            ]
        """
        v = {(r, c): self.zero_matrix().__copy__()
             for r in range(self.__nrows)
             for c in range(self.__ncols)}
        one = self.base_ring().one()
        keys = []
        for r in range(self.__nrows):
            for c in range(self.__ncols):
                keys.append((r, c))
                v[r, c][r, c] = one
                v[r, c].set_immutable()
        from sage.sets.family import Family
        return Family(keys, v.__getitem__)

    def dimension(self):
        r"""
        Return (m rows) \* (n cols) of ``self`` as ``Integer``.

        EXAMPLES::

            sage: MS = MatrixSpace(ZZ,4,6)
            sage: u = MS.dimension()
            sage: u - 24 == 0
            True
        """
        return self.__nrows * self.__ncols

    def dims(self):
        """
        Return (m row, n col) representation of ``self`` dimension.

        EXAMPLES::

            sage: MS = MatrixSpace(ZZ,4,6)
            sage: MS.dims()
            (4, 6)
        """
        return (self.__nrows, self.__ncols)

    def submodule(self, gens, check=True, already_echelonized=False,
                  unitriangular=False, support_order=None, category=None,
                  *args, **opts):
        r"""
        The submodule spanned by a finite set of matrices.

        INPUT:

        - ``gens`` -- list or family of elements of ``self``

        - ``check`` -- boolean (default: ``True``); whether to verify that the
          elements of ``gens`` are in ``self``

        - ``already_echelonized`` -- boolean (default: ``False``); whether
          the elements of ``gens`` are already in (not necessarily
          reduced) echelon form

        - ``unitriangular`` -- boolean (default: ``False``); whether
          the lift morphism is unitriangular

        - ``support_order`` -- (optional) either something that can
          be converted into a tuple or a key function

        If ``already_echelonized`` is ``False``, then the
        generators are put in reduced echelon form using
        :meth:`echelonize`, and reindexed by `0, 1, \ldots`.

        .. WARNING::

            At this point, this method only works for finite
            dimensional submodules and if matrices can be
            echelonized over the base ring.

        If in addition ``unitriangular`` is ``True``, then
        the generators are made such that the coefficients of
        the pivots are 1, so that lifting map is unitriangular.

        The basis of the submodule uses the same index set as the
        generators, and the lifting map sends `y_i` to `gens[i]`.

        .. SEEALSO::

             :meth:`ModulesWithBasis.ParentMethods.submodule`

        EXAMPLES::

            sage: M = MatrixSpace(QQ, 2)
            sage: mat = M.matrix([[1, 2], [3, 4]])
            sage: X = M.submodule([mat], already_echelonized=True); X
            Free module generated by {0} over Rational Field

            sage: mat2 = M.matrix([[1, 0], [-3, 2]])
            sage: X = M.submodule([mat, mat2])
            sage: [X.lift(b) for b in X.basis()]
            [
            [ 1  0]  [0 1]
            [-3  2], [3 1]
            ]

            sage: A = matrix([[1, 1], [0, -1]])
            sage: B = matrix([[0, 1], [0, 2]])
            sage: X = M.submodule([A, B])
            sage: Xp = M.submodule([A, B], support_order=[(0,1), (1,1), (0,0)])
            sage: [X.lift(b) for b in X.basis()]
            [
            [ 1  0]  [0 1]
            [ 0 -3], [0 2]
            ]
            sage: [Xp.lift(b) for b in Xp.basis()]
            [
            [2/3   1]  [-1/3    0]
            [  0   0], [   0    1]
            ]
        """
        support_order = self._compute_support_order(gens, support_order)
        if not already_echelonized:
            gens = self.echelon_form(gens, unitriangular, order=support_order)
        else:
            from copy import copy
            # We will be making gens immutable, so copy the mutable matrices
            gens = [copy(g) if g.is_mutable() else g for g in gens]

        # We need to make sure the result immutable
        for g in gens:
            g.set_immutable()

        from sage.modules.with_basis.subquotient import SubmoduleWithBasis
        return SubmoduleWithBasis(gens, ambient=self,
                                  support_order=support_order,
                                  unitriangular=unitriangular,
                                  category=category, *args, **opts)

    from sage.misc.cachefunc import cached_method

    @cached_method
    def identity_matrix(self):
        """
        Return the identity matrix in ``self``.

        ``self`` must be a space of square
        matrices. The returned matrix is immutable. Please use ``copy`` if
        you want a modified copy.

        EXAMPLES::

            sage: MS1 = MatrixSpace(ZZ,4)
            sage: MS2 = MatrixSpace(QQ,3,4)
            sage: I = MS1.identity_matrix()
            sage: I
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: Er = MS2.identity_matrix()
            Traceback (most recent call last):
            ...
            TypeError: identity matrix must be square

        TESTS::

            sage: MS1.one()[1,2] = 3
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).

        Check different implementations::

            sage: M1 = MatrixSpace(ZZ, 2, implementation='flint')                       # needs sage.libs.linbox
            sage: M2 = MatrixSpace(ZZ, 2, implementation='generic')

            sage: type(M1.identity_matrix())                                            # needs sage.libs.linbox
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: type(M2.identity_matrix())
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        """
        if self.__nrows != self.__ncols:
            raise TypeError("identity matrix must be square")
        A = self.zero_matrix().__copy__()
        for i in range(self.__nrows):
            A[i, i] = 1
        A.set_immutable()
        return A

    one = identity_matrix

    def diagonal_matrix(self, entries):
        """
        Create a diagonal matrix in ``self`` using the specified elements.

        INPUT:

        - ``entries`` -- the elements to use as the diagonal entries

        ``self`` must be a space of square matrices. The length of
        ``entries`` must be less than or equal to the matrix
        dimensions. If the length of ``entries`` is less than the
        matrix dimensions, ``entries`` is padded with zeroes at the
        end.

        EXAMPLES::

            sage: MS1 = MatrixSpace(ZZ,4)
            sage: MS2 = MatrixSpace(QQ,3,4)
            sage: I = MS1.diagonal_matrix([1, 2, 3, 4])
            sage: I
            [1 0 0 0]
            [0 2 0 0]
            [0 0 3 0]
            [0 0 0 4]
            sage: MS2.diagonal_matrix([1, 2])
            Traceback (most recent call last):
            ...
            TypeError: diagonal matrix must be square
            sage: MS1.diagonal_matrix([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            ValueError: number of diagonal matrix entries (5) exceeds the matrix size (4)
            sage: MS1.diagonal_matrix([1/2, 2, 3, 4])
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        Check different implementations::

            sage: M1 = MatrixSpace(ZZ, 2, implementation='flint')                       # needs sage.libs.linbox
            sage: M2 = MatrixSpace(ZZ, 2, implementation='generic')

            sage: type(M1.diagonal_matrix([1, 2]))                                      # needs sage.libs.linbox
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: type(M2.diagonal_matrix([1, 2]))
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        """
        if self.__nrows != self.__ncols:
            raise TypeError("diagonal matrix must be square")
        if self.__nrows < len(entries):
            raise ValueError('number of diagonal matrix entries (%s) exceeds the matrix size (%s)' % (len(entries), self.__nrows))
        A = self.zero_matrix().__copy__()
        for i in range(len(entries)):
            A[i, i] = entries[i]
        return A

    def is_dense(self):
        """
        Return whether matrices in ``self`` are dense.

        EXAMPLES::

            sage: Mat(RDF,2,3).is_sparse()
            False
            sage: Mat(RR,123456,22,sparse=True).is_sparse()
            True
        """
        return not self.__is_sparse

    def is_sparse(self):
        """
        Return whether matrices in ``self`` are sparse.

        EXAMPLES::

            sage: Mat(GF(2011), 10000).is_sparse()                                      # needs sage.rings.finite_rings
            False
            sage: Mat(GF(2011), 10000, sparse=True).is_sparse()                         # needs sage.rings.finite_rings
            True
        """
        return self.__is_sparse

    def is_finite(self):
        """
        Return whether this matrix space is finite.

        EXAMPLES::

            sage: MatrixSpace(GF(101), 10000).is_finite()
            True
            sage: MatrixSpace(QQ, 2).is_finite()
            False
        """
        return self.base_ring().is_finite()

    def gen(self, n):
        """
        Return the `n`-th generator of this matrix space.

        This does not compute all basis matrices, so it is reasonably
        intelligent.

        EXAMPLES::

            sage: M = Mat(GF(7), 10000, 5); M.ngens()
            50000
            sage: a = M.10
            sage: a[:4]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [1 0 0 0 0]
            [0 0 0 0 0]
        """
        if hasattr(self, '__basis'):
            return self.__basis[n]
        r = n // self.__ncols
        c = n - (r * self.__ncols)
        z = self.zero_matrix().__copy__()
        z[r, c] = 1
        return z

    @cached_method
    def zero_matrix(self):
        """
        Return the zero matrix in ``self``.

        ``self`` must be a space of square matrices. The returned matrix is
        immutable. Please use ``copy`` if you want a modified copy.

        EXAMPLES::

            sage: z = MatrixSpace(GF(7), 2, 4).zero_matrix(); z
            [0 0 0 0]
            [0 0 0 0]
            sage: z.is_mutable()
            False

        TESTS::

            sage: MM = MatrixSpace(RDF,1,1,sparse=False); mat = MM.zero_matrix()
            sage: copy(mat)
            [0.0]
            sage: MM = MatrixSpace(RDF,0,0,sparse=False); mat = MM.zero_matrix()
            sage: copy(mat)
            []
            sage: mat.is_mutable()
            False
            sage: MM.zero().is_mutable()
            False

        Check that :issue:`38221` is fixed::

            sage: # needs sage.groups
            sage: G = CyclicPermutationGroup(7)
            sage: R = GF(2)
            sage: A = G.algebra(R)
            sage: S = MatrixSpace(A, 3, 3)
            sage: S.zero_matrix()
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        res = self.element_class(self, None, False, False)
        res.set_immutable()
        return res

    zero = zero_matrix

    def ngens(self):
        """
        Return the number of generators of this matrix space.

        This is the number of entries in the matrices in this space.

        EXAMPLES::

            sage: M = Mat(GF(7), 100, 200); M.ngens()
            20000
        """
        return self.dimension()

    def matrix(self, x=None, **kwds):
        r"""
        Create a matrix in ``self``.

        INPUT:

        - ``x`` -- data to construct a new matrix from. See :func:`matrix`

        - ``coerce`` -- boolean (default: ``True``); if ``False``, assume
          without checking that the values in ``x`` lie in the base ring

        OUTPUT: a matrix in ``self``

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2)
            sage: M.matrix([[1,0],[0,-1]])
            [ 1  0]
            [ 0 -1]
            sage: M.matrix([1,0,0,-1])
            [ 1  0]
            [ 0 -1]
            sage: M.matrix([1,2,3,4])
            [1 2]
            [3 4]

        Note that the last "flip" cannot be performed if ``x`` is a
        matrix, no matter what is ``rows`` (it used to be possible but
        was fixed by :issue:`10793`)::

            sage: projection = matrix(ZZ,[[1,0,0],[0,1,0]])
            sage: projection
            [1 0 0]
            [0 1 0]
            sage: projection.parent()
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
            sage: M = MatrixSpace(ZZ, 3 , 2)
            sage: M
            Full MatrixSpace of 3 by 2 dense matrices over Integer Ring
            sage: M(projection)
            Traceback (most recent call last):
            ...
            ValueError: inconsistent number of rows: should be 3 but got 2

        If you really want to make from a matrix another matrix of different
        dimensions, use either transpose method or explicit conversion to a
        list::

            sage: M(projection.list())
            [1 0]
            [0 0]
            [1 0]

        TESTS:

        The following corner cases were problematic while working on
        :issue:`10628`::

            sage: MS = MatrixSpace(ZZ,2,1)
            sage: MS([[1],[2]])
            [1]
            [2]

            sage: # needs sage.rings.real_mpfr
            sage: MS = MatrixSpace(CC, 2, 1)
            sage: x = polygen(ZZ, 'x')
            sage: F = NumberField(x^2 + 1, name='x')                                    # needs sage.rings.number_field
            sage: MS([F(1), F(0)])                                                      # needs sage.rings.number_field
            [ 1.00000000000000]
            [0.000000000000000]

        :issue:`10628` allowed to provide the data as lists of matrices, but
        :issue:`13012` prohibited it::

            sage: MS = MatrixSpace(ZZ, 4,2)
            sage: MS0 = MatrixSpace(ZZ, 2)
            sage: MS.matrix([MS0([1,2,3,4]), MS0([5,6,7,8])])
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'> to an integer

        A mixed list of matrices and vectors is prohibited as well::

            sage: MS.matrix( [MS0([1,2,3,4])] + list(MS0([5,6,7,8])) )
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'> to an integer

        Check that :issue:`13302` is fixed::

            sage: MatrixSpace(Qp(3), 1,1)([Qp(3).zero()])                               # needs sage.rings.padics
            [0]
            sage: MatrixSpace(Qp(3), 1,1)([Qp(3)(4/3)])                                 # needs sage.rings.padics
            [3^-1 + 1 + O(3^19)]

        One-rowed matrices over combinatorial free modules used to break
        the constructor (:issue:`17124`). Check that this is fixed::

            sage: # needs sage.combinat
            sage: Sym = SymmetricFunctions(ZZ)
            sage: h = Sym.h()
            sage: MatrixSpace(h, 1,1)([h[1]])
            [h[1]]
            sage: MatrixSpace(h, 2,1)([h[1], h[2]])
            [h[1]]
            [h[2]]

        Converting sparse to dense matrices used to be too slow
        (:issue:`20470`). Check that this is fixed::

            sage: m = identity_matrix(GF(2), 2000, sparse=True)
            sage: MS = MatrixSpace(GF(2), 2000, sparse=False)
            sage: md = MS(m)
            sage: md.parent() is MS
            True
        """
        return self(x, **kwds)

    def matrix_space(self, nrows=None, ncols=None, sparse=False):
        """
        Return the matrix space with given number of rows, columns and
        sparsity over the same base ring as self, and defaults the same as
        ``self``.

        EXAMPLES::

            sage: M = Mat(GF(7), 100, 200)
            sage: M.matrix_space(5000)
            Full MatrixSpace of 5000 by 200 dense matrices over Finite Field of size 7
            sage: M.matrix_space(ncols=5000)
            Full MatrixSpace of 100 by 5000 dense matrices over Finite Field of size 7
            sage: M.matrix_space(sparse=True)
            Full MatrixSpace of 100 by 200 sparse matrices over Finite Field of size 7
        """
        if nrows is None:
            nrows = self.__nrows
        if ncols is None:
            ncols = self.__ncols
        base = self._base
        return MatrixSpace(base, nrows, ncols, sparse=sparse)

    def ncols(self):
        """
        Return the number of columns of matrices in this space.

        EXAMPLES::

            sage: M = Mat(ZZ['x'], 200000, 500000, sparse=True)
            sage: M.ncols()
            500000
        """
        return self.__ncols

    def nrows(self):
        """
        Return the number of rows of matrices in this space.

        EXAMPLES::

            sage: M = Mat(ZZ, 200000, 500000)
            sage: M.nrows()
            200000
        """
        return self.__nrows

    def row_space(self):
        """
        Return the module spanned by all rows of matrices in this matrix
        space. This is a free module of rank the number of rows. It will be
        sparse or dense as this matrix space is sparse or dense.

        EXAMPLES::

            sage: M = Mat(ZZ,20,5,sparse=False); M.row_space()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        try:
            return self.__row_space
        except AttributeError:
            self.__row_space = sage.modules.free_module.FreeModule(self.base_ring(),
                                                self.ncols(), sparse=self.is_sparse())
            return self.__row_space

    def column_space(self):
        """
        Return the module spanned by all columns of matrices in this matrix
        space. This is a free module of rank the number of columns. It will
        be sparse or dense as this matrix space is sparse or dense.

        EXAMPLES::

            sage: M = Mat(GF(9,'a'), 20, 5, sparse=True); M.column_space()              # needs sage.rings.finite_rings
            Sparse vector space of dimension 20 over Finite Field in a of size 3^2
        """
        try:
            return self.__column_space
        except AttributeError:
            self.__column_space = sage.modules.free_module.FreeModule(self.base_ring(), self.nrows(),
                                                                   sparse=self.is_sparse())
            return self.__column_space

    def random_element(self, density=None, *args, **kwds):
        """
        Return a random element from this matrix space.

        INPUT:

        - ``density`` -- ``float`` or ``None`` (default: ``None``);  rough
          measure of the proportion of nonzero entries in the random matrix;
          if set to ``None``, all entries of the matrix are randomized,
          allowing for any element of the underlying ring, but if set to
          a ``float``, a proportion of entries is selected and randomized to
          nonzero elements of the ring

        - ``*args, **kwds`` -- remaining parameters, which may be passed to
          the random_element function of the base ring. ("may be", since this
          function calls the ``randomize`` function on the zero matrix, which
          need not call the ``random_element`` function of the base ring at
          all in general.)

        OUTPUT: Matrix

        .. NOTE::

            This method will randomize a proportion of roughly ``density`` entries
            in a newly allocated zero matrix.

            By default, if the user sets the value of ``density`` explicitly, this
            method will enforce that these entries are set to nonzero values.
            However, if the test for equality with zero in the base ring is too
            expensive, the user can override this behaviour by passing the
            argument ``nonzero=False`` to this method.

            Otherwise, if the user does not set the value of ``density``, the
            default value is taken to be 1, and the option ``nonzero=False`` is
            passed to the ``randomize`` method.

        EXAMPLES::

            sage: M = Mat(ZZ, 2, 5).random_element()
            sage: TestSuite(M).run()

            sage: M = Mat(QQ, 2, 5).random_element(density=0.5)
            sage: TestSuite(M).run()

            sage: M = Mat(QQ, 3, sparse=True).random_element()
            sage: TestSuite(M).run()                                                    # needs sage.libs.pari

            sage: M = Mat(GF(9,'a'), 3, sparse=True).random_element()                   # needs sage.rings.finite_rings
            sage: TestSuite(M).run()                                                    # needs sage.rings.finite_rings
        """
        Z = self.zero_matrix().__copy__()
        if density is None:
            Z.randomize(density=float(1), nonzero=kwds.pop('nonzero', False),
                *args, **kwds)
        else:
            Z.randomize(density=density, nonzero=kwds.pop('nonzero', True),
                *args, **kwds)
        return Z

    def _an_element_(self):
        """
        Create a typical element of this matrix space.

        This uses ``some_elements`` of the base ring.

        EXAMPLES::

            sage: MatrixSpace(QQ, 3, 3).an_element()  # indirect doctest
            [ 1/2 -1/2    2]
            [  -2    0    1]
            [  -1   42  2/3]

        TESTS::

            sage: MatrixSpace(ZZ, 0, 0).an_element()
            []

        Check that this works for large matrices and that it returns a
        matrix which is not too trivial::

            sage: M = MatrixSpace(GF(2), 100, 100).an_element()
            sage: M.rank() >= 2
            True

        Check that this works for sparse matrices::

            sage: M = MatrixSpace(ZZ, 1000, 1000, sparse=True).an_element()
            sage: 96 <= M.density() * 10^6 <= 99
            True
        """
        from .args import MatrixArgs
        dim = self.dimension()
        if dim > 100 and self.is_sparse():
            # Sparse case: add 100 elements
            D = {}
            nr = self.nrows()
            nc = self.ncols()
            from random import randrange
            n = 0
            while True:
                for el in self.base().some_elements():
                    if n == 100:
                        ma = MatrixArgs(D, space=self)
                        del D
                        return ma.matrix()
                    D[randrange(nr), randrange(nc)] = el
                    n += 1
                assert D
        else:
            # Dense case
            # Keep appending to L until we have enough elements
            L = []
            while True:
                for el in self.base().some_elements():
                    if len(L) == dim:
                        ma = MatrixArgs(L, space=self)
                        del L  # for efficiency: this may avoid a copy of L
                        return ma.matrix()
                    L.append(el)
                assert L

    def some_elements(self):
        r"""
        Return some elements of this matrix space.

        See :class:`TestSuite` for a typical use case.

        OUTPUT: an iterator

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, 2)
            sage: tuple(M.some_elements())
            (
            [ 0  1]  [1 0]  [0 1]  [0 0]  [0 0]
            [-1  2], [0 0], [0 0], [1 0], [0 1]
            )
            sage: M = MatrixSpace(QQ, 2, 3)
            sage: tuple(M.some_elements())
            (
            [ 1/2 -1/2    2]  [1 0 0]  [0 1 0]  [0 0 1]  [0 0 0]  [0 0 0]  [0 0 0]
            [  -2    0    1], [0 0 0], [0 0 0], [0 0 0], [1 0 0], [0 1 0], [0 0 1]
            )
            sage: M = MatrixSpace(SR, 2, 2)                                             # needs sage.symbolic
            sage: tuple(M.some_elements())                                              # needs sage.symbolic
            (
            [some_variable some_variable]  [1 0]  [0 1]  [0 0]  [0 0]
            [some_variable some_variable], [0 0], [0 0], [1 0], [0 1]
            )
        """
        yield self.an_element()
        yield from self.gens()

    def _magma_init_(self, magma):
        r"""
        EXAMPLES: We first coerce a square matrix.

        ::

            sage: magma(MatrixSpace(QQ, 3))                     # optional - magma
            Full Matrix Algebra of degree 3 over Rational Field

        ::

            sage: magma(MatrixSpace(Integers(8), 2, 3))         # optional - magma
            Full RMatrixSpace of 2 by 3 matrices over IntegerRing(8)
        """
        K = magma(self.base_ring())
        if self.__nrows == self.__ncols:
            return 'MatrixAlgebra(%s,%s)' % (K.name(), self.__nrows)
        return 'RMatrixSpace(%s,%s,%s)' % (K.name(), self.__nrows, self.__ncols)

    def _polymake_init_(self):
        r"""
        Return the polymake representation of the matrix space.

        EXAMPLES::

            sage: polymake(MatrixSpace(QQ, 3))                  # optional - jupymake
            Matrix<Rational>
            sage: polymake(MatrixSpace(QuadraticField(5), 3))   # optional - jupymake, needs sage.rings.number_field
            Matrix<QuadraticExtension>
        """
        from sage.interfaces.polymake import polymake
        K = polymake(self.base_ring())
        return '"Matrix<{}>"'.format(K)

    def _random_nonzero_element(self, *args, **kwds):
        """
        Return a random nonzero matrix.

        This function repeatedly calls ``random_element`` until a nonzero
        matrix is obtained.

        INPUT:

        - ``*args``, ``**kwds`` -- parameters that can be forwarded to the
          ``random_element`` method

        OUTPUT: random nonzero matrix

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 4)
            sage: A = M._random_nonzero_element()
            sage: A.is_zero()
            False

            sage: M = MatrixSpace(ZZ, 0)
            sage: A = M._random_nonzero_element()
            Traceback (most recent call last):
            ...
            ValueError: Full MatrixSpace of 0 by 0 dense matrices over Integer Ring only has zero elements
        """
        if 0 in self.dims():
            raise ValueError("{} only has zero elements".format(self))
        rand_matrix = self.random_element(*args, **kwds)
        while rand_matrix.is_zero():
            rand_matrix = self.random_element(*args, **kwds)
        return rand_matrix

    def from_vector(self, vector, order=None, coerce=True):
        r"""
        Build an element of ``self`` from a vector.

        EXAMPLES::

            sage: A = matrix([[1,2,3], [4,5,6]])
            sage: v = vector(A); v
            (1, 2, 3, 4, 5, 6)
            sage: MS = A.parent()
            sage: MS.from_vector(v)
            [1 2 3]
            [4 5 6]
            sage: order = [(1,2), (1,0), (0,1), (0,2), (0,0), (1,1)]
            sage: MS.from_vector(v, order=order)
            [5 3 4]
            [2 6 1]
        """
        if order is None:
            if self.is_dense():
                return self.element_class(self, vector, coerce=coerce)
            else:
                nc = self.ncols()
                d = {(k // nc, k % nc): c for k, c in vector.dict().items()}
                return self.element_class(self, d, coerce=coerce)
        return super().from_vector(vector, order=order, coerce=coerce)

    def _from_dict(self, d, coerce=True, remove_zeros=True):
        r"""
        Construct an element of ``self`` from the dictionary ``d``.

        INPUT:

        - ``coerce`` -- boolean; coerce the coefficients to the base ring
        - ``remove_zeros`` -- ignored; for compatibility

        EXAMPLES::

            sage: MS = MatrixSpace(QQ['x'], 10, 5)
            sage: A = MS.random_element()
            sage: MS._from_dict(A.monomial_coefficients()) == A
            True
        """
        return self.element_class(self, d, coerce=coerce)


def dict_to_list(entries, nrows, ncols):
    r"""
    Given a dictionary of coordinate tuples, return the list given by
    reading off the nrows\*ncols matrix in row order.

    EXAMPLES::

        sage: from sage.matrix.matrix_space import dict_to_list
        sage: d = {}
        sage: d[(0,0)] = 1
        sage: d[(1,1)] = 2
        sage: dict_to_list(d, 2, 2)
        [1, 0, 0, 2]
        sage: dict_to_list(d, 2, 3)
        [1, 0, 0, 0, 2, 0]
    """
    v = [0] * (nrows * ncols)
    for ij, y in entries.items():
        i, j = ij
        v[i * ncols + j] = y
    return v


def _test_trivial_matrices_inverse(ring, sparse=True, implementation=None, checkrank=True):
    """
    Test inversion, determinant and is_invertible for trivial matrices.

    This function is a helper to check that the inversion of trivial matrices
    (of size 0x0, nx0, 0xn or 1x1) is handled consistently by the various
    implementation of matrices. The coherency is checked through a bunch of
    assertions. If an inconsistency is found, an :exc:`AssertionError` is
    raised which should make clear what is the problem.

    INPUT:

    - ``ring`` -- a ring
    - ``sparse`` -- boolean
    - ``checkrank`` -- boolean

    OUTPUT: nothing if everything is correct, otherwise raise an AssertionError

    The methods determinant, is_invertible, rank and inverse are checked for
     - the 0x0 empty identity matrix
     - the 0x3 and 3x0 matrices
     - the 1x1 null matrix [0]
     - the 1x1 identity matrix [1]

    If ``checkrank`` is ``False`` then the rank is not checked. This is used
    the check matrix over ring where echelon form is not implemented.

    .. TODO::

        This must be adapted to category check framework when ready
        (see :issue:`5274`).

    TESTS::

        sage: from sage.matrix.matrix_space import _test_trivial_matrices_inverse as tinv
        sage: tinv(ZZ, sparse=True)
        sage: tinv(ZZ, sparse=False, implementation='flint')                            # needs sage.libs.linbox
        sage: tinv(ZZ, sparse=False, implementation='generic')
        sage: tinv(QQ, sparse=True)
        sage: tinv(QQ, sparse=False, implementation='flint')                            # needs sage.libs.linbox
        sage: tinv(QQ, sparse=False, implementation='generic')
        sage: tinv(GF(11), sparse=True)
        sage: tinv(GF(11), sparse=False)
        sage: tinv(GF(2), sparse=True)
        sage: tinv(GF(2), sparse=False)
        sage: tinv(SR, sparse=True)                                                     # needs sage.symbolic
        sage: tinv(SR, sparse=False)                                                    # needs sage.symbolic

        sage: # needs scipy
        sage: tinv(RDF, sparse=True)
        sage: tinv(RDF, sparse=False)
        sage: tinv(CDF, sparse=True)                                                    # needs sage.rings.complex_double
        sage: tinv(CDF, sparse=False)                                                   # needs sage.rings.complex_double
        sage: tinv(CyclotomicField(7), sparse=True)                                     # needs sage.rings.number_field
        sage: tinv(CyclotomicField(7), sparse=False)                                    # needs sage.rings.number_field
        sage: tinv(QQ['x,y'], sparse=True)
        sage: tinv(QQ['x,y'], sparse=False)
    """
    # Check that the empty 0x0 matrix is it's own inverse with det=1.
    ms00 = MatrixSpace(ring, 0, 0, sparse=sparse)
    m00 = ms00(0)
    assert m00.determinant() == ring(1)
    assert m00.is_invertible()
    assert m00.inverse() == m00
    if checkrank:
        assert m00.rank() == 0

    # Check that the empty 0x3 and 3x0 matrices are not invertible and that
    # computing the determinant raise the proper exception.
    for ms0 in [MatrixSpace(ring, 0, 3, sparse=sparse),
                MatrixSpace(ring, 3, 0, sparse=sparse)]:
        mn0 = ms0(0)
        assert not mn0.is_invertible()
        try:
            d = mn0.determinant()
            print(d)
            res = False
        except ValueError:
            res = True
        assert res
        try:
            mn0.inverse()
            res = False
        except ArithmeticError:
            res = True
        assert res
        if checkrank:
            assert mn0.rank() == 0

    # Check that the null 1x1 matrix is not invertible and that det=0
    ms1 = MatrixSpace(ring, 1, 1, sparse=sparse)
    m0 = ms1(0)
    assert not m0.is_invertible()
    assert m0.determinant() == ring(0)
    try:
        m0.inverse()
        res = False
    except (ZeroDivisionError, RuntimeError):
        # FIXME: Make pynac throw a ZeroDivisionError on division by
        # zero instead of a runtime Error
        res = True
    assert res
    if checkrank:
        assert m0.rank() == 0

    # Check that the identity 1x1 matrix is its own inverse with det=1
    m1 = ms1(1)
    assert m1.is_invertible()
    assert m1.determinant() == ring(1)
    inv = m1.inverse()
    assert inv == m1
    if checkrank:
        assert m1.rank() == 1


test_trivial_matrices_inverse = deprecated_function_alias(33612, _test_trivial_matrices_inverse)


# Fix unpickling Matrix_modn_dense and Matrix_integer_2x2
lazy_import('sage.matrix.matrix_modn_dense_double', 'Matrix_modn_dense_double')
lazy_import('sage.matrix.matrix_integer_dense', 'Matrix_integer_dense')


def _MatrixSpace_ZZ_2x2():
    from sage.rings.integer_ring import ZZ
    return MatrixSpace(ZZ, 2)


register_unpickle_override('sage.matrix.matrix_modn_dense',
    'Matrix_modn_dense', Matrix_modn_dense_double)
register_unpickle_override('sage.matrix.matrix_integer_2x2',
    'Matrix_integer_2x2', Matrix_integer_dense)
register_unpickle_override('sage.matrix.matrix_integer_2x2',
    'MatrixSpace_ZZ_2x2_class', MatrixSpace)
register_unpickle_override('sage.matrix.matrix_integer_2x2',
    'MatrixSpace_ZZ_2x2', _MatrixSpace_ZZ_2x2)
lazy_import('sage.matrix.matrix_gf2e_dense', 'unpickle_matrix_gf2e_dense_v0')
register_unpickle_override('sage.matrix.matrix_mod2e_dense',
    'unpickle_matrix_mod2e_dense_v0', unpickle_matrix_gf2e_dense_v0)
