# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.modules
r"""
Cubic Hecke matrix representations

This module contains the class :class:`CubicHeckeMatrixRep` which is used to
treat the matrix representations of the elements of the cubic Hecke algebra
(:class:`~sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeAlgebra`)
together with its parent class :class:`CubicHeckeMatrixSpace`. Furthermore,
it contains enums for their types (:class:`RepresentationType`) and names
(:class:`AbsIrreducibeRep`).


AUTHORS:

- Sebastian Oehms May 2020: initial version
"""
# ###########################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ###########################################################################
from enum import Enum

from sage.misc.cachefunc import cached_method
from sage.misc.verbose import verbose
from sage.rings.integer import Integer
from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.matrix.special import block_diagonal_matrix
from sage.databases.cubic_hecke_db import CubicHeckeDataSection as sc


# -------------------------------------------
# Enum for the generators sign (of exponent)
# -------------------------------------------
class GenSign(Enum):
    r"""
    Enum class to select the braid generators sign.

    EXAMPLES::

        sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
        sage: chmr.GenSign.pos
        <GenSign.pos: 1>
        sage: chmr.GenSign.neg
        <GenSign.neg: -1>
    """
    pos = 1
    neg = -1


# -------------------------------------------
# Enum for typ of matrix representation
# -------------------------------------------
class RepresentationType(Enum):
    r"""
    Enum class to select a representation type for the cubic Hecke algebra.

    - ``RegularLeft`` -- left regular representations
    - ``RegularRight`` -- right regular representations
    - ``SplitIrredMarin`` -- split irreducible representations obtained from
      Ivan Marin's data
    - ``SplitIrredChevie`` -- the split irreducible representations obtained
      from CHEVIE via the ``GAP3`` interface

    EXAMPLES::

        sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
        sage: chmr.RepresentationType.RegularLeft.is_regular()
        True
    """
    def is_split(self):
        r"""
        Return ``True`` if this representation type is absolutely split,
        ``False`` else-wise.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chevie = chmr.RepresentationType.SplitIrredChevie
            sage: chevie.is_split()
            True
        """
        return self.value['split']

    def is_regular(self):
        r"""
        Return ``True`` if this representation type is regular, ``False``
        else-wise.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: reg_left = chmr.RepresentationType.RegularLeft
            sage: reg_left.is_regular()
            True
        """
        return self.value['regular']

    def data_section(self):
        r"""
        Return the name of the data file. For more information see
        :class:`~sage.databases.cubic_hecke_db.CubicHeckeDataBase`.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: reg_left = chmr.RepresentationType.RegularLeft
            sage: reg_left.data_section()
            <CubicHeckeDataSection.regular_left: 'regular_left'>
        """
        return self.value['data']

    def number_of_representations(self, nstrands):
        r"""
        Return the number of representations existing to that type.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.RepresentationType.SplitIrredChevie.number_of_representations(4)
            24
            sage: chmr.RepresentationType.SplitIrredMarin.number_of_representations(4)
            24
        """
        if self.value['data'] is None:
            if nstrands < 1 or nstrands > 5:
                raise ValueError("nstrands must be between 1 and 5")
        elif nstrands < 1 or nstrands > 4:
            raise ValueError("nstrands must be between 1 and 4")
        return self.value['num_rep'][nstrands - 1]

    RegularLeft = {'split': False, 'regular': True, 'data': sc.regular_left, 'num_rep': [1, 1, 1, 1]}
    RegularRight = {'split': False, 'regular': True, 'data': sc.regular_right, 'num_rep': [1, 1, 1, 1]}
    SplitIrredMarin = {'split': True, 'regular': False, 'data': sc.split_irred, 'num_rep': [1, 3, 7, 24]}
    SplitIrredChevie = {'split': True, 'regular': False, 'data': None, 'num_rep': [1, 3, 7, 24, 30]}


# ---------------------------------------------
# Enum for absolute irreducible representations
# ---------------------------------------------
class AbsIrreducibeRep(Enum):
    r"""
    Enum class to select an absolutely irreducible representation for the cubic
    Hecke algebra (``CHAn``) on `n`-strands.

    The names are build as follows: Take the determinant of one of the
    generators of the ``CHAn``. This is a monomial in the generic extension
    ring (``GER``) of ``CHA``, say ``a^ib^jc^k`` where ``a, b`` and ``c`` are
    the generators of ``GER``. This does not depend on the choice of the
    generator of ``CHA``, since these are conjugated to each other. This
    monomial might be looked as the weight of the representation. Therefore we
    use it as a name:

        ``Wn_ijk``

    The only ambiguity among the available irreducible representations occurs for the two nine-dimensional modules, which
    are conjugated to each other and distinguished by these names:

        ``W4_333`` and ``W4_333bar``

    Examples of names:

    - ``W2_100`` -- one dimensional representation of the cubic Hecke algebra on 2 strands corresponding to the first root
      of the cubic equation
    - ``W3_111`` -- three dimensional irreducible representation of the cubic Hecke algebra on 3 strands
    - ``W4_242`` -- eight dimensional irreducible representation of the cubic Hecke algebra on 4 strands having the second
      root of the cubic equation as weight of dimension 4

    Alternative names are taken from [MW2012]_ and can be shown by
    :meth:`alternative_name`.

    EXAMPLES::

        sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
        sage: [irr.name for irr in chmr.AbsIrreducibeRep]
        ['W2_100', 'W2_001', 'W2_010', 'W3_100', 'W3_001', 'W3_010', 'W3_011', 'W3_110',
         'W3_101', 'W3_111', 'W4_100', 'W4_001', 'W4_010', 'W4_011', 'W4_110', 'W4_101',
         'W4_111', 'W4_120', 'W4_201', 'W4_012', 'W4_102', 'W4_210', 'W4_021', 'W4_213',
         'W4_132', 'W4_321', 'W4_231', 'W4_123', 'W4_312', 'W4_422', 'W4_224', 'W4_242',
         'W4_333', 'W4_333bar', 'W5_100', 'W5_001', 'W5_010', 'W5_013', 'W5_130', 'W5_301',
         'W5_031', 'W5_103', 'W5_310', 'W5_203', 'W5_032', 'W5_320', 'W5_230', 'W5_023',
         'W5_302', 'W5_033', 'W5_330', 'W5_303', 'W5_163', 'W5_631', 'W5_316', 'W5_136',
         'W5_613', 'W5_361', 'W5_366', 'W5_663', 'W5_636', 'W5_933', 'W5_339', 'W5_393']

    REFERENCES:

    - [MW2012]_
    """
    def alternative_name(self):
        r"""
        Return the name of the split irreducible representation for cubic Hecke
        algebras for up to four strands as given in [MW2012]_.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.AbsIrreducibeRep.W3_011.alternative_name()
            'Tbc'
        """
        return self.value['alt_name']

    def dimension(self):
        r"""
        Return the dimension of the representation.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.AbsIrreducibeRep.W3_111.dimension()
            3
        """
        return self.value['dim']

    def number_gens(self):
        r"""
        Return the number of generators of the underlying cubic Hecke algebra.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.AbsIrreducibeRep.W3_001.number_gens()
            2
            sage: chmr.AbsIrreducibeRep.W4_001.number_gens()
            3
        """
        return self.value['ngens']

    def length_orbit(self):
        r"""
        Return the length of the orbit of this representation under the action
        of the Galois group of the cubic equation.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.AbsIrreducibeRep.W3_001.length_orbit()
            3
            sage: chmr.AbsIrreducibeRep.W3_111.length_orbit()
            1
        """
        return self.value['len_orbit']

    def gap_index(self):
        r"""
        Return the array index of this representation for the access
        to the ``GAP3`` package ``CHEVIE``.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.AbsIrreducibeRep.W3_111.gap_index()
            6
        """
        return self.value['gap_ind']

    def internal_index(self):
        r"""
        Return the array index of this representation for the internal access.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: chmr.AbsIrreducibeRep.W3_111.internal_index()
            6
        """
        return self.value['intern_ind']

    # -------------------------------------------------------------------------------------------------
    # absolutely irreducible representations corresponding to braids on 2 strands
    # -------------------------------------------------------------------------------------------------
    W2_100 = {'alt_name': 'Sa', 'dim': 1, 'ngens': 1, 'len_orbit': 3, 'gap_ind': 0, 'intern_ind': 0}
    W2_001 = {'alt_name': 'Sc', 'dim': 1, 'ngens': 1, 'len_orbit': 3, 'gap_ind': 1, 'intern_ind': 1}
    W2_010 = {'alt_name': 'Sb', 'dim': 1, 'ngens': 1, 'len_orbit': 3, 'gap_ind': 2, 'intern_ind': 2}

    # -------------------------------------------------------------------------------------------------
    # absolutely irreducible representations corresponding to braids on 3 strands
    # -------------------------------------------------------------------------------------------------
    W3_100 = {'alt_name': 'Sa', 'dim': 1, 'ngens': 2, 'len_orbit': 3, 'gap_ind': 0, 'intern_ind': 0}
    W3_001 = {'alt_name': 'Sc', 'dim': 1, 'ngens': 2, 'len_orbit': 3, 'gap_ind': 1, 'intern_ind': 1}
    W3_010 = {'alt_name': 'Sb', 'dim': 1, 'ngens': 2, 'len_orbit': 3, 'gap_ind': 2, 'intern_ind': 2}

    W3_011 = {'alt_name': 'Tbc', 'dim': 2, 'ngens': 2, 'len_orbit': 3, 'gap_ind': 3, 'intern_ind': 3}
    W3_110 = {'alt_name': 'Tab', 'dim': 2, 'ngens': 2, 'len_orbit': 3, 'gap_ind': 4, 'intern_ind': 4}
    W3_101 = {'alt_name': 'Tac', 'dim': 2, 'ngens': 2, 'len_orbit': 3, 'gap_ind': 5, 'intern_ind': 5}

    W3_111 = {'alt_name': 'V', 'dim': 3, 'ngens': 2, 'len_orbit': 1, 'gap_ind': 6, 'intern_ind': 6}

    # -------------------------------------------------------------------------------------------------
    # absolutely irreducible representations corresponding to braids on 4 strands
    # -------------------------------------------------------------------------------------------------
    W4_100 = {'alt_name': 'Sa', 'dim': 1, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 0, 'intern_ind': 0}
    W4_001 = {'alt_name': 'Sc', 'dim': 1, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 1, 'intern_ind': 1}
    W4_010 = {'alt_name': 'Sb', 'dim': 1, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 2, 'intern_ind': 2}

    W4_011 = {'alt_name': 'Tbc', 'dim': 2, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 3, 'intern_ind': 3}
    W4_110 = {'alt_name': 'Tab', 'dim': 2, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 4, 'intern_ind': 4}
    W4_101 = {'alt_name': 'Tac', 'dim': 2, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 5, 'intern_ind': 5}

    W4_111 = {'alt_name': 'V', 'dim': 3, 'ngens': 3, 'len_orbit': 1, 'gap_ind': 6, 'intern_ind': 6}

    W4_120 = {'alt_name': 'Uba', 'dim': 3, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 7, 'intern_ind': 7}
    W4_201 = {'alt_name': 'Uac', 'dim': 3, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 8, 'intern_ind': 8}
    W4_012 = {'alt_name': 'Ucb', 'dim': 3, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 9, 'intern_ind': 9}
    W4_102 = {'alt_name': 'Uca', 'dim': 3, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 10, 'intern_ind': 10}
    W4_210 = {'alt_name': 'Uab', 'dim': 3, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 11, 'intern_ind': 11}
    W4_021 = {'alt_name': 'Ubc', 'dim': 3, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 12, 'intern_ind': 12}

    W4_213 = {'alt_name': 'Vcab', 'dim': 6, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 13, 'intern_ind': 13}
    W4_132 = {'alt_name': 'Vbca', 'dim': 6, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 14, 'intern_ind': 14}
    W4_321 = {'alt_name': 'Vabc', 'dim': 6, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 15, 'intern_ind': 15}
    W4_231 = {'alt_name': 'Vbac', 'dim': 6, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 16, 'intern_ind': 16}
    W4_123 = {'alt_name': 'Vcba', 'dim': 6, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 17, 'intern_ind': 17}
    W4_312 = {'alt_name': 'Vacb', 'dim': 6, 'ngens': 3, 'len_orbit': 6, 'gap_ind': 18, 'intern_ind': 18}

    W4_422 = {'alt_name': 'Wa', 'dim': 8, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 19, 'intern_ind': 19}
    W4_224 = {'alt_name': 'Wc', 'dim': 8, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 20, 'intern_ind': 20}
    W4_242 = {'alt_name': 'Wb', 'dim': 8, 'ngens': 3, 'len_orbit': 3, 'gap_ind': 21, 'intern_ind': 21}

    W4_333 = {'alt_name': 'X', 'dim': 9, 'ngens': 3, 'len_orbit': 2, 'gap_ind': 22, 'intern_ind': 22}
    W4_333bar = {'alt_name': 'Xbar', 'dim': 9, 'ngens': 3, 'len_orbit': 2, 'gap_ind': 23, 'intern_ind': 23}

    # -------------------------------------------------------------------------------------------------
    # absolutely irreducible representations corresponding to braids on 5 strands
    # -------------------------------------------------------------------------------------------------
    W5_100 = {'alt_name': None, 'dim': 1, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 0, 'intern_ind': 0}
    W5_001 = {'alt_name': None, 'dim': 1, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 1, 'intern_ind': 1}
    W5_010 = {'alt_name': None, 'dim': 1, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 2, 'intern_ind': 2}

    W5_013 = {'alt_name': None, 'dim': 4, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 3, 'intern_ind': 3}
    W5_130 = {'alt_name': None, 'dim': 4, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 4, 'intern_ind': 4}
    W5_301 = {'alt_name': None, 'dim': 4, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 5, 'intern_ind': 5}
    W5_031 = {'alt_name': None, 'dim': 4, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 6, 'intern_ind': 6}
    W5_103 = {'alt_name': None, 'dim': 4, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 7, 'intern_ind': 7}
    W5_310 = {'alt_name': None, 'dim': 4, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 8, 'intern_ind': 8}

    W5_203 = {'alt_name': None, 'dim': 5, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 9, 'intern_ind': 9}
    W5_032 = {'alt_name': None, 'dim': 5, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 10, 'intern_ind': 10}
    W5_320 = {'alt_name': None, 'dim': 5, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 11, 'intern_ind': 11}
    W5_230 = {'alt_name': None, 'dim': 5, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 12, 'intern_ind': 12}
    W5_023 = {'alt_name': None, 'dim': 5, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 13, 'intern_ind': 13}
    W5_302 = {'alt_name': None, 'dim': 5, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 14, 'intern_ind': 14}

    W5_033 = {'alt_name': None, 'dim': 6, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 15, 'intern_ind': 15}
    W5_330 = {'alt_name': None, 'dim': 6, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 16, 'intern_ind': 16}
    W5_303 = {'alt_name': None, 'dim': 6, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 17, 'intern_ind': 17}

    W5_163 = {'alt_name': None, 'dim': 10, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 18, 'intern_ind': 18}
    W5_631 = {'alt_name': None, 'dim': 10, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 19, 'intern_ind': 19}
    W5_316 = {'alt_name': None, 'dim': 10, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 20, 'intern_ind': 20}
    W5_136 = {'alt_name': None, 'dim': 10, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 21, 'intern_ind': 21}
    W5_613 = {'alt_name': None, 'dim': 10, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 22, 'intern_ind': 22}
    W5_361 = {'alt_name': None, 'dim': 10, 'ngens': 4, 'len_orbit': 6, 'gap_ind': 23, 'intern_ind': 23}

    W5_366 = {'alt_name': None, 'dim': 15, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 24, 'intern_ind': 24}
    W5_663 = {'alt_name': None, 'dim': 15, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 26, 'intern_ind': 25}
    W5_636 = {'alt_name': None, 'dim': 15, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 27, 'intern_ind': 26}

    W5_933 = {'alt_name': None, 'dim': 15, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 25, 'intern_ind': 27}
    W5_339 = {'alt_name': None, 'dim': 15, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 28, 'intern_ind': 28}
    W5_393 = {'alt_name': None, 'dim': 15, 'ngens': 4, 'len_orbit': 3, 'gap_ind': 29, 'intern_ind': 29}


# ------------------------------------------------------------------------------------------------------------------
# Definition of CubicHeckeMatrixRep
# --------------------------------------------------------------------------------------------------------
class CubicHeckeMatrixRep(Matrix_generic_dense):
    r"""
    Class to supervise the diagonal block matrix structure arising from
    cubic Hecke algebra-representations.

    EXAMPLES::

        sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
        sage: CHA2.<c1> = algebras.CubicHecke(2)
        sage: MS = chmr.CubicHeckeMatrixSpace(CHA2)
        sage: m1 = MS(c1); m1
        [         a          0          0]
        [         0          b          0]
        [         0          0 -b - a + u]
        sage: type(m1)
        <class 'sage.algebras.hecke_algebras.cubic_hecke_matrix_rep.CubicHeckeMatrixSpace_with_category.element_class'>
        sage: m1.block_diagonal_list()
        [[a], [b], [-b - a + u]]

        sage: MSo = chmr.CubicHeckeMatrixSpace(CHA2, original=True)
        sage: MSo(c1)
        [a 0 0]
        [0 b 0]
        [0 0 c]

        sage: reg_left = chmr.RepresentationType.RegularLeft
        sage: MSreg = chmr.CubicHeckeMatrixSpace(CHA2, representation_type=reg_left)
        sage: MSreg(c1)
        [ 0 -v  1]
        [ 1  u  0]
        [ 0  w  0]
        sage: len(_.block_diagonal_list())
        1

    TESTS:

    The minpoly does not work over more generic rings::

        sage: TestSuite(m1).run(skip='_test_minpoly')
    """

    @cached_method
    def _get_block(self, ind):
        r"""
        Return the ``ind``-th sub-matrix block of ``self`` considered
        as block diagonal matrix.

        INPUT:

        - ``ind`` -- integer specifying the list index according to
          :meth:`internal_index` respectively :meth:`gap_index`

        OUTPUT:

        An instance of :class:`Matrix_generic_dense` representing
        the specified block of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: c1.matrix()._get_block(0)      # indirect doctest
            [a]
        """
        representation_type = self.parent()._representation_type
        if not representation_type.is_split():
            return matrix(self)
        n = self.parent()._cubic_hecke_algebra.ngens()
        s = sum(irr_rep.dimension() for irr_rep in AbsIrreducibeRep if irr_rep.number_gens() == n and irr_rep.internal_index() < ind)
        for irr_rep in AbsIrreducibeRep:
            if irr_rep.number_gens() == n and irr_rep.internal_index() == ind:
                d = irr_rep.dimension()
                return matrix(self.submatrix(s, s, d, d))
        raise ValueError('no irreducible representation for this index')

    @cached_method
    def _irr_to_ind(self, irr):
        r"""
        Return the index if the given split irreducible representation
        of ``self``.

        INPUT:

        - ``irr`` -- an instance of :class:`AbsIrreducibeRep` specifying an
          absolute irreducible representation of the cubic Hecke algebra

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: m1 = c1.matrix()
            sage: m1._irr_to_ind(CHA2.irred_repr.W2_001)
            1
            sage: m1._irr_to_ind(CHA2.irred_repr.W3_001)
            Traceback (most recent call last):
            ...
            TypeError: representation must have 1 generators
        """
        representation_type = self.parent()._representation_type
        if not representation_type.is_split():
            raise TypeError('representation type is non split')

        ch_algebra = self.parent()._cubic_hecke_algebra
        if ch_algebra.strands() != irr.number_gens() + 1:
            raise TypeError('representation must have %s generators' % (ch_algebra.strands() - 1))

        ind = irr.gap_index()
        if representation_type == RepresentationType.SplitIrredMarin:
            ind = irr.internal_index()
        return ind

    @cached_method
    def __getitem__(self, item):
        r"""
        Return the sub-matrix block of ``self`` considered as block diagonal
        matrix specified by `item`.

        Overloading builtin-method to select a list-item.

        INPUT:

        - ``item`` -- an :class:`AbsIrreducibeRep` specifying an
          absolute irreducible representation of the cubic Hecke algebra;
          alternatively, it can be specified by list index
          (see :meth:`internal_index` respectively :meth:`gap_index`)

        OUTPUT:

        An instance of :class:`Matrix_generic_dense` representing
        the specified block of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: m1 = c1.matrix()
            sage: m1[0]                       # indirect doctest
            [a]
            sage: m1[CHA2.irred_repr.W2_001]  # indirect doctest
            [b]
        """
        if isinstance(item, AbsIrreducibeRep):
            return self._get_block(self._irr_to_ind(item))
        elif isinstance(item, (Integer, int)):
            return self._get_block(item)

        return super().__getitem__(item)

    @cached_method
    def block_diagonal_list(self):
        r"""
        Return the list of sub-matrix blocks of ``self`` considered
        as block diagonal matrix.

        OUTPUT:

        A list of instances of :class:`Matrix_generic_dense` each of
        which represents a diagonal block of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: c1.matrix().block_diagonal_list()
            [[a], [b], [-b - a + u]]
        """
        representation_type = self.parent()._representation_type
        n = self.parent()._cubic_hecke_algebra.strands()
        m = representation_type.number_of_representations(n)
        return [self._get_block(i) for i in range(m)]

    @cached_method
    def reduce_to_irr_block(self, irr):
        r"""
        Return a copy of ``self`` with zeroes outside the block corresponding to
        ``irr`` but the block according to the input identical to that of ``self``.

        INPUT:

        - ``irr`` -- an :class:`AbsIrreducibeRep` specifying an
          absolute irreducible representation of the cubic Hecke algebra;
          alternatively, it can be specified by list index (see
          :meth:`internal_index` respectively :meth:`gap_index`)

        OUTPUT:

        An instance of :class:`Matrix_generic_dense` with exactly one nonzero block
        according to ``irr``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: m1 = c1.matrix()
            sage: m1.reduce_to_irr_block(0)
            [a 0 0]
            [0 0 0]
            [0 0 0]
            sage: m1.reduce_to_irr_block(CHA2.irred_repr.W2_001)
            [0 0 0]
            [0 b 0]
            [0 0 0]
        """
        if isinstance(irr, AbsIrreducibeRep):
            ind = self._irr_to_ind(irr)
        else:
            ind = Integer(irr)
        from copy import copy
        mat_list = copy(self.parent().zero().block_diagonal_list())
        mat_list[ind] = self[ind]
        return block_diagonal_matrix(mat_list, subdivide=self.parent()._subdivide, sparse=True)


# ------------------------------------------------------------------------------------------------------------------
# Definition of CubicHeckeMatrixSpace
# --------------------------------------------------------------------------------------------------------
class CubicHeckeMatrixSpace(MatrixSpace):
    r"""
    The matrix space of cubic Hecke algebra representations.

    INPUT:

    - ``cubic_hecke_algebra`` -- (optional)
      :class:`~sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeAlgebra`
      must be given if ``element`` fails to be an instance of its element class
    - ``representation_type`` -- (default: ``RepresentationType.SplitIrredChevie``)
      :class:`RepresentationType` specifying the type of the representation
    - ``subdivide`` -- boolean (default: ``False``); whether or not to subdivide
      the resulting matrices

    - ``original`` -- boolean (default: ``False``); if ``True``, the matrix
      will have coefficients in the generic base / extension ring

    EXAMPLES::

        sage: CHA2.<c1> = algebras.CubicHecke(2)
        sage: c1.matrix()      # indirect doctest
        [         a          0          0]
        [         0          b          0]
        [         0          0 -b - a + u]
        sage: c1.matrix(original=True)
        [a 0 0]
        [0 b 0]
        [0 0 c]
        sage: c1.matrix(representation_type = CHA2.repr_type.RegularLeft)   # indirect doctest
        [ 0 -v  1]
        [ 1  u  0]
        [ 0  w  0]
    """
    @staticmethod
    def __classcall_private__(cls, cubic_hecke_algebra, representation_type=None, subdivide=False, original=False):
        r"""
        Normalize the arguments to call the ``__init__`` constructor.

        See the documentation in ``__init__``.

        TESTS::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA2)
            sage: MS2 = chmr.CubicHeckeMatrixSpace(CHA2, representation_type=CHA2.repr_type.SplitIrredMarin, subdivide=False)
            sage: MS is MS2
            True
        """
        from sage.algebras.hecke_algebras.cubic_hecke_algebra import CubicHeckeAlgebra

        if not isinstance(cubic_hecke_algebra, CubicHeckeAlgebra):
            raise TypeError('cubic_hecke_algebra must be an instance of CubicHeckeAlgebra')

        if representation_type is None:
            representation_type = RepresentationType.SplitIrredMarin

        if representation_type == RepresentationType.SplitIrredChevie:
            from sage.combinat.root_system.reflection_group_real import is_chevie_available
            if not is_chevie_available():
                raise ValueError('CHEVIE is not available')

        base_ring = cubic_hecke_algebra.base_ring(generic=original)
        dimension = cubic_hecke_algebra.dimension()
        if representation_type.is_split():
            dimension = cubic_hecke_algebra._dim_irr_rep
            base_ring = cubic_hecke_algebra.extension_ring(generic=original)
        # Bypass the MatrixSpace.__classcall__
        return super(MatrixSpace, cls).__classcall__(cls, base_ring, int(dimension),
                                                     cubic_hecke_algebra=cubic_hecke_algebra,
                                                     representation_type=representation_type,
                                                     subdivide=subdivide)

    def __init__(self, base_ring,
                 dimension,
                 cubic_hecke_algebra,
                 representation_type,
                 subdivide):
        r"""
        Initialize ``self``.

        TESTS::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA3, original=True)

        The minpoly does not work over more generic rings::

            sage: TestSuite(MS).run(skip='_test_elements')     # long time
        """
        from sage.algebras.hecke_algebras.cubic_hecke_algebra import CubicHeckeAlgebra

        if not isinstance(cubic_hecke_algebra, CubicHeckeAlgebra):
            raise TypeError('cubic_hecke_algebra must be an instance of CubicHeckeAlgebra')

        # -------------------------------------------------------------------------------------------------
        # saving input parameters
        # -------------------------------------------------------------------------------------------------
        self._cubic_hecke_algebra = cubic_hecke_algebra
        self._representation_type = representation_type
        self._subdivide = subdivide

        original_base_ring = cubic_hecke_algebra.base_ring(generic=True)

        if representation_type.is_split():
            original_base_ring = cubic_hecke_algebra.extension_ring(generic=True)
            specialize = cubic_hecke_algebra._generic_extension_ring_map
        else:
            specialize = cubic_hecke_algebra._ring_of_definition_map

        verbose("original_base_ring %s base_ring %s" % (original_base_ring, base_ring), level=2)

        self._original_base_ring = original_base_ring
        self._specialize = specialize

        super().__init__(base_ring, dimension, dimension, sparse=True, implementation=CubicHeckeMatrixRep)

    def construction(self):
        r"""
        Return ``None`` since this construction is not functorial.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: MS = c1.matrix().parent()
            sage: MS._test_category()   # indirect doctest
        """
        return None

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: MS = c1.matrix().parent()
            sage: loads(dumps(MS)) == MS      # indirect doctest
            True
        """
        original = self.base_ring() == self._original_base_ring
        return CubicHeckeMatrixSpace, (self._cubic_hecke_algebra, self._representation_type, self._subdivide, original)

    def _element_constructor_(self, x):
        r"""
        INPUT:

        - ``x`` -- an element of a
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeAlgebra`
          or an element whose parent is a :class:`MatrixSpace`

        EXAMPlES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA3, original=True)
            sage: m1 = MS._element_constructor_(c1)
            sage: isinstance(m1, MS.element_class)
            True
            sage: isinstance(MS._element_constructor_(m1), MS.element_class)
            True

            sage: m = matrix(MS.base_ring(), 12, 12, lambda i, j: 1)
            sage: MS._element_constructor_(m)
            Traceback (most recent call last):
            ...
            TypeError: incompatible block structure
        """
        # -------------------------------------------------------------------------------------------------
        # checking input and setting the self._cubic_hecke_algebra
        # -------------------------------------------------------------------------------------------------
        ch_algebra = self._cubic_hecke_algebra
        ele_parent = x.parent()
        ori_base_ring = self._original_base_ring
        if isinstance(ele_parent, MatrixSpace):
            # TODO: Find preimage in cubic hecke algebra
            d1, d2 = x.dimensions()
            if d1 != self.ncols() or d2 != self.nrows():
                raise ValueError('incompatible dimensions!')

            if ele_parent.base_ring() == ori_base_ring:
                x = self._specialize_matrix(x)
            elif ele_parent.base_ring() != self.base_ring():
                raise ValueError('incompatible base ring!')
            x_in_self = self.element_class(self, x)
            matrix_list = x_in_self.block_diagonal_list()
            matrix = block_diagonal_matrix(matrix_list, subdivide=self._subdivide, sparse=True)
            if matrix != x:
                raise TypeError('incompatible block structure')
            return self.element_class(self, matrix)

        if ele_parent == ch_algebra:
            mat = ch_algebra._apply_module_morphism(x, self._image_on_basis)
            return self(mat)

        raise TypeError('element must be an instance of CubicHeckeElement or a matrix')

    @cached_method
    def __call__(self, entries=None, coerce=True, copy=None):
        r"""
        Construct an element of ``self``.

        This method needs to be overloaded here since
        :class:`MatrixSpace` has an own implementation of it.

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA2)
            sage: MS(c1)
            [         a          0          0]
            [         0          b          0]
            [         0          0 -b - a + u]
        """
        from sage.algebras.hecke_algebras.cubic_hecke_algebra import CubicHeckeAlgebra
        if entries is None:
            return super().__call__(entries=entries, coerce=coerce, copy=copy)
        if not hasattr(entries, 'parent'):
            return super().__call__(entries=entries, coerce=coerce, copy=copy)
        ele_parent = entries.parent()
        if not isinstance(ele_parent, (CubicHeckeAlgebra, MatrixSpace)):
            return super().__call__(entries=entries, coerce=coerce, copy=copy)
        return self._element_constructor_(entries)

    @cached_method
    def _specialize_matrix(self, mat):
        r"""
        Return the given matrix specializing the original coefficients
        from data import to the base ring of ``self``.

        INPUT:

        - ``mat`` -- matrix over the original base ring

        OUTPUT: matrix over the base ring of ``self``

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA2)
            sage: B = MS._original_base_ring
            sage: a, b, c = B.gens()
            sage: mat = matrix(B, [[a, b], [0, c]])
            sage: MS._specialize_matrix(mat)
            [         a          b]
            [         0 -b - a + u]
        """
        base_ring = self.base_ring()
        original_base_ring = self._original_base_ring
        specialize = self._specialize

        if base_ring == original_base_ring:
            return mat

        mat_dict = {k: specialize(original_base_ring(v)) for k, v in mat.dict().items()}
        return matrix(base_ring, mat_dict)

    @cached_method
    def _image_on_gen(self, gen_ind):
        r"""
        Return the matrix list corresponding to the generator given by
        ``(gen_ind,)`` in Tietze form under the representation_type of
        ``self`` from the data-file or via the ``GAP3`` interface

        INPUT:

        - ``gen_ind`` -- integer; index of a generator of the cubic Hecke
          algebra attached to ``self + 1``; negative values correspond to
          the according inverses

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA3)
            sage: MS._image_on_gen(1)
            [
                           [  b   0]  [  a   0]  [  a   0]
            [a], [c], [b], [b*c   c], [a*b   b], [a*c   c],
            <BLANKLINE>
            [        c         0         0]
            [b^2 + a*c         b         0]
            [        b         1         a]
            ]

            sage: CHA2 = CHA3.cubic_hecke_subalgebra()
            sage: MSreg = chmr.CubicHeckeMatrixSpace(CHA2, representation_type=CHA2.repr_type.RegularRight)
            sage: MSreg._image_on_gen(-1)
            [
            [     0      1 (-u)/w]
            [     0      0    1/w]
            [     1      0    v/w]
            ]
        """
        representation_type = self._representation_type
        original_base_ring = self._original_base_ring
        ch_algebra = self._cubic_hecke_algebra
        n = ch_algebra.strands()

        def invert_gen(matr):
            r"""
            Return the inverse matrix of generators.
            """
            cfs = ch_algebra.cubic_equation(as_coefficients=True, generic=True)
            fac = - 1 / cfs[0]
            cf0, cf1, cf2, cf3 = (original_base_ring(cf * fac) for cf in cfs)

            matri = cf1 * matr.parent().one()
            matri += cf2 * matr
            matri += cf3 * matr**2
            d1, d2 = matr.dimensions()
            matrI = matrix(original_base_ring, d1, d2, lambda i, j: original_base_ring(matri[i, j]))
            return matrI

        if n == 2:
            if representation_type.is_split():
                # Split representations for n == 2 are missing in CHEVIE and data files
                a, b, c = original_base_ring.gens()
                matrix_list = [matrix(1, 1, [a]), matrix(1, 1, [b]), matrix(1, 1, [c])]
                if gen_ind < 0:
                    matrix_list = [invert_gen(mat) for mat in matrix_list]
                return matrix_list

        num_rep = representation_type.number_of_representations(n)

        if representation_type == RepresentationType.SplitIrredChevie:
            rep_list = [ch_algebra._fetch_matrix_list_from_chevie(i + 1) for i in range(num_rep)]
            if gen_ind > 0:
                matrix_list = [rep[gen_ind - 1] for rep in rep_list]
            else:
                matrix_list = [invert_gen(rep[-gen_ind - 1]) for rep in rep_list]
        else:
            database = ch_algebra._database
            matrix_list = database.read_matrix_representation(representation_type, gen_ind, n, original_base_ring)
        return matrix_list

    @cached_method
    def _image_on_basis(self, basis_element):
        r"""
        Return the image of the given basis element of the cubic Hecke algebra
        in ``self``.

        INPUT:

        - ``basis_element`` -- a
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeElement`
          that is a monomial

        EXAMPLES::

            sage: import sage.algebras.hecke_algebras.cubic_hecke_matrix_rep as chmr
            sage: CHA3.<c1, c2> = algebras.CubicHecke(3)
            sage: MS = chmr.CubicHeckeMatrixSpace(CHA3, original=True)
            sage: MS._image_on_basis(c1)
            [        a         0         0         0         0         0         0         0         0         0         0         0]
            [        0         c         0         0         0         0         0         0         0         0         0         0]
            [        0         0         b         0         0         0         0         0         0         0         0         0]
            [        0         0         0         b         0         0         0         0         0         0         0         0]
            [        0         0         0       b*c         c         0         0         0         0         0         0         0]
            [        0         0         0         0         0         a         0         0         0         0         0         0]
            [        0         0         0         0         0       a*b         b         0         0         0         0         0]
            [        0         0         0         0         0         0         0         a         0         0         0         0]
            [        0         0         0         0         0         0         0       a*c         c         0         0         0]
            [        0         0         0         0         0         0         0         0         0         c         0         0]
            [        0         0         0         0         0         0         0         0         0 b^2 + a*c         b         0]
            [        0         0         0         0         0         0         0         0         0         b         1         a]
        """
        representation_type = self._representation_type
        ch_algebra = self._cubic_hecke_algebra
        filecache = ch_algebra._filecache

        original_base_ring = self._original_base_ring

        ele_Tietze = basis_element.Tietze()
        matrix_list = filecache.read_matrix_representation(representation_type, ele_Tietze, original_base_ring)
        if matrix_list is None:
            verbose('not in  memory %s (Tietze %s)' % (basis_element, ele_Tietze), level=2)
            if len(ele_Tietze) == 0:
                matrix_list = ch_algebra._create_matrix_list_for_one(representation_type)
            else:
                for gen_ind in ele_Tietze:
                    gen_matrix_list = self._image_on_gen(gen_ind)
                    if matrix_list is None:
                        matrix_list = list(gen_matrix_list)
                    else:
                        for i in range(len(matrix_list)):
                            matrix_list[i] *= gen_matrix_list[i]

            filecache.write_matrix_representation(representation_type, ele_Tietze, matrix_list)
            verbose('%s saved to memory' % basis_element, level=2)

        mat = block_diagonal_matrix(matrix_list, subdivide=self._subdivide, sparse=True)
        return self._specialize_matrix(mat)

    @cached_method
    def zero(self):
        r"""
        Return the zero element of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: m1   = c1.matrix()
            sage: m1rl = c1.matrix(representation_type = CHA2.repr_type.RegularLeft)
            sage: z   = m1.parent().zero()
            sage: zrl = m1rl.parent().zero()
            sage: matrix(z) == matrix(zrl), z.is_zero(), zrl.is_zero()
            (True, True, True)
            sage: z.block_diagonal_list()
            [[0], [0], [0]]
            sage: zrl.block_diagonal_list()
            [
            [0 0 0]
            [0 0 0]
            [0 0 0]
            ]
        """
        z = self.element_class(self, super().zero())
        z._cubic_hecke_element = self._cubic_hecke_algebra.zero()
        z.set_immutable()
        return z

    @cached_method
    def one(self):
        r"""
        Return the one element of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2)
            sage: m1   = c1.matrix()
            sage: m1rl = c1.matrix(representation_type = CHA2.repr_type.RegularLeft)
            sage: o   = m1.parent().one()
            sage: orl = m1rl.parent().one()
            sage: matrix(o) == matrix(orl), o.is_one(), orl.is_one()
            (True, True, True)
            sage: o.block_diagonal_list()
            [[1], [1], [1]]
            sage: orl.block_diagonal_list()
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]
        """
        o = self.element_class(self, super().one())
        o._cubic_hecke_element = self._cubic_hecke_algebra.one()
        o.set_immutable()
        return o

    @cached_method
    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2, cubic_equation_roots=(2, 3, 5))
            sage: c1.matrix()
            [2 0 0]
            [0 3 0]
            [0 0 5]
            sage: _.parent()._an_element_()
            [ 94/3     0     0]
            [    0 187/3     0]
            [    0     0 373/3]
        """
        x = self._cubic_hecke_algebra.an_element()
        return self(x)

    @cached_method
    def some_elements(self):
        r"""
        Return a generator of elements of ``self``.

        EXAMPLES::

            sage: CHA2.<c1> = algebras.CubicHecke(2, cubic_equation_roots=(2, 3, 5))
            sage: M = c1.matrix(); M
            [2 0 0]
            [0 3 0]
            [0 0 5]
            sage: MS = M.parent()
            sage: MS.some_elements()
            (
            [ 94/3     0     0]
            [    0 187/3     0]
            [    0     0 373/3]
            )
            sage: MS.some_elements() == tuple(MS(x) for x in CHA2.some_elements())
            True
        """
        return tuple([self(x) for x in self._cubic_hecke_algebra.some_elements()])
