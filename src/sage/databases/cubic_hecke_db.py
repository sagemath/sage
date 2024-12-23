r"""
Cubic Hecke database

This module contains the class :class:`CubicHeckeDataBase` which serves as an
interface to `Ivan Marin's data files
<http://www.lamfa.u-picardie.fr/marin/representationH4-en.html>`__
with respect to the cubic Hecke algebras. The data is available via a Python
wrapper as a pip installable package `database_cubic_hecke
<https://pypi.org/project/database-cubic-hecke/>`__.
For installation hints please see the documentation there.
All data needed for the cubic Hecke algebras on less than four strands is
included in this module for demonstration purpose (see for example
:func:`read_basis`, :func:`read_irr` , ... generated with the help of
:func:`create_demo_data`).

In addition to Ivan Marin's data the package contains a function
:func:`read_markov` to obtain the coefficients of Markov traces on the
cubic Hecke algebras. This data has been precomputed with the help of
``create_markov_trace_data.py`` in the `database_cubic_hecke repository
<https://github.com/soehms/database_cubic_hecke>`__.
Again, for less than four strands, this data is includes here for
demonstration purposes.

Furthermore, this module contains the class :class:`CubicHeckeFileCache`
that enables
:class:`~sage.algebras.hecke_algebras.cubic_hecke_algebras.CubicHeckeAlgebra`
to keep intermediate results of calculations in the file system.

The enum :class:`MarkovTraceModuleBasis` serves as basis for the submodule
of linear forms on the cubic Hecke algebra on at most four strands
satisfying the Markov trace condition for its cubic Hecke subalgebras.

To use the database, you need to install the optional
:ref:`database_cubic_hecke <spkg_database_cubic_hecke>` package by the Sage
command ::

    sage -i database_cubic_hecke

EXAMPLES::

    sage: # optional - database_cubic_hecke
    sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
    sage: cha_db = CubicHeckeDataBase()
    sage: cha_db
    <sage.databases.cubic_hecke_db.CubicHeckeDataBase object at ...>

AUTHORS:

- Sebastian Oehms (2020-05): initial version
- Sebastian Oehms (2022-03): PyPi version and Markov trace functionality
"""

##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

import os
from enum import Enum

from sage.structure.sage_object import SageObject
from sage.misc.persist import _base_dumps, load
from sage.misc.temporary_file import atomic_write
from sage.misc.verbose import verbose
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.algebras.hecke_algebras.cubic_hecke_base_ring import CubicHeckeExtensionRing


# ------------------------------------------------------------------------------
# functions to convert matrices and ring elements to and from flat python
# dictionaries in order to save matrices avoiding compatibility problems with
# older or newer sage versions and to save disc space
# ------------------------------------------------------------------------------
def simplify(mat):
    r"""
    Convert a matrix to a dictionary consisting of flat Python objects.

    INPUT:

    - ``mat`` -- matrix to be converted into a ``dict``

    OUTPUT:

    A ``dict`` from which ``mat`` can be reconstructed via
    element construction. The values of the dictionary may be
    dictionaries of tuples of integers or strings.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import simplify
        sage: import sage.algebras.hecke_algebras.cubic_hecke_base_ring as chbr
        sage: ER.<a, b, c> = chbr.CubicHeckeExtensionRing()
        sage: mat = matrix(ER, [[2*a, -3], [c, 4*b*~c]]); mat
        [     2*a       -3]
        [       c 4*b*c^-1]
        sage: simplify(mat)
        {(0, 0): {(1, 0, 0): {0: 2}},
         (0, 1): {(0, 0, 0): {0: -3}},
         (1, 0): {(0, 0, 1): {0: 1}},
         (1, 1): {(0, 1, -1): {0: 4}}}
        sage: mat == matrix(ER, _)
        True
        sage: F = ER.fraction_field()
        sage: matf = mat.change_ring(F)
        sage: simplify(matf)
        {(0, 0): '2*a', (0, 1): '-3', (1, 0): 'c', (1, 1): '4*b/c'}
        sage: matf == matrix(F, _)
        True
    """
    B = mat.base_ring()
    d = mat.dict()
    if isinstance(B, CubicHeckeExtensionRing):
        # Laurent polynomial cannot be reconstructed from string
        res = {k: {tuple(j): u.monomial_coefficients()
                   for j, u in v.monomial_coefficients().items()}
               for k, v in d.items()}
    else:
        res = {k: str(v) for k, v in d.items()}
    return res


class CubicHeckeDataSection(Enum):
    r"""
    Enum for the different sections of the database.

    The following choices are possible:

    - ``basis`` -- list of basis elements
    - ``reg_left_reprs`` -- data for the left regular representation
    - ``reg_right_reprs`` -- data for the right regular representation
    - ``irr_reprs`` -- data for the split irreducible representations
    - ``markov_tr_cfs`` -- data for the coefficients of the formal Markov traces

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
        sage: cha_db = CubicHeckeDataBase()
        sage: cha_db.section
        <enum 'CubicHeckeDataSection'>
    """
    basis = 'basis'
    regular_left = 'regular_left'
    regular_right = 'regular_right'
    split_irred = 'split_irred'
    markov_tr_cfs = 'markov_tr_cfs'


# -------------------------------------------------------------------------------
# Class to supply data for the basis and matrix representation for the cubic
# Hecke algebra
# -------------------------------------------------------------------------------
class CubicHeckeDataBase(SageObject):
    r"""
    Database interface for
    :class:`~sage.algebras.hecke_algebras.cubic_hecke_algebras.CubicHeckeAlgebra`

    The original data are obtained from `Ivan Marin's web page
    <http://www.lamfa.u-picardie.fr/marin/representationH4-en.html>`__

    The data needed to work with the cubic Hecke algebras on less than 4 strands
    is completely contained in this module. Data needed for the larger algebras
    can be installed as an optional Sage package which comes as a ``pip`` installable
    `Python wrapper <https://pypi.org/project/database-cubic-hecke/>`__ of
    Ivan Marin's data. For more information see the `corresponding repository
    <https://github.com/soehms/database_cubic_hecke>`__.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
        sage: cha_db = CubicHeckeDataBase()
        sage: cha_db._feature
        Feature('database_cubic_hecke')
    """
    section = CubicHeckeDataSection

    def __init__(self):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db._data_library
            {}
        """
        from sage.features.databases import DatabaseCubicHecke
        self._feature = DatabaseCubicHecke()
        self._data_library = {}
        self._demo = None

    def version(self):
        r"""
        Return the current version of the database.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.version() > '2022.1.1'   # optional - database_cubic_hecke
            True
        """
        self._feature.require()
        from database_cubic_hecke import version
        return version()

    def demo_version(self):
        r"""
        Return whether the cubic Hecke database is installed completely or
        just the demo version is used.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.demo_version()       # optional - database_knotinfo
            False
        """
        if self._demo is None:
            if self._feature.is_present():
                self._demo = False
            else:
                self._demo = True
        return self._demo

    # --------------------------------------------------------------------------
    # read from an sobj-file obtained from Ivan Marin's database
    # --------------------------------------------------------------------------
    def read(self, section, variables=None, nstrands=4):
        r"""
        Access various static data libraries.

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeDataSection`
          to select the data to be read in

        OUTPUT: a dictionary containing the data corresponding to the section

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: basis = cha_db.read(cha_db.section.basis, nstrands=3)
            sage: len(basis)
            24
        """
        if not isinstance(section, CubicHeckeDataSection):
            raise TypeError('section must be an instance of enum %s' % CubicHeckeDataBase.section)

        data_lib = self._data_library

        nstrands = int(nstrands)
        if (section, nstrands) in data_lib.keys():
            return data_lib[(section, nstrands)]

        verbose('loading data library %s for %s strands ...' % (section.value, nstrands))

        from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import GenSign

        if self.demo_version():
            if nstrands >= 4:
                self._feature.require()
            from .cubic_hecke_db import read_basis, read_irr, read_regl, read_regr, read_markov
        else:
            from database_cubic_hecke import read_basis, read_irr, read_reg
            from database_cubic_hecke.markov_trace_coeffs import read_markov

            def read_regl(variables, num_strands):
                return read_reg(variables, num_strands=num_strands)

            def read_regr(variables, num_strands):
                return read_reg(variables, right=True, num_strands=num_strands)

        if section == CubicHeckeDataSection.basis:
            data_lib[(section, nstrands)] = read_basis(nstrands)
        elif section == CubicHeckeDataSection.markov_tr_cfs:
            keys = [k for k in MarkovTraceModuleBasis if k.strands() <= nstrands]
            res = {k: read_markov(k.name, variables, num_strands=nstrands) for k in keys}
            data_lib[(section, nstrands)] = res
        elif section == CubicHeckeDataSection.split_irred:
            dim_list, repr_list, repr_list_inv = read_irr(variables, nstrands)
            data_lib[(section, nstrands)] = {GenSign.pos: repr_list, GenSign.neg: repr_list_inv}
        else:
            if section == CubicHeckeDataSection.regular_right:
                dim_list, repr_list, repr_list_inv = read_regr(variables, nstrands)
            else:
                dim_list, repr_list, repr_list_inv = read_regl(variables, nstrands)
            data_lib[(section, nstrands)] = {GenSign.pos: repr_list, GenSign.neg: repr_list_inv}

        verbose('... finished!')
        return data_lib[(section, nstrands)]

    # --------------------------------------------------------------------------
    # matrix_reprs_from_file_cache_
    # --------------------------------------------------------------------------
    def read_matrix_representation(self, representation_type, gen_ind, nstrands, ring_of_definition):
        r"""
        Return the matrix representations from the database.

        INPUT:

        - ``representation_type`` -- an element of
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: CHA3 = algebras.CubicHecke(2)
            sage: GER = CHA3.extension_ring(generic=True)
            sage: cha_db = CHA3._database
            sage: rt = CHA3.repr_type
            sage: m1 =cha_db.read_matrix_representation(rt.SplitIrredMarin, 1, 3, GER)
            sage: len(m1)
            7
            sage: GBR = CHA3.base_ring(generic=True)
            sage: m1rl = cha_db.read_matrix_representation(rt.RegularLeft, 1, 3, GBR)
            sage: m1rl[0].dimensions()
            (24, 24)
        """
        from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import RepresentationType, GenSign
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' % RepresentationType)

        td = ring_of_definition.gens_dict_recursive()
        if 'e3' in td.keys():
            td['j'] = td['e3']
            td.pop('e3')
        v = tuple(td.values())

        num_rep = representation_type.number_of_representations(nstrands)
        rep_list = self.read(representation_type.data_section(), variables=v, nstrands=nstrands)
        if gen_ind > 0:
            rep_list = [rep_list[GenSign.pos][i] for i in range(num_rep)]
            matrix_list = [matrix(ring_of_definition, rep[gen_ind-1], sparse=True) for rep in rep_list]
        else:
            # data of inverse of generators is stored under negative strand-index
            rep_list = [rep_list[GenSign.neg][i] for i in range(num_rep)]
            matrix_list = [matrix(ring_of_definition, rep[-gen_ind-1], sparse=True) for rep in rep_list]
        for m in matrix_list:
            m.set_immutable()
        return matrix_list


class MarkovTraceModuleBasis(Enum):
    r"""
    Enum for the basis elements for the Markov trace module.

    The choice of the basis elements doesn't have a systematically background
    apart from generating the submodule of maximal rank in the module of linear
    forms on the cubic Hecke algebra for which the Markov trace condition with
    respect to its cubic Hecke subalgebras hold. The number of crossings in
    the corresponding links is chosen as minimal as possible.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
        sage: MarkovTraceModuleBasis.K92.description()
        'knot 9_34'
    """
    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U2    # indirect doctest
            U2
        """
        return self.name

    def __gt__(self, other):
        r"""
        Implement comparison of different items in order to have ``sorted`` work.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: sorted(MarkovTraceModuleBasis)
            [U1, U2, U3, K4, U4, K4U, K6, K7, K91, K92]
        """
        if self.__class__ is other.__class__:
            tups = (self.strands(), len(self.braid_tietze()), self.name)
            tupo = (other.strands(), len(other.braid_tietze()), other.name)
            return tups > tupo
        return NotImplemented

    def strands(self):
        r"""
        Return the number of strands of the minimal braid representative
        of the link corresponding to ``self``.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.K7.strands()
            4
        """
        return self.value[1]

    def braid_tietze(self, strands_embed=None):
        r"""
        Return the Tietze representation of the braid corresponding to this basis
        element.

        INPUT:

        - ``strands_embed`` -- (optional) the number of strands of the braid
          if strands should be added

        OUTPUT: a tuple representing the braid in Tietze form

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U2.braid_tietze()
            ()
            sage: MarkovTraceModuleBasis.U2.braid_tietze(strands_embed=4)
            (2, 3)
        """
        if not strands_embed:
            strands_embed = self.strands()

        if strands_embed > self.strands():
            last_gen = strands_embed-1
            return self.braid_tietze(strands_embed=last_gen) + (last_gen,)

        return self.value[2]

    def writhe(self):
        r"""
        Return the writhe of the link corresponding to this basis element.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.K4.writhe()
            0
            sage: MarkovTraceModuleBasis.K6.writhe()
            1
        """
        from sage.functions.generalized import sign
        return sum(sign(t) for t in self.braid_tietze())

    def description(self):
        r"""
        Return a description of the link corresponding to this basis element.

        In the case of knots it refers to the naming according to
        `KnotInfo <https://knotinfo.math.indiana.edu/>`__.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U3.description()
            'three unlinks'
        """
        return self.value[0]

    def link(self):
        r"""
        Return the :class:`Link` that represents this basis element.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U1.link()
            Link with 1 component represented by 0 crossings
            sage: MarkovTraceModuleBasis.K4.link()
            Link with 1 component represented by 4 crossings
        """
        from sage.knots.link import Link
        pd_code = self.value[3]
        if pd_code is not None:
            # since :class:`Link` does not construct disjoint union of unlinks
            # from the braid representation, we need a pd_code here
            return Link(pd_code)
        else:
            from sage.groups.braid import BraidGroup
            B = BraidGroup(self.strands())
            return Link(B(self.braid_tietze()))

    def regular_homfly_polynomial(self):
        r"""
        Return the regular variant of the HOMFLY-PT polynomial of the link that
        represents this basis element.

        This is the HOMFLY-PT polynomial renormalized by the writhe factor
        such that it is an invariant of regular isotopy.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U1.regular_homfly_polynomial()
            1
            sage: u2 = MarkovTraceModuleBasis.U2.regular_homfly_polynomial(); u2
            -L*M^-1 - L^-1*M^-1
            sage: u2**2 == MarkovTraceModuleBasis.U3.regular_homfly_polynomial()
            True
            sage: u2**3 == MarkovTraceModuleBasis.U4.regular_homfly_polynomial()
            True
        """
        H = self.link().homfly_polynomial()
        L, M = H.parent().gens()
        return H * L**self.writhe()

    def regular_kauffman_polynomial(self):
        r"""
        Return the regular variant of the Kauffman polynomial of the link that
        represents this basis element.

        This is the Kauffman polynomial renormalized by the writhe factor
        such that it is an invariant of regular isotopy.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U1.regular_homfly_polynomial()
            1
            sage: u2 = MarkovTraceModuleBasis.U2.regular_kauffman_polynomial(); u2
            a*z^-1 - 1 + a^-1*z^-1
            sage: u2**2 == MarkovTraceModuleBasis.U3.regular_kauffman_polynomial()
            True
            sage: u2**3 == MarkovTraceModuleBasis.U4.regular_kauffman_polynomial()
            True
        """
        from sage.knots.knotinfo import KnotInfo
        K = KnotInfo.L2a1_1.kauffman_polynomial().parent()
        a, z = K.gens()
        d = kauffman[self.name]
        if d:
            return K(d)*a**self.writhe()
        U2rkp = MarkovTraceModuleBasis.U2.regular_kauffman_polynomial()
        if self.name == 'K4U':
            K4rkp = MarkovTraceModuleBasis.K4.regular_kauffman_polynomial()
            return K4rkp * U2rkp
        exp = self.strands() - 1
        return U2rkp**exp

    def links_gould_polynomial(self):
        r"""
        Return the Links-Gould polynomial of the link that represents this
        basis element.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import MarkovTraceModuleBasis
            sage: MarkovTraceModuleBasis.U1.links_gould_polynomial()
            1
            sage: MarkovTraceModuleBasis.U2.links_gould_polynomial()
            0
            sage: MarkovTraceModuleBasis.K4.links_gould_polynomial()
            2*t0*t1 - 3*t0 - 3*t1 + t0*t1^-1 + 7 + t0^-1*t1
            - 3*t1^-1 - 3*t0^-1 + 2*t0^-1*t1^-1
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        R = LaurentPolynomialRing(ZZ, 't0, t1')
        return R(links_gould[self.name])

    U1 = ['one unlink',    1, (), []]
    U2 = ['two unlinks',   2, (), [[3, 1, 4, 2], [4, 1, 3, 2]]]
    U3 = ['three unlinks', 3, (), [[3, 7, 4, 8], [4, 7, 5, 8],
                                   [5, 1, 6, 2], [6, 1, 3, 2]]]
    U4 = ['four unlinks',  4, (), [[3, 9, 4, 10], [4, 9, 5, 10], [5, 11, 6, 12],
                                   [6, 11, 7, 12], [7, 1, 8, 2], [8, 1, 3, 2]]]
    K4U = ['knot 4_1 plus one unlink', 4, (1, -2, 1, -2),
           [[3, 8, 4, 9], [9, 7, 10, 6], [7, 4, 8, 5], [5, 11, 6, 10],
           [11, 1, 12, 2], [12, 1, 3, 2]]]
    K4 = ['knot 4_1',  3, (1, -2, 1, -2), None]
    K6 = ['knot 6_1',  4, (1, 1, 2, -1, -3, 2, -3), None]
    K7 = ['knot 7_4',  4, (1, 1, 2, -1, 2, 2, 3, -2, 3), None]
    K91 = ['knot 9_29', 4, (1, -2, -2, 3, -2, 1, -2, 3, -2), None]
    K92 = ['knot 9_34', 4, (-1, 2, -1, 2, -3, 2, -1, 2, -3), None]


kauffman = {
 'U1': 1,
 'U2': {(1, -1): 1, (0, 0): -1, (-1, -1): 1},
 'U3': None,
 'U4': None,
 'K4U': None,
 'K4': {(2, 2): 1, (1, 3): 1, (2, 0): -1, (1, 1): -1, (0, 2): 2, (-1, 3): 1,
        (0, 0): -1, (-1, 1): -1, (-2, 2): 1, (-2, 0): -1},
 'K6': {(2, 2): 1, (1, 3): 1, (0, 4): 1, (-1, 5): 1, (2, 0): -1, (-1, 3): -2,
        (-2, 4): 2, (-3, 5): 1, (-1, 1): 2, (-2, 2): -4, (-3, 3): -3, (-4, 4): 1,
        (-2, 0): 1, (-3, 1): 2, (-4, 2): -3, (-4, 0): 1},
 'K7': {(-2, 2): 1, (-3, 3): 2, (-4, 4): 3, (-5, 5): 2, (-6, 6): 1, (-4, 2): -4,
        (-5, 3): -2, (-7, 5): 3, (-8, 6): 1, (-4, 0): 2, (-6, 2): -3, (-7, 3): -8,
        (-8, 4): -3, (-9, 5): 1, (-7, 1): 4, (-8, 2): 2, (-9, 3): -4, (-8, 0): -1,
        (-9, 1): 4},
 'K91': {(7, 3): 1, (6, 4): 3, (5, 5): 6, (4, 6): 8, (3, 7): 6, (2, 8): 2,
         (5, 3): -5, (4, 4): -13, (3, 5): -8, (2, 6): 6, (1, 7): 9, (0, 8): 2,
         (5, 1): 2, (4, 2): 8, (3, 3): -1, (2, 4): -24, (1, 5): -24, (0, 6): -1,
         (-1, 7): 3, (4, 0): -2, (3, 1): 2, (2, 2): 17, (1, 3): 14, (0, 4): -11,
         (-1, 5): -10, (-2, 6): 1, (2, 0): -5, (1, 1): -1, (0, 2): 12, (-1, 3): 9,
         (-2, 4): -3, (0, 0): -3, (-1, 1): -1, (-2, 2): 3, (-2, 0): -1},
 'K92': {(5, 5): 1, (4, 6): 4, (3, 7): 6, (2, 8): 3, (5, 3): -1, (4, 4): -7,
         (3, 5): -11, (2, 6): 5, (1, 7): 14, (0, 8): 3, (4, 2): 3, (3, 3): 5,
         (2, 4): -19, (1, 5): -26, (0, 6): 9, (-1, 7): 8, (2, 2): 10, (1, 3): 12,
         (0, 4): -23, (-1, 5): -10, (-2, 6): 8, (2, 0): -1, (1, 1): -1,
         (0, 2): 11, (-1, 3): 4, (-2, 4): -10, (-3, 5): 4, (0, 0): -1,
         (-1, 1): -1, (-2, 2): 4, (-3, 3): -2, (-4, 4): 1, (-2, 0): -1}}


links_gould = {
 'U1': 1,
 'U2': 0,
 'U3': 0,
 'U4': 0,
 'K4U': 0,
 'K4': {(1, 1): 2, (1, 0): -3, (0, 1): -3, (1, -1): 1, (0, 0): 7, (-1, 1): 1,
        (0, -1): -3, (-1, 0): -3, (-1, -1): 2},
 'K6': {(2, 2): 2, (2, 1): -3, (1, 2): -3, (2, 0): 1, (1, 1): 10, (0, 2): 1,
        (1, 0): -10, (0, 1): -10, (1, -1): 3, (0, 0): 17, (-1, 1): 3, (0, -1): -7,
        (-1, 0): -7, (-1, -1): 4},
 'K7': {(4, 3): -1, (3, 4): -1, (4, 2): 1, (3, 3): 6, (2, 4): 1, (3, 2): -11,
        (2, 3): -11, (3, 1): 6, (2, 2): 28, (1, 3): 6, (2, 1): -27, (1, 2): -27,
        (2, 0): 9, (1, 1): 38, (0, 2): 9, (1, 0): -17, (0, 1): -17, (0, 0): 9},
 'K91': {(2, 2): 6, (2, 1): -20, (1, 2): -20, (2, 0): 29, (1, 1): 76, (0, 2): 29,
         (2, -1): -25, (1, 0): -123, (0, 1): -123, (-1, 2): -25, (2, -2): 14,
         (1, -1): 116, (0, 0): 217, (-1, 1): 116, (-2, 2): 14, (2, -3): -5,
         (1, -2): -71, (0, -1): -216, (-1, 0): -216, (-2, 1): -71, (-3, 2): -5,
         (2, -4): 1, (1, -3): 27, (0, -2): 136, (-1, -1): 214, (-2, 0): 136,
         (-3, 1): 27, (-4, 2): 1, (1, -4): -5, (0, -3): -50, (-1, -2): -122,
         (-2, -1): -122, (-3, 0): -50, (-4, 1): -5, (0, -4): 8, (-1, -3): 37,
         (-2, -2): 52, (-3, -1): 37, (-4, 0): 8, (-1, -4): -4, (-2, -3): -9,
         (-3, -2): -9, (-4, -1): -4},
 'K92': {(3, 1): 6, (2, 2): 12, (1, 3): 6, (3, 0): -15, (2, 1): -63, (1, 2): -63,
         (0, 3): -15, (3, -1): 14, (2, 0): 112, (1, 1): 216, (0, 2): 112,
         (-1, 3): 14, (3, -2): -6, (2, -1): -92, (1, 0): -334, (0, 1): -334,
         (-1, 2): -92, (-2, 3): -6, (3, -3): 1, (2, -2): 37, (1, -1): 262,
         (0, 0): 503, (-1, 1): 262, (-2, 2): 37, (-3, 3): 1, (2, -3): -6,
         (1, -2): -104, (0, -1): -400, (-1, 0): -400, (-2, 1): -104, (-3, 2): -6,
         (1, -3): 17, (0, -2): 162, (-1, -1): 330, (-2, 0): 162, (-3, 1): 17,
         (0, -3): -27, (-1, -2): -136, (-2, -1): -136, (-3, 0): -27, (-1, -3): 22,
         (-2, -2): 54, (-3, -1): 22, (-2, -3): -7, (-3, -2): -7}}


class CubicHeckeFileCache(SageObject):
    """
    A class to cache calculations of
    :class:`~sage.algebras.hecke_algebras.cubic_hecke_algebras.CubicHeckeAlgebra`
    in the local file system.
    """

    class section(Enum):
        r"""
        Enum for the different sections of file cache. The following choices are
        possible:

        - ``matrix_representations`` -- file cache for representation matrices
          of basis elements
        - ``braid_images`` -- file cache for images of braids
        - ``basis_extensions`` -- file cache for a dynamical growing basis used
          in the case of cubic Hecke algebras on more than 4 strands
        - ``markov_trace`` -- file cache for intermediate results of long
          calculations in order to recover the results already obtained by
          previous attemps of calculation until the corresponding intermediate
          step

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: cha_fc = CubicHeckeFileCache(CHA2)
            sage: cha_fc.section
            <enum 'section'>
        """
        def filename(self, nstrands=None):
            r"""
            Return the file name under which the data of this file cache section
            is stored as an sobj-file.

            INPUT:

            - ``nstrands`` -- (optional) :class:`Integer`; number of strands of
              the underlying braid group if the data file depends on it

            EXAMPLES::

                sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
                sage: CHA2 = algebras.CubicHecke(2)
                sage: cha_fc = CubicHeckeFileCache(CHA2)
                sage: cha_fc.section.matrix_representations.filename(2)
                'matrix_representations_2.sobj'
                sage: cha_fc.section.braid_images.filename(2)
                'braid_images_2.sobj'
            """
            if nstrands is None:
                return '%s.sobj' % self.value
            else:
                return '%s_%s.sobj' % (self.value, nstrands)

        matrix_representations = 'matrix_representations'
        braid_images = 'braid_images'
        basis_extensions = 'basis_extensions'
        markov_trace = 'markov_trace'

    def __init__(self, num_strands):
        r"""
        Initialize ``self``.

        INPUT:

        - ``num_strands`` -- integer giving the number of strands of the
          corresponding cubic Hecke algebra

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: cha_fc._file_cache_path.endswith('cubic_hecke')
            True
        """
        self._nstrands = num_strands

        from sage.env import DOT_SAGE
        self._file_cache_path = os.path.join(DOT_SAGE, 'cubic_hecke')
        self._data_library = {}
        os.makedirs(self._file_cache_path, exist_ok=True)

    def _warn_incompatibility(self, fname):
        """
        Warn the user that he has an incomaptible file cache under ``Sage_DOT``
        and move it away to another file (marked with timestamp).

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: path = cha2_fc._file_cache_path
            sage: fname = os.path.join(path, 'test')
            sage: os.system('touch %s' % fname)
            0
            sage: new_fname = cha2_fc._warn_incompatibility(fname)
            doctest:...: UserWarning: incompatible file cache ...test has been saved to ...test_...
            sage: os.remove(new_fname)
        """
        from warnings import warn
        from datetime import date
        today = date.today()
        new_fname = '%s_%s' % (fname, today)
        os.rename(fname, new_fname)
        warn('incompatible file cache %s has been saved to %s' % (fname, new_fname))
        return new_fname

    def reset_library(self, section=None):
        r"""
        Reset the file cache corresponding to the specified ``section``.

        INPUT:

        - ``section`` -- an element of :class:`CubicHeckeFileCache.section`
          to select the section of the file cache or ``None`` (default)
          meaning all sections

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.read(cha2_fc.section.braid_images)
            {}
            sage: cha2_fc.reset_library(cha2_fc.section.matrix_representations)
            sage: data_mat = cha2_fc.read(cha2_fc.section.matrix_representations)
            sage: len(data_mat.keys())
            4
        """
        if section is None:
            for sec in self.section:
                self.reset_library(section=sec)
            return

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' % CubicHeckeFileCache.section)

        from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import RepresentationType
        data_lib = self._data_library
        empty_dict = {}
        if section == self.section.matrix_representations:
            for rep_type in RepresentationType:
                new_dict = {}
                empty_dict.update({rep_type.name: new_dict})
        elif section == self.section.basis_extensions:
            empty_dict = []
        data_lib.update({section: empty_dict})

    def is_empty(self, section=None):
        r"""
        Return ``True`` if the cache of the given ``section`` is empty.

        INPUT:

        - ``section`` -- an element of :class:`CubicHeckeFileCache.section`
          to select the section of the file cache or ``None`` (default)
          meaning all sections

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library()
            sage: cha2_fc.is_empty()
            True
        """
        if section is None:
            return all(self.is_empty(section=sec) for sec in self.section)

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' % CubicHeckeFileCache.section)

        self.read(section)
        data_lib = self._data_library[section]
        from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import RepresentationType
        if section == self.section.matrix_representations:
            for rep_type in RepresentationType:
                if len(data_lib[rep_type.name]) > 0:
                    return False
            return True

        if section == self.section.basis_extensions and self._nstrands > 4:
            # the new generators and their inverses are not counted
            # since they are added during initialization
            return len(data_lib) <= 2*(self._nstrands - 4)
        return not data_lib

    # --------------------------------------------------------------------------
    # save data file system
    # --------------------------------------------------------------------------
    def write(self, section=None):
        r"""
        Write data from memory to the file system.

        INPUT:

        - ``section`` -- an element of :class:`CubicHeckeFileCache.section`
          specifying the section where the corresponding cached data belong to;
          if omitted, the data of all sections is written to the file system

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.write(cha2_fc.section.braid_images)
        """
        data_lib = self._data_library
        lib_path = self._file_cache_path

        if section is None:
            for sec in self.section:
                if sec in data_lib.keys():
                    self.write(section=sec)
            return

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' % CubicHeckeFileCache.section)

        if section not in data_lib.keys():
            raise ValueError("No data for file %s in memory" % section)

        verbose('saving file cache %s ...' % section)
        fname = os.path.join(lib_path, section.filename(self._nstrands))
        with atomic_write(fname, binary=True) as f:
            f.write(_base_dumps(data_lib[section]))
            f.close()

    # --------------------------------------------------------------------------
    # read from file system
    # --------------------------------------------------------------------------
    def read(self, section):
        r"""
        Read data into memory from the file system.

        INPUT:

        - ``section`` -- an element of :class:`CubicHeckeFileCache.section`
          specifying the section where the corresponding cached data belong to

        OUTPUT:

        Dictionary containing the data library corresponding to the section
        of file cache

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.read(cha2_fc.section.braid_images)
            {}
        """
        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' % CubicHeckeFileCache.section)

        data_lib = self._data_library
        lib_path = self._file_cache_path

        if section in data_lib.keys():
            return data_lib[section]

        verbose('loading file cache %s ...' % section)
        fname = os.path.join(lib_path, section.filename(self._nstrands))
        try:
            data_lib[section] = load(fname)
            verbose('... finished!')
        except OSError:
            self.reset_library(section)
            verbose('... not found!')
        except (ImportError, ModuleNotFoundError):
            self._warn_incompatibility(fname)
            self.reset_library(section)

        return data_lib[section]

    # --------------------------------------------------------------------------
    # read matrix representation from file cache
    # --------------------------------------------------------------------------
    def read_matrix_representation(self, representation_type, monomial_tietze, ring_of_definition):
        r"""
        Return the matrix representations of the given monomial (in Tietze form)
        if it has been stored in the file cache before.
        INPUT:

        - ``representation_type`` -- an element of
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation
        - ``monomial_tietze`` -- tuple representing the braid in Tietze form
        - ``ring_of_definition`` -- instance of
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_base_ring.CubicHeckeRingOfDefinition`
          respectively
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_base_ring.CubicHeckeExtensionRing`
          depending on whether ``representation_type`` is split or not

        OUTPUT:

        Dictionary containing all matrix representations of ``self`` of the
        given ``representation_type`` which have been stored in the file cache.
        Otherwise ``None`` is returned.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: R = CHA2.base_ring(generic=True)
            sage: cha_fc = CHA2._filecache
            sage: g, = CHA2.gens(); gt = g.Tietze()
            sage: rt = CHA2.repr_type
            sage: g.matrix(representation_type=rt.RegularLeft)
            [ 0 -v  1]
            [ 1  u  0]
            [ 0  w  0]
            sage: [_] == cha_fc.read_matrix_representation(rt.RegularLeft, gt, R)
            True
            sage: cha_fc.reset_library(cha_fc.section.matrix_representations)
            sage: cha_fc.write(cha_fc.section.matrix_representations)
            sage: cha_fc.read_matrix_representation(rt.RegularLeft, gt, R) == None
            True
        """
        from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import RepresentationType
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' % RepresentationType)

        matrix_representations = self.read(self.section.matrix_representations)[representation_type.name]
        if monomial_tietze in matrix_representations.keys():
            matrix_list_dict = matrix_representations[monomial_tietze]
            matrix_list = [matrix(ring_of_definition, mat_dict, sparse=True) for mat_dict in matrix_list_dict]
            for m in matrix_list:
                m.set_immutable()
            return matrix_list
        return None

    # --------------------------------------------------------------------------
    # matrix_representation to file cache
    # --------------------------------------------------------------------------
    def write_matrix_representation(self, representation_type, monomial_tietze, matrix_list):
        r"""
        Write the matrix representation of a monomial to the file cache.

        INPUT:

        - ``representation_type`` -- an element of
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation
        - ``monomial_tietze`` -- tuple representing the braid in Tietze form
        - ``matrix_list`` -- list of matrices corresponding to the irreducible
          representations

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: R = CHA2.base_ring(generic=True)
            sage: cha_fc = CHA2._filecache
            sage: g, = CHA2.gens(); gi = ~g; git = gi.Tietze()
            sage: rt = CHA2.repr_type
            sage: m = gi.matrix(representation_type=rt.RegularRight)
            sage: cha_fc.read_matrix_representation(rt.RegularRight, git, R)
            [
            [     0      1 (-u)/w]
            [     0      0    1/w]
            [     1      0    v/w]
            ]
            sage: CHA2.reset_filecache(cha_fc.section.matrix_representations)
            sage: cha_fc.read_matrix_representation(rt.RegularLeft, git, R) == None
            True
            sage: cha_fc.write_matrix_representation(rt.RegularRight, git, [m])
            sage: [m] == cha_fc.read_matrix_representation(rt.RegularRight, git, R)
            True
        """
        from sage.algebras.hecke_algebras.cubic_hecke_matrix_rep import RepresentationType
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' % RepresentationType)

        sec = self.section.matrix_representations
        all_matrix_representations = self.read(sec)
        if representation_type.name not in all_matrix_representations.keys():
            # old file-cache is not compatible with current dictionary keys.
            fname = os.path.join(self._file_cache_path, sec.filename(self._nstrands))
            self._warn_incompatibility(fname)
            all_matrix_representations = self.read(sec)

        matrix_representations = all_matrix_representations[representation_type.name]

        if monomial_tietze in matrix_representations.keys():
            # entry already registered
            return

        matrix_representation_dict = [simplify(mat) for mat in list(matrix_list)]
        matrix_representations[monomial_tietze] = matrix_representation_dict

        self.write(sec)
        return

    # --------------------------------------------------------------------------
    # read braid images from file cache
    # --------------------------------------------------------------------------
    def read_braid_image(self, braid_tietze, ring_of_definition):
        r"""
        Return the list of pre calculated braid images from file cache.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid in Tietze form
        - ``ring_of_definition`` -- a
          :class:`~sage.algebras.hecke_algebras.cubic_hecke_base_ring.CubicHeckeRingOfDefinition`

        OUTPUT:

        A dictionary containing the pre calculated braid image of the given
        braid.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: ring_of_definition = CHA2.base_ring(generic=True)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens(); b2 = b**2
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.braid_images)
            True
            sage: b2_img = CHA2(b2); b2_img
            w*c^-1 + u*c - v
            sage: cha_fc.write_braid_image(b2.Tietze(), b2_img.to_vector())
            sage: cha_fc.read_braid_image(b2.Tietze(), ring_of_definition)
            (-v, u, w)
        """
        braid_images = self.read(self.section.braid_images)
        if braid_tietze in braid_images.keys():
            braid_image = braid_images[braid_tietze]
            result_list = [ring_of_definition(cf) for cf in list(braid_image)]
            from sage.modules.free_module_element import vector
            return vector(ring_of_definition, result_list)
        return None

    # --------------------------------------------------------------------------
    # braid image to_file cache
    # --------------------------------------------------------------------------
    def write_braid_image(self, braid_tietze, braid_image_vect):
        r"""
        Write the braid image of the given braid to the file cache.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid in Tietze form
        - ``braid_image_vect`` -- image of the given braid as a vector with
          respect to the basis of the cubic Hecke algebra

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: ring_of_definition = CHA2.base_ring(generic=True)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens(); b3 = b**3
            sage: b3_img = CHA2(b3); b3_img
            u*w*c^-1 + (u^2-v)*c - (u*v-w)
            sage: cha_fc.write_braid_image(b3.Tietze(), b3_img.to_vector())
            sage: cha_fc.read_braid_image(b3.Tietze(), ring_of_definition)
            (-u*v + w, u^2 - v, u*w)
            sage: cha_fc.reset_library(CubicHeckeFileCache.section.braid_images)
            sage: cha_fc.write(CubicHeckeFileCache.section.braid_images)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.braid_images)
            True
        """
        braid_images = self.read(self.section.braid_images)

        if braid_tietze in braid_images.keys():
            # entry already registered
            return

        braid_image_dict = [str(cf) for cf in list(braid_image_vect)]
        braid_images[braid_tietze] = braid_image_dict

        self.write(self.section.braid_images)
        return

    # --------------------------------------------------------------------------
    # basis to file cache
    # --------------------------------------------------------------------------
    def update_basis_extensions(self, new_basis_extensions):
        r"""
        Update the file cache for basis extensions for cubic Hecke algebras on
        more than 4 strands according to the given ``new_basis_extensions``.

        INPUT:

        - ``new_basis_extensions`` -- list of additional (to the static basis)
          basis elements which should replace the former such list in the file

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.basis_extensions)
            True
            sage: cha_fc.update_basis_extensions([(1,), (-1,)])
            sage: cha_fc.read(CubicHeckeFileCache.section.basis_extensions)
            [(1,), (-1,)]
            sage: cha_fc.reset_library(CubicHeckeFileCache.section.basis_extensions)
            sage: cha_fc.write(CubicHeckeFileCache.section.basis_extensions)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.basis_extensions)
            True
        """
        self._data_library.update({self.section.basis_extensions: new_basis_extensions})
        self.write(self.section.basis_extensions)
        return


# -----------------------------------------------------------------------------
# Demo data section
# -----------------------------------------------------------------------------

func_name = 'read_%s'

var_decl = "\n    %s = variables"
var_doc_input = "\n    - ``variables`` -- tuple containing the indeterminates of the representation"
var_doc_decl = "\n        sage: L.<%s> = LaurentPolynomialRing(ZZ)"

template = """def %s(%snum_strands=3):
    %s%s
    data = {}
    data[2] = %s
    data[3] = %s
    return data[num_strands]
"""


doc = r"""{}
    Return precomputed data of Ivan Marin.

    This code was generated by :func:`create_demo_data`, please do not edit.

    INPUT:
%s
    - ``num_strands`` -- integer; number of strands of the cubic Hecke algebra

    EXAMPLES::

        {}age: from sage.databases.cubic_hecke_db import %s%s
        {}age: %s(%s2)
    {}
""".format('r"""', 's', 's', '"""')  # s in the middle to hide these lines from _test_enough_doctests


def create_demo_data(filename='demo_data.py'):
    r"""
    Generate code for the functions inserted below to access the
    small cases of less than 4 strands as demonstration cases.

    The code is written to a file with the given name from where
    it can be copied into this source file (in case a change is needed).

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import create_demo_data
        sage: create_demo_data()   # not tested
    """
    # ---------------------------------------------------------------
    # preparations
    # ---------------------------------------------------------------
    def create_repr_func(name, variables, data2, data3):
        fname = func_name % name
        if variables:
            v = str(variables)
            vars2 = v + ', '
            decl = var_decl % v
            doc_dec = var_doc_decl % v[1: len(v)-1]
            doc_str = doc % (var_doc_input, fname,  doc_dec, fname, vars2)
            res = template % (fname, 'variables, ', doc_str, decl, data2,  data3)
        else:
            doc_str = doc % ('', fname, '', fname, '')
            res = template % (fname, '', doc_str, '', data2,  data3)
        return res

    from textwrap import fill

    def fl(s):
        return fill(str(s), subsequent_indent='        ')

    from sympy import var
    from database_cubic_hecke import read_basis, read_irr, read_reg

    vari = var('a, b, c, j')
    varr = var('u, v, w')

    # ---------------------------------------------------------------
    # read data
    # ---------------------------------------------------------------
    bas2 = fl(read_basis(num_strands=2))
    bas3 = fl(read_basis(num_strands=3))

    irr2 = fl(read_irr(variables=vari, num_strands=2))
    irr3 = fl(read_irr(variables=vari, num_strands=3))

    regl2 = fl(read_reg(variables=varr, num_strands=2))
    regl3 = fl(read_reg(variables=varr, num_strands=3))

    regr2 = fl(read_reg(variables=varr, num_strands=2, right=True))
    regr3 = fl(read_reg(variables=varr, num_strands=3, right=True))

    # ---------------------------------------------------------------
    # create functions and write them to file
    # ---------------------------------------------------------------
    bas = create_repr_func('basis',  '',   bas2,  bas3)
    irr = create_repr_func('irr',    vari, irr2,  irr3)
    regl = create_repr_func('regl',   varr, regl2, regl3)
    regr = create_repr_func('regr',   varr, regr2, regr3)

    with open(filename, 'w') as f:
        f.write(bas)
        f.write(irr)
        f.write(regl)
        f.write(regr)


def read_basis(num_strands=3):
    r"""
    Return precomputed data of Ivan Marin.

    This code was generated by :func:`create_demo_data`, please do not edit.

    INPUT:

    - ``num_strands`` -- integer; number of strands of the cubic Hecke algebra

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_basis
        sage: read_basis(2)
        [[], [1], [-1]]
    """
    data = {}
    data[2] = [[], [1], [-1]]
    data[3] = [[], [1], [-1], [2], [-2], [1, 2], [1, -2], [-1, 2], [-1, -2], [1, 2,
        1], [1, 2, -1], [-1, 2, 1], [-1, 2, -1], [1, -2, 1], [-1, -2,
        1], [2, 1], [-2, 1], [2, -1], [-2, -1], [1, -2, -1], [-1, -2,
        -1], [2, -1, 2], [1, 2, -1, 2], [-1, 2, -1, 2]]
    return data[num_strands]


def read_irr(variables, num_strands=3):
    r"""
    Return precomputed data of Ivan Marin.

    This code was generated by :func:`create_demo_data`, please do not edit.

    INPUT:

    - ``variables`` -- tuple containing the indeterminates of the representation
    - ``num_strands`` -- integer; number of strands of the cubic Hecke algebra

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_irr
        sage: L.<a, b, c, j> = LaurentPolynomialRing(ZZ)
        sage: read_irr((a, b, c, j), 2)
        ([1, 1, 1],
        [[{(0, 0): a}], [{(0, 0): c}], [{(0, 0): b}]],
        [[{(0, 0): a^-1}], [{(0, 0): c^-1}], [{(0, 0): b^-1}]])
    """
    (a, b, c, j) = variables
    data = {}
    data[2] = ([1, 1, 1], [[{(0, 0): a}], [{(0, 0): c}], [{(0, 0): b}]], [[{(0, 0):
        1/a}], [{(0, 0): 1/c}], [{(0, 0): 1/b}]])
    data[3] = ([1, 1, 1, 2, 2, 2, 3], [[{(0, 0): a}, {(0, 0): a}], [{(0, 0): c},
        {(0, 0): c}], [{(0, 0): b}, {(0, 0): b}], [{(0, 0): b, (1, 0):
        b*c, (1, 1): c}, {(0, 0): c, (0, 1): -1, (1, 1): b}], [{(0,
        0): a, (1, 0): a*b, (1, 1): b}, {(0, 0): b, (0, 1): -1, (1,
        1): a}], [{(0, 0): a, (1, 0): a*c, (1, 1): c}, {(0, 0): c, (0,
        1): -1, (1, 1): a}], [{(0, 0): c, (1, 0): a*c + b**2, (1, 1):
        b, (2, 0): b, (2, 1): 1, (2, 2): a}, {(0, 0): a, (0, 1): -1,
        (0, 2): b, (1, 1): b, (1, 2): -a*c - b**2, (2, 2): c}]],
        [[{(0, 0): 1/a}, {(0, 0): 1/a}], [{(0, 0): 1/c}, {(0, 0):
        1/c}], [{(0, 0): 1/b}, {(0, 0): 1/b}], [{(0, 0): 1/b, (1, 0):
        -1, (1, 1): 1/c}, {(0, 0): 1/c, (0, 1): 1/(b*c), (1, 1):
        1/b}], [{(0, 0): 1/a, (1, 0): -1, (1, 1): 1/b}, {(0, 0): 1/b,
        (0, 1): 1/(a*b), (1, 1): 1/a}], [{(0, 0): 1/a, (1, 0): -1, (1,
        1): 1/c}, {(0, 0): 1/c, (0, 1): 1/(a*c), (1, 1): 1/a}], [{(0,
        0): 1/c, (1, 0): -a/b - b/c, (1, 1): 1/b, (2, 0): 1/b, (2, 1):
        -1/(a*b), (2, 2): 1/a}, {(0, 0): 1/a, (0, 1): 1/(a*b), (0, 2):
        1/b, (1, 1): 1/b, (1, 2): a/b + b/c, (2, 2): 1/c}]])
    return data[num_strands]


def read_regl(variables, num_strands=3):
    r"""
    Return precomputed data of Ivan Marin.

    This code was generated by :func:`create_demo_data`, please do not edit.

    INPUT:

    - ``variables`` -- tuple containing the indeterminates of the representation
    - ``num_strands`` -- integer; number of strands of the cubic Hecke algebra

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_regl
        sage: L.<u, v, w> = LaurentPolynomialRing(ZZ)
        sage: read_regl((u, v, w), 2)
        ([3],
        [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w}]],
        [[{(0, 1): 1, (0, 2): -u*w^-1, (1, 2): w^-1, (2, 0): 1, (2, 2): v*w^-1}]])
    """
    (u, v, w) = variables
    data = {}
    data[2] = ([3], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w}]],
        [[{(0, 1): 1, (0, 2): -u/w, (1, 2): 1/w, (2, 0): 1, (2, 2):
        v/w}]])
    data[3] = ([24], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w, (3,
        5): -v, (3, 7): 1, (4, 6): -v, (4, 8): 1, (5, 3): 1, (5, 5):
        u, (6, 4): 1, (6, 6): u, (7, 5): w, (8, 6): w, (9, 9): u, (9,
        15): 1, (10, 10): u, (10, 17): 1, (11, 9): w, (12, 10): w,
        (13, 13): u, (13, 16): 1, (14, 13): w, (15, 9): -v, (15, 11):
        1, (16, 13): -v, (16, 14): 1, (17, 10): -v, (17, 12): 1, (18,
        19): -v, (18, 20): 1, (19, 18): 1, (19, 19): u, (20, 19): w,
        (21, 22): -v, (21, 23): 1, (22, 21): 1, (22, 22): u, (23, 22):
        w}, {(0, 3): -v, (0, 4): 1, (1, 15): -v, (1, 16): 1, (1, 22):
        -v, (1, 23): u*v/w, (2, 17): -v, (2, 18): 1, (2, 23): v*(u*v -
        w)/w, (3, 0): 1, (3, 3): u, (4, 3): w, (5, 9): -v, (5, 10): 1,
        (5, 12): -u/w, (5, 22): u, (5, 23): -u**2/w, (6, 11): -v, (6,
        22): w, (6, 23): -u, (7, 12): -u*v/w, (7, 13): -v, (7, 19): 1,
        (7, 21): -v, (7, 23): -u**2*v/w, (8, 12): -v, (8, 14): -v, (8,
        20): 1, (8, 23): -u*v, (9, 5): 1, (9, 9): u, (9, 23): u/w,
        (10, 9): w, (10, 11): -u, (11, 6): 1, (11, 11): u, (11, 13):
        u, (12, 13): w, (12, 23): -v, (13, 23): 1, (14, 8): 1, (14,
        14): u, (14, 23): v, (15, 1): 1, (15, 12): u/w, (15, 15): u,
        (16, 11): v, (16, 15): w, (16, 23): -u, (17, 2): 1, (17, 12):
        u*v/w, (17, 17): u, (18, 12): v, (18, 17): w, (19, 21): w,
        (19, 23): v, (20, 14): w, (21, 7): 1, (21, 21): u, (21, 23):
        u*v/w, (22, 11): 1, (23, 12): 1, (23, 23): u}]], [[{(0, 1): 1,
        (0, 2): -u/w, (1, 2): 1/w, (2, 0): 1, (2, 2): v/w, (3, 5): 1,
        (3, 7): -u/w, (4, 6): 1, (4, 8): -u/w, (5, 7): 1/w, (6, 8):
        1/w, (7, 3): 1, (7, 7): v/w, (8, 4): 1, (8, 8): v/w, (9, 11):
        1/w, (10, 12): 1/w, (11, 11): v/w, (11, 15): 1, (12, 12): v/w,
        (12, 17): 1, (13, 14): 1/w, (14, 14): v/w, (14, 16): 1, (15,
        9): 1, (15, 11): -u/w, (16, 13): 1, (16, 14): -u/w, (17, 10):
        1, (17, 12): -u/w, (18, 19): 1, (18, 20): -u/w, (19, 20): 1/w,
        (20, 18): 1, (20, 20): v/w, (21, 22): 1, (21, 23): -u/w, (22,
        23): 1/w, (23, 21): 1, (23, 23): v/w}, {(0, 3): 1, (0, 4):
        -u/w, (1, 15): 1, (1, 16): -u/w, (1, 22): u*v/w, (1, 23):
        -u/w, (2, 17): 1, (2, 18): -u/w, (3, 4): 1/w, (4, 0): 1, (4,
        4): v/w, (5, 9): 1, (5, 10): -u/w, (5, 13): -u/w, (5, 22):
        -u**2/w, (6, 11): 1, (6, 12): -u/w, (6, 13): -u*v/w, (6, 22):
        -u, (7, 19): -u/w, (7, 21): 1, (8, 13): -v, (8, 14): 1, (8,
        20): -u/w, (9, 10): 1/w, (9, 22): u/w, (10, 5): 1, (10, 6):
        -u/w, (10, 10): v/w, (10, 13): -u**2/w, (10, 23): u/w, (11,
        22): 1, (12, 13): -u, (12, 23): 1, (13, 12): 1/w, (13, 13):
        v/w, (14, 20): 1/w, (15, 13): u/w, (15, 16): 1/w, (15, 22):
        -v/w, (16, 1): 1, (16, 6): v/w, (16, 13): u*v/w, (16, 16):
        v/w, (17, 13): u*v/w, (17, 18): 1/w, (17, 23): -v/w, (18, 2):
        1, (18, 13): v, (18, 18): v/w, (18, 23): -v**2/w, (19, 7): 1,
        (19, 12): v/w, (19, 19): v/w, (19, 23): u*v/w, (20, 8): 1,
        (20, 20): v/w, (20, 23): v, (21, 13): -v/w, (21, 19): 1/w,
        (22, 6): 1/w, (22, 13): u/w, (22, 22): v/w, (23, 13): 1}]])
    return data[num_strands]


def read_regr(variables, num_strands=3):
    r"""
    Return precomputed data of Ivan Marin.

    This code was generated by :func:`create_demo_data`, please do not edit.

    INPUT:

    - ``variables`` -- tuple containing the indeterminates of the representation
    - ``num_strands`` -- integer; number of strands of the cubic Hecke algebra

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_regr
        sage: L.<u, v, w> = LaurentPolynomialRing(ZZ)
        sage: read_regr((u, v, w), 2)
        ([3],
        [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w}]],
        [[{(0, 1): 1, (0, 2): -u*w^-1, (1, 2): w^-1, (2, 0): 1, (2, 2): v*w^-1}]])
    """
    (u, v, w) = variables
    data = {}
    data[2] = ([3], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w}]],
        [[{(0, 1): 1, (0, 2): -u/w, (1, 2): 1/w, (2, 0): 1, (2, 2):
        v/w}]])
    data[3] = ([24], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w, (3,
        15): -v, (3, 17): 1, (4, 16): -v, (4, 18): 1, (4, 22): v**2,
        (4, 23): -v, (5, 9): -v, (5, 10): 1, (6, 13): -v, (6, 19): 1,
        (6, 21): -v, (6, 22): -u*v, (7, 11): -v, (7, 12): 1, (8, 14):
        -v, (8, 20): 1, (8, 22): -v*w, (9, 5): 1, (9, 9): u, (9, 23):
        u/w, (10, 9): w, (10, 21): -u, (10, 22): -u**2, (11, 7): 1,
        (11, 11): u, (11, 21): u, (11, 23): u*v/w, (12, 11): w, (12,
        22): -u*w, (13, 6): 1, (13, 13): u, (13, 22): v, (14, 8): 1,
        (14, 14): u, (14, 23): v, (15, 3): 1, (15, 15): u, (15, 22):
        u, (15, 23): -u**2/w, (16, 4): 1, (16, 16): u, (16, 21): v,
        (17, 15): w, (17, 22): u*v, (17, 23): -u, (18, 16): w, (19,
        13): w, (20, 14): w, (21, 22): -v, (21, 23): 1, (22, 21): 1,
        (22, 22): u, (23, 22): w}, {(0, 3): -v, (0, 4): 1, (1, 5): -v,
        (1, 6): 1, (2, 7): -v, (2, 8): 1, (3, 0): 1, (3, 3): u, (4,
        3): w, (5, 1): 1, (5, 5): u, (6, 5): w, (7, 2): 1, (7, 7): u,
        (8, 7): w, (9, 9): u, (9, 15): 1, (10, 13): u, (10, 16): 1,
        (10, 22): -v, (11, 9): w, (12, 13): w, (12, 23): -v, (13, 23):
        1, (14, 21): w, (14, 23): v, (15, 9): -v, (15, 11): 1, (16,
        22): w, (16, 23): -u, (17, 13): -v, (17, 14): 1, (17, 21): -v,
        (18, 19): -v, (18, 20): 1, (19, 18): 1, (19, 19): u, (20, 19):
        w, (21, 17): 1, (21, 21): u, (22, 10): 1, (22, 22): u, (23,
        12): 1, (23, 23): u}]], [[{(0, 1): 1, (0, 2): -u/w, (1, 2):
        1/w, (2, 0): 1, (2, 2): v/w, (3, 15): 1, (3, 17): -u/w, (3,
        23): u*(u*v - w)/w**2, (4, 16): 1, (4, 18): -u/w, (4, 22): -v,
        (4, 23): u*v/w, (5, 9): 1, (5, 10): -u/w, (5, 21): -u/w, (5,
        22): -u**2/w, (5, 23): -u*v/w**2, (6, 13): 1, (6, 19): -u/w,
        (6, 23): -v/w, (7, 11): 1, (7, 12): -u/w, (7, 21): -u*v/w, (7,
        22): -u, (7, 23): -u*v**2/w**2, (8, 14): 1, (8, 20): -u/w, (8,
        21): -v, (8, 23): -v**2/w, (9, 10): 1/w, (9, 22): u/w, (10,
        5): 1, (10, 10): v/w, (10, 22): u*v/w, (11, 12): 1/w, (11,
        23): u/w, (12, 7): 1, (12, 12): v/w, (12, 23): u*v/w, (13,
        19): 1/w, (14, 20): 1/w, (15, 17): 1/w, (15, 21): u/w, (16,
        18): 1/w, (17, 3): 1, (17, 17): v/w, (17, 21): u*v/w, (18, 4):
        1, (18, 18): v/w, (18, 21): v, (19, 6): 1, (19, 19): v/w, (19,
        22): v, (20, 8): 1, (20, 20): v/w, (20, 23): v, (21, 22): 1,
        (21, 23): -u/w, (22, 23): 1/w, (23, 21): 1, (23, 23): v/w},
        {(0, 3): 1, (0, 4): -u/w, (1, 5): 1, (1, 6): -u/w, (2, 7): 1,
        (2, 8): -u/w, (3, 4): 1/w, (4, 0): 1, (4, 4): v/w, (5, 6):
        1/w, (6, 1): 1, (6, 6): v/w, (7, 8): 1/w, (8, 2): 1, (8, 8):
        v/w, (9, 11): 1/w, (10, 13): -u**2/w, (10, 16): -u/w, (10,
        22): 1, (11, 11): v/w, (11, 15): 1, (12, 13): -u, (12, 23): 1,
        (13, 12): 1/w, (13, 13): v/w, (14, 12): v/w, (14, 14): v/w,
        (14, 17): 1, (15, 9): 1, (15, 11): -u/w, (16, 10): 1, (16,
        12): -u/w, (16, 16): v/w, (17, 13): u*v/w, (17, 14): -u/w,
        (17, 21): 1, (18, 19): 1, (18, 20): -u/w, (19, 20): 1/w, (20,
        18): 1, (20, 20): v/w, (21, 13): -v/w, (21, 14): 1/w, (22,
        13): u/w, (22, 16): 1/w, (23, 13): 1}]])
    return data[num_strands]


# generated data function for cubic Hecke algebra
def read_markov(bas_ele, variables, num_strands=4):
    r"""
    Return precomputed Markov trace coefficients.

    This code was generated by ``create_markov_trace_data.py`` (from
    the ``database_cubic_hecke`` repository), please do not edit.

    INPUT:

    - ``bas_ele`` -- an element of :class:`MarkovTraceModuleBasis`
    - ``variables`` -- tuple consisting of the variables used in
      the coefficients
    - ``num_strands`` -- integer (default: 4); the number of strands

    OUTPUT:

    A list of the coefficients. The i'th member corresponds to the i'th
    basis element.

    EXAMPLES::

        sage: # needs sympy
        sage: from sage.databases.cubic_hecke_db import read_markov
        sage: from sympy import var
        sage: u, v, w, s = var('u, v, w, s')
        sage: variables = (u, v, w, s)
        sage: read_markov('U2', variables, num_strands=3)
        [0, s, 1/s, s, 1/s, 0, 0, 0, 0, -s*v, s, s, -s*u/w, -v/s, 1/s,
         0, 0, 0, 0, 1/s, -u/(s*w), -v/s, 0, 0]
    """
    u, v, w, s = variables
    data = {}
    data[2] = {'U1': [0, s, 1/s], 'U2': [1, 0, 0]}
    data[3] = {'U1': [0, 0, 0, 0, 0, s**2, 1, 1, 1/s**2, u*s**2 + w, 0, 0, (s**2 +
        v)/w, (u*s**2 + w)/s**2, 0, s**2, 1, 1, 1/s**2, 0, (s**2 +
        v)/(w*s**2), (u*s**2 + w)/s**2, s**2, 0], 'U2': [0, s, 1/s, s,
        1/s, 0, 0, 0, 0, -v*s, s, s, (-u*s)/w, (-v)/s, 1/s, 0, 0, 0,
        0, 1/s, (-u)/(w*s), (-v)/s, 0, 0], 'U3': [1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'K4': [0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1]}

    return data[num_strands][bas_ele]
