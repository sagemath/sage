r"""
Library interface to Kenzo

Kenzo is a set of lisp functions to compute homology and
homotopy groups of topological spaces.

AUTHORS:

- Miguel Marco, Ana Romero (2019-01): Initial version


For this interface, Kenzo is loaded into ECL which is itself loaded
as a C library in Sage. Kenzo objects in this interface are nothing
but wrappers around ECL objects.
"""
# ****************************************************************************
#       Copyright (C) 2019 Miguel Marco <mmarco@unizar.es>
#                      and Ana Romero <ana.romero@unirioja.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.sage_object import SageObject
from sage.homology.homology_group import HomologyGroup
from sage.rings.integer_ring import ZZ
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup

from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex
from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet

from sage.libs.ecl import EclObject, ecl_eval, EclListIterator
from sage.features.kenzo import Kenzo

# defining the auxiliary functions as wrappers over the kenzo ones
kenzo_names = ['add',
               'array-dimensions',
               'basis_aux1',
               'basis_aux1',
               'bicomplex-spectral-sequence',
               'build-finite-ss2',
               'build-mrph-aux',
               'change-sorc-trgt-aux',
               'chcm-mat',
               'chcm-mat2',
               'classifying-space',
               'cmps',
               'convertmatrice',
               'crts-prdc',
               'degr-aux',
               'dffr-aux',
               'dffr_aux1',
               'dgop',
               'dgop-int-ext',
               'dstr-change-sorc-trgt-aux',
               'echcm',
               'eilenberg-moore-spectral-sequence',
               'evaluation-aux1',
               'gmsm',
               'homologie',
               'homotopy-list',
               'idnt-mrph',
               'join',
               'k-z',
               'k-z2',
               'k-zp',
               'kabstractsimplex_aux1',
               'kchaincomplex_aux1',
               'kmorphismchaincomplex_aux1',
               'loop-space',
               'make-array-from-lists',
               'make-array-to-lists',
               'moore',
               'ncol',
               'nlig',
               'nreverse',
               'nth',
               'opps',
               'orgn_aux1',
               'sbtr',
               'serre-spectral-sequence-product',
               'serre-whitehead-spectral-sequence',
               'sfinitesimplicialset_aux1',
               'smash-product',
               'sorc-aux',
               'spectral-sequence-differential-matrix',
               'spectral-sequence-group',
               'sphere',
               'suspension',
               'tnsr-prdc',
               'trgt-aux',
               'wedge',
               'zero-mrph']


# Now initialize Kenzo. For each string s in kenzo_names, the
# following defines __s__, a wrapper for a Kenzo function. For
# example __sphere__ is defined as EclObject("sphere"). Hyphens
# are replaced with underscores to get valid Python identifiers.
if Kenzo().is_present():
    from sage.env import KENZO_FAS
    if KENZO_FAS:
        ecl_eval("(require :kenzo \"{}\")".format(KENZO_FAS))
    else:
        ecl_eval("(require :kenzo)")

    ecl_eval("(in-package :cat)")
    ecl_eval("(setf *HOMOLOGY-VERBOSE* nil)")
    for s in kenzo_names:
        name = '__{}__'.format(s.replace('-', '_'))
        exec('{} = EclObject("{}")'.format(name, s))


def Sphere(n):
    r"""
    Return the ``n`` dimensional sphere as a Kenzo simplicial set.

    INPUT:

    - ``n`` -- the dimension of the sphere

    OUTPUT:

    - A :class:`KenzoSimplicialSet`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
        sage: s2 = Sphere(2)                           # optional - kenzo
        sage: s2                                       # optional - kenzo
        [K1 Simplicial-Set]
        sage: [s2.homology(i) for i in range(8)]       # optional - kenzo
        [Z, 0, Z, 0, 0, 0, 0, 0]
    """
    kenzosphere = __sphere__(n)
    return KenzoSimplicialSet(kenzosphere)


def MooreSpace(m, n):
    r"""
    Return the Moore space ``M(m, n)`` as a Kenzo simplicial set.

    The Moore space ``M(m, n)`` is the space whose n'th homology group
    is isomorphic to the cyclic group of order ``m``, and the rest of the
    homology groups are trivial.

    INPUT:

    - ``m`` -- A positive integer. The order of the nontrivial homology group.

    - ``n`` -- The dimension in which the homology is not trivial

    OUTPUT:

    - A KenzoSimplicialSet

    EXAMPLES::

        sage: from sage.interfaces.kenzo import MooreSpace   # optional - kenzo
        sage: m24 = MooreSpace(2,4)                          # optional - kenzo
        sage: m24                                            # optional - kenzo
        [K10 Simplicial-Set]
        sage: [m24.homology(i) for i in range(8)]            # optional - kenzo
        [Z, 0, 0, 0, C2, 0, 0, 0]
    """
    kenzomoore = __moore__(m, n)
    return KenzoSimplicialSet(kenzomoore)


def EilenbergMacLaneSpace(G, n):
    r"""
    Return the Eilenberg-MacLane space ``K(G, n)`` as a Kenzo simplicial group.

    The Eilenberg-MacLane space ``K(G, n)`` is the space whose has n'th homotopy
    group isomorphic to ``G``, and the rest of the homotopy groups are trivial.

    INPUT:

    - ``G`` -- group. Currently only ``ZZ`` and the additive group of two
      elements are supported.

    - ``n`` -- the dimension in which the homotopy is not trivial

    OUTPUT:

    - A :class:`KenzoSimplicialGroup`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
        sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
        sage: [e3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, Z, 0, C2, 0, C3]
        sage: f3 = EilenbergMacLaneSpace(AdditiveAbelianGroup([2]), 3)   # optional - kenzo
        sage: [f3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, C2, 0, C2, C2, C2]
    """
    if G == ZZ:
        kenzospace = __k_z__(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G == AdditiveAbelianGroup([2]):
        kenzospace = __k_z2__(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G in CommutativeAdditiveGroups() and G.is_cyclic():
        kenzospace = __k_zp__(G.cardinality(), n)
        return KenzoSimplicialGroup(kenzospace)
    else:
        raise NotImplementedError("Eilenberg-MacLane spaces are only supported over ZZ and ZZ_n")


class KenzoObject(SageObject):
    r"""
    Wrapper to Kenzo objects

    INPUT:

    - ``kenzo_object`` -- a wrapper around a Kenzo object
      (which is an ecl object).
    """

    def __init__(self, kenzo_object):
        r"""
        Construct the chain complex.

        TESTS::

            sage: from sage.interfaces.kenzo import KenzoObject  # optional -kenzo
            sage: from sage.interfaces.kenzo import __sphere__   # optional -kenzo
            sage: ks = __sphere__(2)                             # optional -kenzo
            sage: ks                                             # optional -kenzo
            <ECL: [K1 Simplicial-Set]>
            sage: s2 = KenzoObject(ks)                           # optional -kenzo
            sage: s2                                             # optional -kenzo
            [K1 Simplicial-Set]
            sage: TestSuite(s2).run(skip='_test_pickling')       # optional -kenzo

        """
        self._kenzo = kenzo_object

    def _repr_(self):
        r"""
        Represent the object.

        It just uses the ecl representation, removing the
        ecl decoration.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace    # optional - kenzo
            sage: m2 = MooreSpace(2,4)                            # optional - kenzo
            sage: m2._repr_()                                     # optional - kenzo
            '[K10 Simplicial-Set]'
        """
        kenzo_string = repr(self._kenzo)
        return kenzo_string[6:-1]


class KenzoSpectralSequence(KenzoObject):
    r"""
    Wrapper around Kenzo spectral sequences
    """

    def group(self, p, i, j):
        r"""
        Return the ``i,j``'th group of the ``p`` page.

        INPUT:

        - ``p`` -- the page to take the group from.

        - ``i`` -- the column where the group is taken from.

        - ``j`` -- the row where the group is taken from.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
            sage: S2 = Sphere(2)                           # optional - kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional - kenzo
            sage: EMS.group(0, -1, 2)                      # optional - kenzo
            Additive abelian group isomorphic to Z
            sage: EMS.group(0, -1, 3)                      # optional - kenzo
            Trivial group
        """
        invs = __spectral_sequence_group__(self._kenzo, p, i, j).python()
        if not invs:
            invs = []
        return AdditiveAbelianGroup(invs)

    def matrix(self, p, i, j):
        r"""
        Return the matrix that determines the differential from the
        ``i,j``'th group of the ``p``'th page.

        INPUT:

        - ``p`` -- the page.

        - ``i`` -- the column of the differential domain.

        - ``j`` -- the row of the differential domain.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
            sage: S3 = Sphere(3)                             # optional - kenzo
            sage: L = S3.loop_space()                        # optional - kenzo
            sage: EMS = L.em_spectral_sequence()             # optional - kenzo
            sage: EMS.table(1, -5, -2, 5, 8)                 # optional - kenzo
              0   Z   Z + Z + Z   Z + Z + Z
              0   0   0           0
              0   0   Z           Z + Z
              0   0   0           0
            sage: EMS.matrix(1, -2 ,8)                       # optional - kenzo
            [ 3 -2  0]
            [ 3  0 -3]
            [ 0  2 -3]
        """
        klist = __spectral_sequence_differential_matrix__(self._kenzo, p, i, j)
        plist = klist.python()
        if plist is None or plist == [None]:
            i = len(self.group(p, i, j).invariants())
            j = len(self.group(p, i - p, j + p - 1).invariants())
            return matrix(i, j)
        return matrix(plist)

    def differential(self, p, i, j):
        r"""
        Return the ``(p, i, j)`` differential morphism of the spectral sequence.

        INPUT:

        - ``p`` -- the page.

        - ``i`` -- the column of the differential domain.

        - ``j`` -- the row of the differential domain.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
            sage: S3 = Sphere(3)                             # optional - kenzo
            sage: L = S3.loop_space()                        # optional - kenzo
            sage: EMS = L.em_spectral_sequence()             # optional - kenzo
            sage: EMS.table(1,-5,-2,5,8)                     # optional - kenzo
              0   Z   Z + Z + Z   Z + Z + Z
              0   0   0           0
              0   0   Z           Z + Z
              0   0   0           0
            sage: EMS.matrix(1, -3, 8)                       # optional - kenzo
            [ 2 -2  2]
            sage: EMS.differential(1, -3, 8)                 # optional - kenzo
            Morphism from module over Integer Ring with invariants (0, 0, 0) to module with invariants (0,) that sends the generators to [(2), (-2), (2)]
        """
        domain = self.group(p, i, j)
        codomain = self.group(p, i - p, j + p - 1)
        M = self.matrix(p, i, j)
        images = [codomain(r) for r in M.columns()]
        return domain.hom(images, codomain=codomain)

    def table(self, p, i1, i2, j1, j2):
        r"""
        Return a table printing the groups in the ``p`` page.

        INPUT:

        - ``p`` -- the page to print.

        -- ``i1`` -- the first column to print.

        -- ``i2`` -- the last column to print.

        -- ``j1`` -- the first row to print.

        -- ``j2`` -- the last row to print.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
            sage: S2 = Sphere(2)                           # optional - kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional - kenzo
            sage: EMS.table(0, -2, 2, -2, 2)               # optional - kenzo
              0   Z   0   0   0
              0   0   0   0   0
              0   0   Z   0   0
              0   0   0   0   0
              0   0   0   0   0
        """
        from sage.misc.table import table
        groups = []
        for j in range(j2 - j1 + 1):
            row = []
            for i in range(i1, i2 + 1):
                group = self.group(p, i, j2 - j)
                if group.invariants():
                    row.append(group.short_name())
                else:
                    row.append('0')
            groups.append(row)
        return table(groups)


class KenzoChainComplex(KenzoObject):
    r"""
    Wrapper to Kenzo chain complexes. Kenzo simplicial sets are a particular case
    of Kenzo chain complexes.
    """
    def homology(self, n):
        r"""
        Return the ``n``'th homology group of the chain complex associated to this
        kenzo object.

        INPUT:

        - ``n`` -- the dimension in which compute the homology

        OUTPUT:

        - An homology group.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
            sage: s2 = Sphere(2)                             # optional - kenzo
            sage: s2                                         # optional - kenzo
            [K1 Simplicial-Set]
            sage: s2.homology(2)                             # optional - kenzo
            Z
        """
        echcm1 = __echcm__(self._kenzo)
        m1 = __chcm_mat__(echcm1, n)
        m2 = __chcm_mat__(echcm1, n + 1)
        homology = __homologie__(m1, m2)
        lhomomology = [i for i in EclListIterator(homology)]
        res = []
        for component in lhomomology:
            pair = [i for i in EclListIterator(component)]
            res.append(pair[0].python())
        return HomologyGroup(len(res), ZZ, res)

    def tensor_product(self, other):
        r"""
        Return the tensor product of ``self`` and ``other``.

        INPUT:

        - ``other`` --  The Kenzo object with which to compute the tensor product

        OUTPUT:

        - A :class:`KenzoChainComplex`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: p = s2.tensor_product(s3)                   # optional - kenzo
            sage: type(p)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoChainComplex'>
            sage: [p.homology(i) for i in range(8)]           # optional - kenzo
            [Z, 0, Z, Z, 0, Z, 0, 0]
        """
        return KenzoChainComplex(__tnsr_prdc__(self._kenzo, other._kenzo))

    def basis(self, dim):
        r"""
        Return the list of generators of the chain complex associated to the kenzo
        object ``self`` in dimension ``dim``.

        INPUT:

        - ``dim`` -- An integer number

        OUTPUT:

        - A list of the form ['G"dim"G0', 'G"dim"G1', 'G"dim"G2', ...].

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex   # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)         # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                                # optional - kenzo
            sage: kenzo_chcm                                                           # optional - kenzo
            [K... Chain-Complex]
            sage: for i in range(6):                                                   # optional - kenzo
            ....:     print("Basis in dimension %i: %s" % (i, kenzo_chcm.basis(i)))    # optional - kenzo
            Basis in dimension 0: ['G0G0', 'G0G1', 'G0G2']
            Basis in dimension 1: ['G1G0', 'G1G1']
            Basis in dimension 2: None
            Basis in dimension 3: ['G3G0', 'G3G1']
            Basis in dimension 4: ['G4G0', 'G4G1']
            Basis in dimension 5: ['G5G0', 'G5G1', 'G5G2']

        """
        return __basis_aux1__(self._kenzo, dim).python()

    def identity_morphism(self):
        r"""
        Return the identity morphism (degree 0) between ``self`` and itself.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                   # optional - kenzo
            sage: s2 = Sphere(2)                                             # optional - kenzo
            sage: tp = s2.tensor_product(s2)                                 # optional - kenzo
            sage: idnt = tp.identity_morphism()                              # optional - kenzo
            sage: type(idnt)                                                 # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoChainComplexMorphism'>
        """
        return KenzoChainComplexMorphism(__idnt_mrph__(self._kenzo))

    def null_morphism(self, target=None, degree=None):
        r"""
        Return the null morphism between the chain complexes ``self`` and ``target``
        of degree ``degree``.

        INPUT:

        - ``target`` -- A KenzoChainComplex or None (default).
        - ``degree`` -- An integer number or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism` representing the null morphism between
          ``self`` and ``target`` of degree ``degree``. If ``target`` takes None value,
          ``self`` is assumed as the target chain complex; if ``degree`` takes None value,
          0 is assumed as the degree of the null morphism.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                   # optional - kenzo
            sage: s2 = Sphere(2)                                             # optional - kenzo
            sage: s3 = Sphere(3)                                             # optional - kenzo
            sage: tp22 = s2.tensor_product(s2)                               # optional - kenzo
            sage: tp22                                                       # optional - kenzo
            [K... Chain-Complex]
            sage: tp23 = s2.tensor_product(s3)                               # optional - kenzo
            sage: tp23                                                       # optional - kenzo
            [K... Chain-Complex]
            sage: null1 = tp22.null_morphism()                               # optional - kenzo
            sage: null1                                                      # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null2 = tp22.null_morphism(target = tp23, degree = -3)     # optional - kenzo
            sage: null2                                                      # optional - kenzo
            [K... Morphism (degree -3): K... -> K...]
        """
        if target is None:
            target = self
        if degree is None:
            degree = 0
        if not isinstance(target, KenzoChainComplex):
            raise ValueError("'target' parameter must be a KenzoChainComplex instance")
        elif (not degree == 0) and (not degree.is_integer()):
            raise ValueError("'degree' parameter must be an Integer number")
        else:
            return KenzoChainComplexMorphism(__zero_mrph__(self._kenzo, target._kenzo, degree))

    def differential(self, dim=None, comb=None):
        r"""
        Return the differential of a combination.

        INPUT:

        - ``dim`` -- An integer number or None (default)

        - ``comb`` -- A list representing a formal sum of generators in the module
          of dimension ``dim`` or None (default). For example, to represent
          G7G12 + 3*G7G0 - 5*G7G3 we use the list [3, 'G7G0', -5, 'G7G3', 1, 'G7G12'].
          Note that the generators must be in ascending order respect to the number
          after the second G in their representation; the parameter
          ``comb`` = [1, 'G7G12', 3, 'G7G0', -5, 'G7G3'] will produce an error in
          Kenzo.

        OUTPUT:

        - If ``dim`` and ``comb`` are not None, it returns a Kenzo combination
          representing the differential of the formal combination represented by
          ``comb`` in the chain complex ``self`` in dimension ``dim``. On the other
          hand, if `dim`` or ``comb`` (or both) take None value, the differential
          :class:`KenzoMorphismChainComplex` of ``self`` is returned.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)
            sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
            sage: kenzo_chcm                                                      # optional - kenzo
            [K... Chain-Complex]
            sage: kenzo_chcm.basis(4)                                             # optional - kenzo
            ['G4G0', 'G4G1']
            sage: kenzo_chcm.differential(4, [1, 'G4G0'])                         # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 3}
            <1 * G3G0>
            <3 * G3G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: kenzo_chcm.basis(5)                                             # optional - kenzo
            ['G5G0', 'G5G1', 'G5G2']
            sage: kenzo_chcm.differential(5, [1, 'G5G0', 2, 'G5G2'])              # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <6 * G4G0>
            <-3 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
        """
        if dim is not None and comb is not None:
            cmbn_list = pairing(comb)
            return KenzoObject(__dffr_aux1__(self._kenzo, dim, cmbn_list))
        else:
            return KenzoChainComplexMorphism(__dffr_aux__(self._kenzo))

    def orgn(self):
        r"""
        Return the :orgn slot of Kenzo, which stores as a list the origin of the object

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                # optional - kenzo
            sage: s2 = Sphere(2)                                          # optional - kenzo
            sage: l2 = s2.loop_space()                                    # optional - kenzo
            sage: l2.orgn()                                               # optional - kenzo
            '(LOOP-SPACE [K... Simplicial-Set])'
            sage: A = l2.cartesian_product(s2)                            # optional - kenzo
            sage: A.orgn()                                                # optional - kenzo
            '(CRTS-PRDC [K... Simplicial-Group] [K... Simplicial-Set])'
        """
        return str(__orgn_aux1__(self._kenzo))


class KenzoSimplicialSet(KenzoChainComplex):
    r"""
    Wrapper to Kenzo simplicial sets.

    In Kenzo, the homology of a simplicial set in computed from its associated
    chain complex. Hence, this class inherits from `KenzoChainComplex`.
    """

    def loop_space(self, n=1):
        r"""
        Return the ``n`` th iterated loop space.

        INPUT:

        - ``n`` -- (default: 1) the number of times to iterate the loop space
          construction

        OUTPUT:

        - A :class:`KenzoSimplicialGroup`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: type(l2)                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: [l2.homology(i) for i in range(8)]          # optional - kenzo
            [Z, Z, Z, Z, Z, Z, Z, Z]
        """
        return KenzoSimplicialGroup(__loop_space__(self._kenzo, n))

    def cartesian_product(self, other):
        r"""
        Return the cartesian product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the product is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: p = s2.cartesian_product(s3)                # optional - kenzo
            sage: type(p)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [p.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, Z, Z, 0, Z]
        """
        prod_kenzo = __crts_prdc__(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(prod_kenzo)

    def suspension(self):
        r"""
        Return the suspension of the simplicial set.

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
            sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
            sage: s = e3.suspension()                                        # optional - kenzo
            sage: type(s)                                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [s.homology(i) for i in range(6)]                          # optional - kenzo
            [Z, 0, 0, 0, Z, 0]
        """
        return KenzoSimplicialSet(__suspension__(self._kenzo))

    def homotopy_group(self, n):
        """
        Return the n'th homotopy group of ``self``

        INPUT:

        - ``n`` -- the dimension of the homotopy group to be computed

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere      # optional - kenzo
            sage: s2 = Sphere(2)                                # optional - kenzo
            sage: p = s2.cartesian_product(s2)                  # optional - kenzo
            sage: p.homotopy_group(3)                           # optional - kenzo
            Multiplicative Abelian group isomorphic to Z x Z


        .. WARNING::

            This method assumes that the underlying space is simply connected.
            You might get wrong answers if it is not.
        """
        if n not in ZZ or n < 2:
            raise ValueError("""homotopy groups can only be computed
                for dimensions greater than 1""")
        lgens = __homotopy_list__(self._kenzo, n).python()
        if lgens is not None:
            trgens = [0 if i == 1 else i for i in sorted(lgens)]
            return AbelianGroup(trgens)
        else:
            return AbelianGroup([])

    def em_spectral_sequence(self):
        r"""
        Return the Eilenberg-Moore spectral sequence of ``self``.

        OUTPUT:

        - A :class:`KenzoSpectralSequence`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
            sage: S2 = Sphere(2)                           # optional - kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional - kenzo
            sage: EMS.table(0, -2, 2, -2, 2)               # optional - kenzo
              0   Z   0   0   0
              0   0   0   0   0
              0   0   Z   0   0
              0   0   0   0   0
              0   0   0   0   0


        .. WARNING::

            This method assumes that the underlying space is simply connected.
            You might get wrong answers if it is not.
        """
        if self.homology(1).invariants():
            raise ValueError("""Eilenberg-Moore spectral sequence implemented
                only for 1-reduced simplicial sets""")
        return KenzoSpectralSequence(__eilenberg_moore_spectral_sequence__(self._kenzo))

    def sw_spectral_sequence(self):
        r"""
        Return the Serre sequence of the first step of the Whitehead tower.

        OUTPUT:

        - A :class:`KenzoSpectralSequence`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere  # optional - kenzo
            sage: S3 = Sphere(3)                            # optional - kenzo
            sage: E = S3.sw_spectral_sequence()             # optional - kenzo
            sage: T = E.table(0, 0, 4, 0, 4)                # optional - kenzo
            sage: T                                         # optional - kenzo
              Z   0   0   Z   0
              0   0   0   0   0
              Z   0   0   Z   0
              0   0   0   0   0
              Z   0   0   Z   0
        """
        if self.homology(1).invariants():
            raise ValueError("""Eilenberg-Moore spectral sequence implemented
                only for 1-reduced simplicial sets""")
        return KenzoSpectralSequence(__serre_whitehead_spectral_sequence__(self._kenzo))

    def serre_spectral_sequence(self):
        r"""
        Return the spectral sequence of ``self``.

        The object self must be created as a cartesian product (twisted or not).

        OUTPUT:

        - A :class:`KenzoSpectralSequence`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere
            sage: S2 = Sphere(2)                            # optional - kenzo
            sage: S3 = Sphere(3)                            # optional - kenzo
            sage: P = S2.cartesian_product(S3)              # optional - kenzo
            sage: E = P.serre_spectral_sequence()           # optional - kenzo
            sage: E.table(0, 0, 2, 0, 3)                    # optional - kenzo
              Z   0   Z
              0   0   0
              0   0   0
              Z   0   Z

        .. WARNING::

            This method assumes that the underlying space is simply connected.
            You might get wrong answers if it is not.
        """
        if self.homology(1).invariants():
            raise ValueError("""Eilenberg-Moore spectral sequence implemented
                only for 1-reduced simplicial sets""")
        return KenzoSpectralSequence(__serre_spectral_sequence_product__(self._kenzo))

    def wedge(self, other):
        r"""
        Return the wedge of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the wedge is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: w = s2.wedge(s3)                            # optional - kenzo
            sage: type(w)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [w.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, Z, Z, 0, 0]
        """
        wedge_kenzo = __wedge__(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(wedge_kenzo)

    def join(self, other):
        r"""
        Return the join of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the join is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: j = s2.join(s3)                             # optional - kenzo
            sage: type(j)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [j.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, 0, 0, 0, 0]
        """
        join_kenzo = __join__(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(join_kenzo)

    def smash_product(self, other):
        r"""
        Return the smash product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the smash product is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: s = s2.smash_product(s3)                    # optional - kenzo
            sage: type(s)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [s.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, 0, 0, 0, Z]
        """
        smash_kenzo = __smash_product__(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(smash_kenzo)


class KenzoSimplicialGroup(KenzoSimplicialSet):
    r"""
    Wrapper around Kenzo simplicial groups.
    """

    def classifying_space(self):
        r"""
        Return the classifying space.

        OUTPUT:

        - A :class:`KenzoSimplicialGroup`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace      # optional - kenzo
            sage: m2 = MooreSpace(2,4)                              # optional - kenzo
            sage: l2 = m2.loop_space()                              # optional - kenzo
            sage: c = l2.classifying_space()                        # optional - kenzo
            sage: type(c)                                           # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: [c.homology(i) for i in range(8)]                 # optional - kenzo
            [Z, 0, 0, 0, C2, 0, 0, 0]
        """
        return KenzoSimplicialGroup(__classifying_space__(self._kenzo))


def k2s_matrix(kmatrix):
    r"""
    Convert an array of ECL to a matrix of Sage.

    INPUT:

    - ``kmatrix`` -- An array in ECL

    EXAMPLES::

        sage: from sage.interfaces.kenzo import k2s_matrix         # optional - kenzo
        sage: from sage.libs.ecl import EclObject
        sage: M = EclObject("#2A((1 2 3) (3 2 1) (1 1 1))")
        sage: k2s_matrix(M)                                        # optional - kenzo
        [1 2 3]
        [3 2 1]
        [1 1 1]
    """
    dimensions = __array_dimensions__(kmatrix).python()
    kmatrix_list = __make_array_to_lists__(kmatrix).python()
    return matrix(dimensions[0], dimensions[1], kmatrix_list)


def s2k_matrix(smatrix):
    r"""
    Convert a matrix of Sage to an array of ECL.

    INPUT:

    - ``smatrix`` -- A matrix in Sage

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import s2k_matrix      # optional - kenzo
        sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
        sage: s2k_matrix(A)                                     # optional - kenzo
        <ECL: #2A((1 2 3) (3 2 1) (1 1 1))>
    """
    initcontents = []
    dimensions = smatrix.dimensions()
    for i in smatrix.rows():
        initcontents.append(i.list())
    return __make_array_from_lists__(dimensions[0], dimensions[1], initcontents)


def s2k_dictmat(sdictmat):
    r"""
    Convert a dictionary in Sage, whose values are matrices, to an assoc list
    in ECL.

    INPUT:

    - ``sdictmat`` -- A dictionary in Sage

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import s2k_dictmat   # optional - kenzo
        sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
        sage: B = Matrix([[1,2],[2,1],[1,1]])
        sage: d = {1 : A, 2 : B}
        sage: s2k_dictmat(d)                                  # optional - kenzo
        <ECL: ((2 . #2A((1 2) (2 1) (1 1))) (1 . #2A((1 2 3) (3 2 1) (1 1 1))))>

    """
    rslt = EclObject([])
    for k in sdictmat.keys():
        rslt = EclObject(k).cons(s2k_matrix(sdictmat[k])).cons(rslt)
    return rslt


def pairing(slist):
    r"""
    Convert a list of Sage (which has an even length) to an assoc list in ECL.

    INPUT:

    - ``slist`` -- A list in Sage

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import pairing   # optional - kenzo
        sage: l = [1,2,3]
        sage: pairing(l)                                  # optional - kenzo
        <ECL: ((2 . 3))>
    """
    rslt = EclObject([])
    for k in range(len(slist) - 1, 0, -2):
        rslt = EclObject(slist[k - 1]).cons(EclObject(slist[k])).cons(rslt)
    return rslt


def KChainComplex(chain_complex):
    r"""
    Construct a KenzoChainComplex from a ChainComplex of degree = -1 in
    Sage.

    INPUT:

    - ``chain_complex`` -- A ChainComplex of degree = -1

    OUTPUT:

    - A KenzoChainComplex

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
        sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
        sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
        sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
        sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)    # optional - kenzo
        sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
        sage: kenzo_chcm                                                      # optional - kenzo
        [K... Chain-Complex]
        sage: kenzo_chcm.homology(5)                                          # optional - kenzo
        Z x Z
    """
    d = chain_complex.differential()
    chcm = s2k_dictmat(d)
    str_orgn = str(d)[1:-1].replace(":", " ").replace(" ", ".").replace("\n", "").replace(",", "")
    return KenzoChainComplex(__kchaincomplex_aux1__(chcm, str_orgn))


def SChainComplex(kchaincomplex, start=0, end=15):
    r"""
    Convert the KenzoChainComplex ``kchcm`` (between dimensions ``start`` and
    ``end``) to a ChainComplex.

    INPUT:

    - ``kchaincomplex`` -- A KenzoChainComplex

    - ``start`` -- An integer number (optional, default 0)

    - ``end`` -- An integer number greater than or equal to ``start`` (optional, default 15)

    OUTPUT:

    - A ChainComplex

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplex, SChainComplex   # optional - kenzo
        sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
        sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
        sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
        sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)     # optional - kenzo
        sage: SChainComplex(KChainComplex(sage_chcm)) == sage_chcm             # optional - kenzo
        True

    ::

        sage: from sage.interfaces.kenzo import SChainComplex, Sphere     # optional - kenzo
        sage: S4 = Sphere(4)                       # optional - kenzo
        sage: C = SChainComplex(S4)                # optional - kenzo
        sage: C                                    # optional - kenzo
        Chain complex with at most 3 nonzero terms over Integer Ring
        sage: C._ascii_art_()                      # optional - kenzo
        0 <-- C_4 <-- 0  ...  0 <-- C_0 <-- 0
        sage: [C.homology(i) for i in range(6)]    # optional - kenzo
        [Z, 0, 0, 0, Z, 0]
    """
    matrices = {}
    for i in range(start, end):
        dffr_i = __chcm_mat2__(kchaincomplex._kenzo, i)
        nlig = __nlig__(dffr_i).python()
        ncol = __ncol__(dffr_i).python()
        if ((nlig != 0) and (ncol != 0)):
            matrices[i] = k2s_matrix(__convertmatrice__(dffr_i))
        else:
            matrices[i] = matrix(nlig, ncol)
    return ChainComplex(matrices, degree=-1)


def SAbstractSimplex(simplex, dim):
    r"""
    Convert an abstract simplex of Kenzo to an AbstractSimplex.

    INPUT:

    - ``simplex`` -- An abstract simplex of Kenzo.

    - ``dim``-- The dimension of ``simplex``.

    OUTPUT:

    - An AbstractSimplex.

    EXAMPLES::

        sage: from sage.libs.ecl import EclObject, ecl_eval
        sage: from sage.interfaces.kenzo import KenzoObject,\
        ....: SAbstractSimplex                                  # optional - kenzo
        sage: KAbSm = KenzoObject(ecl_eval("(ABSM 15 'K)"))     # optional - kenzo
        sage: SAbSm1 = SAbstractSimplex(KAbSm, 2)               # optional - kenzo
        sage: SAbSm2 = SAbstractSimplex(KAbSm, 7)               # optional - kenzo
        sage: SAbSm1.degeneracies()                             # optional - kenzo
        [3, 2, 1, 0]
        sage: SAbSm1.dimension()                                # optional - kenzo
        6
        sage: SAbSm2.dimension()                                # optional - kenzo
        11
    """
    degeneracies = __dgop_int_ext__(__dgop__(simplex._kenzo)).python()
    if degeneracies is None:
        degeneracies = []
    else:
        degeneracies = tuple(degeneracies)
    name = __gmsm__(simplex._kenzo).python()
    return AbstractSimplex(dim, degeneracies, name=name)


def KAbstractSimplex(simplex):
    r"""
    Convert an AbstractSimplex in Sage to an abstract simplex of Kenzo.

    INPUT:

    - ``simplex`` -- An AbstractSimplex.

    OUTPUT:

    - An abstract simplex of Kenzo.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import AbstractSimplex
        sage: from sage.interfaces.kenzo import KAbstractSimplex,\
        ....: SAbstractSimplex                                          # optional - kenzo
        sage: SAbSm = AbstractSimplex(1, (2,0,3,2,1), name = 'SAbSm')   # optional - kenzo
        sage: KAbSm = KAbstractSimplex(SAbSm)                           # optional - kenzo
        sage: SAbSm2 = SAbstractSimplex(KAbSm, 1)                       # optional - kenzo
        sage: SAbSm.degeneracies() == SAbSm2.degeneracies()             # optional - kenzo
        True
        sage: SAbSm.dimension() == SAbSm2.dimension()                   # optional - kenzo
        True
    """
    return KenzoObject(__kabstractsimplex_aux1__(simplex.degeneracies(),
                                             's' + str(hash(simplex))))


def KFiniteSimplicialSet(sset):
    r"""
    Convert a finite SimplicialSet in Sage to a finite simplicial set of Kenzo.

    INPUT:

    - ``sset`` -- A finite SimplicialSet.

    OUTPUT:

    - A finite simplicial set of Kenzo.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: from sage.interfaces.kenzo import KFiniteSimplicialSet    # optional - kenzo
        sage: s0 = AbstractSimplex(0, name='s0')
        sage: s1 = AbstractSimplex(0, name='s1')
        sage: s2 = AbstractSimplex(0, name='s2')
        sage: s01 = AbstractSimplex(1, name='s01')
        sage: s02 = AbstractSimplex(1, name='s02')
        sage: s12 = AbstractSimplex(1, name='s12')
        sage: s012 = AbstractSimplex(2, name='s012')
        sage: Triangle = SimplicialSet({s01: (s1, s0),\
        ....: s02: (s2, s0), s12: (s2, s1)}, base_point = s0)
        sage: KTriangle = KFiniteSimplicialSet(Triangle)                # optional - kenzo
        sage: KTriangle.homology(1)                                     # optional - kenzo
        Z
        sage: KTriangle.basis(1)                                        # optional - kenzo
        ['CELL_1_0', 'CELL_1_1', 'CELL_1_2']
        sage: S1 = simplicial_sets.Sphere(1)
        sage: S3 = simplicial_sets.Sphere(3)
        sage: KS1vS3 = KFiniteSimplicialSet(S1.wedge(S3))               # optional - kenzo
        sage: KS1vS3.homology(3)                                        # optional - kenzo
        Z
    """
    from sage.topology.simplicial_set_constructions import ProductOfSimplicialSets
    if isinstance(sset, ProductOfSimplicialSets):
        f0 = KFiniteSimplicialSet(sset.factor(0))
        for f1 in sset.factors()[1:]:
            f0 = f0.cartesian_product(KFiniteSimplicialSet(f1))
        return f0
    else:
        allcells = sset.cells()
        namecells = {c: 'cell_{}_{}'.format(d, allcells[d].index(c))
                     for d in allcells for c in allcells[d]}
        dim = sset.dimension()
        list_rslt = [namecells[i] for i in sset.n_cells(0)]
        if (dim > 0):
            for k in range(1, dim + 1):
                k_cells = sset.n_cells(k)
                if k_cells:
                    list_rslt.append(k)
                    for x in k_cells:
                        list_rslt.append(namecells[x])
                        auxiliar_list = []
                        for z in sset.faces(x):
                            degen_z = z.degeneracies()
                            name = namecells[z.nondegenerate()]
                            degen_z.append(name)
                            auxiliar_list.append(degen_z)
                        list_rslt.append(auxiliar_list)
        return KenzoSimplicialSet(__build_finite_ss2__(list_rslt))


def SFiniteSimplicialSet(ksimpset, limit):
    r"""
    Convert the ``limit``-skeleton of a finite simplicial set in Kenzo to a
    finite SimplicialSet in Sage.

    INPUT:

    - ``ksimpset`` -- A finite simplicial set in Kenzo.

    - ``limit`` -- A natural number.

    OUTPUT:

    - A finite SimplicialSet.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import SimplicialSet
        sage: from sage.interfaces.kenzo import AbstractSimplex,\
        ....:  KFiniteSimplicialSet, SFiniteSimplicialSet, Sphere   # optional - kenzo
        sage: s0 = AbstractSimplex(0, name='s0')                    # optional - kenzo
        sage: s1 = AbstractSimplex(0, name='s1')                    # optional - kenzo
        sage: s2 = AbstractSimplex(0, name='s2')                    # optional - kenzo
        sage: s01 = AbstractSimplex(1, name='s01')                  # optional - kenzo
        sage: s02 = AbstractSimplex(1, name='s02')                  # optional - kenzo
        sage: s12 = AbstractSimplex(1, name='s12')                  # optional - kenzo
        sage: s012 = AbstractSimplex(2, name='s012')                # optional - kenzo
        sage: Triangle = SimplicialSet({s01: (s1, s0),\
        ....: s02: (s2, s0), s12: (s2, s1)}, base_point = s0)       # optional - kenzo
        sage: KTriangle = KFiniteSimplicialSet(Triangle)            # optional - kenzo
        sage: STriangle = SFiniteSimplicialSet(KTriangle, 1)        # optional - kenzo
        sage: STriangle.homology()                                  # optional - kenzo
        {0: 0, 1: Z}
        sage: S1 = simplicial_sets.Sphere(1)                        # optional - kenzo
        sage: S3 = simplicial_sets.Sphere(3)                        # optional - kenzo
        sage: KS1vS3 = KFiniteSimplicialSet(S1.wedge(S3))           # optional - kenzo
        sage: SS1vS3 = SFiniteSimplicialSet(KS1vS3, 3)              # optional - kenzo
        sage: SS1vS3.homology()                                     # optional - kenzo
        {0: 0, 1: Z, 2: 0, 3: Z}
    """
    list_orgn = __orgn_aux1__(ksimpset._kenzo).python()
    if __nth__(0, list_orgn).python()[0] == 'CRTS-PRDC':
        return SFiniteSimplicialSet(
            KenzoSimplicialSet(__nth__(1, list_orgn)), limit).cartesian_product(
                SFiniteSimplicialSet(KenzoSimplicialSet(__nth__(2, list_orgn)), limit))
    rslt = {}
    simplices = []
    faces = []
    bases = []
    names = []
    for k in range(limit + 1):
        basis_k = __basis_aux1__(ksimpset._kenzo, k)
        names_k = ksimpset.basis(k)
        lbasis_k = [AbstractSimplex(k, name=i) for i in EclListIterator(basis_k)]
        bases.append(lbasis_k)
        names.append(names_k)
    all_simplices = __sfinitesimplicialset_aux1__(ksimpset._kenzo, limit)
    lall_simplices = [i for i in EclListIterator(all_simplices)]
    dim = 1
    for Kdim in lall_simplices:
        for simp in Kdim:
            index1 = names[dim].index(str(simp.car()))
            lKdim_cdr = []
            for i in EclListIterator(simp.cdr()):
                degenop = __dgop_int_ext__(__dgop__(i)).python()
                if degenop is None:
                    degenop = []
                index2 = names[dim - len(degenop) - 1].index(str(__gmsm__(i)))
                lKdim_cdr.append(bases[dim - len(degenop) - 1][index2].apply_degeneracies(*degenop))
            simplices.append(bases[dim][index1])
            faces.append(tuple(lKdim_cdr))
        dim += 1
    for i in range(len(simplices)):
        rslt[simplices[i]] = faces[i]
    return SimplicialSet(rslt)


class KenzoChainComplexMorphism(KenzoObject):
    r"""
    Wrapper to Kenzo morphisms between chain complexes.
    """

    def source_complex(self):
        r"""
        Return the source chain complex of the morphism.

        OUTPUT:

        - A :class:`KenzoChainComplex`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()              # optional - kenzo
            sage: differential_morphism                                          # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: differential_morphism.source_complex()                         # optional - kenzo
            [K... Chain-Complex]
        """
        return KenzoChainComplex(__sorc_aux__(self._kenzo))

    def target_complex(self):
        r"""
        Return the target chain complex of the morphism.

        OUTPUT:

        - A :class:`KenzoChainComplex`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()              # optional - kenzo
            sage: differential_morphism                                          # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: differential_morphism.target_complex()                         # optional - kenzo
            [K... Chain-Complex]
        """
        return KenzoChainComplex(__trgt_aux__(self._kenzo))

    def degree(self):
        r"""
        Return the degree of the morphism.

        OUTPUT:

        - An integer number, the degree of the morphism.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()              # optional - kenzo
            sage: differential_morphism                                          # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: differential_morphism.degree()                                 # optional - kenzo
            -1
            sage: differential_morphism.composite(differential_morphism).degree() # optional - kenzo
            -2
            sage: kenzo_chcm.null_morphism().degree()                            # optional - kenzo
            0
        """
        return __degr_aux__(self._kenzo).python()

    def evaluation(self, dim, comb):
        r"""
        Apply the morphism on a combination ``comb`` of dimension ``dim``.

        INPUT:

        - ``dim`` -- An integer number

        - ``comb`` -- A list representing a formal sum of generators in the module
          of dimension ``dim``. For example, to represent G7G12 + 3*G7G0 - 5*G7G3
          we use the list [3, 'G7G0', -5, 'G7G3', 1, 'G7G12']. Note that the
          generators must be in ascending order respect to the number after the
          second G in their representation; the parameter
          ``comb`` = [1, 'G7G12', 3, 'G7G0', -5, 'G7G3'] will produce an error in
          Kenzo.

        OUTPUT:

        - A Kenzo combination representing the result of applying the morphism on the formal
          combination represented by ``comb`` in the chain complex ``self`` in
          dimension ``dim``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)
            sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
            sage: kenzo_chcm                                                      # optional - kenzo
            [K... Chain-Complex]
            sage: differential_morphism = kenzo_chcm.differential()               # optional - kenzo
            sage: differential_morphism                                           # optional - kenzo
            [K... Morphism (degree -1): K... -> K...]
            sage: dif_squared = differential_morphism.composite(differential_morphism)  # optional - kenzo
            sage: dif_squared                                                     # optional - kenzo
            [K... Morphism (degree -2): K... -> K...]
            sage: kenzo_chcm.basis(5)                                             # optional - kenzo
            ['G5G0', 'G5G1', 'G5G2']
            sage: kenzo_chcm.differential(5, [1, 'G5G0', 2, 'G5G2'])              # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <6 * G4G0>
            <-3 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: differential_morphism.evaluation(5, [1, 'G5G0', 2, 'G5G2'])     # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <6 * G4G0>
            <-3 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: dif_squared.evaluation(5, [1, 'G5G0', 2, 'G5G2'])               # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 3}
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: idnt = kenzo_chcm.identity_morphism()                             # optional - kenzo
            sage: idx2 = idnt.sum(idnt)                                             # optional - kenzo
            sage: idnt.evaluation(5, [1, 'G5G0', 2, 'G5G2'])                        # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 5}
            <1 * G5G0>
            <2 * G5G2>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: idx2.evaluation(5, [1, 'G5G0', 2, 'G5G2'])                      # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 5}
            <2 * G5G0>
            <4 * G5G2>
            ------------------------------------------------------------------------------
            <BLANKLINE>
        """
        if dim.is_integer():
            if isinstance(comb, list):
                cmbn_list = pairing(comb)
                return KenzoObject(__evaluation_aux1__(self._kenzo, dim, cmbn_list))
            else:
                raise ValueError("'comb' parameter must be a list")
        else:
            raise ValueError("'dim' parameter must be an integer number")

    def opposite(self):
        r"""
        Return the opposite morphism of ``self``, i.e., -1 x ``self``.

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: idnt = kenzo_chcm.identity_morphism()                          # optional - kenzo
            sage: idnt                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: opps_id = idnt.opposite()                                      # optional - kenzo
            sage: opps_id                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: idnt.evaluation(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <2 * G4G0>
            <-5 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: opps_id.evaluation(4, [2, 'G4G0', -5, 'G4G1'])                 # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <-2 * G4G0>
            <5 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
        """
        return KenzoChainComplexMorphism(__opps__(self._kenzo))

    def composite(self, object=None):
        r"""
        Return the composite of ``self`` and the morphism(s) given by the parameter ``object``.

        INPUT:

        - ``object`` -- A KenzoChainComplexMorphism instance, a KenzoChainComplex instance, a tuple
          of KenzoChainComplexMorphism and KenzoChainComplex instances, or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`: if ``object`` is a KenzoChainComplexMorphism, the
          composite of ``self`` and ``object`` is returned; if ``object`` is a KenzoChainComplex,
          the composite of ``self`` and the differential morphism of ``object`` is returned; if
          ``object`` is a tuple, the composite of ``self`` and the morphisms or the differential
          morphisms of the given chain complexes in ``object`` is returned (if ``object`` is
          None, ``self`` morphism is returned).

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere                   # optional - kenzo
            sage: s2 = Sphere(2)                                             # optional - kenzo
            sage: s3 = Sphere(3)                                             # optional - kenzo
            sage: tp22 = s2.tensor_product(s2)                               # optional - kenzo
            sage: tp23 = s2.tensor_product(s3)                               # optional - kenzo
            sage: idnt = tp22.identity_morphism()                            # optional - kenzo
            sage: idnt                                                       # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null = tp23.null_morphism(target = tp22, degree = 4)       # optional - kenzo
            sage: null                                                       # optional - kenzo
            [K... Morphism (degree 4): K... -> K...]
            sage: idnt.composite((tp22, null))                               # optional - kenzo
            [K... Morphism (degree 3): K... -> K...]
        """
        if object is None:
            return self
        if isinstance(object, KenzoChainComplexMorphism):
            return KenzoChainComplexMorphism(__cmps__(self._kenzo, object._kenzo))
        elif isinstance(object, KenzoChainComplex):
            return KenzoChainComplexMorphism(__cmps__(self._kenzo, __dffr_aux__(object._kenzo)))
        elif isinstance(object, tuple):
            rslt = self._kenzo
            for mrph in object:
                rslt = __cmps__(rslt, mrph._kenzo)
            return KenzoChainComplexMorphism(rslt)

    def sum(self, object=None):
        r"""
        Return a morphism, sum of the morphism ``self`` and the morphism(s) given
        by the parameter ``object``.

        INPUT:

        - ``object`` -- A KenzoChainComplexMorphism instance, a tuple of KenzoChainComplexMorphism
          instances or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`, sum of the morphism ``self`` and the morphism(s)
          given by ``object`` (if ``object`` is None, ``self`` morphism is returned).

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: idnt = kenzo_chcm.identity_morphism()                          # optional - kenzo
            sage: idnt                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: opps_id = idnt.opposite()                                      # optional - kenzo
            sage: opps_id                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null = kenzo_chcm.null_morphism()                              # optional - kenzo
            sage: null                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: idx2 = idnt.sum(idnt)                                          # optional - kenzo
            sage: idx5 = idx2.sum(\
            ....: (opps_id, idnt, idnt, null, idx2.sum(idnt), opps_id))          # optional - kenzo
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: idx2.evaluation(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <4 * G4G0>
            <-10 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: idx5.evaluation(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <10 * G4G0>
            <-25 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
        """
        if object is None:
            return self
        if isinstance(object, KenzoChainComplexMorphism):
            return KenzoChainComplexMorphism(__add__(self._kenzo, object._kenzo))
        elif isinstance(object, tuple):
            rslt = self._kenzo
            for mrph in object:
                rslt = __add__(rslt, mrph._kenzo)
            return KenzoChainComplexMorphism(rslt)

    def substract(self, object=None):
        r"""
        Return a morphism, difference of the morphism ``self`` and the morphism(s) given by the
        parameter ``object``.

        INPUT:

        - ``object`` -- A KenzoChainComplexMorphism instance, a tuple of KenzoChainComplexMorphism
          instances or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`, difference of the morphism ``self`` and the
          morphism(s) given by ``object`` (if ``object`` is None, ``self`` morphism is returned).
          For example, if ``object`` = (mrph1, mrph2, mrph3) the result is
          ``self`` - mrph1 - mrph2 - mrph3.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)   # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                          # optional - kenzo
            sage: kenzo_chcm                                                     # optional - kenzo
            [K... Chain-Complex]
            sage: idnt = kenzo_chcm.identity_morphism()                          # optional - kenzo
            sage: idnt                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: opps_id = idnt.opposite()                                      # optional - kenzo
            sage: opps_id                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null = kenzo_chcm.null_morphism()                              # optional - kenzo
            sage: null                                                           # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: idx2 = idnt.substract(opps_id)                                 # optional - kenzo
            sage: opps_idx2 = idx2.substract\
            ....: ((opps_id, idnt, idnt, null, idx2.substract(opps_id)))         # optional - kenzo
            sage: kenzo_chcm.basis(4)                                            # optional - kenzo
            ['G4G0', 'G4G1']
            sage: idx2.evaluation(4, [2, 'G4G0', -5, 'G4G1'])                    # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <4 * G4G0>
            <-10 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: opps_idx2.evaluation(4, [2, 'G4G0', -5, 'G4G1'])               # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <-4 * G4G0>
            <10 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
        """
        if object is None:
            return self
        if isinstance(object, KenzoChainComplexMorphism):
            return KenzoChainComplexMorphism(__sbtr__(self._kenzo, object._kenzo))
        elif isinstance(object, tuple):
            rslt = self._kenzo
            for mrph in object:
                rslt = __sbtr__(rslt, mrph._kenzo)
            return KenzoChainComplexMorphism(rslt)

    def change_source_target_complex(self, source=None, target=None):
        r"""
        Build, from the morphism ``self``, a new morphism with ``source``
        and ``target`` as source and target Kenzo chain complexes, respectively.

        INPUT:

        - ``source`` -- A KenzoChainComplex instance or None (default).

        - ``target`` -- A KenzoChainComplex instance or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism` inheriting from ``self`` the
          degree (:degr slot in Kenzo), the algorithm (:intr slot in Kenzo)
          and the strategy (:strt slot in Kenzo). The source and target slots
          of this new morphism are given by the parameters ``source`` and
          ``target`` respectively; if any parameter is ommited, the corresponding
          slot is inherited from ``self``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere, KenzoChainComplex # optional - kenzo
            sage: from sage.libs.ecl import ecl_eval
            sage: ZCC = KenzoChainComplex(ecl_eval("(z-chcm)"))               # optional - kenzo
            sage: ZCC                                                         # optional - kenzo
            [K... Chain-Complex]
            sage: s2 = Sphere(2)                                              # optional - kenzo
            sage: s3 = Sphere(3)                                              # optional - kenzo
            sage: tp = s2.tensor_product(s3)                                  # optional - kenzo
            sage: tp                                                          # optional - kenzo
            [K... Filtered-Chain-Complex]
            sage: null = ZCC.null_morphism(tp)                                # optional - kenzo
            sage: null                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null.source_complex()                                       # optional - kenzo
            [K... Chain-Complex]
            sage: null2 = null.change_source_target_complex(source = tp)      # optional - kenzo
            sage: null2                                                       # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null2.source_complex()                                      # optional - kenzo
            [K... Filtered-Chain-Complex]
        """
        source = source or self.source_complex()
        target = target or self.target_complex()
        return KenzoChainComplexMorphism(
            __change_sorc_trgt_aux__(self._kenzo, source._kenzo, target._kenzo))

    def destructive_change_source_target_complex(self, source=None, target=None):
        r"""
        Modify destructively the morphism ``self`` taking ``source`` and ``target`` as source and
        target Kenzo chain complexes of ``self``, respectively.

        INPUT:

        - ``source`` -- A KenzoChainComplex instance or None (default).

        - ``target`` -- A KenzoChainComplex instance or None (default).

        OUTPUT:

        - A :class:`KenzoChainComplexMorphism`. The source and target slots of ``self`` are replaced
          respectively by the parameters ``source`` and ``target``; if any parameter is ommited, the
          corresponding slot is inherited from ``self``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere, KenzoChainComplex # optional - kenzo
            sage: from sage.libs.ecl import ecl_eval
            sage: ZCC = KenzoChainComplex(ecl_eval("(z-chcm)"))               # optional - kenzo
            sage: ZCC                                                         # optional - kenzo
            [K... Chain-Complex]
            sage: s2 = Sphere(2)                                              # optional - kenzo
            sage: s3 = Sphere(3)                                              # optional - kenzo
            sage: tp = s2.tensor_product(s3)                                  # optional - kenzo
            sage: tp                                                          # optional - kenzo
            [K... Filtered-Chain-Complex]
            sage: null = ZCC.null_morphism(tp)                                # optional - kenzo
            sage: null                                                        # optional - kenzo
            [K... Morphism (degree 0): K... -> K...]
            sage: null.target_complex()                                       # optional - kenzo
            [K... Filtered-Chain-Complex]
            sage: null.destructive_change_source_target_complex(target = ZCC) # optional - kenzo
            [K... Cohomology-Class on K... of degree 0]
            sage: null.target_complex()                                       # optional - kenzo
            [K... Chain-Complex]
        """
        source = source or self.source_complex()
        target = target or self.target_complex()
        return KenzoChainComplexMorphism(
            __dstr_change_sorc_trgt_aux__(self._kenzo, source._kenzo, target._kenzo))


def build_morphism(source_complex, target_complex, degree, algorithm, strategy, orgn):
    r"""
    Build a morphism of chain complexes by means of the corresponding build-mrph Kenzo
    function.

    INPUT:

    - ``source_complex`` -- The source object as a KenzoChainComplex instance

    - ``target_complex`` -- The target object as a KenzoChainComplex instance

    - ``degree`` -- An integer number representing the degree of the morphism

    - ``algorithm`` -- A Lisp function defining the mapping (:intr slot in Kenzo)

    - ``strategy`` -- The strategy (:strt slot in Kenzo), which must be one of
      the two strings ``gnrt`` or ``cmbn``, depending if the ``algorithm`` (a Lisp
      function) uses as arguments a degree and a generator or a combination,
      respectively.

    - ``orgn`` -- A list containing a description about the origin of the morphism

    OUTPUT:

    - A :class:`KenzoChainComplexMorphism`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KenzoChainComplex,\
        ....: build_morphism                                            # optional - kenzo
        sage: from sage.libs.ecl import ecl_eval
        sage: ZCC = KenzoChainComplex(ecl_eval("(z-chcm)"))             # optional - kenzo
        sage: A = build_morphism(ZCC, ZCC, -1,\
        ....: ecl_eval("#'(lambda (comb) (cmbn (1- (degr comb))))"),\
        ....: "cmbn", ["zero morphism on ZCC"])                         # optional - kenzo
        sage: A.target_complex()                                        # optional - kenzo
        [K... Chain-Complex]
        sage: A.degree()                                                # optional - kenzo
        -1
        sage: type(A)                                                   # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoChainComplexMorphism'>
    """
    return KenzoChainComplexMorphism(
        __build_mrph_aux__(source_complex._kenzo, target_complex._kenzo,
                       degree, algorithm, ":"+strategy, orgn))


def morphism_dictmat(morphism):
    r"""
    Computes a list of matrices in ECL associated to a morphism in Sage.

    INPUT:

    - ``morphism`` -- A morphism of chain complexes

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import morphism_dictmat    # optional - kenzo
        sage: X = simplicial_complexes.Simplex(1)
        sage: Y = simplicial_complexes.Simplex(0)
        sage: g = Hom(X,Y)({0:0, 1:0})
        sage: f = g.associated_chain_complex_morphism()
        sage: morphism_dictmat(f)                                   # optional - kenzo
        <ECL: ((2 . #2A()) (1 . #2A()) (0 . #2A((1 1))))>
    """
    rslt = EclObject([])
    source = morphism.domain()
    d = source.differential()
    for k in d.keys():
        rslt = EclObject(k).cons(s2k_matrix(morphism.in_degree(k))).cons(rslt)
    return rslt


def KChainComplexMorphism(morphism):
    r"""
    Construct a KenzoChainComplexMorphism from a ChainComplexMorphism in Sage.

    INPUT:

    - ``morphism`` -- A morphism of chain complexes

    OUTPUT:

    - A :class:`KenzoChainComplexMorphism`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplexMorphism           # optional - kenzo
        sage: C = ChainComplex({0: identity_matrix(ZZ, 1)})
        sage: D = ChainComplex({0: zero_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
        sage: f = Hom(C,D)({0: identity_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
        sage: g = KChainComplexMorphism(f)                                      # optional - kenzo
        sage: g                                                                 # optional - kenzo
        [K... Morphism (degree 0): K... -> K...]
        sage: g.source_complex()                                                # optional - kenzo
        [K... Chain-Complex]
        sage: g.target_complex()                                                # optional - kenzo
        [K... Chain-Complex]
    """
    source = KChainComplex(morphism.domain())
    target = KChainComplex(morphism.codomain())
    matrix_list = morphism_dictmat(morphism)
    return KenzoChainComplexMorphism(
        __kmorphismchaincomplex_aux1__(matrix_list, source._kenzo, target._kenzo))


def s2k_listofmorphisms(l):
    r"""
    Computes a list of morphisms of chain complexes in Kenzo from a list of morphisms in Sage.

    INPUT:

    - ``l`` -- A list of morphisms of chain complexes

    OUTPUT:

    - A :class:`EclObject`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import s2k_listofmorphisms  # optional - kenzo
        sage: C1 = ChainComplex({1: matrix(ZZ, 0, 2, [])}, degree_of_differential=-1)
        sage: C2 = ChainComplex({1: matrix(ZZ, 1, 2, [1, 0])},degree_of_differential=-1)
        sage: C3 = ChainComplex({0: matrix(ZZ, 0,2 , [])},degree_of_differential=-1)
        sage: M1 = Hom(C2,C1)({1: matrix(ZZ, 2, 2, [2, 0, 0, 2])})
        sage: M2 = Hom(C3,C2)({0: matrix(ZZ, 1, 2, [2, 0])})
        sage: l = [M1, M2]
        sage: s2k_listofmorphisms(l)                                 # optional - kenzo
        <ECL: ([K... Morphism (degree 0): K... -> K...] [K... Morphism (degree 0): K... -> K...])>
    """
    rslt = EclObject([])
    for m in l:
        rslt = EclObject(KChainComplexMorphism(m)._kenzo).cons(rslt)
    return __nreverse__(rslt)


def BicomplexSpectralSequence(l):
    r"""
    Construct the spectral sequence associated to the bicomplex given by a list of morphisms.

    INPUT:

    - ``l`` -- A list of morphisms of chain complexes

    OUTPUT:

    - A :class:`KenzoSpectralSequence`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import BicomplexSpectralSequence # optional - kenzo
        sage: C1 = ChainComplex({1: matrix(ZZ, 0, 2, [])}, degree_of_differential=-1)
        sage: C2 = ChainComplex({1: matrix(ZZ, 1, 2, [1, 0])},degree_of_differential=-1)
        sage: C3 = ChainComplex({0: matrix(ZZ, 0,2 , [])},degree_of_differential=-1)
        sage: M1 = Hom(C2,C1)({1: matrix(ZZ, 2, 2, [2, 0, 0, 2])})
        sage: M2 = Hom(C3,C2)({0: matrix(ZZ, 1, 2, [2, 0])})
        sage: l = [M1, M2]
        sage: E = BicomplexSpectralSequence(l)                        # optional - kenzo
        sage: E.group(2,0,1)                                          # optional - kenzo
        Additive abelian group isomorphic to Z/2 + Z
        sage: E.table(3,0,2,0,2)                                      # optional - kenzo
        0           0   0
        Z/2 + Z/4   0   0
        0           0   Z
        sage: E.matrix(2,2,0)                                         # optional - kenzo
        [ 0  0]
        [-4  0]
    """
    return KenzoSpectralSequence(__bicomplex_spectral_sequence__(s2k_listofmorphisms(l)))
