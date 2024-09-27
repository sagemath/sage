# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of Python modules in the Sage library

All of these features are present in a monolithic installation of the Sage library,
such as the one made by the SageMath distribution.

The features are defined for the purpose of separately testing modularized
distributions such as :ref:`sagemath-categories <spkg_sagemath_categories>`
and :ref:`sagemath-repl <spkg_sagemath_repl>`.

Often, doctests in a module of the Sage library illustrate the
interplay with a range of different objects; this is a form of integration testing.
These objects may come from modules shipped in
other distributions. For example, :mod:`sage.structure.element`
(shipped by :ref:`sagemath-objects <spkg_sagemath_objects>`,
one of the most fundamental distributions) contains the
doctest::

    sage: G = SymmetricGroup(4)                                                         # needs sage.groups
    sage: g = G([2, 3, 4, 1])                                                           # needs sage.groups
    sage: g.powers(4)                                                                   # needs sage.groups
    [(), (1,2,3,4), (1,3)(2,4), (1,4,3,2)]

This test cannot pass when the distribution :ref:`sagemath-objects <spkg_sagemath_objects>`
is tested separately (in a virtual environment): In this situation,
:class:`SymmetricGroup` is not defined anywhere (and thus not present
in the top-level namespace).
Hence, we conditionalize this doctest on the presence of the feature
:class:`sage.groups <sage__groups>`.
"""

# *****************************************************************************
#       Copyright (C) 2021-2023 Matthias Koeppe
#                     2021      Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule, StaticFile
from .join_feature import JoinFeature


class SAGE_SRC(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the
    monolithic source tree of the Sage library.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import SAGE_SRC
            sage: isinstance(SAGE_SRC(), SAGE_SRC)
            True
        """
        from sage.env import SAGE_SRC
        # We check the file bin/sage-src-env-config.in, which by design is:
        # - never installed,
        # - not included in the sagemath-standard sdist,
        # - included only in one modularized sdist, of pkgs/sage-conf_pypi,
        #   where it appears in a subdirectory (sage_root/src/bin/)
        StaticFile.__init__(self, 'SAGE_SRC',
                            filename='bin/sage-src-env-config.in',
                            search_path=(SAGE_SRC,) if SAGE_SRC else ())


class sagemath_doc_html(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of the documentation
    of the Sage library in HTML format.

    Developers often use ``make build`` instead of ``make`` to avoid the
    long time it takes to compile the documentation. Although commands
    such as ``make ptest`` build the documentation before testing, other
    test commands such as ``make ptestlong-nodoc`` or ``./sage -t --all``
    do not.

    All doctests that refer to the built documentation need to be marked
    ``# needs sagemath_doc_html``.

    TESTS::

        sage: from sage.features.sagemath import sagemath_doc_html
        sage: sagemath_doc_html().is_present()                                          # needs sagemath_doc_html
        FeatureTestResult('sagemath_doc_html', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sagemath_doc_html
            sage: isinstance(sagemath_doc_html(), sagemath_doc_html)
            True
        """
        from sage.env import SAGE_DOC
        StaticFile.__init__(self, 'sagemath_doc_html',
                            filename='html',
                            search_path=(SAGE_DOC,),
                            spkg='sagemath_doc_html',
                            type='standard')


class sage__combinat(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.combinat`.

    EXAMPLES:

    Python modules that provide elementary combinatorial objects such as :mod:`sage.combinat.subset`,
    :mod:`sage.combinat.composition`, :mod:`sage.combinat.permutation` are always available;
    there is no need for an ``# optional/needs`` tag::

        sage: Permutation([1,2,3]).is_even()
        True
        sage: Permutation([6,1,4,5,2,3]).bruhat_inversions()
        [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [3, 4], [3, 5]]

    Use ``# needs sage.combinat`` for doctests that use any other Python modules
    from :mod:`sage.combinat`, for example :mod:`sage.combinat.tableau_tuple`::

        sage: TableauTuple([[[7,8,9]],[],[[1,2,3],[4,5],[6]]]).shape()                  # needs sage.combinat
        ([3], [], [3, 2, 1])

    Doctests that use Python modules from :mod:`sage.combinat` that involve trees,
    graphs, hypergraphs, posets, quivers, combinatorial designs,
    finite state machines etc. should be marked ``# needs sage.combinat sage.graphs``::

        sage: L = Poset({0: [1], 1: [2], 2:[3], 3:[4]})                                 # needs sage.combinat sage.graphs
        sage: L.is_chain()                                                              # needs sage.combinat sage.graphs
        True

    Doctests that use combinatorial modules/algebras, or root systems should use the tag
    ``# needs sage.combinat sage.modules``::

        sage: # needs sage.combinat sage.modules
        sage: A = SchurAlgebra(QQ, 2, 3)
        sage: a = A.an_element(); a
        2*S((1, 1, 1), (1, 1, 1)) + 2*S((1, 1, 1), (1, 1, 2))
         + 3*S((1, 1, 1), (1, 2, 2))
        sage: L = RootSystem(['A',3,1]).root_lattice()
        sage: PIR = L.positive_imaginary_roots(); PIR
        Positive imaginary roots of type ['A', 3, 1]

    Doctests that use lattices, semilattices, or Dynkin diagrams should use the tag
    ``# needs sage.combinat sage.graphs sage.modules``::

        sage: L = LatticePoset({0: [1,2], 1: [3], 2: [3,4], 3: [5], 4: [5]})            # needs sage.combinat sage.graphs sage.modules
        sage: L.meet_irreducibles()                                                     # needs sage.combinat sage.graphs sage.modules
        [1, 3, 4]

    TESTS::

        sage: from sage.features.sagemath import sage__combinat
        sage: sage__combinat().is_present()                                             # needs sage.combinat
        FeatureTestResult('sage.combinat', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__combinat
            sage: isinstance(sage__combinat(), sage__combinat)
            True
        """
        # sage.combinat will be a namespace package.
        # Testing whether sage.combinat itself can be imported is meaningless.
        # Some modules providing basic combinatorics are already included in sagemath-categories.
        # Hence, we test a Python module within the package.
        JoinFeature.__init__(self, 'sage.combinat',
                             [PythonModule('sage.combinat'),                        # namespace package
                              PythonModule('sage.combinat.tableau'),                # representative
                             ],
                             spkg='sagemath_combinat', type='standard')


class sage__geometry__polyhedron(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.geometry.polyhedron`.

    EXAMPLES:

    Doctests that use polyhedra, cones, geometric complexes, triangulations, etc. should use
    the tag ``# needs sage.geometry.polyhedron``::

        sage: co = polytopes.truncated_tetrahedron()                                    # needs sage.geometry.polyhedron
        sage: co.volume()                                                               # needs sage.geometry.polyhedron
        184/3

    Some constructions of polyhedra require additional tags::

        sage: # needs sage.combinat sage.geometry.polyhedron sage.rings.number_field
        sage: perm_a3_reg_nf = polytopes.generalized_permutahedron(
        ....:    ['A',3], regular=True, backend='number_field'); perm_a3_reg_nf
        A 3-dimensional polyhedron in AA^3 defined as the convex hull of 24 vertices

    TESTS::

        sage: from sage.features.sagemath import sage__geometry__polyhedron
        sage: sage__geometry__polyhedron().is_present()                                 # needs sage.geometry.polyhedron
        FeatureTestResult('sage.geometry.polyhedron', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__geometry__polyhedron
            sage: isinstance(sage__geometry__polyhedron(), sage__geometry__polyhedron)
            True
        """
        JoinFeature.__init__(self, 'sage.geometry.polyhedron',
                             [PythonModule('sage.geometry'),                        # namespace package
                              PythonModule('sage.geometry.polyhedron'),             # representative
                              PythonModule('sage.schemes.toric'),                   # namespace package
                              PythonModule('sage.schemes.toric.variety'),           # representative
                             ],
                             spkg='sagemath_polyhedra', type='standard')


class sage__graphs(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.graphs`.

    EXAMPLES:

    Doctests that use anything from :mod:`sage.graphs` (:class:`Graph`, :class:`DiGraph`, ...)
    should be marked ``# needs sage.graphs``. The same applies to any doctest that
    uses a :class:`~sage.combinat.posets.posets.Poset`, cluster algebra quiver, finite
    state machines, abelian sandpiles, or Dynkin diagrams::

        sage: g = graphs.PetersenGraph()                                                # needs sage.graphs
        sage: r, s = g.is_weakly_chordal(certificate=True); r                           # needs sage.graphs
        False

    Also any use of tree classes defined in :mod:`sage.combinat` (:class:`BinaryTree`,
    :class:`RootedTree`, ...) in doctests should be marked the same.

    By way of generalization, any use of :class:`SimplicialComplex` or other abstract complexes from
    :mod:`sage.topology`, hypergraphs, and combinatorial designs, should be marked
    ``# needs sage.graphs`` as well::

        sage: X = SimplicialComplex([[0,1,2], [1,2,3]])                                 # needs sage.graphs
        sage: X.link(Simplex([0]))                                                      # needs sage.graphs
        Simplicial complex with vertex set (1, 2) and facets {(1, 2)}

        sage: IncidenceStructure([[1,2,3],[1,4]]).degrees(2)                            # needs sage.graphs
        {(1, 2): 1, (1, 3): 1, (1, 4): 1, (2, 3): 1, (2, 4): 0, (3, 4): 0}

    On the other hand, matroids are not implemented as posets in Sage but are instead
    closely tied to linear algebra over fields; hence use ``# needs sage.modules`` instead::

        sage: # needs sage.modules
        sage: M = Matroid(Matrix(QQ, [[1, 0, 0, 0, 1, 1, 1],
        ....:                         [0, 1, 0, 1, 0, 1, 1],
        ....:                         [0, 0, 1, 1, 1, 0, 1]]))
        sage: N = (M / [2]).delete([3, 4])
        sage: sorted(N.groundset())
        [0, 1, 5, 6]

    However, many constructions (and some methods) of matroids do involve graphs::

        sage: # needs sage.modules
        sage: W = matroids.Wheel(3)     # despite the name, not created via graphs
        sage: W.is_isomorphic(N)           # goes through a graph isomorphism test      # needs sage.graphs
        False
        sage: K4 = matroids.CompleteGraphic(4)    # this one is created via graphs      # needs sage.graphs
        sage: K4.is_isomorphic(W)                                                       # needs sage.graphs
        True

    TESTS::

        sage: from sage.features.sagemath import sage__graphs
        sage: sage__graphs().is_present()                                               # needs sage.graphs
        FeatureTestResult('sage.graphs', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__graphs
            sage: isinstance(sage__graphs(), sage__graphs)
            True
        """
        JoinFeature.__init__(self, 'sage.graphs',
                             # These lists of modules are an (incomplete) duplication
                             # of information in the distribution's MANIFEST.
                             # But at least as long as the monolithic Sage library is
                             # around, we need this information here for use by
                             # sage-fixdoctests.
                             [PythonModule('sage.graphs'),                          # namespace package
                              PythonModule('sage.graphs.graph'),                    # representative
                              PythonModule('sage.combinat.designs'),                # namespace package
                              PythonModule('sage.combinat.designs.block_design'),   # representative
                              PythonModule('sage.combinat.posets'),                 # namespace package
                              PythonModule('sage.combinat.posets.posets'),          # representative
                              PythonModule('sage.topology'),                        # namespace package
                              PythonModule('sage.topology.simplicial_complex'),     # representative
                             ],
                             spkg='sagemath_graphs', type='standard')


class sage__groups(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``sage.groups``.

    EXAMPLES:

    Permutations and sets of permutations are always available, but permutation groups are
    implemented in Sage using the :ref:`GAP <spkg_gap>` system and require the tag
    ``# needs sage.groups``::

        sage: p = Permutation([2,1,4,3])
        sage: p.to_permutation_group_element()                                          # needs sage.groups
        (1,2)(3,4)

    TESTS::

        sage: from sage.features.sagemath import sage__groups
        sage: sage__groups().is_present()                                               # needs sage.groups
        FeatureTestResult('sage.groups', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__groups
            sage: isinstance(sage__groups(), sage__groups)
            True
        """
        JoinFeature.__init__(self, 'sage.groups',
                             [PythonModule('sage.groups.perm_gps.permgroup')],
                             spkg='sagemath_groups', type='standard')


class sage__libs__braiding(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.libs.braiding`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__libs__braiding
        sage: sage__libs__braiding().is_present()                                            # needs sage.libs.braiding
        FeatureTestResult('sage.libs.braiding', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__braiding
            sage: isinstance(sage__libs__braiding(), sage__libs__braiding)
            True
        """
        PythonModule.__init__(self, 'sage.libs.braiding',
                              spkg='sagemath_libbraiding', type='standard')


class sage__libs__ecl(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.libs.ecl`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__libs__ecl
        sage: sage__libs__ecl().is_present()                        # optional - sage.libs.ecl
        FeatureTestResult('sage.libs.ecl', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__ecl
            sage: isinstance(sage__libs__ecl(), sage__libs__ecl)
            True
        """
        PythonModule.__init__(self, 'sage.libs.ecl',
                              spkg='sagemath_symbolics', type='standard')


class sage__libs__flint(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.flint`
    and other modules depending on FLINT.

    In addition to the modularization purposes that this tag serves, it also provides attribution
    to the upstream project.

    TESTS::

        sage: from sage.features.sagemath import sage__libs__flint
        sage: sage__libs__flint().is_present()                                          # needs sage.libs.flint
        FeatureTestResult('sage.libs.flint', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__flint
            sage: isinstance(sage__libs__flint(), sage__libs__flint)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.flint',
                             [PythonModule('sage.libs.flint.arith_sage'),
                              PythonModule('sage.libs.flint.flint_sage')],
                             spkg='sagemath_flint', type='standard')


class sage__libs__gap(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.gap`
    (the library interface to :ref:`GAP <spkg_gap>`) and :mod:`sage.interfaces.gap` (the pexpect
    interface to GAP). By design, we do not distinguish between these two, in order
    to facilitate the conversion of code from the pexpect interface to the library
    interface.

    .. SEEALSO::

        :class:`Features for GAP packages <~sage.features.gap.GapPackage>`

    TESTS::

        sage: from sage.features.gap import sage__libs__gap
        sage: sage__libs__gap().is_present()                                            # needs sage.libs.gap
        FeatureTestResult('sage.libs.gap', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.gap import sage__libs__gap
            sage: isinstance(sage__libs__gap(), sage__libs__gap)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.gap',
                             [PythonModule('sage.libs.gap.libgap'),
                              PythonModule('sage.interfaces.gap'),
                              PythonModule('sage.groups.matrix_gps.finitely_generated_gap'),
                              PythonModule('sage.groups.matrix_gps.group_element_gap'),
                              PythonModule('sage.groups.matrix_gps.heisenberg'),
                              PythonModule('sage.groups.matrix_gps.isometries'),
                              PythonModule('sage.groups.matrix_gps.linear_gap'),
                              PythonModule('sage.groups.matrix_gps.matrix_group_gap'),
                              PythonModule('sage.groups.matrix_gps.named_group_gap'),
                              PythonModule('sage.groups.matrix_gps.orthogonal_gap'),
                              PythonModule('sage.groups.matrix_gps.symplectic_gap'),
                              PythonModule('sage.groups.matrix_gps.unitary_gap'),
                              PythonModule('sage.matrix.matrix_gap'),
                              PythonModule('sage.rings.universal_cyclotomic_field')])


class sage__libs__linbox(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.linbox`
    and other modules depending on Givaro, FFLAS-FFPACK, LinBox.

    In addition to the modularization purposes that this tag serves, it also provides attribution
    to the upstream project.

    TESTS::

        sage: from sage.features.sagemath import sage__libs__linbox
        sage: sage__libs__linbox().is_present()                                         # needs sage.libs.linbox
        FeatureTestResult('sage.libs.linbox', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__linbox
            sage: isinstance(sage__libs__linbox(), sage__libs__linbox)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.linbox',
                             [PythonModule('sage.rings.finite_rings.element_givaro'),
                              PythonModule('sage.matrix.matrix_modn_dense_float'),
                              PythonModule('sage.matrix.matrix_modn_dense_double')],
                             spkg='sagemath_linbox', type='standard')


class sage__libs__m4ri(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of Cython modules
    depending on the M4RI and/or M4RIe libraries.

    In addition to the modularization purposes that this tag serves,
    it also provides attribution to the upstream project.

    TESTS::

        sage: from sage.features.sagemath import sage__libs__m4ri
        sage: sage__libs__m4ri().is_present()                                           # needs sage.libs.m4ri
        FeatureTestResult('sage.libs.m4ri', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__m4ri
            sage: isinstance(sage__libs__m4ri(), sage__libs__m4ri)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.m4ri',
                             [PythonModule('sage.matrix.matrix_gf2e_dense'),
                              PythonModule('sage.matrix.matrix_mod2_dense')],
                             spkg='sagemath_m4ri', type='standard')


class sage__libs__ntl(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.ntl`
    and other modules depending on NTL.

    In addition to the modularization purposes that this tag serves,
    it also provides attribution to the upstream project.

    TESTS::

        sage: from sage.features.sagemath import sage__libs__ntl
        sage: sage__libs__ntl().is_present()                                            # needs sage.libs.ntl
        FeatureTestResult('sage.libs.ntl', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__ntl
            sage: isinstance(sage__libs__ntl(), sage__libs__ntl)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.ntl',
                             [PythonModule('sage.libs.ntl.convert')],
                             spkg='sagemath_ntl', type='standard')


class sage__libs__homfly(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.homfly`.

    In addition to the modularization purposes that this tag serves,
    it also provides attribution to the upstream project.

    TESTS::

        sage: from sage.features.sagemath import sage__libs__homfly
        sage: sage__libs__homfly().is_present()                                         # needs sage.libs.homfly
        FeatureTestResult('sage.libs.homfly', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__homfly
            sage: isinstance(sage__libs__homfly(), sage__libs__homfly)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.homfly',
                             [PythonModule('sage.libs.homfly')],
                             spkg='sagemath_homfly', type='standard')


class sage__libs__pari(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.libs.pari`.

    SageMath uses the :ref:`PARI <spkg_pari>` library (via :ref:`cypari2
    <spkg_cypari>`) for numerous purposes.  Doctests that involves such features
    should be marked ``# needs sage.libs.pari``.

    In addition to the modularization purposes that this tag serves, it also
    provides attribution to the upstream project.

    EXAMPLES::

        sage: R.<a> = QQ[]
        sage: S.<x> = R[]
        sage: f = x^2 + a; g = x^3 + a
        sage: r = f.resultant(g); r                                                     # needs sage.libs.pari
        a^3 + a^2

    TESTS::

        sage: from sage.features.sagemath import sage__libs__pari
        sage: sage__libs__pari().is_present()                                           # needs sage.libs.pari
        FeatureTestResult('sage.libs.pari', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__pari
            sage: isinstance(sage__libs__pari(), sage__libs__pari)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.pari',
                             [PythonModule('sage.libs.pari.convert_sage')],
                             spkg='sagemath_pari', type='standard')


class sage__libs__singular(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.singular`
    (the library interface to Singular) and :mod:`sage.interfaces.singular` (the pexpect
    interface to Singular). By design, we do not distinguish between these two, in order
    to facilitate the conversion of code from the pexpect interface to the library
    interface.

    .. SEEALSO::

        :class:`Feature singular <~sage.features.singular.Singular>`

    TESTS::

        sage: from sage.features.singular import sage__libs__singular
        sage: sage__libs__singular().is_present()                                       # needs sage.libs.singular
        FeatureTestResult('sage.libs.singular', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.singular import sage__libs__singular
            sage: isinstance(sage__libs__singular(), sage__libs__singular)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.singular',
                             [PythonModule('sage.libs.singular.singular'),
                              PythonModule('sage.interfaces.singular')])


class sage__modular(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.modular`.

    TESTS::

        sage: from sage.features.sagemath import sage__modular
        sage: sage__modular().is_present()                                              # needs sage.modular
        FeatureTestResult('sage.modular', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__modular
            sage: isinstance(sage__modular(), sage__modular)
            True
        """
        JoinFeature.__init__(self, 'sage.modular',
                             [PythonModule('sage.modular.modform.eisenstein_submodule')],
                             spkg='sagemath_schemes', type='standard')


class sage__modules(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.modules`.

    EXAMPLES:

    All uses of implementations of vector spaces / free modules in SageMath, whether
    :class:`sage.modules.free_module.FreeModule`,
    :class:`sage.combinat.free_module.CombinatorialFreeModule`,
    :class:`sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`, or
    additive abelian groups, should be marked ``# needs sage.modules``.

    The same holds for matrices, tensors, algebras, quadratic forms,
    point lattices, root systems, matrix/affine/Weyl/Coxeter groups, matroids,
    and ring derivations.

    Likewise, all uses of :mod:`sage.coding`, :mod:`sage.crypto`, and :mod:`sage.homology`
    in doctests should be marked ``# needs sage.modules``.

    TESTS::

        sage: from sage.features.sagemath import sage__modules
        sage: sage__modules().is_present()                                              # needs sage.modules
        FeatureTestResult('sage.modules', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__modules
            sage: isinstance(sage__modules(), sage__modules)
            True
        """
        JoinFeature.__init__(self, 'sage.modules',
                             [PythonModule('sage.modules'),                         # namespace package
                              PythonModule('sage.modules.free_module'),             # representative
                              PythonModule('sage.matrix'),                          # namespace package
                              PythonModule('sage.matrix.matrix2'),                  # representative
                              PythonModule('sage.combinat.free_module'),
                              PythonModule('sage.quadratic_forms'),                 # namespace package
                              PythonModule('sage.quadratic_forms.quadratic_form'),  # representative
                              PythonModule('sage.groups.additive_abelian'),         # namespace package
                              PythonModule('sage.groups.additive_abelian.qmodnz'),  # representative
                              PythonModule('sage.groups.affine_gps'),               # namespace package
                              PythonModule('sage.groups.affine_gps.affine_group'),  # representative
                              PythonModule('sage.groups.matrix_gps'),               # namespace package
                              PythonModule('sage.groups.matrix_gps.named_group'),   # representative
                              PythonModule('sage.homology'),                        # namespace package
                              PythonModule('sage.homology.chain_complex'),          # representative
                              PythonModule('sage.matroids'),                        # namespace package
                              PythonModule('sage.matroids.matroid'),                # representative
                             ],
                             spkg='sagemath_modules', type='standard')


class sage__numerical__mip(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.numerical.mip`.

    TESTS::

        sage: from sage.features.sagemath import sage__numerical__mip
        sage: sage__numerical__mip().is_present()                                       # needs sage.numerical.mip
        FeatureTestResult('sage.numerical.mip', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__numerical__mip
            sage: isinstance(sage__numerical__mip(), sage__numerical__mip)
            True
        """
        PythonModule.__init__(self, 'sage.numerical.mip',
                              spkg='sagemath_polyhedra')


class sage__plot(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.plot`.

    TESTS::

        sage: from sage.features.sagemath import sage__plot
        sage: sage__plot().is_present()                                                 # needs sage.plot
        FeatureTestResult('sage.plot', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__plot
            sage: isinstance(sage__plot(), sage__plot)
            True
        """
        JoinFeature.__init__(self, 'sage.plot',
                             [PythonModule('sage.plot.plot')],
                             spkg='sagemath_plot', type='standard')


class sage__rings__complex_double(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.complex_double`.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__complex_double
        sage: sage__rings__complex_double().is_present()                                # needs sage.rings.complex_double
        FeatureTestResult('sage.rings.complex_double', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__complex_double
            sage: isinstance(sage__rings__complex_double(), sage__rings__complex_double)
            True
        """
        PythonModule.__init__(self, 'sage.rings.complex_double',
                              spkg='sagemath_modules', type='standard')


class sage__rings__finite_rings(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.finite_rings`;
    specifically, the element implementations using the :ref:`PARI <spkg_pari>` library.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__finite_rings
        sage: sage__rings__finite_rings().is_present()                                  # needs sage.rings.finite_rings
        FeatureTestResult('sage.rings.finite_rings', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__finite_rings
            sage: isinstance(sage__rings__finite_rings(), sage__rings__finite_rings)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.finite_rings',
                             [PythonModule('sage.rings.finite_rings.element_pari_ffelt'),
                              PythonModule('sage.rings.algebraic_closure_finite_field'),
                              sage__libs__pari()],
                             type='standard')


class sage__rings__function_field(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.function_field`.

    EXAMPLES:

    Rational function fields are always available::

        sage: K.<x> = FunctionField(QQ)
        sage: K.maximal_order()
        Maximal order of Rational function field in x over Rational Field

    Use the tag ``# needs sage.rings.function_field`` whenever extensions
    of function fields (by adjoining a root of a univariate polynomial) come into play::

        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L                         # needs sage.rings.function_field
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    Such extensions of function fields are implemented using Gröbner bases of polynomial rings;
    Sage makes essential use of the :ref:`Singular <spkg_singular>` system for this.
    (It is not necessary to use the tag ``# needs sage.libs.singular``; it is
    implied by ``# needs sage.rings.function_field``.)

    TESTS::

        sage: from sage.features.sagemath import sage__rings__function_field
        sage: sage__rings__function_field().is_present()                                # needs sage.rings.function_field
        FeatureTestResult('sage.rings.function_field', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__function_field
            sage: isinstance(sage__rings__function_field(), sage__rings__function_field)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.function_field',
                             [PythonModule('sage.rings.function_field.function_field_polymod'),
                              sage__libs__singular()],
                             type='standard')


class sage__rings__number_field(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.number_field`.

    Number fields are implemented in Sage using a complicated mixture of various libraries,
    including :ref:`FLINT <spkg_flint>`, :ref:`GAP <spkg_gap>`,
    :ref:`MPFI <spkg_mpfi>`, :ref:`NTL <spkg_ntl>`, and :ref:`PARI <spkg_pari>`.

    EXAMPLES:

    Rational numbers are, of course, always available::

        sage: QQ in NumberFields()
        True

    Doctests that construct algebraic number fields should be marked ``# needs sage.rings.number_field``::

        sage: # needs sage.rings.number_field
        sage: K.<cuberoot2> = NumberField(x^3 - 2)
        sage: L.<cuberoot3> = K.extension(x^3 - 3)
        sage: S.<sqrt2> = L.extension(x^2 - 2); S
        Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field

        sage: # needs sage.rings.number_field
        sage: K.<zeta> = CyclotomicField(15)
        sage: CC(zeta)
        0.913545457642601 + 0.406736643075800*I

    Doctests that make use of the algebraic field ``QQbar`` or the algebraic real field ``AA``
    should be marked likewise::

        sage: # needs sage.rings.number_field
        sage: AA(-1)^(1/3)
        -1
        sage: QQbar(-1)^(1/3)
        0.500000000000000? + 0.866025403784439?*I

    Use of the universal cyclotomic field should be marked
    ``# needs sage.libs.gap sage.rings.number_field``.

        sage: # needs sage.libs.gap sage.rings.number_field
        sage: UCF = UniversalCyclotomicField(); UCF
        Universal Cyclotomic Field
        sage: E = UCF.gen
        sage: f = E(2) + E(3); f
        2*E(3) + E(3)^2
        sage: f.galois_conjugates()
        [2*E(3) + E(3)^2, E(3) + 2*E(3)^2]

    TESTS::

        sage: from sage.features.sagemath import sage__rings__number_field
        sage: sage__rings__number_field().is_present()                                  # needs sage.rings.number_field
        FeatureTestResult('sage.rings.number_field', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__number_field
            sage: isinstance(sage__rings__number_field(), sage__rings__number_field)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.number_field',
                             [PythonModule('sage.rings.number_field.number_field_element'),
                              PythonModule('sage.rings.universal_cyclotomic_field'),
                              PythonModule('sage.rings.qqbar')],
                             type='standard')


class sage__rings__padics(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``sage.rings.padics``.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__padics
        sage: sage__rings__padics().is_present()                                        # needs sage.rings.padics
        FeatureTestResult('sage.rings.padics', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__padics
            sage: isinstance(sage__rings__padics(), sage__rings__padics)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.padics',
                             [PythonModule('sage.rings.padics.factory')],
                             type='standard')


class sage__rings__polynomial__pbori(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.rings.polynomial.pbori`.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__polynomial__pbori
        sage: sage__rings__polynomial__pbori().is_present()                             # needs sage.rings.polynomial.pbori
        FeatureTestResult('sage.rings.polynomial.pbori', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__polynomial__pbori
            sage: isinstance(sage__rings__polynomial__pbori(), sage__rings__polynomial__pbori)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.polynomial.pbori',
                             [PythonModule('sage.rings.polynomial.pbori.pbori')],
                             spkg='sagemath_brial', type='standard')


class sage__rings__real_double(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.real_double`.

    EXAMPLES:

    The Real Double Field is basically always available, and no ``# optional/needs`` tag is needed::

        sage: RDF.characteristic()
        0

    The feature exists for use in doctests of Python modules that are shipped by the
    most fundamental distributions.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__real_double
        sage: sage__rings__real_double().is_present()                                   # needs sage.rings.real_double
        FeatureTestResult('sage.rings.real_double', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__real_double
            sage: isinstance(sage__rings__real_double(), sage__rings__real_double)
            True
        """
        PythonModule.__init__(self, 'sage.rings.real_double', type='standard')


class sage__rings__real_mpfr(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.real_mpfr`.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__real_mpfr
        sage: sage__rings__real_mpfr().is_present()                                     # needs sage.rings.real_mpfr
        FeatureTestResult('sage.rings.real_mpfr', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__real_mpfr
            sage: isinstance(sage__rings__real_mpfr(), sage__rings__real_mpfr)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.real_mpfr',
                             [PythonModule('sage.rings.real_mpfr'),
                              PythonModule('sage.rings.complex_mpfr'),
                             ],
                             spkg='sagemath_modules', type='standard')


class sage__sat(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.sat`.

    TESTS::

        sage: from sage.features.sagemath import sage__sat
        sage: sage__sat().is_present()                                                  # needs sage.sat
        FeatureTestResult('sage.sat', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__sat
            sage: isinstance(sage__sat(), sage__sat)
            True
        """
        JoinFeature.__init__(self, 'sage.sat',
                             [PythonModule('sage.sat.expression')],
                             spkg='sagemath_combinat', type='standard')


class sage__schemes(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.schemes`.

    TESTS::

        sage: from sage.features.sagemath import sage__schemes
        sage: sage__schemes().is_present()                                              # needs sage.schemes
        FeatureTestResult('sage.schemes', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__schemes
            sage: isinstance(sage__schemes(), sage__schemes)
            True
        """
        JoinFeature.__init__(self, 'sage.schemes',
                             [PythonModule('sage.schemes.elliptic_curves.ell_generic')],
                             spkg='sagemath_schemes', type='standard')


class sage__symbolic(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.symbolic`.

    EXAMPLES:

    The symbolics subsystem of Sage will be provided by the distribution
    sagemath-symbolics, in preparation at :issue:`35095`. If it is not installed,
    Sage will be able to provide installation advice::

        sage: from sage.features.sagemath import sage__symbolic
        sage: print(sage__symbolic().resolution())                                      # optional - sage_spkg, not tested
        ...To install sagemath_symbolics...you can try to run...
        pip install sagemath-symbolics
        ...

    TESTS::

        sage: from sage.features.sagemath import sage__symbolic
        sage: sage__symbolic().is_present()                                             # needs sage.symbolic
        FeatureTestResult('sage.symbolic', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__symbolic
            sage: isinstance(sage__symbolic(), sage__symbolic)
            True
        """
        JoinFeature.__init__(self, 'sage.symbolic',
                             [PythonModule('sage.symbolic.expression'),
                              PythonModule('sage.manifolds'),
                              PythonModule('sage.calculus.calculus'),
                              PythonModule('sage.calculus.desolvers'),
                              PythonModule('sage.calculus.predefined'),
                              PythonModule('sage.calculus.tests'),
                              PythonModule('sage.calculus.var'),
                              PythonModule('sage.geometry.riemannian_manifolds'),
                              PythonModule('sage.geometry.hyperbolic_space'),
                              PythonModule('sage.dynamics.complex_dynamics'),
                              PythonModule('sage.libs.ecl'),
                              PythonModule('sage.interfaces.fricas'),
                              PythonModule('sage.interfaces.giac'),
                              PythonModule('sage.interfaces.magma'),
                              PythonModule('sage.interfaces.magma_free'),
                              PythonModule('sage.interfaces.maple'),
                              PythonModule('sage.interfaces.mathematica'),
                              PythonModule('sage.interfaces.mathics'),
                              PythonModule('sage.interfaces.maxima'),
                              PythonModule('sage.interfaces.maxima_abstract'),
                              PythonModule('sage.interfaces.maxima_lib'),
                              PythonModule('sage.interfaces.qepcad'),
                              PythonModule('sage.interfaces.sympy'),
                              PythonModule('sage.interfaces.sympy_wrapper'),
                             ], spkg='sagemath_symbolics', type='standard')


def all_features():
    r"""
    Return features corresponding to parts of the Sage library.

    These features are named after Python packages/modules (e.g., :mod:`sage.symbolic`),
    not distribution packages (**sagemath-symbolics**).

    This design is motivated by a separation of concerns: The author of a module that depends
    on some functionality provided by a Python module usually already knows the
    name of the Python module, so we do not want to force the author to also
    know about the distribution package that provides the Python module.

    Instead, we associate distribution packages to Python modules in
    :mod:`sage.features.sagemath` via the ``spkg`` parameter of
    :class:`~sage.features.Feature`.

    EXAMPLES::

        sage: from sage.features.sagemath import all_features
        sage: list(all_features())
        [...Feature('sage.combinat'), ...]
    """
    return [SAGE_SRC(),
            sagemath_doc_html(),
            sage__combinat(),
            sage__geometry__polyhedron(),
            sage__graphs(),
            sage__groups(),
            sage__libs__braiding(),
            sage__libs__ecl(),
            sage__libs__flint(),
            sage__libs__gap(),
            sage__libs__homfly(),
            sage__libs__linbox(),
            sage__libs__m4ri(),
            sage__libs__ntl(),
            sage__libs__pari(),
            sage__libs__singular(),
            sage__modular(),
            sage__modules(),
            sage__numerical__mip(),
            sage__plot(),
            sage__rings__complex_double(),
            sage__rings__finite_rings(),
            sage__rings__function_field(),
            sage__rings__number_field(),
            sage__rings__padics(),
            sage__rings__polynomial__pbori(),
            sage__rings__real_double(),
            sage__rings__real_mpfr(),
            sage__sat(),
            sage__schemes(),
            sage__symbolic()]
