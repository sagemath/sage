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

    sage: G = SymmetricGroup(4)                                                         # optional - sage.groups
    sage: g = G([2, 3, 4, 1])                                                           # optional - sage.groups
    sage: g.powers(4)                                                                   # optional - sage.groups
    [(), (1,2,3,4), (1,3)(2,4), (1,4,3,2)]

This test cannot pass when the distribution :ref:`sagemath-objects <sagemath-objects>`
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
from .singular import sage__libs__singular


class sagemath_doc_html(StaticFile):
    r"""
    A :class:`Feature` which describes the presence of the documentation
    of the Sage library in HTML format.

    Developers often use ``make build`` instead of ``make`` to avoid the
    long time it takes to compile the documentation. Although commands
    such as ``make ptest`` build the documentation before testing, other
    test commands such as ``make ptestlong-nodoc`` or ``./sage -t --all``
    do not.

    All doctests that refer to the built documentation need to be marked
    ``# optional - sagemath_doc_html``.

    TESTS::

        sage: from sage.features.sagemath import sagemath_doc_html
        sage: sagemath_doc_html().is_present()                                          # optional - sagemath_doc_html
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
                            spkg='sagemath_doc_html')


class sage__combinat(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.combinat`.

    EXAMPLES:

    Python modules that provide elementary combinatorial objects such as :mod:`sage.combinat.subset`,
    :mod:`sage.combinat.composition`, :mod:`sage.combinat.permutation` are always available;
    there is no need for an ``# optional`` annotation::

        sage: Permutation([1,2,3]).is_even()
        True
        sage: Permutation([6,1,4,5,2,3]).bruhat_inversions()
        [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [3, 4], [3, 5]]

    Use ``# optional - sage.combinat`` for doctests that use any other Python modules
    from :mod:`sage.combinat`, for example :mod:`sage.combinat.tableau_tuple`::

        sage: TableauTuple([[[7,8,9]],[],[[1,2,3],[4,5],[6]]]).shape()                  # optional - sage.combinat
        ([3], [], [3, 2, 1])

    Doctests that use Python modules from :mod:`sage.combinat` that involve trees,
    graphs, hypergraphs, posets, quivers, combinatorial designs,
    finite state machines etc. should be marked ``# optional - sage.combinat sage.graphs``::

        sage: L = Poset({0: [1], 1: [2], 2:[3], 3:[4]})                                 # optional - sage.combinat sage.graphs
        sage: L.is_chain()                                                              # optional - sage.combinat sage.graphs
        True

    Doctests that use combinatorial modules/algebras, or root systems should use the annotation
    ``# optional - sage.combinat sage.modules``::

        sage: A = SchurAlgebra(QQ, 2, 3)                                                # optional - sage.combinat sage.modules
        sage: a = A.an_element(); a                                                     # optional - sage.combinat sage.modules
        2*S((1, 1, 1), (1, 1, 1)) + 2*S((1, 1, 1), (1, 1, 2))
         + 3*S((1, 1, 1), (1, 2, 2))
        sage: L = RootSystem(['A',3,1]).root_lattice()                                  # optional - sage.combinat sage.modules
        sage: PIR = L.positive_imaginary_roots(); PIR                                   # optional - sage.combinat sage.modules
        Positive imaginary roots of type ['A', 3, 1]

    Doctests that use lattices, semilattices, or Dynkin diagrams should use the annotation
    ``# optional - sage.combinat sage.graphs sage.modules``::

        sage: L = LatticePoset({0: [1,2], 1: [3], 2: [3,4], 3: [5], 4: [5]})            # optional - sage.combinat sage.graphs sage.modules
        sage: L.meet_irreducibles()                                                     # optional - sage.combinat sage.graphs sage.modules
        [1, 3, 4]

    TESTS::

        sage: from sage.features.sagemath import sage__combinat
        sage: sage__combinat().is_present()                                             # optional - sage.combinat
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
        # Hence, we test a Python module within the package.
        JoinFeature.__init__(self, 'sage.combinat',
                             [PythonModule('sage.combinat.combination')])


class sage__geometry__polyhedron(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.geometry.polyhedron`.

    EXAMPLES:

    Doctests that use polyhedra, cones, geometric complexes, triangulations, etc. should use
    the annotations ``# optional - sage.geometry.polyhedron``::

        sage: co = polytopes.truncated_tetrahedron()                                    # optional - sage.geometry.polyhedron
        sage: co.volume()                                                               # optional - sage.geometry.polyhedron
        184/3

    Some constructions of polyhedra require additional annotations::

        sage: perm_a3_reg_nf = polytopes.generalized_permutahedron(                     # optional - sage.combinat sage.geometry.polyhedron sage.rings.number_field
        ....:    ['A',3], regular=True, backend='number_field'); perm_a3_reg_nf
        A 3-dimensional polyhedron in AA^3 defined as the convex hull of 24 vertices

    TESTS::

        sage: from sage.features.sagemath import sage__geometry__polyhedron
        sage: sage__geometry__polyhedron().is_present()                                 # optional - sage.geometry.polyhedron
        FeatureTestResult('sage.geometry.polyhedron', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__geometry__polyhedron
            sage: isinstance(sage__geometry__polyhedron(), sage__geometry__polyhedron)
            True
        """
        PythonModule.__init__(self, 'sage.geometry.polyhedron')


class sage__graphs(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.graphs`.

    EXAMPLES:


    Matroids are not implemented as posets in Sage but are instead closely tied to
    linear algebra over fields; hence use ``# optional - sage.modules``::

        sage: M = Matroid(Matrix(QQ, [[1, 0, 0, 0, 1, 1, 1],                            # optional - sage.modules
        ....:                         [0, 1, 0, 1, 0, 1, 1],
        ....:                         [0, 0, 1, 1, 1, 0, 1]]))
        sage: N = M / [2] \ [3, 4]                                                      # optional - sage.modules
        sage: sorted(N.groundset())                                                     # optional - sage.modules
        [0, 1, 5, 6]

    However, many constructions (and some methods) of matroids do involve graphs::

        sage: W = matroids.Wheel(3)     # despite the name, not created via graphs      # optional - sage.modules
        sage: W.is_isomorphic(N)        # goes through a graph isomorphism test         # optional - sage.graphs sage.modules
        False
        sage: K4 = matroids.CompleteGraphic(4)    # this one is created via graphs      # optional - sage.graphs sage.modules
        sage: K4.is_isomorphic(W)                                                       # optional - sage.graphs sage.modules
        True

    TESTS::

        sage: from sage.features.sagemath import sage__graphs
        sage: sage__graphs().is_present()                                               # optional - sage.graphs
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
                             [PythonModule('sage.graphs.graph')])


class sage__groups(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.groups``.

    EXAMPLES:

    Permutations and sets of permutations are always available, but permutation groups are
    implemented in Sage using the :ref:`GAP <spkg_gap>` system and require the annotation
    ``# optional - sage.groups``::

        sage: p = Permutation([2,1,4,3])
        sage: p.to_permutation_group_element()                                          # optional - sage.groups
        (1,2)(3,4)

    TESTS::

        sage: from sage.features.sagemath import sage__groups
        sage: sage__groups().is_present()                                               # optional - sage.groups
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
                             [PythonModule('sage.groups.perm_gps.permgroup')])


class sage__libs__pari(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.pari`.

    SageMath uses the :ref:`PARI <spkg_pari>` library (via :ref:`cypari2 <spkg_cypari>`) for numerous purposes.
    Doctests that involves such features should be marked ``# optional - sage.libs.pari``.

    In addition to the modularization purposes that this annotation serves, it also provides attribution
    to the upstream project.

    EXAMPLES::

        sage: R.<a> = QQ[]
        sage: S.<x> = R[]
        sage: f = x^2 + a; g = x^3 + a
        sage: r = f.resultant(g); r                                                     # optional - sage.libs.pari
        a^3 + a^2

    TESTS::

        sage: from sage.features.sagemath import sage__libs__pari
        sage: sage__libs__pari().is_present()                                           # optional - sage.libs.pari
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
                             [PythonModule('sage.libs.pari.convert_sage')])


class sage__modules(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.modules`.

    TESTS::

        sage: from sage.features.sagemath import sage__modules
        sage: sage__modules().is_present()  # optional - sage.modules
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
                             [PythonModule('sage.modules.free_module')])


class sage__plot(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.plot`.

    TESTS::

        sage: from sage.features.sagemath import sage__plot
        sage: sage__plot().is_present()  # optional - sage.plot
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
                             [PythonModule('sage.plot.plot')])


class sage__rings__finite_rings(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.finite_rings`;
    specifically, the element implementations using the :ref:`PARI <spkg_pari>` library.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__finite_rings
        sage: sage__rings__finite_rings().is_present()                                  # optional - sage.rings.finite_rings
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
                             [PythonModule('sage.rings.finite_rings.element_pari_ffelt')])


class sage__rings__function_field(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.function_field`.

    EXAMPLES:

    Rational function fields are always available::

        sage: K.<x> = FunctionField(QQ)
        sage: K.maximal_order()
        Maximal order of Rational function field in x over Rational Field

    Use the annotation ``# optional - sage.rings.function_field`` whenever extensions
    of function fields (by adjoining a root of a univariate polynomial) come into play::

        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L                         # optional - sage.rings.function_field
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    Such extensions of function fields are implemented using Gr\"obner bases of polynomial rings;
    Sage makes essential use of the :ref:`Singular <spkg_singular>` system for this.
    (It is not necessary to use the annotation ``# optional - sage.libs.singular``; it is
    implied by ``# optional - sage.rings.function_field``.)

    TESTS::

        sage: from sage.features.sagemath import sage__rings__function_field
        sage: sage__rings__function_field().is_present()                                # optional - sage.rings.function_field
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
                              sage__libs__singular()])


class sage__rings__number_field(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.number_field`.

    Number fields are implemented in Sage using a complicated mixture of various libraries,
    including :ref:`arb <spkg_arb>`, :ref:`FLINT <spkg_flint>`, :ref:`GAP <spkg_gap>`,
    :ref:`MPFI <spkg_mpfi>`, :ref:`NTL <spkg_ntl>`, and :ref:`PARI <spkg_pari>`.

    EXAMPLES:

    Rational numbers are, of course, always available::

        sage: QQ in NumberFields()
        True

    Doctests that construct algebraic number fields should be marked ``# optional - sage.rings.number_field``::

        sage: K.<cuberoot2> = NumberField(x^3 - 2)                                      # optional - sage.rings.number_field
        sage: L.<cuberoot3> = K.extension(x^3 - 3)                                      # optional - sage.rings.number_field
        sage: S.<sqrt2> = L.extension(x^2 - 2); S                                       # optional - sage.rings.number_field
        Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field

        sage: K.<zeta> = CyclotomicField(15)                                            # optional - sage.rings.number_field
        sage: CC(zeta)                                                                  # optional - sage.rings.number_field
        0.913545457642601 + 0.406736643075800*I

    Doctests that make use of the Algebraic Field ``QQbar``, the Algebraic Real Field ``AA``,
    or the Universal Cyclotomic Field should be marked likewise::

        sage: AA(-1)^(1/3)                                                              # optional - sage.rings.number_field
        -1
        sage: QQbar(-1)^(1/3)                                                           # optional - sage.rings.number_field
        0.500000000000000? + 0.866025403784439?*I

        sage: UCF = UniversalCyclotomicField(); UCF                                     # optional - sage.rings.number_field
        Universal Cyclotomic Field
        sage: E = UCF.gen                                                               # optional - sage.rings.number_field
        sage: f = E(2) + E(3); f                                                        # optional - sage.rings.number_field
        2*E(3) + E(3)^2
        sage: f.galois_conjugates()                                                     # optional - sage.rings.number_field
        [2*E(3) + E(3)^2, E(3) + 2*E(3)^2]

    TESTS::

        sage: from sage.features.sagemath import sage__rings__number_field
        sage: sage__rings__number_field().is_present()                                  # optional - sage.rings.number_field
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
                             [PythonModule('sage.rings.number_field.number_field_element')])


class sage__rings__padics(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.rings.padics``.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__padics
        sage: sage__rings__padics().is_present()  # optional - sage.rings.padics
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
                             [PythonModule('sage.rings.padics.factory')])


class sage__rings__real_double(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.real_double`.

    EXAMPLES:

    The Real Double Field is basically always available, and no ``# optional`` annotation is needed::

        sage: RDF.characteristic()
        0

    The feature exists for use in doctests of Python modules that are shipped by the
    most fundamental distributions.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__real_double
        sage: sage__rings__real_double().is_present()  # optional - sage.rings.real_double
        FeatureTestResult('sage.rings.real_double', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__real_double
            sage: isinstance(sage__rings__real_double(), sage__rings__real_double)
            True
        """
        PythonModule.__init__(self, 'sage.rings.real_double')


class sage__rings__real_mpfr(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.real_mpfr`.

    TESTS::

        sage: from sage.features.sagemath import sage__rings__real_mpfr
        sage: sage__rings__real_mpfr().is_present()  # optional - sage.rings.real_mpfr
        FeatureTestResult('sage.rings.real_mpfr', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__real_mpfr
            sage: isinstance(sage__rings__real_mpfr(), sage__rings__real_mpfr)
            True
        """
        PythonModule.__init__(self, 'sage.rings.real_mpfr')


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
        sage: sage__symbolic().is_present()                                             # optional - sage.symbolic
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
                             [PythonModule('sage.symbolic.expression')])


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
    return [sagemath_doc_html(),
            sage__combinat(),
            sage__geometry__polyhedron(),
            sage__graphs(),
            sage__groups(),
            sage__libs__pari(),
            sage__modules(),
            sage__plot(),
            sage__rings__finite_rings(),
            sage__rings__function_field(),
            sage__rings__number_field(),
            sage__rings__padics(),
            sage__rings__real_double(),
            sage__rings__real_mpfr(),
            sage__symbolic()]
