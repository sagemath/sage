r"""
Features for testing the presence of Python modules in the Sage library
"""

# *****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#                     2021 Kwankyu Lee
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

    EXAMPLES::

        sage: from sage.features.sagemath import sagemath_doc_html
        sage: sagemath_doc_html().is_present()  # optional - sagemath_doc_html
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

    EXAMPLES::

        sage: from sage.features.sagemath import sage__combinat
        sage: sage__combinat().is_present()  # optional - sage.combinat
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
                             [PythonModule('sage.combinat.tableau')],
                             spkg='sagemath_combinat', type="standard")


class sage__geometry__polyhedron(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.geometry.polyhedron`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__geometry__polyhedron
        sage: sage__geometry__polyhedron().is_present()  # optional - sage.geometry.polyhedron
        FeatureTestResult('sage.geometry.polyhedron', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__geometry__polyhedron
            sage: isinstance(sage__geometry__polyhedron(), sage__geometry__polyhedron)
            True
        """
        PythonModule.__init__(self, 'sage.geometry.polyhedron',
                              spkg='sagemath_polyhedra', type="standard")


class sage__graphs(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.graphs`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__graphs
        sage: sage__graphs().is_present()  # optional - sage.graphs
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
                             [PythonModule('sage.graphs.graph')],
                             spkg='sagemath_graphs', type="standard")


class sage__modular(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.modular`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__modular
        sage: sage__modular().is_present()                                              # optional - sage.modular
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


class sage__groups(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.groups``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__groups
        sage: sage__groups().is_present()  # optional - sage.groups
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


class sage__libs__cmr(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.cmr`
    and other modules depending on CMR, the Combinatorial Matrix Recognition library.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__libs__cmr
        sage: sage__libs__cmr().is_present()                                            # optional - sage.libs.cmr
        FeatureTestResult('sage.libs.cmr', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__libs__cmr
            sage: isinstance(sage__libs__cmr(), sage__libs__cmr)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.cmr',
                             [PythonModule('sage.matrix.matrix_cmr_sparse')],
                             spkg='sagemath_cmr')


class sage__libs__flint(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.flint`
    and other modules depending on FLINT and arb.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__libs__flint
        sage: sage__libs__flint().is_present()                                          # optional - sage.libs.flint
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
                             [PythonModule('sage.libs.flint.flint'),
                              PythonModule('sage.libs.arb.arith')],
                             spkg='sagemath_flint', type='standard')


class sage__libs__ntl(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.ntl`
    and other modules depending on NTL and arb.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__libs__ntl
        sage: sage__libs__ntl().is_present()                                            # optional - sage.libs.ntl
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


class sage__libs__pari(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.pari`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__libs__pari
        sage: sage__libs__pari().is_present()                       # optional - sage.libs.pari
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


class sage__modules(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.modules`.

    EXAMPLES::

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
                             [PythonModule('sage.modules.free_module')],
                             spkg='sagemath_modules', type='standard')


class sage__plot(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.plot`.

    EXAMPLES::

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
                             [PythonModule('sage.plot.plot')],
                             spkg='sagemath_symbolics', type='standard')


class sage__rings__finite_rings(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.finite_rings`;
    specifically, the element implementations using PARI.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__rings__finite_rings
        sage: sage__rings__finite_rings().is_present()  # optional - sage.rings.finite_rings
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
                             [PythonModule('sage.rings.finite_rings.element_pari_ffelt')],
                             type='standard')


class sage__rings__function_field(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.rings.function_field`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__rings__function_field
        sage: sage__rings__function_field().is_present()  # optional - sage.rings.function_field
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

    EXAMPLES::

        sage: from sage.features.sagemath import sage__rings__number_field
        sage: sage__rings__number_field().is_present()  # optional - sage.rings.number_field
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
                             [PythonModule('sage.rings.number_field.number_field_element')],
                             type='standard')


class sage__rings__padics(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.rings.padics``.

    EXAMPLES::

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
                             [PythonModule('sage.rings.padics.factory')],
                             type='standard')


class sage__rings__polynomial__pbori(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.rings.polynomial.pbori`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__rings__polynomial__pbori
        sage: sage__rings__polynomial__pbori().is_present()                             # optional - sage.rings.polynomial.pbori
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

    EXAMPLES::

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

    EXAMPLES::

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
        PythonModule.__init__(self, 'sage.rings.real_mpfr',
                              spkg='sagemath_modules')


class sage__schemes(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.schemes`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__schemes
        sage: sage__schemes().is_present()                                              # optional - sage.schemes
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
                             spkg="sagemath_schemes")


class sage__symbolic(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`sage.symbolic`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__symbolic
        sage: sage__symbolic().is_present()  # optional - sage.symbolic
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
                             [PythonModule('sage.symbolic.expression')],
                             spkg="sagemath_symbolics")


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
            sage__libs__cmr(),
            sage__libs__flint(),
            sage__libs__ntl(),
            sage__libs__pari(),
            sage__modular(),
            sage__modules(),
            sage__plot(),
            sage__rings__finite_rings(),
            sage__rings__function_field(),
            sage__rings__number_field(),
            sage__rings__padics(),
            sage__rings__polynomial__pbori(),
            sage__rings__real_double(),
            sage__rings__real_mpfr(),
            sage__schemes(),
            sage__symbolic()]
