# sage_setup: distribution = sagemath-polyhedra
# sage.doctest: needs sage.rings.real_double
"""
Base class for polyhedra over ``RDF``
"""

from sage.rings.real_double import RDF
from .base import Polyhedron_base


class Polyhedron_RDF(Polyhedron_base):
    """
    Base class for polyhedra over ``RDF``.

    TESTS::

        sage: p = Polyhedron([(0,0)], base_ring=RDF);  p
        A 0-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run()
    """
    # 1e-6 is the cddf+ default fuzzy zero cutoff

    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=RDF)
            sage: p._is_zero(0)
            True
            sage: p._is_zero(1e-3)
            False

        This is a fuzzy zero for floating-point numbers::

            sage: p._is_zero(1e-10)
            True
        """
        return abs(x) <= 1e-6

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=RDF)
            sage: p._is_nonneg(1)
            True
            sage: p._is_nonneg(-1e-3)
            False

        This is a fuzzy zero for floating-point numbers::

            sage: p._is_nonneg(-1e-10)
            True
        """
        return x >= -1e-6

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=RDF)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            True

        This is a fuzzy zero for floating-point numbers::

            sage: p._is_positive(-1e-10)
            True
        """
        return x >= -1e-6

    _base_ring = RDF
