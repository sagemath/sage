r"""
Strata of quadratic differentials on Riemann surfaces

.. WARNING::

    This module is deprecated. You are advised to install and use the
    surface_dynamics package instead available at
    https://pypi.python.org/pypi/surface_dynamics/
"""
from __future__ import print_function
from six import iteritems

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer


class QuadraticStratum(SageObject):
    r"""
    Stratum of quadratic differentials.
    """
    def __init__(self, *l):
        """
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            doctest:warning
            ...
            DeprecationWarning: QuadraticStratum is deprecated and will be removed from Sage.
            You are advised to install the surface_dynamics package via:
            sage -pip install surface_dynamics
            If you do not have write access to the Sage installation you can
            alternatively do
            sage -pip install surface_dynamics --user
            The package surface_dynamics subsumes all flat surface related
            computation that are currently available in Sage. See more
            information at
            http://www.labri.fr/perso/vdelecro/surface-dynamics/latest/
            See http://trac.sagemath.org/20695 for details.
            sage: loads(dumps(a)) == a
            True
            sage: QuadraticStratum([])
            Traceback (most recent call last):
            ...
            ValueError: the list must be non empty !
        """
        from sage.dynamics.surface_dynamics_deprecation import surface_dynamics_deprecation
        surface_dynamics_deprecation("QuadraticStratum")

        if isinstance(l[0], list) or isinstance(l[0], tuple):
            if not l[0]:
                raise ValueError("the list must be non empty !")
            self._zeroes = []
            for (i, j) in iteritems(l):
                i = Integer(i)
                j = Integer(j)
                self._zeroes += [i]*j
        else:
            for i in l:
                i = Integer(i)
            self._zeroes = sorted(list(l), reverse=True)


        self._genus = sum(l)/4 + 1
        self._genus = Integer(self._genus)

    def __repr__(self):
        r"""
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            sage: print(a)
            Q(-1, -1, -1, -1)
        """
        return "Q(" + str(self._zeroes)[1:-1] + ")"

    def __str__(self):
        r"""
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            sage: print(a)
            Q(-1, -1, -1, -1)
        """
        return "Q(" + str(self._zeroes)[1:-1] + ")"

    def __eq__(self, other):
        r"""
        TESTS::

            sage: QuadraticStratum(0) == QuadraticStratum(0)
            True
            sage: QuadraticStratum(4) == QuadraticStratum(0)
            False
        """
        return type(self) is type(other) and self._zeroes == other._zeroes

    def __ne__(self, other):
        r"""
        TESTS::

            sage: QuadraticStratum(0) != QuadraticStratum(0)
            False
            sage: QuadraticStratum(4) != QuadraticStratum(0)
            True
        """
        return type(self) is not type(other) or self._zeroes != other._zeroes

    def genus(self):
        r"""
        Returns the genus.

        EXAMPLES:

        ::

            sage: QuadraticStratum(-1,-1,-1,-1).genus()
            0
        """
        return self._genus
