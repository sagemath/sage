r"""
Algebra of differential forms

Algebra of differential forms defined on a CoordinatePatch (an open subset of
Euclidian space, see ``CoordinatePatch`` for details).

AUTHORS:

 - Joris Vankerschaver (2010-05-26)

.. TODO::

    - Allow for forms with values in a vector space

    - Incorporate Kahler differentials

REFERENCES:

- R. Abraham, J. E. Marsden, and T. S. Ratiu: Manifolds, tensor analysis,
  and applications.  Springer-Verlag 1988, texts in Applied Mathematical
  Sciences, volume 75, 2nd edition.

- :wikipedia:`Differential_form`

"""

#*****************************************************************************
#    Copyright (C) 2010 Joris Vankerschaver (joris.vankerschaver@gmail.com)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import Algebra
from sage.tensor.coordinate_patch import CoordinatePatch
from sage.tensor.differential_form_element import DifferentialForm
from sage.symbolic.ring import SR, var



class DifferentialForms(Algebra):
    """
    The algebra of all differential forms on an open subset of Euclidian space
    of arbitrary dimension.

    EXAMPLES:

    To define an algebra of differential forms, first create a coordinate
    patch::

        sage: p, q = var('p, q')
        sage: U = CoordinatePatch((p, q)); U
        doctest:...: DeprecationWarning: Use Manifold instead.
        See http://trac.sagemath.org/24444 for details.
        Open subset of R^2 with coordinates p, q
        sage: F = DifferentialForms(U); F
        doctest:...:  DeprecationWarning: For the set of differential forms of
         degree p, use U.diff_form_module(p), where U is the base manifold
         (type U.diff_form_module? for details).
        See http://trac.sagemath.org/24444 for details.
        Algebra of differential forms in the variables p, q

    If no coordinate patch is supplied, a default one (using the variables
    x, y, z) will be used::

        sage: F = DifferentialForms(); F
        doctest:...: DeprecationWarning: Use Manifold instead.
        See http://trac.sagemath.org/24444 for details.
        Algebra of differential forms in the variables x, y, z

    """

    Element = DifferentialForm

    def __init__(self, coordinate_patch = None):
        """
        Construct the algebra of differential forms on a given coordinate patch.

        See ``DifferentialForms`` for details.

        INPUT:

        - ``coordinate_patch`` -- Coordinate patch where the algebra lives.

        If no coordinate patch is given, a default coordinate patch with
        coordinates (x, y, z) is used.

        EXAMPLES::

            sage: p, q = var('p, q')
            sage: U = CoordinatePatch((p, q)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^2 with coordinates p, q
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables p, q
        """
        from sage.categories.graded_algebras_with_basis \
            import GradedAlgebrasWithBasis
        from sage.structure.parent_gens import ParentWithGens
        from sage.misc.superseded import deprecation
        deprecation(24444, 'For the set of differential forms of degree p, ' +
                    'use U.diff_form_module(p), where U is the base ' +
                    'manifold (type U.diff_form_module? for details).')

        if not coordinate_patch:
            x, y, z = var('x, y, z')
            coordinate_patch = CoordinatePatch((x, y, z))

        if not isinstance(coordinate_patch, CoordinatePatch):
            raise TypeError("%s not a valid Coordinate Patch" % coordinate_patch)
        self._patch = coordinate_patch

        ParentWithGens.__init__(self, SR, \
                                category = GradedAlgebrasWithBasis(SR))


    def __eq__(self, other):
        """
        Return True if self is equal to other.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: p, q = var('p, q')
            sage: V = CoordinatePatch((p, q)); V
            Open subset of R^2 with coordinates p, q
            sage: G = DifferentialForms(V); G
            Algebra of differential forms in the variables p, q
            sage: H = DifferentialForms(U); H
            Algebra of differential forms in the variables x, y, z
            sage: F == G
            False
            sage: F == H
            True
        """

        if type(other) is type(self):
            return self._patch == other._patch
        else:
            return False


    def __ne__(self, other):
        """
        Return True if self is not equal to other.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: p, q = var('p, q')
            sage: V = CoordinatePatch((p, q)); V
            Open subset of R^2 with coordinates p, q
            sage: G = DifferentialForms(V); G
            Algebra of differential forms in the variables p, q
            sage: F != G
            True
        """
        return not self == other

    def ngens(self):
        """
        Return the number of generators of this algebra.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F.ngens()
            3
        """
        return len(self._patch.coordinates())


    def gen(self, i=0):
        """
        Return the `i^{th}` generator of ``self``.  This is a one-form,
        more precisely the exterior derivative of the i-th coordinate.

        INPUT:

        - ``i`` - integer (optional, default 0)


        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F.gen(0)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            dx
            sage: F.gen(1)
            dy
            sage: F.gen(2)
            dz

        """

        form = DifferentialForm(self, 0, self._patch.coordinate(i))
        return form.diff()

    def gens(self):
        """
        Return a list of the generators of ``self``.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F.gens()
            (dx, dy, dz)
        """
        return tuple(self.gen(n) for n in range(self._patch.dim()))

    def base_space(self):
        """
        Return the coordinate patch on which this algebra is defined.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F.base_space()
            Open subset of R^3 with coordinates x, y, z
        """
        return self._patch


    def _element_constructor_(self, fun):
        """
        Coerce a given function (element of the symbolic ring)
        into a differential form of degree zero.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F(sin(x*y))    # indirect doctest
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sin(x*y)

        """


        fun = SR(fun)
        if fun not in self:
            raise ValueError("Function not an element of this algebra of differential forms.")

        return DifferentialForm(self, 0, fun)


    def __contains__(self, element):
        """
        Check if a given element belongs to this algebra of differential forms.

        EXAMPLES::

            sage: x, y, p, q = var('x, y, p, q')
            sage: U = CoordinatePatch((x, y)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^2 with coordinates x, y
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y
            sage: x in F
            True
            sage: sin(y) in F
            True
            sage: p in F
            False
            sage: cos(q) in F
            False
        """

        parent = None
        try:
            parent = element.parent()
        except AttributeError:
            pass

        if parent == self:
            return True

        if parent == SR:
            for coordinate in element.variables():
                if coordinate not in self._patch.coordinates():
                    return False
            return True

        return False


    def _coerce_map_from_(self, S):
        """
        Only the symbolic ring coerces into the algebra of differential forms.

        EXAMPLES::

            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F._coerce_map_from_(SR)
            True
            sage: F._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(CC)
            False
            sage: F._coerce_map_from_(RR)
            False

        """
        return S is SR or S is self


    def _repr_(self):
        r"""
        String representation of this algebra of differential forms.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: F._repr_()
            'Algebra of differential forms in the variables x, y, z'
        """

        return "Algebra of differential forms in the variables " + \
            ', '.join(str(var) for var in self._patch.coordinates())


    def _latex_(self):
        r"""
        Latex representation of this algebra of differential forms.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: F = DifferentialForms(U); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: latex(F)
            \Omega^\ast(\mathbb{\RR}^3)
            sage: latex(F) == F._latex_()
            True
        """
        return "\\Omega^\\ast(\\mathbb{\\RR}^%s)" % self._patch.dim()
