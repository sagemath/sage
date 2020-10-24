r"""
Elements of the algebra of differential forms

AUTHORS:

- Joris Vankerschaver (2010-07-25)

"""

# ****************************************************************************
#  Copyright (C) 2010 Joris Vankerschaver <joris.vankerschaver@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.symbolic.ring import SR
from sage.structure.element import RingElement, AlgebraElement
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutation


def sort_subscript(subscript):
    """
    A subscript is a range of integers.  This function sorts a subscript
    in the sense of arranging it in ascending order.  The return values
    are the sign of the subscript and the sorted subscript, where the
    sign is defined as follows:

    #. sign == 0 if two or more entries in the subscript were equal.
    #. sign == +1, -1 if a positive (resp. negative) permutation was used to sort the subscript.


    INPUT:

    - ``subscript`` -- a subscript, i.e. a range of not necessarily
      distinct integers


    OUTPUT:

    - Sign of the permutation used to arrange the subscript, where 0
        means that the original subscript had two or more entries that
        were the same

    - Sorted subscript.


    EXAMPLES::

        sage: from sage.tensor.differential_form_element import sort_subscript
        sage: sort_subscript((1, 3, 2))
        (-1, (1, 2, 3))
        sage: sort_subscript((1, 3))
        (1, (1, 3))
        sage: sort_subscript((4, 2, 7, 9, 8))
        (1, (2, 4, 7, 8, 9))

    """
    if len(subscript) == 0:
        return 1, ()

    sub_list = sorted(subscript)
    offsets = [subscript.index(x)+1 for x in sub_list]

    # Check that offsets is a true permutation of 1..n
    n = len(offsets)
    if sum(offsets) != n*(n+1)//2:
        sign = 0
    else:
        sign = Permutation(offsets).signature()

    return sign, tuple(sub_list)



class DifferentialFormFormatter:
    r"""
    This class contains all the functionality to print a differential form in a
    graphically pleasing way.  This class is called by the ``_latex_`` and
    ``_repr_`` methods of the DifferentialForm class.

    In a nutshell (see the documentation of ``DifferentialForm`` for more
    details), differential forms are represented internally as a dictionary,
    where the keys are tuples representing the non-zero components of the
    form and the values are the component functions.  The methods of this
    class create string and latex representations out of the specification
    of a subscript and a component function.


    EXAMPLES::

        sage: from sage.tensor.differential_form_element import DifferentialFormFormatter
        sage: x, y, z = var('x, y, z')
        sage: U = CoordinatePatch((x, y, z))
        doctest:...: DeprecationWarning: Use Manifold instead.
        See http://trac.sagemath.org/24444 for details.
        sage: D = DifferentialFormFormatter(U)
        sage: D.repr((0, 2), sin(x*y))
        'sin(x*y)*dx/\\dz'
        sage: D.latex((0, 2), sin(x*y))
        '\\sin\\left(x y\\right) d x \\wedge d z'
        sage: D.latex((1, 2), exp(z))
        'e^{z} d y \\wedge d z'

    """
    def __init__(self, space):
        r"""
        Construct a differential form formatter.  See
        ``DifferentialFormFormatter`` for more information.

        INPUT:

        - space -- CoordinatePatch where the differential forms live.

        EXAMPLES::

            sage: from sage.tensor.differential_form_element import DifferentialFormFormatter
            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: D = DifferentialFormFormatter(U)
            sage: D.repr((0, 2), sin(x*y))
            'sin(x*y)*dx/\\dz'

        """
        self._space = space


    def repr(self, comp, fun):
        r"""
        String representation of a primitive differential form, i.e. a function
        times a wedge product of d's of the coordinate functions.

        INPUT:

        - ``comp`` -- a subscript of a differential form.

        - ``fun`` -- the component function of this form.


        EXAMPLES::

            sage: from sage.tensor.differential_form_element import DifferentialFormFormatter
            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: D = DifferentialFormFormatter(U)
            sage: D.repr((0, 1), z^3)
            'z^3*dx/\\dy'

        """

        str = "/\\".join( \
            [('d%r' % self._space.coordinate(c)) for c in comp])

        if fun == 1 and len(comp) > 0:
            # We have a non-trivial form whose component function is 1,
            # so we just return the formatted form part and ignore the 1.
            return str
        else:
            funstr = fun._repr_()

            if not self._is_atomic(funstr):
                funstr = '(' + funstr + ')'

            if len(str) > 0:
                return funstr + "*" + str
            else:
                return funstr


    def latex(self, comp, fun):
        r"""
        Latex representation of a primitive differential form, i.e. a function
        times a wedge product of d's of the coordinate functions.

        INPUT:

        - ``comp`` -- a subscript of a differential form.

        - ``fun`` -- the component function of this form.

        EXAMPLES::

            sage: from sage.tensor.differential_form_element import DifferentialFormFormatter
            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: D = DifferentialFormFormatter(U)
            sage: D.latex((0, 1), z^3)
            'z^{3} d x \\wedge d y'
            sage: D.latex((), 1)
            '1'
            sage: D.latex((), z^3)
            'z^{3}'
            sage: D.latex((0,), 1)
            'd x'
        """

        from sage.misc.latex import latex

        s = " \\wedge ".join( \
            [('d %s' % latex(self._space.coordinate(c))) for c in comp])

        # Make sure this is a string and not a LatexExpr
        s = str(s)

        # Add coefficient except if it's 1
        if fun == 1:
            if s:
                return s
            else:
                return "1"

        funstr = fun._latex_()
        if not self._is_atomic(funstr):
            funstr = '(' + funstr + ')'

        if s:
            s = " " + s
        return funstr + s


    def _is_atomic(self, str):
        r"""
        Helper function to check whether a given string expression
        is atomic.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: U = CoordinatePatch((x, y, z))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: from sage.tensor.differential_form_element import DifferentialFormFormatter
            sage: D = DifferentialFormFormatter(U)
            sage: D._is_atomic('a + b')
            False
            sage: D._is_atomic('(a + b)')
            True
        """
        level = 0
        for n, c in enumerate(str):
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1

            if c == '+' or c == '-':
                if level == 0 and n > 0:
                    return False
        return True



class DifferentialForm(AlgebraElement):
    r"""
    Differential form class.

    EXAMPLES:

    In order to instantiate differential forms of various degree, we begin
    by specifying the CoordinatePatch on which they live, as well as their
    parent DifferentialForms algebra.

    ::

        sage: x, y, z = var('x, y, z')
        sage: U = CoordinatePatch((x, y, z))
        doctest:...: DeprecationWarning: Use Manifold instead.
        See http://trac.sagemath.org/24444 for details.
        sage: F = DifferentialForms(U)
        doctest:...:  DeprecationWarning: For the set of differential forms of
         degree p, use U.diff_form_module(p), where U is the base manifold
         (type U.diff_form_module? for details).
        See http://trac.sagemath.org/24444 for details.
        sage: form1 = DifferentialForm(F, 0, sin(x*y)); form1
        doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
         where U is the base manifold (type U.diff_form? for details).
        See http://trac.sagemath.org/24444 for details.
        sin(x*y)

    In the previous example, we created a zero-form from a given function.
    To create forms of higher degree, we can use the subscript operator
    access the various components::

        sage: form2 = DifferentialForm(F, 1); form2
        0
        sage: form2[0] = 1
        sage: form2[1] = exp(cos(x))
        sage: form2[2] = 1/ln(y)
        sage: form2
        dx + e^cos(x)*dy + 1/log(y)*dz

    We may calculate the exterior derivative of a form, and observe that
    applying the exterior derivative twice always yields zero::

        sage: dform = form1.diff(); dform
        doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
         where U is the base manifold (type U.diff_form? for details).
        See http://trac.sagemath.org/24444 for details.
        y*cos(x*y)*dx + x*cos(x*y)*dy
        sage: dform.diff()
        0

    As can be seen from the previous example, the exterior derivative
    increases the degree of a form by one::

        sage: form2.degree()
        1
        sage: form2.diff().degree()
        2

    The ``d`` function provides a convenient shorthand for applying the
    diff member function.  Since d appears in other areas of mathematics
    as well, this function is not imported in the global namespace
    automatically::

        sage: from sage.tensor.differential_form_element import d
        sage: form2
        dx + e^cos(x)*dy + 1/log(y)*dz
        sage: d(form2)
        -e^cos(x)*sin(x)*dx/\dy + -1/(y*log(y)^2)*dy/\dz
        sage: form2.diff()
        -e^cos(x)*sin(x)*dx/\dy + -1/(y*log(y)^2)*dy/\dz
        sage: d(form1) == form1.diff()
        True

    The wedge product of two forms can be computed by means of the wedge
    member function::

        sage: form1 = DifferentialForm(F, 2)
        sage: form1[0, 1] = exp(z); form1
        e^z*dx/\dy
        sage: form2 = DifferentialForm(F, 1)
        sage: form2[2] = exp(-z)
        sage: form1.wedge(form2)
        doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
         where U is the base manifold (type U.diff_form? for details).
        See http://trac.sagemath.org/24444 for details.
        dx/\dy/\dz

    For this member function, there exists again a procedural function
    which is completely equivalent::

        sage: from sage.tensor.differential_form_element import wedge
        sage: form1.wedge(form2)
        dx/\dy/\dz
        sage: wedge(form1, form2)
        dx/\dy/\dz
        sage: form1.wedge(form2) == wedge(form1, form2)
        True


    .. NOTE::

        Differential forms are stored behind the screens as dictionaries,
        where the keys are the subscripts of the non-zero components, and
        the values are those components.

        For example, on a
        space with coordinates x, y, z, the form

            f = sin(x*y) dx /\\ dy + exp(z) dy /\\ dz

        would be represented as the dictionary

            {(0, 1): sin(x*y), (1, 2): exp(z)}.

        Most differential forms are ''sparse'' in the sense that most of
        their components are zero, so that this representation is more
        efficient than storing all of the components in a vector.

    """

    def __init__(self, parent, degree, fun = None):
        r"""
        Construct a differential form.

        INPUT:

        - ``parent`` -- Parent algebra of differential forms.

        - ``degree`` -- Degree of the differential form.

        - ``fun`` (default: None) -- Initialize this differential form with the given function.  If the degree is not zero, this argument is silently ignored.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: f = DifferentialForm(F, 0, sin(z)); f
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sin(z)

        """

        from sage.tensor.differential_forms import DifferentialForms
        from sage.misc.superseded import deprecation
        deprecation(24444, 'Use U.diff_form(degree) instead, where U is ' +
                    'the base manifold (type U.diff_form? for details).')
        if not isinstance(parent, DifferentialForms):
            raise TypeError("Parent not an algebra of differential forms.")

        RingElement.__init__(self, parent)

        self._degree = degree
        self._components = {}

        if degree == 0 and fun is not None:
            self[[]] = fun

    def __getitem__(self, subscript):
        r"""
        Return a given component of the differential form.

        INPUT:

        - ``subscript`` -- subscript of the component.  Must be an integer
          or a list of integers.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: f = DifferentialForm(F, 0, sin(x*y)); f
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sin(x*y)
            sage: f[()]
            sin(x*y)

            sage: df = f.diff(); df
            y*cos(x*y)*dx + x*cos(x*y)*dy
            sage: df[0]
            y*cos(x*y)
            sage: df[1]
            x*cos(x*y)
            sage: df[2]
            0
        """
        if isinstance(subscript, (Integer, int)):
            subscript = (subscript, )
        else:
            subscript = tuple(subscript)

        dim = self.parent().base_space().dim()
        if any(s >= dim for s in subscript):
            raise ValueError("Index out of bounds.")

        if len(subscript) != self._degree:
            raise TypeError("%s is not a subscript of degree %s" %
                            (subscript, self._degree))

        sign, subscript = sort_subscript(subscript)

        if subscript in self._components:
            return sign*self._components[subscript]
        else:
            return 0

    def __setitem__(self, subscript, fun):
        r"""
        Modify a given component of the differential form.

        INPUT:

        - ``subscript`` -- subscript of the component.  Must be an integer
          or a list of integers.

        EXAMPLES::

            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: f = DifferentialForm(F, 2)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1, 2] = x; f
            x*dy/\dz
        """

        if isinstance(subscript, (Integer, int)):
            subscript = (subscript, )
        else:
            subscript = tuple(subscript)

        dim = self.parent().base_space().dim()
        if any(s >= dim for s in subscript):
            raise ValueError("Index out of bounds.")

        if len(subscript) != self._degree:
            raise TypeError("%s is not a subscript of degree %s" %
                            (subscript, self._degree))

        sign, subscript = sort_subscript(subscript)
        self._components[subscript] = sign*SR(fun)


    def is_zero(self):
        r"""
        Return True if ``self`` is the zero form.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1); f
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            0
            sage: f.is_zero()
            True
            sage: f[1] = 1
            sage: f.is_zero()
            False
            sage: f.diff()
            0
            sage: f.diff().is_zero()
            True
        """

        self._cleanup()
        return len(self._components) == 0


    def degree(self):
        r"""
        Return the degree of self.

        EXAMPLES::

            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: f = DifferentialForm(F, 2)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1, 2] = x; f
            x*dy/\dz
            sage: f.degree()
            2

        The exterior differential increases the degree of forms by one::

            sage: g = f.diff(); g
            dx/\dy/\dz
            sage: g.degree()
            3
        """
        return self._degree


    def __eq__(self, other):
        r"""
        Test whether two differential forms are equal.

        EXAMPLES::

            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: f = DifferentialForm(F, 2)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1,2] = x; f
            x*dy/\dz
            sage: f == f
            True

            sage: g = DifferentialForm(F, 3)
            sage: g[0, 1, 2] = 1; g
            dx/\dy/\dz
            sage: f == g
            False
            sage: f.diff() == g
            True
        """
        if type(other) is type(self):
            if self._degree != other._degree:
                return False
            else:
                # TODO: the following two lines are where most of the
                # execution time is spent.
                self._cleanup()
                other._cleanup()

                if len(self._components) != len(other._components):
                    return False


                # We compare the component dictionary of both differential
                # forms, keeping in mind that the set of keys is
                # lexicographically ordered, so that we can simply iterate
                # over both dictionaries in one go and compare (key, value)
                # pairs as we go along.

                for (key1, val1), (key2, val2) in zip(self._components.items(),
                                                      other._components.items()):
                    if key1 != key2 or str(val1) != str(val2):
                        return False
                return True
        else:
            return False


    def __ne__(self, other):
        r"""
        Test whether two differential forms are not equal.

        EXAMPLES::

            sage: F = DifferentialForms(); F
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            Algebra of differential forms in the variables x, y, z
            sage: f = DifferentialForm(F, 2)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1,2] = x; f
            x*dy/\dz
            sage: g = DifferentialForm(F, 3)
            sage: g[0, 1, 2] = 1; g
            dx/\dy/\dz
            sage: f != g
            True
        """
        return not self == other

    def _neg_(self):
        r"""
        Return the negative of self.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[0] = y
            sage: f[1] = -x
            sage: f
            y*dx + -x*dy
            sage: -f
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            -y*dx + x*dy
            sage: -f == f._neg_()
            True
        """

        neg = DifferentialForm(self.parent(), self._degree)
        for comp in self._components:
            neg._components[comp] = -self._components[comp]
        return neg


    def _add_(self, other):
        r"""
        Add self and other

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: g = DifferentialForm(F, 1)
            sage: f[0] = exp(x); f
            e^x*dx
            sage: g[1] = sin(y); g
            sin(y)*dy
            sage: f + g
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            e^x*dx + sin(y)*dy
            sage: f + g == f._add_(g)
            True

        Forms must have the same degree to be added::

            sage: h = DifferentialForm(F, 2)
            sage: h[1, 2] = x; h
            x*dy/\dz
            sage: f + h
            Traceback (most recent call last):
            ...
            TypeError: Cannot add forms of degree 1 and 2

        Subtraction is implemented by adding the negative::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            sage: f = DifferentialForm(F, 1)
            sage: g = DifferentialForm(F, 1)
            sage: f[0] = exp(x); f
            e^x*dx
            sage: g[1] = sin(y); g
            sin(y)*dy
            sage: f - g
            e^x*dx + -sin(y)*dy
            sage: f - g == f._sub_(g)
            True

        Forms must have the same degree to be subtracted::

            sage: h = DifferentialForm(F, 2)
            sage: h[1, 2] = x; h
            x*dy/\dz
            sage: f - h
            Traceback (most recent call last):
            ...
            TypeError: Cannot add forms of degree 1 and 2
        """

        if self.is_zero():
            return other
        if other.is_zero():
            return self

        if self._degree != other._degree:
            raise TypeError("Cannot add forms of degree %s and %s" % \
                    (self._degree, other._degree))

        sumform = DifferentialForm(self.parent(), self._degree)
        sumform._components = self._components.copy()
        for comp, fun in other._components.items():
            sumform[comp] += fun

        sumform._cleanup()
        return sumform

    def _cleanup(self):
        r"""
        Helper function to clean up self, i.e. to remove any
        zero components from the dictionary of components.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[0] = 0
            sage: f[1] = 1
            sage: f[2] = 0
            sage: f._dump_all()
            {(0,): 0, (1,): 1, (2,): 0}
            sage: f._cleanup()
            sage: f._dump_all()
            {(1,): 1}

        """
        zeros = []

        for comp in self._components:
            if self._components[comp].is_zero():
                zeros.append(comp)

        for comp in zeros:
            del self._components[comp]


    def _dump_all(self):
        r"""
        Helper function to dump the internal dictionary of form components.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1] = exp(cos(x))
            sage: f[2] = sin(ln(y))
            sage: f
            e^cos(x)*dy + sin(log(y))*dz
            sage: f._dump_all()
            {(1,): e^cos(x), (2,): sin(log(y))}
            sage: g = DifferentialForm(F, 2)
            sage: g[1, 2] = x+y+z
            sage: g
            (x + y + z)*dy/\dz
            sage: g._dump_all()
            {(1, 2): x + y + z}

        """
        from pprint import pprint
        pprint(self._components)

    def diff(self):
        r"""
        Compute the exterior differential of ``self``.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 0, sin(x*y)); f
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sin(x*y)
            sage: f.diff()
            y*cos(x*y)*dx + x*cos(x*y)*dy
            sage: g = DifferentialForm(F, 1)
            sage: g[0] = y/2
            sage: g[1] = -x/2
            sage: g
            1/2*y*dx + -1/2*x*dy
            sage: g.diff()
            -1*dx/\dy
            sage: h = DifferentialForm(F, 2)
            sage: h[0, 1] = exp(z)
            sage: h.diff()
            e^z*dx/\dy/\dz

        The square of the exterior differential operator is
        identically zero::

            sage: f
            sin(x*y)
            sage: f.diff()
            y*cos(x*y)*dx + x*cos(x*y)*dy
            sage: f.diff().diff()
            0

            sage: g.diff().diff()
            0

        The exterior differential operator is a derivation of degree one
        on the space of differential forms.  In this example we import the
        operator d() as a short-hand for having to call the diff()
        member function.

        ::

            sage: from sage.tensor.differential_form_element import d
            sage: d(f)
            y*cos(x*y)*dx + x*cos(x*y)*dy

            sage: d(f).wedge(g) + f.wedge(d(g))
            (-x*y*cos(x*y) - sin(x*y))*dx/\dy
            sage: d(f.wedge(g))
            (-x*y*cos(x*y) - sin(x*y))*dx/\dy

            sage: d(f.wedge(g)) == d(f).wedge(g) + f.wedge(d(g))
            True
        """

        diff_form = DifferentialForm(self.parent(), self._degree + 1)

        for comp in self._components:
            fun = self._components[comp]
            for n, coord in enumerate(self.parent().base_space().coordinates()):
                diff_form[(n, ) + comp] += fun.differentiate(coord)

        diff_form._cleanup()
        return diff_form


    def derivative(self, *args, **kwargs):
        r"""
        Compute the exterior derivative of ``self``.  This is the same as
        calling the ``diff`` member function.


        EXAMPLES::

            sage: x, y = var('x, y')
            sage: U = CoordinatePatch((x, y))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: F = DifferentialForms(U)
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: q = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: q[0] = -y/2
            sage: q[1] =  x/2
            sage: q.diff()
            dx/\dy
            sage: q.derivative()
            dx/\dy

        Invoking ``diff`` on a differential form has the same effect as
        calling this member function::

            sage: diff(q)
            dx/\dy
            sage: diff(q) == q.derivative()
            True

        When additional arguments are supplied to ``diff``, an error is raised,
        since only the exterior derivative has intrinsic meaning while
        derivatives with respect to the coordinate variables (in whichever
        way) are coordinate dependent, and hence not intrinsic.

        ::

            sage: diff(q, x)
            Traceback (most recent call last):
            ...
            ValueError: Differentiation of a form does not take any arguments.
        """

        if len(args) > 0 or len(kwargs) > 0:
            raise ValueError("Differentiation of a form does not take any arguments.")
        return self.diff()


    def wedge(self, other):
        r"""
        Returns the wedge product of ``self`` and other.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[0] = x^2
            sage: f[1] = y
            sage: f
            x^2*dx + y*dy
            sage: g = DifferentialForm(F, 1)
            sage: g[2] = z^3
            sage: g
            z^3*dz
            sage: f.wedge(g)
            x^2*z^3*dx/\dz + y*z^3*dy/\dz

        The wedge product is graded commutative::

            sage: f.wedge(g)
            x^2*z^3*dx/\dz + y*z^3*dy/\dz
            sage: g.wedge(f)
            -x^2*z^3*dx/\dz + -y*z^3*dy/\dz
            sage: f.wedge(f)
            0

        When the wedge product of forms belonging to different algebras
        is computed, an error is raised::

            sage: x, y, p, q = var('x, y, p, q')
            sage: F = DifferentialForms(CoordinatePatch((x, y)))
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            sage: G = DifferentialForms(CoordinatePatch((p, q)))
            sage: f = DifferentialForm(F, 0, 1); f
            1
            sage: g = DifferentialForm(G, 0, x); g
            x
            sage: f.parent()
            Algebra of differential forms in the variables x, y
            sage: g.parent()
            Algebra of differential forms in the variables p, q
            sage: f.wedge(g)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parents for wedge: 'Algebra of differential forms in the variables x, y' and  'Algebra of differential forms in the variables p, q'

        """

        if self.parent() != other.parent():
            raise TypeError("unsupported operand parents for wedge: " +\
                "\'%s\' and  \'%s\'" % (self.parent(), other.parent()))

        output = DifferentialForm(self.parent(), self._degree + other._degree)
        if self._degree + other._degree > self.parent().ngens():
            return output

        for lcomp, lfun in self._components.items():
            for rcomp, rfun in other._components.items():
                output[lcomp + rcomp] += lfun*rfun

        output._cleanup()
        return output



    def _mul_(self, other):
        r"""
        Multiply self and other.  This is identical to the wedge operator.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = F.gen(0); f
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            dx
            sage: g = F.gen(1); g
            dy
            sage: f*g
            dx/\dy
            sage: f.wedge(g)
            dx/\dy
            sage: f*g == f.wedge(g)
            True
            sage: f*g == f._mul_(g)
            True

        """
        return self.wedge(other)


    def _latex_(self):
        r"""
        Return a latex representation of self.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1] = exp(z); f
            e^z*dy
            sage: latex(f)
            e^{z} d y
            sage: g = f.diff(); g
            -e^z*dy/\dz
            sage: latex(g)
            -e^{z} d y \wedge d z
            sage: latex(g) == g._latex_()
            True

        """
        if len(self._components) == 0:
            return '0'

        format = DifferentialFormFormatter(self.parent().base_space())
        output = [format.latex(comp, fun) \
                      for (comp, fun) in sorted(self._components.items())]
        return ' + '.join(output)


    def _repr_(self):
        r"""
        Return string representation of self.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f[1] = exp(z); f
            e^z*dy
            sage: print(f)
            e^z*dy
            sage: f._repr_()
            'e^z*dy'
        """
        if len(self._components) == 0:
            return '0'

        format = DifferentialFormFormatter(self.parent().base_space())
        output = [format.repr(comp, fun) \
                      for (comp, fun) in sorted(self._components.items())]
        return ' + '.join(output)



    # Unsupported methods

    def abs(self):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.abs()
            Traceback (most recent call last):
            ...
            NotImplementedError: Absolute value not defined for differential forms.

        """
        raise NotImplementedError("Absolute value not defined for differential forms.")


    def leading_coefficient(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.leading_coefficient()
            Traceback (most recent call last):
            ...
            NotImplementedError: leading_coefficient not defined for differential forms.

        """
        raise NotImplementedError("leading_coefficient not defined for differential forms.")


    def leading_item(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.leading_item()
            Traceback (most recent call last):
            ...
            NotImplementedError: leading_item not defined for differential forms.

        """
        raise NotImplementedError("leading_item not defined for differential forms.")


    def leading_monomial(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.leading_monomial()
            Traceback (most recent call last):
            ...
            NotImplementedError: leading_monomial not defined for differential forms.

        """
        raise NotImplementedError("leading_monomial not defined for differential forms.")


    def leading_support(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.leading_support()
            Traceback (most recent call last):
            ...
            NotImplementedError: leading_support not defined for differential forms.

        """
        raise NotImplementedError("leading_support not defined for differential forms.")


    def leading_term(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.leading_term()
            Traceback (most recent call last):
            ...
            NotImplementedError: leading_term not defined for differential forms.

        """
        raise NotImplementedError("leading_term not defined for differential forms.")


    def trailing_coefficient(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.trailing_coefficient()
            Traceback (most recent call last):
            ...
            NotImplementedError: trailing_coefficient not defined for differential forms.

        """
        raise NotImplementedError("trailing_coefficient not defined for differential forms.")


    def trailing_item(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.trailing_item()
            Traceback (most recent call last):
            ...
            NotImplementedError: leading_coefficient not defined for differential forms.

        """
        raise NotImplementedError("leading_coefficient not defined for differential forms.")


    def trailing_monomial(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.trailing_monomial()
            Traceback (most recent call last):
            ...
            NotImplementedError: trailing_monomial not defined for differential forms.

        """
        raise NotImplementedError("trailing_monomial not defined for differential forms.")


    def trailing_support(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.trailing_support()
            Traceback (most recent call last):
            ...
            NotImplementedError: trailing_support not defined for differential forms.

        """
        raise NotImplementedError("trailing_support not defined for differential forms.")


    def trailing_term(self, key=None):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.trailing_term()
            Traceback (most recent call last):
            ...
            NotImplementedError: trailing_term not defined for differential forms.

        """
        raise NotImplementedError("trailing_term not defined for differential forms.")


    def map_coefficients(self, f):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.map_coefficients(lambda x: x)
            Traceback (most recent call last):
            ...
            NotImplementedError: map_coefficients not defined for differential forms.

        """
        raise NotImplementedError("map_coefficients not defined for differential forms.")


    def map_item(self, f):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.map_item(lambda x: x)
            Traceback (most recent call last):
            ...
            NotImplementedError: map_item not defined for differential forms.

        """
        raise NotImplementedError("map_item not defined for differential forms.")


    def map_support(self, f):
        """
        Method not defined for differential forms.

        EXAMPLES::

            sage: F = DifferentialForms()
            doctest:...:  DeprecationWarning: For the set of differential forms of
             degree p, use U.diff_form_module(p), where U is the base manifold
             (type U.diff_form_module? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f = DifferentialForm(F, 1)
            doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
             where U is the base manifold (type U.diff_form? for details).
            See http://trac.sagemath.org/24444 for details.
            sage: f.map_support(lambda x: x)
            Traceback (most recent call last):
            ...
            NotImplementedError: map_support not defined for differential forms.

        """
        raise NotImplementedError("map_support not defined for differential forms.")




def d(form):
    r"""
    Returns the exterior derivative of a given form, i.e. calls the diff()
    member function.

    EXAMPLES::

        sage: from sage.tensor.differential_form_element import d
        sage: x, y, z = var('x, y, z')
        sage: F = DifferentialForms()
        doctest:...:  DeprecationWarning: For the set of differential forms of
         degree p, use U.diff_form_module(p), where U is the base manifold
         (type U.diff_form_module? for details).
        See http://trac.sagemath.org/24444 for details.
        sage: f = DifferentialForm(F, 1)
        doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
         where U is the base manifold (type U.diff_form? for details).
        See http://trac.sagemath.org/24444 for details.
        sage: f[2] = cos(x); f
        cos(x)*dz
        sage: d(f)
        -sin(x)*dx/\dz
        sage: f.diff()
        -sin(x)*dx/\dz
        sage: d(f) == f.diff()
        True
    """
    return form.diff()


def wedge(left, right):
    r"""
    Computes the wedge product of two forms, i.e. calls the wedge()
    member function.

    EXAMPLES::

        sage: from sage.tensor.differential_form_element import wedge
        sage: x, y, z = var('x, y, z')
        sage: F = DifferentialForms()
        doctest:...:  DeprecationWarning: For the set of differential forms of
         degree p, use U.diff_form_module(p), where U is the base manifold
         (type U.diff_form_module? for details).
        See http://trac.sagemath.org/24444 for details.
        sage: f = DifferentialForm(F, 1)
        doctest:...: DeprecationWarning: Use U.diff_form(degree) instead,
         where U is the base manifold (type U.diff_form? for details).
        See http://trac.sagemath.org/24444 for details.
        sage: f[2] = cos(x); f
        cos(x)*dz
        sage: g = DifferentialForm(F, 1)
        sage: g[1] = sin(y); g
        sin(y)*dy
        sage: wedge(f, g)
        -cos(x)*sin(y)*dy/\dz
        sage: f.wedge(g)
        -cos(x)*sin(y)*dy/\dz
        sage: wedge(f, g) == f.wedge(g)
        True
    """
    return left.wedge(right)
