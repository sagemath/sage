r"""
Calculus method control

The class :class:`CalculusMethod` governs the calculus methods (symbolic and
numerical) to be used for coordinate computations on manifolds.

AUTHORS:

- Marco Mancini, (2017): initial version

"""

#******************************************************************************
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR
from sage.manifolds.utilities import (simplify_chain_real,
                                      simplify_chain_generic,
                                      simplify_chain_real_sympy,
                                      simplify_chain_generic_sympy,)
from sage.misc.latex import latex
import sympy

# Conversion functions
def _SR_to_Sympy(expression):
    r"""
    Convert an expression from ``SR`` to ``sympy``.

    In the case that the expression is already a ``sympy`` expression, then it
    is returned

    INPUT:

    - ``expression`` -- ``SR`` or ``sympy`` symbolic expression

    OUTPUT:

    - ``expression`` -- ``sympy`` symbolic expression

    EXAMPLES::

        sage: from sage.manifolds.calculus_method import _SR_to_Sympy
        sage: a = x^2+sin(x)^2; a
        x^2 + sin(x)^2
        sage: type(a)
        <type 'sage.symbolic.expression.Expression'>
        sage: b = _SR_to_Sympy(a); b
        x**2 + sin(x)**2
        sage: type(b)
        <class 'sympy.core.add.Add'>

    The conversion of a ``sympy`` expression is not converted::

        sage: _SR_to_Sympy(b) == b
        True

    """
    # Nothing to do if expression is already a SymPy object:
    if type(expression) in sympy.core.all_classes:
        return expression
    return SR(expression)._sympy_()

def _Sympy_to_SR(expression):
    r"""
    Convert an expression from ``sympy`` to ``SR``.

    In the case that the expression is already a ``SR`` expression, then it
    is returned.

    INPUT:

    - ``expression`` -- ``sympy`` symbolic expression

    OUTPUT:

    - ``expression`` -- SR or ``sympy`` symbolic expression


    EXAMPLES::

        sage: from sage.manifolds.calculus_method import _Sympy_to_SR, _SR_to_Sympy
        sage: a = x^2+sin(x)^2
        sage: type(a)
        <type 'sage.symbolic.expression.Expression'>
        sage: b = _SR_to_Sympy(a); b
        x**2 + sin(x)**2
        sage: type(b)
        <class 'sympy.core.add.Add'>
        sage: _Sympy_to_SR(b) == a
        x^2 + sin(x)^2 == x^2 + sin(x)^2
        sage: bool(_Sympy_to_SR(b) == a)
        True

    """
    try:
        return SR(expression)
    except:
        # If SR cannot transform a sympy expression is because it is a
        # sympy abstract function.
        a = expression._sage_()
        # As all sage objects have a ._sage_ operator, they have to be
        # catched
        if type(a) == type(expression):
            raise TypeError
        return a

class CalculusMethod(SageObject):
    r"""
    Control of calculus methods used on coordinate charts of manifolds.

    This class stores the possible calculus methods and permits to
    select some basic operations working on them.
    For the moment, only two calculus methods are implemented:

    - Sage's symbolic engine (Pynac/Maxima), implemented via the
      symbolic ring ``SR``
    - SymPy engine, denoted ``sympy`` hereafter

    INPUT:

    - ``current`` -- (default: ``None``) current symbolic method

    - ``base_field_type`` -- (default: ``'real'``) base field type of the
      manifold (cf.
      :meth:`~sage.manifolds.manifold.TopologicalManifold.base_field_type`)

    EXAMPLES::

        sage: from sage.manifolds.calculus_method import CalculusMethod
        sage: calc_meth = CalculusMethod()

    In the display, the current method is pointed out by `*`::

        sage: calc_meth
        Possible calculus methods:
         - SR (*) (default)
         - sympy

    The current method is changed by :meth:`set`::

        sage: calc_meth.set('sympy')
        sage: calc_meth
        Possible calculus methods:
         - SR (default)
         - sympy (*)
        sage: calc_meth.reset()
        sage: calc_meth
        Possible calculus methods:
         - SR (*) (default)
         - sympy

    """
    _default = 'SR'  # default symbolic method
    _methods = ('SR','sympy') # implemented methods
    _tranf = {'SR':  _Sympy_to_SR,'sympy' : _SR_to_Sympy} # translators

    def __init__(self, current=None, base_field_type='real'):
        r"""
        Initializes ``self``.

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='complex')
            sage: calc_meth
            Possible calculus methods:
             - SR (*) (default)
             - sympy

        """
        self._current = self._default if current is None else current
        self._simplify_dict = {}
        # initialization of the dict for simplifications
        if base_field_type == 'real':
            self._simplify_dict['sympy'] = simplify_chain_real_sympy
            self._simplify_dict['SR'] = simplify_chain_real
        else:
            self._simplify_dict['sympy'] = simplify_chain_generic_sympy
            self._simplify_dict['SR'] = simplify_chain_generic
        self._base_field_type = base_field_type
        self._latex_dict = {'sympy': sympy.latex, 'SR': latex}

    def simplify(self, expression, method=None):
        r"""
        Apply some simplification chain to a given symbolic expression.

        INPUT:

        - ``expression`` -- symbolic expression to simplify

        - ``method`` -- (default: ``None``) string defining the calculus method
          to use; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the current calculus method of ``self`` is used.

        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y')
            sage: f = x^2+sin(x)^2+cos(x)^2
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='real')
            sage: calc_meth.simplify(f)
            x^2 + 1

        Methods cannot be mixed::

            sage: calc_meth.set('sympy')
            sage: calc_meth.simplify(f)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.symbolic.expression.Expression' object has no attribute 'combsimp'

        In the present case, one should either transform ``f`` to a SymPy
        object::

            sage: calc_meth.simplify(f._sympy_())
            x**2 + 1

        or force the calculus method to be ``'SR'``::

            sage: calc_meth.simplify(f, method='SR')
            x^2 + 1

        """
        if method is None:
            method = self._current
        return self._simplify_dict[method](expression)

    def is_trivial_zero(self, expression, method=None):
        r"""
        Check if an expression is trivially equal to zero without any
        simplification.

        INPUT:

        - ``expression`` -- expression

        - ``method`` -- (default: ``None``) string defining the calculus method
          to use; if ``None`` the current calculus method of ``self`` is used.

        OUTPUT:

        - ``True`` is expression is trivially zero, ``False`` elsewhere.

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='real')
            sage: f = sin(x) - sin(x)
            sage: calc_meth.is_trivial_zero(f)
            True
            sage: calc_meth.is_trivial_zero(f._sympy_(), method='sympy')
            True

        ::

            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: calc_meth.is_trivial_zero(f)
            False
            sage: calc_meth.is_trivial_zero(f._sympy_(), method='sympy')
            False

        """
        if method is None:
            method = self._current
        if method == 'SR':
            return expression.is_trivial_zero()
        elif method == 'sympy':
            # we have to test SymPy's is_zero because it could be 'NoneType'
            if expression.is_zero:
                return True
            return False

    def set(self, method):
        r"""
        Set the current calculus method.

        - ``method`` -- string defining the calculus method

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='complex')
            sage: calc_meth.set('sympy')
            sage: calc_meth
            Possible calculus methods:
             - SR (default)
             - sympy (*)
            sage: calc_meth.reset()

        """
        self._test(method)
        self._current = method

    def _test(self, method):
        r"""
        Test if a calculus method is available.

        INPUT:

        - ``method`` -- string defining the calculus method

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='complex')
            sage: calc_meth._test('lala')
            Traceback (most recent call last):
            ...
            NotImplementedError: method lala not implemented yet

        """
        if method not in self._methods:
            raise NotImplementedError("method " + method + " not " +
                                      "implemented yet")

    def reset(self):
        r"""
        Set the current symbolic method to default one.

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='complex')
            sage: calc_meth.set('sympy')
            sage: calc_meth
            Possible calculus methods:
             - SR (default)
             - sympy (*)
            sage: calc_meth.reset()
            sage: calc_meth
            Possible calculus methods:
             - SR (*) (default)
             - sympy

        """
        self._current = self._default

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_meth = CalculusMethod(base_field_type='complex')
            sage: calc_meth._repr_()
            'Possible calculus methods:\n - SR (*) (default)\n - sympy'

        """
        resu = 'Possible calculus methods:\n'
        for method in self._methods:
            resu += ' - {}'.format(method)
            if method == self._current: resu += ' (*)'
            if method == self._default: resu += ' (default)'
            resu += '\n'
        return resu[:-1]
