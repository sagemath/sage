# sage_setup: distribution = sagemath-symbolics
r"""
Coordinate calculus methods

The class :class:`CalculusMethod` governs the calculus methods (symbolic and
numerical) used for coordinate computations on manifolds.

AUTHORS:

- Marco Mancini (2017): initial version
- Eric Gourgoulhon (2019): add :meth:`~CalculusMethod.set_simplify_function`
  and various accessors
"""

# *****************************************************************************
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR
from sage.manifolds.utilities import (simplify_chain_real,
                                      simplify_chain_generic,
                                      simplify_chain_real_sympy,
                                      simplify_chain_generic_sympy,)
from sage.misc.latex import latex

try:
    import sympy
    from sympy import latex as sympy_latex
except ImportError:
    sympy_latex = None


# Conversion functions
def _SR_to_Sympy(expression):
    r"""
    Convert an expression from ``SR`` to ``sympy``.

    In the case where the expression is already ``sympy``, it is simply
    returned.

    INPUT:

    - ``expression`` -- ``SR`` or ``sympy`` symbolic expression

    OUTPUT: ``expression`` -- ``sympy`` symbolic expression

    EXAMPLES::

        sage: from sage.manifolds.calculus_method import _SR_to_Sympy
        sage: a = x^2 + sin(x)^2; a
        x^2 + sin(x)^2
        sage: type(a)
        <class 'sage.symbolic.expression.Expression'>
        sage: b = _SR_to_Sympy(a); b
        x**2 + sin(x)**2
        sage: type(b)
        <class 'sympy.core.add.Add'>

    An already ``sympy`` expression is unchanged::

        sage: _SR_to_Sympy(b) is b
        True
    """
    # Nothing to do if expression is already a SymPy object:
    if isinstance(expression, sympy.Basic):
        return expression
    return SR(expression)._sympy_()


def _Sympy_to_SR(expression):
    r"""
    Convert an expression from ``sympy`` to ``SR``.

    In the case where the expression is already ``SR``, it is simply returned.

    INPUT:

    - ``expression`` -- ``sympy`` symbolic expression

    OUTPUT: ``expression`` -- ``SR`` or ``sympy`` symbolic expression

    EXAMPLES::

        sage: from sage.manifolds.calculus_method import _Sympy_to_SR, _SR_to_Sympy
        sage: a = x^2 + sin(x)^2
        sage: type(a)
        <class 'sage.symbolic.expression.Expression'>
        sage: b = _SR_to_Sympy(a); b
        x**2 + sin(x)**2
        sage: type(b)
        <class 'sympy.core.add.Add'>
        sage: _Sympy_to_SR(b) == a
        x^2 + sin(x)^2 == x^2 + sin(x)^2
        sage: bool(_)
        True
    """
    try:
        return SR(expression)
    except Exception:
        # If SR cannot transform a sympy expression this is because it is a
        # sympy abstract function
        a = expression._sage_()
        # As all sage objects have a ._sage_ operator, they have to be
        # caught
        if type(a) is type(expression):
            raise TypeError
        return a


class CalculusMethod(SageObject):
    r"""
    Control of calculus backends used on coordinate charts of manifolds.

    This class stores the possible calculus methods and permits to switch
    between them, as well as to change the simplifying functions associated
    with them.
    For the moment, only two calculus backends are implemented:

    - Sage's symbolic engine (Pynac + Maxima), implemented via the
      Symbolic Ring ``SR``
    - SymPy engine, denoted ``sympy`` hereafter

    INPUT:

    - ``current`` -- (default: ``None``) string defining the calculus method
      that will be considered as the active one, until it is changed by
      :meth:`set`; must be one of

      - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
      - ``'sympy'``: SymPy
      - ``None``: the default calculus method (``'SR'``)

    - ``base_field_type`` -- (default: ``'real'``) base field type of the
      manifold (cf.
      :meth:`~sage.manifolds.manifold.TopologicalManifold.base_field_type`)

    EXAMPLES::

        sage: from sage.manifolds.calculus_method import CalculusMethod
        sage: cm = CalculusMethod()

    In the display, the currently active method is pointed out with a star::

        sage: cm
        Available calculus methods (* = current):
         - SR (default) (*)
         - sympy

    It can be changed with :meth:`set`::

        sage: cm.set('sympy')
        sage: cm
        Available calculus methods (* = current):
         - SR (default)
         - sympy (*)

    while :meth:`reset` brings back to the default::

        sage: cm.reset()
        sage: cm
        Available calculus methods (* = current):
         - SR (default) (*)
         - sympy

    See :meth:`simplify_function` for the default simplification algorithms
    associated with each calculus method and :meth:`set_simplify_function` for
    introducing a new simplification algorithm.
    """
    _default = 'SR'  # default calculus method
    _methods = ('SR', 'sympy')  # implemented methods
    _tranf = {'SR': _Sympy_to_SR, 'sympy': _SR_to_Sympy}  # translators

    def __init__(self, current=None, base_field_type='real'):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='complex')
            sage: cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
        """
        self._current = self._default if current is None else current
        # Initialization of the dictionary of simplifying functions:
        self._simplify_dict = {}
        if base_field_type == 'real':
            self._simplify_dict['sympy'] = simplify_chain_real_sympy
            self._simplify_dict['SR'] = simplify_chain_real
        else:
            self._simplify_dict['sympy'] = simplify_chain_generic_sympy
            self._simplify_dict['SR'] = simplify_chain_generic
        # The default simplifying functions are saved:
        self._simplify_dict_default = self._simplify_dict.copy()
        self._latex_dict = {'sympy': sympy_latex, 'SR': latex}

    def simplify(self, expression, method=None):
        r"""
        Apply the simplifying function associated with a given calculus method
        to a symbolic expression.

        INPUT:

        - ``expression`` -- symbolic expression to simplify

        - ``method`` -- (default: ``None``) string defining the calculus method
          to use; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the current calculus method of ``self`` is used

        OUTPUT: the simplified version of ``expression``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x, y> = M.chart()
            sage: f = x^2 + sin(x)^2 + cos(x)^2
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='real')
            sage: cm.simplify(f)
            x^2 + 1

        Using a weaker simplifying function, like
        :func:`~sage.calculus.functional.simplify`, does not succeed in this
        case::

            sage: cm.set_simplify_function(simplify)
            sage: cm.simplify(f)
            x^2 + cos(x)^2 + sin(x)^2

        Back to the default simplifying function
        (:func:`~sage.manifolds.utilities.simplify_chain_real` in the present
        case)::

            sage: cm.set_simplify_function('default')
            sage: cm.simplify(f)
            x^2 + 1

        A ``SR`` expression, such as ``f``, cannot be simplified when the
        current calculus method is ``sympy``::

            sage: cm.set('sympy')
            sage: cm.simplify(f)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.symbolic.expression.Expression' object has no attribute 'combsimp'...

        In the present case, one should either transform ``f`` to a SymPy
        object::

            sage: cm.simplify(f._sympy_())
            x**2 + 1

        or force the calculus method to be ``'SR'``::

            sage: cm.simplify(f, method='SR')
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
          to use; if ``None`` the current calculus method of ``self`` is used

        OUTPUT: ``True`` is expression is trivially zero, ``False`` elsewhere

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='real')
            sage: f = sin(x) - sin(x)
            sage: cm.is_trivial_zero(f)
            True
            sage: cm.is_trivial_zero(f._sympy_(), method='sympy')
            True

        ::

            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: cm.is_trivial_zero(f)
            False
            sage: cm.is_trivial_zero(f._sympy_(), method='sympy')
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
        Set the currently active calculus method.

        - ``method`` -- string defining the calculus method

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='complex')
            sage: cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
            sage: cm.set('sympy')
            sage: cm
            Available calculus methods (* = current):
             - SR (default)
             - sympy (*)
            sage: cm.set('lala')
            Traceback (most recent call last):
            ...
            NotImplementedError: method lala not implemented
        """
        if method not in self._methods:
            raise NotImplementedError("method {} not ".format(method) +
                                      "implemented")
        self._current = method

    def current(self):
        r"""
        Return the active calculus method as a string.

        OUTPUT:

        String defining the calculus method; one of

        - ``'SR'`` -- Sage's default symbolic engine (Symbolic Ring)
        - ``'sympy'`` -- SymPy

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(); cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
            sage: cm.current()
            'SR'
            sage: cm.set('sympy')
            sage: cm.current()
            'sympy'
        """
        return self._current

    def set_simplify_function(self, simplifying_func, method=None):
        r"""
        Set the simplifying function associated to a given calculus method.

        INPUT:

        - ``simplifying_func`` -- either the string ``'default'`` for restoring
          the default simplifying function or a function ``f`` of a single
          argument ``expr`` such that ``f(expr)`` returns an object of the same
          type as ``expr`` (hopefully the simplified version of ``expr``), this
          type being

          - :class:`~sage.symbolic.expression.Expression` if ``method`` = ``'SR'``
          - a SymPy type if ``method`` = ``'sympy'``

        - ``method`` -- (default: ``None``) string defining the calculus method
          for which ``simplifying_func`` is provided; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the currently active calculus method of ``self`` is
            assumed

        EXAMPLES:

        On a real manifold, the default simplifying function is
        :func:`~sage.manifolds.utilities.simplify_chain_real` when the
        calculus method is ``SR``::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='real'); cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
            sage: cm.simplify_function() is \
            ....: sage.manifolds.utilities.simplify_chain_real
            True

        Let us change it to :func:`~sage.calculus.functional.simplify`::

            sage: cm.set_simplify_function(simplify)
            sage: cm.simplify_function() is simplify
            True

        Since ``SR`` is the current calculus method, the above is equivalent
        to::

            sage: cm.set_simplify_function(simplify, method='SR')
            sage: cm.simplify_function(method='SR') is simplify
            True

        We revert to the default simplifying function by::

            sage: cm.set_simplify_function('default')

        Then we are back to::

            sage: cm.simplify_function() is \
            ....: sage.manifolds.utilities.simplify_chain_real
            True
        """
        if method is None:
            method = self._current
        if simplifying_func == 'default':
            # the default function is restored
            self._simplify_dict[method] = self._simplify_dict_default[method]
        else:
            self._simplify_dict[method] = simplifying_func

    def simplify_function(self, method=None):
        r"""
        Return the simplifying function associated to a given calculus method.

        The simplifying function is that used in all computations involved
        with the calculus method.

        INPUT:

        - ``method`` -- (default: ``None``) string defining the calculus method
          for which ``simplifying_func`` is provided; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the currently active calculus method of ``self`` is
            assumed

        OUTPUT: the simplifying function

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod()
            sage: cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
            sage: cm.simplify_function()  # random (memory address)
            <function simplify_chain_real at 0x7f958d5b6758>

        The output stands for the function
        :func:`~sage.manifolds.utilities.simplify_chain_real`::

            sage: cm.simplify_function() is \
            ....: sage.manifolds.utilities.simplify_chain_real
            True

        Since ``SR`` is the default calculus method, we have::

            sage: cm.simplify_function() is cm.simplify_function(method='SR')
            True

        The simplifying function associated with ``sympy`` is
        :func:`~sage.manifolds.utilities.simplify_chain_real_sympy`::

            sage: cm.simplify_function(method='sympy')  # random (memory address)
            <function simplify_chain_real_sympy at 0x7f0b35a578c0>
            sage: cm.simplify_function(method='sympy') is \
            ....: sage.manifolds.utilities.simplify_chain_real_sympy
            True

        On complex manifolds, the simplifying functions are
        :func:`~sage.manifolds.utilities.simplify_chain_generic`
        and :func:`~sage.manifolds.utilities.simplify_chain_generic_sympy`
        for respectively ``SR`` and ``sympy``::

            sage: cmc = CalculusMethod(base_field_type='complex')
            sage: cmc.simplify_function(method='SR') is \
            ....: sage.manifolds.utilities.simplify_chain_generic
            True
            sage: cmc.simplify_function(method='sympy') is \
            ....: sage.manifolds.utilities.simplify_chain_generic_sympy
            True

        Note that the simplifying functions can be customized via
        :meth:`set_simplify_function`.
        """
        if method is None:
            method = self._current
        return self._simplify_dict[method]

    def reset(self):
        r"""
        Set the current calculus method to the default one.

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='complex')
            sage: cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
            sage: cm.set('sympy')
            sage: cm
            Available calculus methods (* = current):
             - SR (default)
             - sympy (*)
            sage: cm.reset()
            sage: cm
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy
        """
        self._current = self._default

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: cm = CalculusMethod(base_field_type='complex')
            sage: cm._repr_()
            'Available calculus methods (* = current):\n - SR (default) (*)\n - sympy'
        """
        resu = 'Available calculus methods (* = current):\n'
        for method in self._methods:
            resu += ' - {}'.format(method)
            if method == self._default:
                resu += ' (default)'
            if method == self._current:
                resu += ' (*)'
            resu += '\n'
        return resu[:-1]
