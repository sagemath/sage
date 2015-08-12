r"""
Formatting utilities

This module defines helper functions that are not class methods.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Joris Vankerschaver (2010): for the function :func:`is_atomic()`

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

import six
from sage.structure.sage_object import SageObject

def is_atomic(expression):
    r"""
    Helper function to check whether some LaTeX expression is atomic.

    Adapted from method
    :meth:`~sage.tensor.differential_form_element.DifferentialFormFormatter._is_atomic`
    of class
    :class:`~sage.tensor.differential_form_element.DifferentialFormFormatter`
    written by Joris Vankerschaver (2010).

    INPUT:

    - ``expression`` -- string representing the expression (e.g. LaTeX string)

    OUTPUT:

    - ``True`` if additive operations are enclosed in parentheses and
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import is_atomic
        sage: is_atomic("2*x")
        True
        sage: is_atomic("2+x")
        False
        sage: is_atomic("(2+x)")
        True

    """
    if not isinstance(expression, six.string_types):
        raise TypeError("The argument must be a string")
    level = 0
    for n, c in enumerate(expression):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        if c == '+' or c == '-':
            if level == 0 and n > 0:
                return False
    return True


def is_atomic_wedge_txt(expression):
    r"""
    Helper function to check whether some text-formatted expression is atomic
    in terms of wedge products.

    Adapted from method
    :meth:`~sage.tensor.differential_form_element.DifferentialFormFormatter._is_atomic`
    of class
    :class:`~sage.tensor.differential_form_element.DifferentialFormFormatter`
    written by Joris Vankerschaver (2010).

    INPUT:

    - ``expression`` -- string representing the text-formatted expression

    OUTPUT:

    - ``True`` if wedge products are enclosed in parentheses and
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import is_atomic_wedge_txt
        sage: is_atomic_wedge_txt("a")
        True
        sage: is_atomic_wedge_txt(r"a/\b")
        False
        sage: is_atomic_wedge_txt(r"(a/\b)")
        True
        sage: is_atomic_wedge_txt(r"(a/\b)/\c")
        False
        sage: is_atomic_wedge_txt(r"(a/\b/\c)")
        True

    """
    if not isinstance(expression, six.string_types):
        raise TypeError("The argument must be a string.")
    level = 0
    for n, c in enumerate(expression):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        if c == '/' and expression[n+1:n+2] == '\\':
            if level == 0 and n > 0:
                return False
    return True


def is_atomic_wedge_latex(expression):
    r"""
    Helper function to check whether LaTeX-formatted expression is atomic in
    terms of wedge products.

    Adapted from method
    :meth:`~sage.tensor.differential_form_element.DifferentialFormFormatter._is_atomic`
    of class
    :class:`~sage.tensor.differential_form_element.DifferentialFormFormatter`
    written by Joris Vankerschaver (2010).

    INPUT:

    - ``expression`` -- string representing the LaTeX expression

    OUTPUT:

    - ``True`` if wedge products are enclosed in parentheses and
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import is_atomic_wedge_latex
        sage: is_atomic_wedge_latex(r"a")
        True
        sage: is_atomic_wedge_latex(r"a\wedge b")
        False
        sage: is_atomic_wedge_latex(r"(a\wedge b)")
        True
        sage: is_atomic_wedge_latex(r"(a\wedge b)\wedge c")
        False
        sage: is_atomic_wedge_latex(r"((a\wedge b)\wedge c)")
        True
        sage: is_atomic_wedge_latex(r"(a\wedge b\wedge c)")
        True
        sage: is_atomic_wedge_latex(r"\omega\wedge\theta")
        False
        sage: is_atomic_wedge_latex(r"(\omega\wedge\theta)")
        True
        sage: is_atomic_wedge_latex(r"\omega\wedge(\theta+a)")
        False

    """
    if not isinstance(expression, six.string_types):
        raise TypeError("The argument must be a string.")
    level = 0
    for n, c in enumerate(expression):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        if c == '\\' and expression[n+1:n+6] == 'wedge':
            if level == 0 and n > 0:
                return False
    return True


def format_mul_txt(name1, operator, name2):
    r"""
    Helper function for text-formatted names of results of multiplication or
    tensor product.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import format_mul_txt
        sage: format_mul_txt('a', '*', 'b')
        'a*b'
        sage: format_mul_txt('a+b', '*', 'c')
        '(a+b)*c'
        sage: format_mul_txt('a', '*', 'b+c')
        'a*(b+c)'
        sage: format_mul_txt('a+b', '*', 'c+d')
        '(a+b)*(c+d)'
        sage: format_mul_txt(None, '*', 'b')
        sage: format_mul_txt('a', '*', None)

    """
    if name1 is None or name2 is None:
        return None
    if not is_atomic(name1) or not is_atomic_wedge_txt(name1):
        name1 = '(' + name1 + ')'
    if not is_atomic(name2) or not is_atomic_wedge_txt(name2):
        name2 = '(' + name2 + ')'
    return name1 + operator + name2


def format_mul_latex(name1, operator, name2):
    r"""
    Helper function for LaTeX names of results of multiplication or tensor
    product.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import format_mul_latex
        sage: format_mul_latex('a', '*', 'b')
        'a*b'
        sage: format_mul_latex('a+b', '*', 'c')
        '\\left(a+b\\right)*c'
        sage: format_mul_latex('a', '*', 'b+c')
        'a*\\left(b+c\\right)'
        sage: format_mul_latex('a+b', '*', 'c+d')
        '\\left(a+b\\right)*\\left(c+d\\right)'
        sage: format_mul_latex(None, '*', 'b')
        sage: format_mul_latex('a', '*', None)

    """
    if name1 is None or name2 is None:
        return None
    if not is_atomic(name1) or not is_atomic_wedge_latex(name1):
        name1 = r'\left(' + name1 + r'\right)'
    if not is_atomic(name2) or not is_atomic_wedge_latex(name2):
        name2 = r'\left(' + name2 + r'\right)'
    return name1 + operator + name2


def format_unop_txt(operator, name):
    r"""
    Helper function for text-formatted names of results of unary operator.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import format_unop_txt
        sage: format_unop_txt('-', 'a')
        '-a'
        sage: format_unop_txt('-', 'a+b')
        '-(a+b)'
        sage: format_unop_txt('-', '(a+b)')
        '-(a+b)'
        sage: format_unop_txt('-', None)

    """
    if name is None:
        return None
    if not is_atomic(name) or not is_atomic_wedge_txt(name):
    #!# is_atomic_otimes_txt should be added
        name = '(' + name + ')'
    return operator + name


def format_unop_latex(operator, name):
    r"""
    Helper function for LaTeX names of results of unary operator.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import format_unop_latex
        sage: format_unop_latex('-', 'a')
        '-a'
        sage: format_unop_latex('-', 'a+b')
        '-\\left(a+b\\right)'
        sage: format_unop_latex('-', '(a+b)')
        '-(a+b)'
        sage: format_unop_latex('-', None)

    """
    if name is None:
        return None
    if not is_atomic(name) or not is_atomic_wedge_latex(name):
    #!# is_atomic_otimes_latex should be added
        name = r'\left(' + name + r'\right)'
    return operator + name


class FormattedExpansion(SageObject):
    r"""
    Helper class for displaying tensor expansions.

    EXAMPLES::

        sage: from sage.tensor.modules.format_utilities import FormattedExpansion
        sage: f = FormattedExpansion('v', r'\tilde v')
        sage: f
        v
        sage: latex(f)
        \tilde v
        sage: f = FormattedExpansion('x/2', r'\frac{x}{2}')
        sage: f
        x/2
        sage: latex(f)
        \frac{x}{2}

    """
    def  __init__(self, txt=None, latex=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.format_utilities import FormattedExpansion
            sage: f = FormattedExpansion('v', r'\tilde v')
            sage: f
            v

        """
        self._txt = txt
        self._latex = latex

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.format_utilities import FormattedExpansion
            sage: f = FormattedExpansion('v', r'\tilde v')
            sage: f._repr_()
            'v'

        """
        return self._txt

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLE::

            sage: from sage.tensor.modules.format_utilities import FormattedExpansion
            sage: f = FormattedExpansion('v', r'\tilde v')
            sage: f._latex_()
            '\\tilde v'

        """
        return self._latex
