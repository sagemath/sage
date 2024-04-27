"""
Default Settings

AUTHORS: William Stein and David Kohel
"""
# ****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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

# default variable name
var_name = 'x'


def variable_names(n, name=None):
    r"""
    Convert a root string into a tuple of variable names by adding
    numbers in sequence.

    INPUT:

    - ``n`` a non-negative Integer; the number of variable names to
       output
    - ``names`` a string (default: ``None``); the root of the variable
      name.

    EXAMPLES::

        sage: from sage.misc.defaults import variable_names
        sage: variable_names(0)
        ()
        sage: variable_names(1)
        ('x',)
        sage: variable_names(1,'alpha')
        ('alpha',)
        sage: variable_names(2,'alpha')
        ('alpha0', 'alpha1')
    """
    if name is None:
        name = var_name
    n = int(n)
    if n == 1:
        return (name,)
    return tuple(['%s%s' % (name, i) for i in range(n)])


def latex_variable_names(n, name=None):
    r"""
    Convert a root string into a tuple of variable names by adding
    numbers in sequence.

    INPUT:

    - ``n`` a non-negative Integer; the number of variable names to
      output
    - ``names`` a string (default: ``None``); the root of the variable
      name.

    EXAMPLES::

        sage: from sage.misc.defaults import latex_variable_names
        sage: latex_variable_names(0)
        ()
        sage: latex_variable_names(1,'a')
        ('a',)
        sage: latex_variable_names(3,beta)
        ('beta_{0}', 'beta_{1}', 'beta_{2}')
        sage: latex_variable_names(3,r'\beta')
        ('\\beta_{0}', '\\beta_{1}', '\\beta_{2}')
    """
    if name is None:
        name = var_name
    n = int(n)
    if n == 1:
        return (name,)
    return tuple(['%s_{%s}' % (name, i) for i in range(n)])


def set_default_variable_name(name, separator=''):
    r"""
    Change the default variable name and separator.
    """
    global var_name, var_sep
    var_name = str(name)
    var_sep = str(separator)


# default series precision
series_prec = 20


def series_precision():
    """
    Return the Sage-wide precision for series (symbolic,
    power series, Laurent series).

    EXAMPLES::

        sage: series_precision()
        20
    """
    return series_prec


def set_series_precision(prec):
    """
    Change the Sage-wide precision for series (symbolic,
    power series, Laurent series).

    EXAMPLES::

        sage: set_series_precision(5)
        sage: series_precision()
        5
        sage: set_series_precision(20)
    """
    global series_prec
    series_prec = prec
