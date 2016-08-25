"""
Default Settings

AUTHORS: William Stein and David Kohel
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# default variable name
var_name = 'x'



def variable_names(n, name=None):
    if name is None:
        name = var_name
    n = int(n)
    if n == 1:
        return [name]
    return tuple(['%s%s'%(name,i) for i in range(n)])

def latex_variable_names(n, name=None):
    if name is None:
        name = var_name
    n = int(n)
    if n == 1:
        return [name]
    v = tuple(['%s_{%s}'%(name,i) for i in range(n)])
    return v

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
