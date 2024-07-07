"""
MathML output support

In order to support MathML formatting, an object should define a special
method _mathml_(self) that returns its MathML representation.
"""

# *****************************************************************************
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
# *****************************************************************************


def list_function(x):
    return 'MATHML version of the list %s' % (x,)


def tuple_function(x):
    return 'MATHML version of the tuple %s' % (x,)


def bool_function(x):
    return 'MATHML version of %s' % (x,)


def str_function(x):
    return 'MATHML version of the string %s' % (x,)


# One can add to the latex_table in order to install latexing
# functionality for other types.

mathml_table = {
    list: list_function,
    tuple: tuple_function,
    bool: bool_function,
    str: str_function,
    float: str,
    int: str
}


class MathML(str):
    def __repr__(self):
        return str(self)


def mathml(x):
    """
    Output x formatted for inclusion in a MathML document.
    """
    try:
        return MathML(x._mathml_())
    except (AttributeError, TypeError):
        for k, f in mathml_table.items():
            if isinstance(x, k):
                return MathML(f(x))

        if x is None:
            return MathML("MATHML version of 'None'")

        return MathML(str_function(str(x)))
