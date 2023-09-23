"""
Support Python's numbers abstract base class

.. SEEALSO:: :pep:`3141` for more information about :class:`numbers`.

TESTS::

    sage: import numbers
    sage: isinstance(5, numbers.Integral)
    True
    sage: isinstance(5, numbers.Number)
    True
    sage: isinstance(5/1, numbers.Integral)
    False
    sage: isinstance(22/7, numbers.Rational)
    True
    sage: isinstance(1.3, numbers.Real)
    True
    sage: isinstance(CC(1.3), numbers.Real)
    False
    sage: isinstance(CC(1.3 + I), numbers.Complex)
    True
    sage: isinstance(RDF(1.3), numbers.Real)
    True
    sage: isinstance(CDF(1.3, 4), numbers.Complex)
    True
    sage: isinstance(AA(sqrt(2)), numbers.Real)                                         # needs sage.rings.number_field sage.symbolic
    True
    sage: isinstance(QQbar(I), numbers.Complex)                                         # needs sage.rings.number_field
    True

This doesn't work with symbolic expressions at all::

    sage: isinstance(pi, numbers.Real)                                                  # needs sage.symbolic
    False
    sage: isinstance(I, numbers.Complex)                                                # needs sage.rings.number_field
    False
    sage: isinstance(sqrt(2), numbers.Real)                                             # needs sage.rings.number_field sage.symbolic
    False

Because we do this, NumPy's ``isscalar()`` recognizes Sage types::

    sage: from numpy import isscalar                                                    # needs numpy
    sage: isscalar(3.141)                                                               # needs numpy
    True
    sage: isscalar(4/17)                                                                # needs numpy
    True
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
