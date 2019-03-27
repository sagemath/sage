# -*- encoding: utf-8 -*-
"""
Python 2 and 3 Compatibility

TESTS:

Test the functionality of ``six.with_metaclass`` and various issues
which came up with it. Sage used to have a custom version of
``with_metaclass``, but this is now fixed upstream. ::

    sage: from six import with_metaclass
    sage: class Meta(type): pass
    sage: class X(with_metaclass(Meta)): pass
    sage: type(X) is Meta
    True
    sage: issubclass(X, object)
    True
    sage: class Base(object): pass
    sage: class X(with_metaclass(Meta, Base)): pass
    sage: type(X) is Meta
    True
    sage: issubclass(X, Base)
    True
    sage: class Base2(object): pass
    sage: class X(with_metaclass(Meta, Base, Base2)): pass
    sage: type(X) is Meta
    True
    sage: issubclass(X, Base)
    True
    sage: issubclass(X, Base2)
    True
    sage: X.__mro__ == (X, Base, Base2, object) or X.__mro__
    True

Check that :trac:`18503` is fixed, i.e. that with_metaclass
works with cdef'ed metaclasses::

    sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
    sage: class X(with_metaclass(ClasscallMetaclass)): pass
    sage: type(X) is ClasscallMetaclass
    True
    sage: X.__mro__ == (X, object) or X.__mro__
    True

Check a fix for :trac:`16074`::

    sage: from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
    sage: from sage.modules.with_basis.morphism import ModuleMorphismByLinearity
    sage: from sage.structure.unique_representation import UniqueRepresentation
    sage: class ExteriorAlgebraDifferential(with_metaclass(
    ....:        InheritComparisonClasscallMetaclass,
    ....:        ModuleMorphismByLinearity, UniqueRepresentation
    ....:     )):
    ....:     pass
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

from six import *


def u(x):
    r"""
    Convert `x` to unicode, assuming UTF-8 encoding.

    Python2 behaviour:

    If input is unicode, returns the input.

    If input is str (assumed to be utf-8 encoded), convert to unicode.

    Python3 behaviour:

    If input is str, returns the input.

    If input is bytes (assumed to be utf-8 encoded), convert to unicode.

    EXAMPLES::

        sage: from sage.misc.six import u
        sage: u("500 €")
        u'500 \u20ac'
        sage: u(u"500 \u20ac")
        u'500 \u20ac'
    """
    if isinstance(x, text_type):  # py2 unicode and py3 str
        return x
    if isinstance(x, bytes):
        return x.decode("utf-8")
    raise TypeError('input has no conversion to unicode')
