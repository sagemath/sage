# sage_setup: distribution = sagemath-symbolics
"""
Substitution Maps

This object wraps Pynac ``exmap`` objects. These encode substitutions
of symbolic expressions. The main use of this module is to hook into
Pynac's ``subs()`` methods and pass a wrapper for the substitution map
back to Python.
"""

# ****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object cimport SageObject


cdef class SubstitutionMap(SageObject):

    cdef GExMap _gmapobj

    cpdef Expression apply_to(self, Expression expr, unsigned options):
        """
        Apply the substitution to a symbolic expression.

        EXAMPLES::

            sage: from sage.symbolic.expression import make_map
            sage: subs = make_map({x:x+1})
            sage: subs.apply_to(x^2, 0)
            (x + 1)^2
        """
        return new_Expression_from_GEx(expr._parent,
                                       expr._gobj.subs_map(self._gmapobj, options))

    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: from sage.symbolic.expression import make_map
            sage: make_map({x:x+1})
            SubsMap
        """
        return 'SubsMap'  # GEx_to_str(&x._gobj)


cdef SubstitutionMap new_SubstitutionMap_from_GExMap(const GExMap& smap):
    """
    Wrap a Pynac object into a Python object.

    INPUT:

    - ``smap`` -- a Pynac ``exmap``

    OUTPUT: a new Python :class:`SubstitutionMap`

    EXAMPLES::

        sage: from sage.symbolic.expression import make_map
        sage: make_map({x:x+1})
        SubsMap
    """
    cdef SubstitutionMap result
    result = <SubstitutionMap>SubstitutionMap.__new__(SubstitutionMap)
    result._gmapobj = smap
    return result


cpdef SubstitutionMap make_map(subs_dict):
    """
    Construct a new substitution map.

    OUTPUT: a new :class:`SubstitutionMap` for doctesting

    EXAMPLES::

        sage: from sage.symbolic.expression import make_map
        sage: make_map({x:x+1})
        SubsMap
    """
    cdef GExMap smap
    for k, v in subs_dict.iteritems():
        smap.insert(make_pair((<Expression>k)._gobj,
                              (<Expression>v)._gobj))
    return new_SubstitutionMap_from_GExMap(smap)
