"""
libSingular: Function Factory

AUTHORS:

- Martin Albrecht (2010-01): initial version

"""
#*****************************************************************************
#       Copyright (C) 2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.singular.function import singular_function, lib, list_of_functions

class SingularFunctionFactory(object):
    """
    A convenient interface to libsingular functions.
    """
    def __getattr__(self, name):
        """
        EXAMPLES::

            sage: import sage.libs.singular.function_factory
            sage: groebner = sage.libs.singular.function_factory.ff.groebner
            sage: groebner
            groebner (singular function)

            sage: import sage.libs.singular.function_factory
            sage: primdecSY = sage.libs.singular.function_factory.ff.primdec__lib.primdecSY
            sage: primdecSY
            primdecSY (singular function)
        """
        if name.startswith("_"):
            raise AttributeError("Singular Function Factory has no attribute '%s'"%name)

        try:
            return singular_function(name)
        except NameError:
            if name.endswith("__lib"):
                name = name[:-5]
                lib(name+".lib")
                return SingularFunctionFactory()
            else:
                raise NameError("function or package '%s' unknown."%(name))

    def trait_names(self):
        """
        EXAMPLES::

            sage: import sage.libs.singular.function_factory
            sage: "groebner" in sage.libs.singular.function_factory.ff.trait_names()
            True
            
        """
        return list_of_functions()

ff = SingularFunctionFactory()
