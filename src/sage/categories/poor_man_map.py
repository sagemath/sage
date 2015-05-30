"""
Poor Man's map
"""
#*****************************************************************************
#       Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import sage.structure.sage_object

class PoorManMap(sage.structure.sage_object.SageObject):
    """
    A class for maps between sets which are not (yet) modeled by parents

    Could possibly disappear when all combinatorial classes / enumerated sets will be parents
    """
    def __init__(self, function, domain = None, codomain = None, name = None):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain = [1,2,3], codomain = [1,2,6])
            sage: f
            A map from [1, 2, 3] to [1, 2, 6]
            sage: f(3)
            6
            sage: TestSuite(f).run()

        """
        self._function = function
        self._domain = domain
        self._codomain = codomain
        self._name = name

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: PoorManMap(lambda x: x+2)  # indirect doctest
            A map
            sage: PoorManMap(lambda x: x+2, domain = [1,2,3])
            A map from [1, 2, 3]
            sage: PoorManMap(lambda x: x+2, domain = (1,2,3))
            A map from (1, 2, 3)
            sage: PoorManMap(lambda x: x+2, codomain = [3,4,5])
            A map to [3, 4, 5]

        """
        return ((self._name if self._name is not None else "A map") +
                (" from %s"%(self._domain,) if self._domain   is not None else ""     ) +
                (" to %s"%(self._codomain,) if self._codomain is not None else ""     ))

    def domain(self):
        """
        Returns the domain of ``self``

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: PoorManMap(lambda x: x+1, domain = [1,2,3], codomain = [2,3,4]).domain()
            [1, 2, 3]
        """
        return self._domain

    def codomain(self):
        """
        Returns the codomain of ``self``

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: PoorManMap(lambda x: x+1, domain = [1,2,3], codomain = [2,3,4]).codomain()
            [2, 3, 4]
        """
        return self._codomain

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain = [1,2,3], codomain = [1,2,6])
            sage: g = PoorManMap(factorial, domain = [1,2,3], codomain = [1,2,6])
            sage: h1 = PoorManMap(factorial, domain = [1,2,3], codomain = [1,2,6,8])
            sage: h2 = PoorManMap(factorial, domain = [1,2,3], codomain = [1,2,6,8])
            sage: h3 = PoorManMap(factorial, domain = [1,2,3,4], codomain = [1,2,6])
            sage: h4 = PoorManMap(lambda x: x, domain = [1,2,3], codomain = [1,2,6])
            sage: f == g, f == h1, f == h2, f == h3, f == h4, f == 1, 1 == f
            (True, False, False, False, False, False, False)
        """
        return self.__class__ is other.__class__ and self.__dict__ == other.__dict__

    def __mul__(self, other):
        """
        Composition

        INPUT:
         - ``self`` -- a map `f`
         - ``other`` -- a map `g`

        Returns the composition map `f\circ g` of `f`` and `g`

        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(lambda x: x+1, domain = [1,2,3], codomain = [2,3,4])
            sage: g = PoorManMap(lambda x: -x,  domain = [2,3,4], codomain = [-2,-3,-4])
            sage: f*g
            A map from [2, 3, 4] to [2, 3, 4]

        """
        return PoorManComposeMap(self, other)

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(lambda x: x+1, domain = [1,2,3], codomain = [2,3,4])
            sage: f(2)
            3
        """
        return self._function(*args)

class PoorManComposeMap(PoorManMap):
    def __init__(self, f, g):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(factorial, domain = [1,2,3], codomain = [1,2,6])
            sage: g = PoorManMap(sqrt,  domain = [1,4,9], codomain = [1,2,3])
            sage: h = f*g
            sage: h.codomain()
            [1, 2, 6]
            sage: h.domain()
            [1, 4, 9]
            sage: TestSuite(h).run()
        """
        self.f = f
        self.g = g
        try:
            self._domain = g.domain()
        except AttributeError:
            self._domain = None
        self._codomain = f.codomain()
        self._name = None

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: from sage.categories.poor_man_map import PoorManMap
            sage: f = PoorManMap(lambda x: x+1, domain = [1,2,3], codomain = [2,3,4])
            sage: g = PoorManMap(lambda x: -x,  domain = [2,3,4], codomain = [-2,-3,-4])
            sage: (f*g)(2)
            -1

        """
        return self.f(self.g(*args))
