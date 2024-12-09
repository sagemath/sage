# sage_setup: distribution = sagemath-objects
"""
Test of the :mod:`~sage.structure.factory` module
"""

# ****************************************************************************
#  Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
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
# *****************************************************************************

from sage.structure.factory import UniqueFactory


class A:
    # something we can weakref
    pass


class UniqueFactoryTester(UniqueFactory):

    def create_key(self, *args, **kwds):
        """
        EXAMPLES::

            sage: from sage.structure.test_factory import UniqueFactoryTester
            sage: test_factory = UniqueFactoryTester('foo')
            sage: test_factory.create_key(1, 2, 3)
            (1, 2, 3)
        """
        return args

    def create_object(self, version, key, **extra_args):
        """
        EXAMPLES::

            sage: from sage.structure.test_factory import UniqueFactoryTester
            sage: test_factory = UniqueFactoryTester('foo')
            sage: test_factory.create_object('version', key=(1, 2, 4))
            Making object (1, 2, 4)
            <sage.structure.test_factory.A object at ...>
        """
        print("Making object", key)
        return A()


test_factory = UniqueFactoryTester('sage.structure.test_factory.test_factory')
