"""
This module is deprecated. Use ``group_algebra`` instead.

TESTS::

    sage: from sage.algebras.group_algebra_new import GroupAlgebra
    sage: GroupAlgebra(GL(3, QQ), ZZ)
    doctest:...: DeprecationWarning: 
    Importing GroupAlgebra from here is deprecated. If you need to use it, please import it directly from sage.algebras.group_algebra
    See http://trac.sagemath.org/17779 for details.
    Group algebra of group "General Linear Group of degree 3 over Rational Field" over base ring Integer Ring
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.group_algebra', '*', deprecation=17779)
