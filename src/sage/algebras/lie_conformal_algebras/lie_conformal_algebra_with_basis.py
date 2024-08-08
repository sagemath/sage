# sage.doctest: needs sage.combinat sage.modules
"""
Lie Conformal Algebras With Basis

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.combinat.free_module import CombinatorialFreeModule


class LieConformalAlgebraWithBasis(CombinatorialFreeModule):
    """
    Abstract base class for a Lie conformal algebra with a
    preferred basis.

    This class provides no functionality, it simply passes the
    arguments to :class:`CombinatorialFreeModule`.

    EXAMPLES::

        sage: R = lie_conformal_algebras.Virasoro(QQbar);R
        The Virasoro Lie conformal algebra over Algebraic Field

    TESTS::

        sage: R = lie_conformal_algebras.Virasoro(QQ)
        sage: R.0
        L
        sage: R._repr_generator(R.0)
        'L'
        sage: R = lie_conformal_algebras.Affine(QQ, 'A1')
        sage: R.0
        B[alpha[1]]
        sage: R._repr_generator(R.0)
        'B[alpha[1]]'
        sage: R = lie_conformal_algebras.Affine(QQ, 'A1', names = ('e', 'h','f'))
        sage: R.0
        e
        sage: R._repr_generator(R.0)
        'e'
    """
    def __init__(self, R, basis_keys=None, element_class=None, category=None,
                 prefix=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: V = lie_conformal_algebras.Affine(QQ,'A1')
            sage: TestSuite(V).run()
        """
        if prefix is None:
            prefix = ''
            kwds['bracket'] = ''
            kwds['string_quotes'] = False

        default_category = LieConformalAlgebras(R).WithBasis()
        try:
            category = default_category.or_subcategory(category)
        except ValueError:
            category = default_category.Super().or_subcategory(category)

        super().__init__(R, basis_keys=basis_keys, element_class=element_class,
                         category=category, prefix=prefix, names=None, **kwds)
