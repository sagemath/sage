"IntegerFactorization objects"

from sage.structure.factorization import Factorization
from sage.rings.integer_ring import ZZ


class IntegerFactorization(Factorization):
    """
    A lightweight class for an ``IntegerFactorization`` object,
    inheriting from the more general ``Factorization`` class.

    In the ``Factorization`` class the user has to create a list
    containing the factorization data, which is then passed to the
    actual ``Factorization`` object upon initialization.

    However, for the typical use of integer factorization via
    the ``Integer.factor()`` method in ``sage.rings.integer``
    this is noticeably too much overhead, slowing down the
    factorization of integers of up to about 40 bits by a factor
    of around 10.  Moreover, the initialization done in the
    ``Factorization`` class is typically unnecessary:  the caller
    can guarantee that the list contains pairs of an ``Integer``
    and an ``int``, as well as that the list is sorted.

    AUTHOR:

    - Sebastian Pancratz (2010-01-10)
    """

    def __init__(self, x, unit=None, cr=False, sort=True, simplify=True,
                 unsafe=False):
        """
        Set ``self`` to the factorization object with list ``x``,
        which must be a sorted list of pairs, where each pair contains
        a factor and an exponent.

        If the flag ``unsafe`` is set to ``False`` this method delegates
        the initialization to the parent class, which means that a rather
        lenient and careful way of initialization is chosen.  For example,
        elements are coerced or converted into the right parents, multiple
        occurrences of the same factor are collected (in the commutative
        case), the list is sorted (unless ``sort`` is ``False``) etc.

        However, if the flag is set to ``True``, no error handling is
        carried out.  The list ``x`` is assumed to list of pairs.  The
        type of the factors is assumed to be constant across all factors:
        either ``Integer`` (the generic case) or ``int`` (as supported
        by the flag ``int_`` of the ``factor()`` method).  The type of
        the exponents is assumed to be ``int``.  The list ``x`` itself
        will be referenced in this factorization object and hence the
        caller is responsible for not changing the list after creating
        the factorization.  The unit is assumed to be either ``None`` or
        of type ``Integer``, taking one of the values `+1` or `-1`.

        EXAMPLES::

            sage: factor(15)
            3 * 5

        We check that :trac:`13139` is fixed::

            sage: from sage.structure.factorization_integer import IntegerFactorization
            sage: IntegerFactorization([(3, 1)], unsafe=True)
            3
        """
        if unsafe:
            if unit is None:
                self._Factorization__unit = ZZ._one_element
            else:
                self._Factorization__unit = unit

            self._Factorization__x        = x
            self._Factorization__universe = ZZ
            self._Factorization__cr       = cr

            if sort:
                self.sort()
            if simplify:
                self.simplify()

        else:
            super(IntegerFactorization, self).__init__(x, unit=unit, cr=cr,
                                                       sort=sort,
                                                       simplify=simplify)

    def __sort__(self, key=None):
        """
        Sort the factors in this factorization.

        INPUT:

        - ``key`` -- (default: ``None``) comparison key

        EXAMPLES::

            sage: F = factor(15)
            sage: F.sort(key=lambda x: -x[0])
            sage: F
            5 * 3
        """
        if key is not None:
            self.__x.sort(key=key)
        else:
            self.__x.sort()
