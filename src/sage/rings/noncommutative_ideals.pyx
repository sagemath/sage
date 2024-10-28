# sage.doctest: needs sage.modules
"""
Ideals of non-commutative rings

Generic implementation of one- and two-sided ideals of non-commutative rings.

AUTHOR:

- Simon King (2011-03-21), <simon.king@uni-jena.de>, :issue:`7797`.

EXAMPLES::

    sage: MS = MatrixSpace(ZZ,2,2)
    sage: MS*MS([0,1,-2,3])
    Left Ideal
    (
      [ 0  1]
      [-2  3]
    )
     of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
    sage: MS([0,1,-2,3])*MS
    Right Ideal
    (
      [ 0  1]
      [-2  3]
    )
     of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
    sage: MS*MS([0,1,-2,3])*MS
    Twosided Ideal
    (
      [ 0  1]
      [-2  3]
    )
     of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

See :mod:`~sage.algebras.letterplace.letterplace_ideal` for a more
elaborate implementation in the special case of ideals in free
algebras.

TESTS::

    sage: # needs sage.combinat
    sage: A = SteenrodAlgebra(2)
    sage: IL = A*[A.1+A.2,A.1^2]; IL
    Left Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
    sage: TestSuite(IL).run(skip=['_test_category'],verbose=True)
    running ._test_eq() . . . pass
    running ._test_new() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_pickling() . . . pass
"""
# ****************************************************************************
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element cimport MonoidElement
from sage.structure.parent cimport Parent
from sage.categories.monoids import Monoids
from sage.rings.ideal_monoid import IdealMonoid_c
from sage.rings.ideal import Ideal_generic

from sage.rings.integer_ring import ZZ


class IdealMonoid_nc(IdealMonoid_c):
    """
    Base class for the monoid of ideals over a non-commutative ring.

    .. NOTE::

        This class is essentially the same as
        :class:`~sage.rings.ideal_monoid.IdealMonoid_c`,
        but does not complain about non-commutative rings.

    EXAMPLES::

        sage: MS = MatrixSpace(ZZ,2,2)
        sage: MS.ideal_monoid()
        Monoid of ideals of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        INPUT:

        - ``R`` -- a ring

        TESTS::

            sage: from sage.rings.noncommutative_ideals import IdealMonoid_nc
            sage: MS = MatrixSpace(ZZ,2,2)
            sage: IdealMonoid_nc(MS)
            Monoid of ideals of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        self._IdealMonoid_c__R = R
        Parent.__init__(self, base=ZZ,
                        category=Monoids())
        self._populate_coercion_lists_()

    def _element_constructor_(self, x):
        r"""
        Create an ideal in this monoid from ``x``.

        INPUT:

        - ``x`` -- an ideal or a list of elements

        TESTS::

            sage: # needs sage.combinat
            sage: A = SteenrodAlgebra(2) # indirect doctest
            sage: IL = A*[A.1+A.2,A.1^2]; IL
            Left Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: IR = [A.1+A.2,A.1^2]*A; IR
            Right Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: IT = A*[A.1+A.2,A.1^2]*A; IT
            Twosided Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: M = IL.parent()
            sage: M([A.0, 1])
            Twosided Ideal (Sq(1), 1) of mod 2 Steenrod algebra, milnor basis

        ::

            sage: # needs sage.combinat
            sage: IL == loads(dumps(IL))
            True
            sage: IR == loads(dumps(IR))
            True
            sage: IT == loads(dumps(IT))
            True
        """
        side = "twosided"
        if isinstance(x, Ideal_nc):
            side = x.side()
            x = x.gens()
        elif isinstance(x, Ideal_generic):
            x = x.gens()
        cdef MonoidElement y = self._IdealMonoid_c__R.ideal(x, side=side)
        y._parent = self
        return y


class Ideal_nc(Ideal_generic):
    """
    Generic non-commutative ideal.

    All fancy stuff such as the computation of Groebner bases must be
    implemented in sub-classes. See :class:`~sage.algebras.letterplace.letterplace_ideal.LetterplaceIdeal`
    for an example.

    EXAMPLES::

        sage: MS = MatrixSpace(QQ,2,2)
        sage: I = MS*[MS.1,MS.2]; I
        Left Ideal
        (
          [0 1]
          [0 0],
        <BLANKLINE>
          [0 0]
          [1 0]
        )
         of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: [MS.1,MS.2]*MS
        Right Ideal
        (
          [0 1]
          [0 0],
        <BLANKLINE>
          [0 0]
          [1 0]
        )
         of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: MS*[MS.1,MS.2]*MS
        Twosided Ideal
        (
          [0 1]
          [0 0],
        <BLANKLINE>
          [0 0]
          [1 0]
        )
         of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
    """
    def __init__(self, ring, gens, coerce=True, side='twosided'):
        """
        Initialize ``self``.

        INPUT:

        - ``ring`` -- a ring

        - ``gens`` -- list or tuple of elements

        - ``coerce`` -- boolean (default: ``True``); first coerce the given
          list of elements into the given ring

        - ``side`` -- string (default ``'twosided'``); must be ``'left'``,
          ``'right'`` or ``'twosided'``. Determines whether the ideal is a
          left, right or twosided ideal.

        TESTS::

            sage: MS = MatrixSpace(ZZ,2,2)
            sage: from sage.rings.noncommutative_ideals import Ideal_nc
            sage: Ideal_nc(MS,[MS.1,MS.2], side='left')
            Left Ideal
            (
              [0 1]
              [0 0],
            <BLANKLINE>
              [0 0]
              [1 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: Ideal_nc(MS,[MS.1,MS.2], side='right')
            Right Ideal
            (
              [0 1]
              [0 0],
            <BLANKLINE>
              [0 0]
              [1 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        if side not in ['left', 'right', 'twosided']:
            raise ValueError("Ideals are left, right or twosided, but not %s" % side)
        self.__side = side
        Ideal_generic.__init__(self, ring, gens, coerce=coerce)

    def __repr__(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: # needs sage.combinat
            sage: A = SteenrodAlgebra(2)
            sage: A*[A.1+A.2,A.1^2]      # indirect doctest
            Left Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: [A.1+A.2,A.1^2]*A
            Right Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: A*[A.1+A.2,A.1^2]*A
            Twosided Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
        """
        return "%s Ideal %s of %s" % (self.__side.capitalize(),
                                      self._repr_short(), self.ring())

    def __eq__(self, right):
        """
        Ideals of different sidedness do not compare equal. Apart from
        that, the generators are compared.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: A = SteenrodAlgebra(2)
            sage: IR = [A.1+A.2,A.1^2]*A
            sage: IL = A*[A.1+A.2,A.1^2]
            sage: IT = A*[A.1+A.2,A.1^2]*A
            sage: IT == IL
            False
            sage: IR == [A.1+A.2,A.1^2]*A
            True
        """
        if not isinstance(right, Ideal_nc):
            return False
        if self.side() != right.side():
            return False
        S = set(self.gens())
        T = set(right.gens())
        if S == T:
            return True
        return False

    def __ne__(self, right):
        """
        Ideals of different sidedness do not compare equal. Apart from
        that, the generators are compared.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: A = SteenrodAlgebra(2)
            sage: IR = [A.1+A.2,A.1^2]*A
            sage: IL = A*[A.1+A.2,A.1^2]
            sage: IT = A*[A.1+A.2,A.1^2]*A
            sage: IT != IL
            True
            sage: IR != [A.1+A.2,A.1^2]*A
            False
        """
        return not self.__eq__(right)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: A = SteenrodAlgebra(2)
            sage: IR = [A.1+A.2,A.1^2]*A
            sage: IL = A*[A.1+A.2,A.1^2]
            sage: IT = A*[A.1+A.2,A.1^2]*A
            sage: hash(IT) == hash(IL)
            False
            sage: hash(IR) == hash([A.1^2,A.1+A.2]*A)
            True
        """
        return hash((self.parent(), self.__side, frozenset(self.gens())))

    def side(self):
        """
        Return a string that describes the sidedness of this ideal.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: A = SteenrodAlgebra(2)
            sage: IL = A*[A.1+A.2,A.1^2]
            sage: IR = [A.1+A.2,A.1^2]*A
            sage: IT = A*[A.1+A.2,A.1^2]*A
            sage: IL.side()
            'left'
            sage: IR.side()
            'right'
            sage: IT.side()
            'twosided'
        """
        return self.__side

    def __mul__(self, other):
        r"""
        Multiply ``self`` with ``other``.

        Multiplication of a one-sided ideal with its ring from the other side
        yields a two-sided ideal.

        Let `L` (resp. `R`) be a left (resp. right) ideal, then the product
        `LR` is a twosided ideal generated by `x y`, where `x` (resp. `y`)
        is a generator of `L` (resp. `R`).

        .. TODO::

            The product of left (resp. right) ideals is a left
            (resp. right) ideal. However, these do not necessarily have
            simple generating sets.

        TESTS::

            sage: MS = MatrixSpace(QQ,2,2)
            sage: IL = MS * [2*MS.0,3*MS.1]; IL
            Left Ideal
            (
              [2 0]
              [0 0],
            <BLANKLINE>
              [0 3]
              [0 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: IR = MS.3*MS; IR
            Right Ideal
            (
              [0 0]
              [0 1]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: IL * MS
            Twosided Ideal
            (
              [2 0]
              [0 0],
            <BLANKLINE>
              [0 3]
              [0 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: IL * IR
            Twosided Ideal
            (
              [0 3]
              [0 0]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            sage: IR * IR
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot multiply non-commutative ideals
        """
        if isinstance(other, Ideal_nc) and self.ring() == other.ring():
            if self.side() == "left" and other.side() == "right":
                gens = [z for z in (x * y for x in self.gens() for y in other.gens()) if z]
                return self.ring().ideal(gens, side='twosided')
            raise NotImplementedError("cannot multiply non-commutative ideals")

        if not isinstance(other, Ideal_nc):
            # Perhaps other is a ring and thus has its own
            # multiplication.
            if other == self.ring():
                if self.side() == 'right':
                    return self
                return self.ring().ideal(self.gens(), side='twosided')
        if not isinstance(self, Ideal_nc):
            # This may happen...
            if self == other.ring():
                if other.side() == 'left':
                    return other
                return other.ring().ideal(other.gens(), side='twosided')

        raise NotImplementedError("cannot multiply non-commutative ideals")
