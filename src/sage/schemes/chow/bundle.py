# -*- coding: utf-8 -*-
r"""
Bundles on a ChowScheme

A bundle on a ChowScheme `X` is represented by its rank `r` and Chern classes
`c_0, \dots, c_r`. Typically, for example in order to define a rank 2 bundles
with Chern classes `0` and `4` on the projective plane::

    sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
    sage: E = Bundle(P2, 2, [1, 0, 4*h^2]); E
    Bundle(P2, 2, [1, 0, 4*h^2])
    sage: E.chern_character()
    -4*h^2 + 2
    sage: Sheaf(P2, ch=E.chern_character())
    Sheaf(P2, 2, [1, 0, 4*h^2])

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""
# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.chow.sheaf import Sheaf
from sage.schemes.chow.scheme import PointChowScheme


class Bundle(Sheaf):
    r"""
    Class for Bundles on ChowSchemes.

    TESTS::

        sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
        sage: E = Bundle(P2, 2, [1, 0, 4*h^2])
        sage: TestSuite(E).run()
    """
    def __init__(self, X, r=None, cc=None, ch=None, name=None, latex_name=None):
        """
        Construct a :class:`Bundle`

        INPUT:

        - ``X`` -- the base ChowScheme
        - ``rank`` -- an integer representing the rank of the bundle
        - ``cc`` -- a list of elements of the chowring of ``X``
          representing the Chern classes.

        OUTPUT:

        - :class:`Bundle`

        TESTS::

            sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
            sage: B = Bundle(P2, 0, [1, h]); B
            Traceback (most recent call last):
            ...
            TypeError: Nontrivial Chern classes > rank.

            sage: B = Bundle(P2, 0, [1, 0]); B
            Bundle(P2, 0, [1])
        """
        if (r is not None) and (cc is not None):
            if len(cc) > r + 1:
                A = X.chowring()
                if any(A(x) != A(0) for x in cc[r + 1:]):
                    raise TypeError("Nontrivial Chern classes > rank.")
        Sheaf.__init__(self, X, r, cc, ch, name, latex_name)

    def _repr_(self):
        r"""
        Return a string representation of this Bundle.

        OUTPUT:

        - a string.

        TESTS::

            sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
            sage: B = Bundle(P2, 2, [1, h, 3*h^2])
            sage: B._repr_()
            'Bundle(P2, 2, [1, h, 3*h^2])'
        """
        X, r = self.chowscheme(), self.rank()
        cc = self.chern_classes()
        result = "Bundle(%s, %d, %s)" % (X, r, cc)
        return result


def is_bundle(x):
    """
    Test whether ``x`` is a bundle.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean. Return ``True`` if ``x`` is a bundle.

    EXAMPLES::

        sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
        sage: B = Bundle(P2, 2, [1, 0, 3*h^2])
        sage: is_bundle(B)
        True
        sage: is_bundle(P2)
        False
    """
    return isinstance(x, Bundle)


def BundleDiffRelations(B, A):
    r"""
    Return the relations given by the difference of two bundles `B` and `A`
    on a ChowScheme `X`.

    INPUT:

    - ``B`` -- a bundle on a ChowScheme `X`

    - ``A`` -- a bundle on a ChowScheme `X`

    OUTPUT:

    A list of elements of the Chow ring of `X`.

    EXAMPLE (get the relations for Grass(6, 4))::

        sage: from sage.schemes.chow.bundle import BundleDiffRelations
        sage: G64 = ChowScheme(8, ['w1', 'w2'], [1, 2])
        sage: O = TrivialBundle(G64, 6)
        sage: S = Bundle(G64, 2, [1, '-w1', 'w1^2-w2'])  # Universal Sub
        sage: rels = BundleDiffRelations(O, S); rels
        [-2*w1^3*w2 + 3*w1*w2^2, w1^6 - 2*w1^4*w2 + w2^3]
        sage: A = ChowRing(['w1', 'w2'], [1, 2], [str(x) for x in rels])
        sage: A.rels()  # Returns the relations in a standard basis.
        [w2^5, w1*w2^4, w1^2*w2^3 - 4/3*w2^4, w1^6 - 3*w1^2*w2^2 + w2^3, w1^3*w2 - 3/2*w1*w2^2]

        sage: Grass(6, 4, 'w').rels()
        [w2^5, w1*w2^4, w1^2*w2^3 - 4/3*w2^4, w1^6 - 3*w1^2*w2^2 + w2^3, w1^3*w2 - 3/2*w1*w2^2]


    """
    if not (isinstance(B, Bundle) and isinstance(A, Bundle)):
        raise TypeError("Expected two bundles %s and %s, got" % (B, A))
    if not (B.chowscheme() is A.chowscheme()):
        raise TypeError("Bundles are not over the same base.")

    D = B - A
    C = Bundle(D.chowscheme(), D.rank(), D.chern_classes()[0: D.rank() + 1])

    b = B.total_chern_class()
    a = A.total_chern_class()
    c = C.total_chern_class()

    return [x for x in (b - a * c).by_degrees() if x]


def TrivialBundle(X, r):
    r"""
    Return the trivial bundle of rank `r` on the ChowScheme `X`.

    INPUT:

    - ``X`` -- a ChowScheme, the base

    - ``r`` -- an integer, the rank.

    OUTPUT:

    The trivial bundle of rank ``r`` on ``X``.

    TESTS:

        sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
        sage: TrivialBundle(P2, 3)
        Bundle(P2, 3, [1])
    """
    return Bundle(X, r, [1] + [0] * min(r, X.dimension()))


PointChowScheme.sheaves["tangent"] = Bundle(PointChowScheme, 0, [1])
