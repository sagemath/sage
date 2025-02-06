# -*- coding: utf-8 -*-
r"""
Sheaves on ChowSchemes

A sheaf on a ChowScheme `X` is an element of the Grothendieck ring `G(X)` of
coherent sheaves on `X` with the rational coefficients.
In view of the ring isomorphism given by the Chern character:

.. math::

    ch: G(X)\otimes_{\ZZ}\QQ\ \xrightarrow{\simeq}\ A(X)^{*}\otimes_{\ZZ}\QQ

a sheaf will be represented by its Chern character.

The class :class:`Sheaf` implements most of the common operations on sheaves.
Note that in view of the above all these operations should be understood
as "virtual". For example, the tensor product of two sheaves `E` and `F`
returns the alternating sum of the `\text{Tor}_i(E,F)` and `\text{SHom}(E,F)`
the alternating sum of the `\text{Ext}^i(E,F)`. Similarly, if `f` is a proper
morphism of ChowSchemes, `f_*` of a sheaf `E` returns the alternating sum
of the higher direct images `R^if_*E`.

A sheaf on a ChowScheme `X` may be specified either by its Chern character
or its rank `r` and Chern classes `c_1, \dots, c_r`. Typically, for example
in order to define a rank `2` sheaf with Chern classes `c_1=0` and `c_2=4`
on `\mathbb{P}^{2}` one writes ::

    sage: P2.<h> = Proj(2, 'h', name='P2')
    sage: E = Sheaf(P2, 2, [1, 0, 4*h^2]); E
    Sheaf(P2, 2, [1, 0, 4*h^2])

Note that if one wants to emphasize that `E` is actually a bundle, one may
use the module :mod:`bundle`::

    sage: E = Bundle(P2, 2, [1, 0, 4*h^2]); E
    Bundle(P2, 2, [1, 0, 4*h^2])

The class :class`bundle.Bundle` derives from :class:`Sheaf` and may have some
additional methods in future versions which apply only for vector bundles.
Also tensor products, etc. of bundles return bundles.

TESTS::


    sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
    sage: O1 = Sheaf(P2, 1, [1, h])
    sage: TestSuite(O1).run()

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

from sage.rings.integer import Integer
from sage.rings.power_series_ring import PowerSeriesRing
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method


class Sheaf(SageObject):
    r"""
    A *sheaf* is an element of the Grothendieck ring of coherent sheaves on a
    ChowScheme (tensored with `\QQ`) and represented by its Chern character.
    This class provides operations on sheaves to be understood as virtual, e.g.
    in the Grothendieck ring.
    """
    def __init__(self, X, r=None, cc=None, ch=None, name=None, latex_name=None):
        """
        Construct a :class:`Sheaf`.

        INPUT:

        - ``X`` -- the base ChowScheme

        - ``r`` -- an (optional) integer representing the rank of the bundle

        - ``cc`` -- an (optional) element or a list of elements of the Chow ring
          of the underlying chowscheme representing the Chern classes.

        - ``ch`` -- an (optional) element of the Chow ring of the underlying
          chowscheme representing the Chern character

        - ``name`` -- an (optional) string representing the name

        - ``latex_name`` -- an (optional) string representing the latex name

        OUTPUT:

        - :class:`Sheaf`.

        EXAMPLES::

            sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
            sage: E = Sheaf(P2, 2, [1, h, -3*h^2]); E
            Sheaf(P2, 2, [1, h, -3*h^2])
            sage: E = Sheaf(P2, 2, 1 + h - 3*h^2); E
            Sheaf(P2, 2, [1, h, -3*h^2])
            sage: E = Sheaf(P2, ch=7/2*h^2 + h + 2); E
            Sheaf(P2, 2, [1, h, -3*h^2])
        """
        # Options
        self._name = name
        self._latex_name = latex_name
        # Check that first argument is a ChowScheme:
        from sage.schemes.chow.scheme import is_chowScheme
        if not is_chowScheme(X):
            raise TypeError("ChowScheme expected, got %s" % X)
        self._chowscheme = X
        # Check rank
        rank = 1  # Default
        if r is not None:
            r = int(r) if isinstance(r, Integer) else r
            if not isinstance(r, int):
                raise TypeError("Expect an integer, got %s" % r)
            rank = r
        # Check Chern classes
        AX, d = X.chowring(), X.dimension()
        chern_classes = [AX(1)] + [AX(0)] * d  # Default
        if cc is not None:
            try:
                cc = AX(cc)
            except TypeError:
                pass
            if cc in AX:
                cc = cc.by_degrees()
            if not type(cc) is list:
                raise TypeError("Expect a list of Chern classes", cc)
            if len(cc) > d + 1:
                raise ValueError("Expect at most %s chern classes:"
                                 % (d + 1), cc)
            # Fill up with 0 if too few Chern classes are given.
            cc = cc + [AX(0)] * (d + 1 - len(cc))
            # Checks and Convert (if necessary)
            for i in range(1, d + 1):
                if AX(cc[i]) != AX(0):
                    if AX(cc[i]).lift().degree() != i:
                        raise ValueError("c_%d not in degree %d" % (i, i),
                                         AX(cc[i]))
                chern_classes[i] = AX(cc[i])
            if chern_classes[0] != AX(1):
                raise ValueError("Expect c_0 = 1.", cc[0])
        total_cc = sum([c for c in chern_classes])
        # Check Chern character
        if ch is None:
            # Expect rank and total_cc (or use defaults r=1 and total_cc=1)
            ch = AX(rank + total_cc._logg())
        else:
            if r is not None or cc is not None:
                err = "Expect either rank and Chern classes or Chern character."
                raise TypeError(err)
            try:
                ch = AX(ch)
                ch = ch.truncate(0, X.dimension())
            except:
                raise TypeError("Expect ch in the Chow ring")

        self._chern_character = ch

    def chowscheme(self):
        r"""
        Return the underlying ChowScheme.

        OUTPUT:

        - a ChowScheme

        EXAMPLE::

            sage: X = ChowScheme(4, 'h', 1, 'h^5')
            sage: S = Sheaf(X, 1, [1])
            sage: S.chowscheme()
            ChowScheme(4, 'h', 1, 'h^5')
        """
        return self._chowscheme

    def chern_character(self):
        r"""
        Return the Chern character of this Sheaf.

        OUTPUT:

        - an element of the underlying Chow ring of this Sheaf.

        EXAMPLE::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: ch = S.chern_character(); ch.by_degrees()
            [3, c1, 1/2*c1^2 - c2, 1/6*c1^3 - 1/2*c1*c2 + 1/2*c3]
        """
        return self._chern_character

    @cached_method
    def rank(self):
        r"""
        Return the rank of this Sheaf.

        OUTPUT:

        - an integer.

        EXAMPLE::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: S.rank()
            3
        """
        return int(self.chern_character().by_degrees()[0])

    @cached_method
    def total_chern_class(self):
        r"""
        Return the total Chern class of this Sheaf

        OUTPUT:

        - an element of the underlying Chow ring of this Sheaf.

        EXAMPLE::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: S.total_chern_class()
            c3 + c2 + c1 + 1
        """
        return self.chern_character()._expp()

    @cached_method
    def chern_classes(self):
        r"""
        Return the Chern classes of this Sheaf

        OUTPUT:

        - a list of elements of the underlying Chow ring of this Sheaf.

        EXAMPLE::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: ch = 3 + c1 + (1/2*c1^2 - c2) + (1/6*c1^3 - 1/2*c1*c2 + 1/2*c3)
            sage: S = Sheaf(X, ch=ch)
            sage: S.chern_classes()
            [1, c1, c2, c3]
        """
        return self.total_chern_class().by_degrees()

    def _repr_(self):
        r"""
        Return a string representation of this Sheaf.

        OUTPUT:

        - a string.

        EXAMPLE::

            sage: X.<c1, c2> = ChowScheme(3, ['c1', 'c2'], [1, 2], name='X')
            sage: S = Sheaf(X, 2, [1, c1, c2]); S._repr_()
            'Sheaf(X, 2, [1, c1, c2])'
        """
        return "Sheaf(%s, %d, %s)" % (self.chowscheme(), self.rank(),
                                      self.chern_classes())

    def _name_(self):
        r"""
        Return the name of this Sheaf. This is an internal function that
        returns the name if defined otherwise None.

        TESTS::

            sage: X.<c1, c2> = ChowScheme(3, ['c1', 'c2'], [1, 2], name='X')
            sage: S = Sheaf(X, 2, [1, c1, c2], name='S'); S._name_()
            'S'
            sage: S = Sheaf(X, 2, [1, c1, c2]); S._name_()
        """
        return self._name

    def __str__(self):
        r"""
        Return a text representation of this Sheaf. This is generally
        the name specified during definition or the string representation
        if a name has not been specified.

        OUTPUT:

        - a string.

        EXAMPLES::

            sage: X.<c1, c2> = ChowScheme(3, ['c1', 'c2'], [1, 2], name='X')
            sage: S = Sheaf(X, 2, [1, c1, c2], name='S'); S
            Sheaf(X, 2, [1, c1, c2])
            sage: str(S)
            'S'
        """
        return self._name if self._name is not None else repr(self)

    def _latex_name_(self):
        r"""
        Return the latex_name of this Sheaf. This is an internal method that
        return the latex_name of this sheaf if defined otherwise None.

        TESTS::

            sage: X.<c1, c2> = ChowScheme(3, ['c1', 'c2'], [1, 2], name='X')
            sage: S = Sheaf(X, 2, [1, c1, c2], latex_name='S'); S._latex_name_()
            'S'
            sage: S = Sheaf(X, 2, [1, c1, c2]); S._latex_name_()
        """
        return self._latex_name

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: X.<c1, c2> = ChowScheme(3, ['c1', 'c2'], [1, 2], name='X')
            sage: S = Sheaf(X, 2, [1, c1, c2], name='S')
            sage: S._latex_()
            '\\mathcal{S}'
            sage: S = Sheaf(X, 2, [1, c1, c2], name='S', latex_name='\\bold(S)')
            sage: S._latex_()
            '\\bold(S)'
            sage: latex(S)
            \bold(S)


        """
        if self._latex_name_() is not None:
            return self._latex_name_()
        return r"\mathcal{%s}" % self._name_() if self._name_() else repr(self)

    @cached_method
    def symbolic_chern_polynomial(self, var_name='t'):
        r"""
        Return the symbolic Chern polynomial of this Sheaf.

        INPUT:

        - ``var_name`` -- a (optional) string representing the formal variable.

        OUTPUT:

        - a symbolic polynomial in ``var_name``.

        EXAMPLE::

            sage: X.<c1, c2> = ChowScheme(3, ['c1', 'c2'], [1, 2], name='X')
            sage: S = Sheaf(X, 2, [1, -c1, c2])
            sage: S.symbolic_chern_polynomial()
            1 - c1*t + c2*t^2
        """
        cc = self.chern_classes()
        ring = cc[0].parent()
        R = PowerSeriesRing(ring, var_name, self.chowscheme().dimension() + 1)
        t = R.gen()
        return sum(cci * t ** i for i, cci in enumerate(cc))

    @cached_method
    def segre_classes(self):
        r"""
        Returns `[s_0, ..., s_r]` where the `s_i` are the Segre classes.

        Nota Bene: We use the convention that `s_t(E)*c_t(E^*) = 1`
        (in Fulton's book, `s_t(E)*c_t(E) = 1` is used).

        OUTPUT:

        - a list of elements of the underlying Chow ring of this Sheaf.

        EXAMPLE::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: S.segre_classes()
            [1, c1, c1^2 - c2, c1^3 - 2*c1*c2 + c3]
        """
        c_t = self.dual().symbolic_chern_polynomial()
        s_t = ~c_t
        return s_t.list()

    @cached_method
    def symbolic_segre_polynomial(self, var_name='t'):
        r"""
        Return the symbolic Segre polynomial `s_t`.

        Nota Bene: We use the convention that `s_t(E)*c_t(E^*) = 1`
        (in Fulton's book, `s_t(E)*c_t(E) = 1` is used).

        INPUT:

        - ``var_name`` -- a (optional) string representing the formal variable.

        OUTPUT:

        - a symbolic polynomial in ``var_name``.

        EXAMPLES::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: S.symbolic_segre_polynomial()
            1 + c1*t + (c1^2 - c2)*t^2 + (c1^3 - 2*c1*c2 + c3)*t^3
        """
        cc = self.segre_classes()
        ring = cc[0].parent()
        R = PowerSeriesRing(ring, var_name, self.chowscheme().dimension() + 1)
        t = R.gen()
        return sum(cc[i] * t ** i for i in range(len(cc)))

    @cached_method
    def total_segre_class(self):
        r"""
        Return the total Segre class of this Sheaf.

        Nota Bene: We use the convention that `s_t(E)*c_t(E^*) = 1`
        (in Fulton's book, `s_t(E)*c_t(E) = 1` is used).

        OUTPUT:

        - an element of the underlying Chow ring of this Sheaf.

        EXAMPLES::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: S.total_segre_class()
            c1^3 - 2*c1*c2 + c3 + c1^2 - c2 + c1 + 1
        """
        sc = self.segre_classes()
        return sum(sc[k] for k in range(len(sc)))

    def __eq__(self, other):
        r"""
        Two sheaves are considered to be "equal" if the underlying varieties
        and their Chern characters are the same.

        OUTPUT:

        - a boolean.

        EXAMPLES::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: T = Sheaf(X, 3, [1, -c1, c2, -c3]).dual()
            sage: S == T
            True
        """
        if not isinstance(other, Sheaf):
            return False
        if not self.chowscheme() == other.chowscheme():
            return False
        if self.chern_character() != other.chern_character():
            return False
        return True

    def __ne__(self, other):
        r"""
        Opposite of equal.

        OUTPUT:

        - a boolean.

        EXAMPLES::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: T = Sheaf(X, 3, [1, -c1, c2, -c3]).dual()
            sage: S != T
            False
        """
        return not self.__eq__(other)

    def __hash__(self):
        r"""
        Return the hash of self.

        EXAMPLES::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3])
            sage: H = hash(S)
        """
        return hash((self.chowscheme(), self.chern_character()))

    def _binary_operation_names(self, other, op, latex_op):
        r"""
        Return the binary operation name and latex_name. Internal to this class.

        TESTS::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Sheaf(X, 3, [1, -c1, c2, -c3], name='T')
            sage: S._binary_operation_names(T, '+', '\\oplus')
            ('S + T', '\\mathcal{S} + \\mathcal{T}')
        """
        E, F = self, other
        name, latex_name = None, None
        # If E and F are bundles use latex_op else op
        from sage.schemes.chow.bundle import is_bundle
        lop = latex_op if (is_bundle(E) and is_bundle(F)) else op
        if E._name_() or F._name_():
            name = '%s %s %s' % (str(E), op, str(F))
            latex_name = '%s %s %s' % (E._latex_(), lop, F._latex_())
        if E._latex_name_() or F._latex_name_():
            latex_name = '%s %s %s' % (E._latex_(), lop, F._latex_())
        return name, latex_name

    def __add__(self, other):
        r"""
        Return the sum of this Sheaf and another Sheaf.

        INPUT:

        - another Sheaf on the underlying ChowScheme.

        OUTPUT:

        - a sheaf (or a bundle if both sheaves are bundles)

        TESTS::

            sage: variables = ['c1', 'c2', 'c3', 'd1', 'd2', 'd3']
            sage: degrees = [1, 2, 3, 1, 2, 3]
            sage: X = ChowScheme(3, variables, degrees, name='X')
            sage: (c1, c2, c3, d1, d2, d3) = X.gens()
            sage: S = Sheaf(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Sheaf(X, 3, [1, d1, d2, d3], name='T')
            sage: R = S + T; R
            Sheaf(X, 6, [1, c1 + d1, c2 + c1*d1 + d2, c3 + c2*d1 + c1*d2 + d3])
            sage: str(R)
            'S + T'
            sage: latex(R)
            \mathcal{S} + \mathcal{T}
            sage: S = Bundle(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Bundle(X, 3, [1, d1, d2, d3], name='T')
            sage: R = S + T; R.chern_classes()
            [1, c1 + d1, c2 + c1*d1 + d2, c3 + c2*d1 + c1*d2 + d3]
            sage: str(R)
            'S + T'
            sage: latex(R)
            \mathcal{S} \oplus \mathcal{T}
        """
        E, F = self, other
        # Checks
        if not E.chowscheme() is F.chowscheme():
            raise ValueError("Base varieties of sheaves differ.")
        # Computation
        ch = E.chern_character() + F.chern_character()
        # Result
        name, latex_name = self._binary_operation_names(other, '+', '\\oplus')
        result = Sheaf(E.chowscheme(), ch=ch, name=name, latex_name=latex_name)
        # If both are called from a derived class like bundle return a bundle
        if E.__class__ == F.__class__:
            result.__class__ = E.__class__
        return result

    def __sub__(self, other):
        r"""
        Return the difference of this Sheaf and another Sheaf.

        INPUT:

        - another Sheaf on the underlying ChowScheme.

        OUTPUT:

        - a sheaf.

        TESTS::

            sage: variables = ['c1', 'c2', 'c3', 'd1', 'd2', 'd3']
            sage: degrees = [1, 2, 3, 1, 2, 3]
            sage: X = ChowScheme(3, variables, degrees, name='X')
            sage: (c1, c2, c3, d1, d2, d3) = X.gens()
            sage: S = Sheaf(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Sheaf(X, 3, [1, d1, d2, d3], name='T')
            sage: R = S - T; R
            Sheaf(X, 0, [1, c1 - d1, c2 - c1*d1 + d1^2 - d2, c3 - c2*d1 + c1*d1^2 - d1^3 - c1*d2 + 2*d1*d2 - d3])
            sage: str(R)
            'S - T'
            sage: latex(R)
            \mathcal{S} - \mathcal{T}
            sage: S = Bundle(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Bundle(X, 3, [1, d1, d2, d3], name='T')
            sage: R = S - T; R
            Sheaf(X, 0, [1, c1 - d1, c2 - c1*d1 + d1^2 - d2, c3 - c2*d1 + c1*d1^2 - d1^3 - c1*d2 + 2*d1*d2 - d3])
            sage: str(R)
            'S - T'
            sage: latex(R)
            \mathcal{S} - \mathcal{T}
        """
        E, F = self, other
        if not E.chowscheme() is F.chowscheme():
            raise ValueError("Base varieties of sheaves differ.")
        ch = E.chern_character() - F.chern_character()
        # Result
        name, latex_name = self._binary_operation_names(other, '-', '-')
        result = Sheaf(E.chowscheme(), ch=ch, name=name, latex_name=latex_name)
        return result

    def __mul__(self, other):
        r"""
        Return the product of this Sheaf and another Sheaf. If both sheaves are
        bundles, this is the tensor product. In general, the result is virtual,
        equal to the alternating sum of the `\text{Tor}_i(E,F)`.

        INPUT:

        - another Sheaf on the underlying ChowScheme.

        OUTPUT:

        - a sheaf (or a bundle if both sheaves are bundles)

        TESTS::

            sage: variables = ['c1', 'c2', 'c3', 'd1', 'd2', 'd3']
            sage: degrees = [1, 2, 3, 1, 2, 3]
            sage: X = ChowScheme(3, variables, degrees, name='X')
            sage: (c1, c2, c3, d1, d2, d3) = X.gens()
            sage: S = Sheaf(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Sheaf(X, 3, [1, d1, d2, d3], name='T')
            sage: R = S * T; R
            Sheaf(X, 9, [1, 3*c1 + 3*d1, 3*c1^2 + 3*c2 + 8*c1*d1 + 3*d1^2 + 3*d2, c1^3 + 6*c1*c2 + 3*c3 + 7*c1^2*d1 + 7*c2*d1 + 7*c1*d1^2 + d1^3 + 7*c1*d2 + 6*d1*d2 + 3*d3])
            sage: str(R)
            'S * T'
            sage: latex(R)
            \mathcal{S} * \mathcal{T}
            sage: S = Bundle(X, 3, [1, c1, c2, c3], name='S')
            sage: T = Bundle(X, 3, [1, d1, d2, d3], name='T')
            sage: R = S * T; R
            Bundle(X, 9, [1, 3*c1 + 3*d1, 3*c1^2 + 3*c2 + 8*c1*d1 + 3*d1^2 + 3*d2, c1^3 + 6*c1*c2 + 3*c3 + 7*c1^2*d1 + 7*c2*d1 + 7*c1*d1^2 + d1^3 + 7*c1*d2 + 6*d1*d2 + 3*d3])
            sage: str(R)
            'S * T'
            sage: latex(R)
            \mathcal{S} \otimes \mathcal{T}
        """
        E, F = self, other
        if not E.chowscheme() is F.chowscheme():
            raise ValueError("Base varieties of sheaves differ.")
        ch = E.chern_character() * F.chern_character()
        # Result
        name, latex_name = self._binary_operation_names(other, '*', '\\otimes')
        result = Sheaf(E.chowscheme(), ch=ch, name=name, latex_name=latex_name)
        # If both are called from a derived class like bundle return a bundle
        if self.__class__ == other.__class__:
            result.__class__ = self.__class__
        return result

    def _unary_op_names(self, op, latex_op):
        r"""
        Return the unary operation name and latex_name. Internal to this class.
        The character string 'x' in ``op`` and ``latex_op`` is replaced by this
        sheaf's name and latex_name respectively.

        EXAMPLES::

            sage: X.<c1, c2, c3> = ChowScheme(3, ['c1', 'c2', 'c3'], [1, 2, 3])
            sage: S = Sheaf(X, 3, [1, c1, c2, c3], name='S')
            sage: S._unary_op_names('x^3', 'x^{3}')
            ('S^3', '\\mathcal{S}^{3}')
        """
        name, latex_name = None, None
        if self._name_() is not None:
            name = op.replace('x', self._name_())
            latex_name = latex_op.replace('x', self._latex_())
        if self._latex_name_() is not None:
            latex_name = latex_op.replace('x', self._latex_name_())
        return name, latex_name

    def __pow__(self, n):
        r"""
        Return the `n`-th power of this Sheaf. If the sheaf is a bundle,
        this is the `n`-th tensor product. In general, the result is virtual.
        For negative `n`, this sheaf is required to be a rank 1 bundle.

        INPUT:

        - an integer

        OUTPUT:

        - a sheaf (or a bundle this sheaf is a bundle)

        TESTS::

            sage: variables = ['c1', 'c2', 'c3', 'd1', 'd2', 'd3']
            sage: degrees = [1, 2, 3, 1, 2, 3]
            sage: X = ChowScheme(3, variables, degrees, name='X')
            sage: (c1, c2, c3, d1, d2, d3) = X.gens()
            sage: S = Sheaf(X, 3, [1, c1, c2, c3], name='S')
            sage: S^3
            Sheaf(X, 27, [1, 27*c1, 342*c1^2 + 27*c2, 2702*c1^3 + 666*c1*c2 + 27*c3])
            sage: S^(-1)
            Traceback (most recent call last):
            ...
            ValueError: Negative powers are only allowed in rank 1.
            sage: T = Bundle(X, 1, [1, c1])
            sage: T^(-3)
            Bundle(X, 1, [1, -3*c1])
        """
        X = self.chowscheme()
        from sage.schemes.chow.bundle import is_bundle
        if n < 0 and not (self.rank() == 1 and is_bundle(self)):
            raise ValueError("Negative powers are only allowed in rank 1.")
        if n == 0:
            result = Sheaf(X, 1, [1, 0])
        elif n == 1:
            result = self
        else:
            p = n if n > 0 else -n
            ch = self.chern_character() ** p
            a, b = self._unary_op_names('x^%s' % n, 'x^{%s}' % n)
            result = Sheaf(X, ch=ch, name=a, latex_name=b)
            if n < 0:
                ch = ch.adams(-1)  # Already checked that rank=1 if n<0
                a, b = self._unary_op_names('x^(%s)' % n, 'x^{%s}' % n)
                result = Sheaf(X, ch=ch, name=a, latex_name=b)
        # If called from a derived class like bundle return a bundle
        result.__class__ = self.__class__
        return result

    @cached_method
    def dual(self):
        r"""
        Return the dual of this Sheaf. If the sheaf is a bundle, this is the
        dual bundle. In general, the result is virtual.

        OUTPUT:

        - a sheaf (or a bundle this sheaf is a bundle)

        TESTS::

            sage: variables = ['c1', 'c2', 'c3', 'c4']
            sage: degrees = [1, 2, 3, 4]
            sage: X = ChowScheme(4, variables, degrees, name='X')
            sage: (c1, c2, c3, c4) = X.gens()
            sage: E = Sheaf(X, 4, [1, c1, c2, c3, c4], name='E')
            sage: ED = E.dual(); ED
            Sheaf(X, 4, [1, -c1, c2, -c3, c4])
            sage: str(ED)
            'E^*'
            sage: latex(ED)
            \mathcal{E}^{*}
            sage: E = Sheaf(X, 1, [1, c1], name='E', latex_name='\\bold{E}')
            sage: latex(E.dual())
            \bold{E}^{*}
        """
        X, ch = self.chowscheme(), self.chern_character()
        ch = ch.adams(-1)
        name, latex_name = self._unary_op_names('x^*', 'x^{*}')
        result = Sheaf(X, ch=ch, name=name, latex_name=latex_name)
        # If called from a derived class like bundle return a bundle
        result.__class__ = self.__class__
        return result

    @cached_method
    def determinant(self):
        r"""
        Return the determinant of this sheaf. Synonym for determinant().

        EXAMPLES::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: E = Sheaf(X, 2, [1, c1, c2], name='E')
            sage: E.det()
            Sheaf(X, 1, [1, c1])

        TESTS::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: E = Sheaf(X, 2, [1, c1, c2], name='E')
            sage: str(E.det())
            'det(E)'
            sage: latex(E.det())
            \det(\mathcal{E})
        """
        name, latex_name = self._unary_op_names('det(x)', r'\det(x)')
        result = Sheaf(self.chowscheme(), 1, [1, self.chern_classes()[1]],
                       name=name, latex_name=latex_name)
        # If called from a derived class like bundle return a bundle
        result.__class__ = self.__class__
        return result

    def det(self):
        r"""
        Return the determinant of this sheaf. Synonym for determinant().

        TESTS::

            sage: P1 = Proj(1, name='P1')
            sage: (P1.o() - P1.o(-1)).det()
            Sheaf(P1, 1, [1, h])
        """
        return self.determinant()

    @cached_method
    def adams(self, k):
        r"""
        Return the k-th Adams operator applied to this sheaf.

        EXAMPLES::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: E = Bundle(X, 1, [1, c1])
            sage: all([(E.adams(n) == E^n) for n in range(-2,3)])
            True
            sage: F = Bundle(X, 2, [1, c1, c2])
            sage: all([E.adams(k) + F.adams(k) == (E + F).adams(k) for k in range(-2, 3)])
            True
            sage: F.adams(2).adams(3) == F.adams(2*3)
            True
        """
        X, ch = self.chowscheme(), self.chern_character()
        ch = ch.adams(k)
        a, b = self._unary_op_names('psi_%s(x)' % k, r'\psi_{%s}(x)' % k)
        result = Sheaf(X, ch=ch, name=a, latex_name=b)
        # If called from a derived class like bundle return a bundle
        result.__class__ = self.__class__
        return result

    @cached_method
    def wedge(self, p):
        r"""
        Return the p-th exterior power of this sheaf.

        EXAMPLES::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: F = Bundle(X, 2, [1, c1, c2])
            sage: F.wedge(2) == F.det()
            True

        TESTS::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: F = Bundle(X, 2, [1, c1, c2])
            sage: F.wedge(0) == TrivialBundle(X, 1)
            True
            sage: F.wedge(1) == F
            True
        """
        X, ch = self.chowscheme(), self.chern_character()
        ch = ch.wedge(p)
        a, b = self._unary_op_names(r'/\^%s(x)' % p, r'\Lambda^{%s}(x)' % p)
        result = Sheaf(X, ch=ch, name=a, latex_name=b)
        # If called from a derived class like bundle return a bundle
        result.__class__ = self.__class__
        return result

    @cached_method
    def symm(self, p):
        r"""
        Return the p-th symmetric power of this sheaf.

        EXAMPLES::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: F = Bundle(X, 2, [1, c1, c2])
            sage: F.adams(2) == (F.symm(2) - F.wedge(2))
            True

        TESTS::

            sage: X.<c1, c2> = ChowScheme(2, ['c1', 'c2'], [1, 2], name='X')
            sage: F = Bundle(X, 2, [1, c1, c2])
            sage: F.symm(0) == TrivialBundle(X, 1)
            True
            sage: F.symm(1) == F
            True
        """
        X, ch = self.chowscheme(), self.chern_character()
        ch = ch.symm(p)
        a, b = self._unary_op_names('Sym^%s(x)' % p, r'Sym^{%s}(x)' % p)
        result = Sheaf(X, ch=ch, name=a, latex_name=b)
        # If called from a derived class like a bundle return a bundle
        result.__class__ = self.__class__
        return result

    @cached_method
    def todd_class(self):
        r"""
        Returns the todd class of this sheaf.

        EXAMPLES::

            sage: P3 = Proj(3)
            sage: P3.todd_class()
            h^3 + 11/6*h^2 + 2*h + 1
            sage: P3.todd_class() == P3.tangent_bundle().todd_class()
            True
        """
        return self.chern_character().todd()

    @cached_method
    def euler_characteristic(self):
        r"""
        Return the euler-characteristic of this sheaf using
        Hirzebruch-Riemann-Roch. Requires the underlying ChowScheme to be
        a point.

        EXAMPLES::

            sage: P3 = Proj(3)
            sage: P3.o(-4).euler_characteristic()
            -1
        """
        if not self.chowscheme().base_chowscheme().is_point():
            raise TypeError('The sheaf is not aver an absolute ChowScheme.')
        ch, td = self.chern_character(), self.chowscheme().todd_class()
        return (ch * td).integral()

    @cached_method
    def sans_denominateurs(self, E):
        r"""
        Returns the expression P(self,E) (see Fulton 296-297, Lemma 15.3).

        EXAMPLES::

            sage: P2 = Proj(2)
            sage: P5 = Proj(5, 'k')
            sage: f = P2.hom(['2*h'], P5)
            sage: N = f.upperstar(P5.tangent_bundle()) - P2.tangent_bundle(); N
            Sheaf(Proj(2, 'h'), 3, [1, 9*h, 30*h^2])
            sage: N.sans_denominateurs(P2.tangent_bundle())
            [4, 36*h, 240*h^2]
        """
        if not isinstance(E, Sheaf):
            raise ValueError("Sheaf expected,", E)
        if not self.chowscheme() is E.chowscheme():
            raise ValueError("Base varieties differs :", E)

        n, A = self.chowscheme().dimension(), self.chowscheme().chowring()
        d, e = self.rank(), E.rank()
        R, P = sans_denominateurs_p(n, d, e)
        cd = self.chern_classes()[1: d + 1]
        ce = E.chern_classes()[1: e + 1]
        f = R.hom(cd + [0] * (d - len(cd)) + ce + [0] * (e - len(ce)), A)
        # Return the list as list of elements of the chowring (eg reduce)
        return [f(p) for p in P]


def is_sheaf(x):
    r"""
        Test whether ``x`` is a sheaf.

        INPUT:

        - ``x`` -- anything.

        OUTPUT:

        Boolean. Return ``True`` if ``x`` is a Sheaf, False otherwise.

        TESTS::

            sage: P5 = Proj(5)
            sage: is_sheaf(P5.o(1))
            True
    """
    return isinstance(x, Sheaf)


def sans_denominateurs_p(n, d, e):
    r"""
    Return the expression `P(D,E)` (see Fulton 296-297, Lemma 15.3) where
    `D` and `E` are bundles of ranks `d` and `e` respectively on a variety
    `V` of dimension `n`.

    Note that we compute *up to dimension* `n + d` as in applications `D`
    is generally the normal bundle of an embedding `V\rightarrow W`.

    INPUT:

    - ``n`` -- the dimension of `V`
    - ``d`` -- the rank of `D`
    - ``e`` -- the rank of `E`

    OUTPUT:

    A tuple `(A, [p_O,\dots,p_{n+d}])` of a ring of coefficients `A` and a list
    of coefficients `p_i` such that `P(D,E)=\sum p_i`.

    EXAMPLES::

        sage: from sage.schemes.chow.sheaf import sans_denominateurs_p
        sage: A, P = sans_denominateurs_p(3, 1, 1); P
        [1, d1 - e1, d1^2 - 2*d1*e1 + e1^2, d1^3 - 3*d1^2*e1 + 3*d1*e1^2 - e1^3]
        sage: A, P = sans_denominateurs_p(3, 1, 4); P
        [4, 10*d1 - e1, 20*d1^2 - 5*d1*e1 + e1^2 - 2*e2, 35*d1^3 - 15*d1^2*e1 + 6*d1*e1^2 - e1^3 - 11*d1*e2 + 3*e1*e2 - 3*e3]
        sage: A, P = sans_denominateurs_p(2, 2, 4); P
        [-4, -4*d1 + 2*e1, -4*d1^2 + 10*d2 + 3*d1*e1 - 3*e1^2 + 6*e2]
        sage: A, P = sans_denominateurs_p(1, 3, 4); P
        [8, 12*d1 - 6*e1]
    """
    # We compute without relations (and in dimension d + rank)
    from sage.schemes.chow.scheme import ChowScheme
    dim = n + d
    dis = ['d%s' % (i + 1) for i in range(d)]  # d_1, ..., d_d
    eis = ['e%s' % (i + 1) for i in range(e)]  # e_1, ..., e_e
    deg = [i + 1 for i in range(d)] + [i + 1 for i in range(e)]
    X = ChowScheme(dim, dis + eis, deg, name='X')
    D = Sheaf(X, d, [str(1)] + [str(dis[i]) for i in range(d)])
    E = Sheaf(X, e, [str(1)] + [str(eis[i]) for i in range(e)])
    ch = D.dual().chern_character()
    q = sum((-1) ** i * ch.wedge(i) for i in range(d + 1))
    lhs = (q * E.chern_character())._expp().by_degrees()
    top = D.chern_classes()[d]
    P = [(t / top) for t in lhs[d: len(lhs) + 1]]
    # Return the list as list of elements of the chowring (eg reduce)
    return X.chowring(), P


def SHom(E, F):
    r"""
    Return the sheaf of homomorphisms from `E` to `F`.

    INPUT:

    - ``E`` - a sheaf
    - ``F`` - a sheaf

    OUTPUT:

    A sheaf.

    EXAMPLES::

        sage: P2 = Proj(2, 'h')
        sage: SHom(P2.sheaves['universal_sub'], P2.sheaves['universal_quotient'])
        Bundle(Proj(2, 'h'), 2, [1, 3*h, 3*h^2])
        sage: P2.tangent_bundle()
        Bundle(Proj(2, 'h'), 2, [1, 3*h, 3*h^2])
    """
    if not (is_sheaf(E) and is_sheaf(F)):
        raise ValueError("SHom is defined between sheaves.")
    if not E.chowscheme() is F.chowscheme():
        raise ValueError("Base varieties of sheaves differ.")
    return E.dual() * F


def SEnd(E):
    r"""
    Return the sheaf of endomorphisms of `E`.

    INPUT:

    - ``E`` - a sheaf

    OUTPUT:

    A sheaf.

    EXAMPLES::

        sage: P2 = Proj(2, 'h')
        sage: SEnd(P2.tangent_bundle())
        Bundle(Proj(2, 'h'), 4, [1, 0, 3*h^2])
    """
    return SHom(E, E)
