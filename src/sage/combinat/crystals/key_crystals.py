r"""
Crystals of Key Tableaux

AUTHORS:

- Travis Scrimshaw: Initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.tensor_product import CrystalOfTableaux
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.rings.integer_ring import ZZ

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from collections.abc import Sequence

class KeyTableau(ClonableArray):
    r"""
    A key tableau.

    For more information, see
    :class:`~sage.combinat.crystals.key_crystals.KeyTableaux`.
    """
    def __init__(self, parent, data, check=True):
        r"""
        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(2)
            sage: mg = Y.module_generators[0]
            sage: TestSuite(mg).run()
        """
        if check:
            data = [tuple(row) for row in data]
        ClonableArray.__init__(self, parent, data, check=check)

    def check(self):
        if list(self) not in self.parent():
            raise ValueError("not a key tableau")

    def _repr_(self):
        r"""
        EXAMPLES::
        """
        return repr([list(row) for row in self])

    def _repr_diagram(self):
        r"""
        Return a string representation of the diagram of ``self``.

        EXAMPLES::
        """
        if not self:
            return '0'
        ret = ""
        width = max(len(repr(val)) for row in self for val in row)
        base = f"{{:^{width}}}"
        return "\n".join("|" + " ".join(base.format(val) for val in row)
                         for row in reversed(self))

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.
        """
        if not self:
            return "{\\emptyset}"
        from sage.combinat.output import tex_from_array
        return tex_from_array(self)

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines())

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.
        """
        from sage.typeset.unicode_art import UnicodeArt
        if not self.data:
            return UnicodeArt(["0"])

        from sage.combinat.output import ascii_art_table
        import unicodedata
        v = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
        vl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND LEFT')
        table = [[None]*(self.cols-len(row)) + row for row in reversed(self)]
        ret = []
        for i,row in enumerate(ascii_art_table(table, use_unicode=True).splitlines()):
            if row[-1] == " ":
                if i % 2 == 0:
                    ret.append(row[:-1] + vl)
                else:
                    ret.append(row[:-1] + v)
            else:
                ret.append(row)
        return UnicodeArt(ret)

    def pp(self):
        r"""
        Pretty print ``self``.
        """
        print(self._repr_diagram())

    def _signature_data(self, i):
        """
        Return data for the signature rule.

        We cancel `-+` pairs (in that order) to compute the ``i``-signature
        and return the relevant data. We replace `i` with `+` and
        `i + 1` with a `-`.

        OUTPUT:

        A tuple consisting of the following:

        - the column of the rightmost unpaired `+` (``None`` if does not exist)
        - the column of the leftmost unpaired `-` (``None`` if does not exist)
        - the number of unpaired `+`
        - the number of unpaired `-`
        """
        P = self.parent()
        num_plus = 0
        pos_plus = []
        num_minus = 0
        pos_minus = None
        for width in range(P._width):
            for ind, row in enumerate(self):
                if len(row) <= width:
                    continue
                if row[width] == i:
                    pos_plus.append((ind, width))
                elif row[width] == i + 1:
                    if pos_plus:
                        pos_plus.pop()
                    else:
                        num_minus += 1
                        pos_minus = (ind, width)
        num_plus = len(pos_plus)
        if pos_plus:
            pos_plus = pos_plus[0]
        else:
            pos_plus = None
        return (pos_plus, pos_minus, num_plus, num_minus)

    def e(self, i, by_tableau=False):
        r"""
        Return the application of the Kashiwara raising operator
        `e_i` on ``self``.

        EXAMPLES::
        """
        P = self.parent()
        n = P._cartan_type.rank()
        if by_tableau:
            ret = self.to_tableau().f(n+1-i)
            if ret is None:
                return None
            return P.from_tableau(ret)

        pm = self._signature_data(i)[1]
        if pm is None:
            return None
        data = [list(row) for row in self]
        P = self.parent()
        ind = pm[0]  # the row index containing the unpaired i+1
        row = data[ind]
        upper = data[ind+1:]
        for j in range(pm[1], len(row)):
            if row[j] != i + 1:
                break
            row[j] = i
            for r2 in upper:
                if len(r2) > j and r2[j] == i:
                    r2[j] = i + 1
        data = tuple(map(tuple, data))
        return P.element_class(P, data, check=False)

# e_1
# 211
# 3322

# 211
# 3322
# `\phi`::
# 1144
# 455
# f_4
# 1144
# 555

# 111
# 3322

# 111
# 3322

# 122
# 3311

    def f(self, i, by_tableau=False):
        r"""
        Return the application of the Kashiwara lowering operator
        `f_i` on ``self``.

        EXAMPLES::
        """
        P = self.parent()
        n = P._cartan_type.rank()
        if by_tableau:
            ret = self.to_tableau().e(n+1-i)
            if ret is None:
                return None
            return P.from_tableau(ret)

        pp = self._signature_data(i)[0]
        if pp is None:
            return None
        data = [list(row) for row in self]
        P = self.parent()
        ind = pp[0]  # the row index containing the unpaired i
        row = data[ind]
        upper = data[ind+1:]
        for j in range(pp[1], -1, -1):
            if row[j] != i:
                break
            row[j] = i + 1
            for r2 in upper:
                if len(r2) > j and r2[j] == i + 1:
                    r2[j] = i
        if not P._check_content(data):
            return None
        data = tuple(map(tuple, data))
        return P.element_class(P, data, check=False)


# f_1
# 1
# 55511
# 422

# 1
# 422
# 55511
# Applying `\phi`::
# 11155
# 244
# 5
# e_4
# 11145
# 244
# 5

# 1
# 422
# 55521

# 1
# 55521
# 422

    def weight(self, root_lattice=False):
        r"""
        Return the weight of ``self``.

        INPUT:

        - ``root_lattice`` -- boolean determining whether weight should appear
          in root lattice or not in extended affine weight lattice.

        EXAMPLES::

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.weight()
            2*Lambda[0] + Lambda[1] - 4*Lambda[2] + Lambda[3] - 2*delta
            sage: x.weight(root_lattice=True)
            -2*alpha[0] - 3*alpha[1] - 5*alpha[2] - 3*alpha[3]
        """
        P = self.parent().weight_lattice_realization()
        wt = [P.base_ring().zero()] * len(L.basis())
        for row in self:
            for val in row:
                wt[val-1] += 1
        return L.from_vector(wt, coerce=False)

    def epsilon(self, i):
        r"""
        Return the number of `i`-colored arrows in the `i`-string above
        ``self`` in the crystal graph.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: y.epsilon(1)
            0
            sage: y.epsilon(2)
            3
            sage: y.epsilon(0)
            0
        """
        n = self.parent()._cartan_type.rank()
        return self.to_tableau().phi(n+1-i)
        return ZZ(self._signature_data(i)[3])  # number of -'s

    def phi(self, i):
        r"""
        Return the value `\varepsilon_i(Y) + \langle h_i,
        \mathrm{wt}(Y)\rangle`, where `h_i` is the `i`-th simple
        coroot and `Y` is ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y.phi(1)
            3
            sage: y.phi(2)
            -1
        """
        n = self.parent()._cartan_type.rank()
        return self.to_tableau().epsilon(n+1-i)
        return ZZ(self._signature_data(i)[2])  # number of +'s

    def to_tableau(self):
        """
        Return ``self`` as an element of the crystal of semistandard
        Young tableaux of type `A_n`.
        """
        P = self.parent()
        C = P._ssyt
        n = P._cartan_type.rank()
        ret = sum((sorted((n + 2 - row[i] for row in self if len(row) > i), reverse=True)
                   for i in range(P._width)), [])
        return C(list=ret)


class KeyTableaux(UniqueRepresentation, Parent):
    r"""
    Key tableaux
    """
    @staticmethod
    def __classcall_private__(cls, shape, n=None, category=None):
        r"""
        Normalize input to ensure a unique representation.

        INPUT:

        - ``shape`` -- the shape
        - ``n`` -- (optional) type `A_n`

        EXAMPLES::

            sage: Yinf = crystals.infinity.GeneralizedYoungWalls(3)
            sage: Yinf2 = crystals.infinity.GeneralizedYoungWalls(int(3))
            sage: Yinf is Yinf2
            True
        """
        shape = list(shape)
        if n is None:
            n = len(shape) - 1
        while shape and not shape[-1]:  # standardize by removing trailing 0's
            shape.pop()
        shape = tuple(shape)
        if n < len(shape) - 1:
            raise ValueError(f"the rank must be at least {len(shape)-1}")
        return super().__classcall__(cls, shape, n, category)

    def __init__(self, shape, n, category):
        r"""
        EXAMPLES::

            sage: Yinf = crystals.infinity.GeneralizedYoungWalls(3)
            sage: TestSuite(Yinf).run()
        """
        self._cartan_type = CartanType(['A', n])
        self._shape = shape
        self._width = max(shape, default=0)
        self._ssyt = CrystalOfTableaux(self._cartan_type, shape=sorted(shape, reverse=True))
        category = HighestWeightCrystals().Finite().or_subcategory(category)
        Parent.__init__(self, category=category)

    Element = KeyTableau

    @lazy_attribute
    def module_generators(self):
        r"""
        Return the highest weight element of ``self``.

        Note that this is not the generator of the corresponding module
        of ``self`` as `U_q(\mathfrak{gl}_n)^+`-module (which is instead
        :meth:`extremal_module_generator()`). This is for implementation
        details in the category of
        :class:`~sage.categories.highest_weight_crystals.HighestWeightCrystals`
        assuming the ``module_generators`` contains all of the highest weight
        elements and those generate ``self`` under `f_i`.
        """
        gen = self.extremal_module_generator()
        return (gen.to_highest_weight()[0],)

    @cached_method
    def extremal_module_generator(self):
        """
        Return the generator of ``self``, considered as an extremal weight
        module.
        """
        data = tuple([(i,)*ell for i, ell in enumerate(self._shape, start=1)])
        return self.element_class(self, data, check=False)

    def _element_constructor_(self, data, check=True):
        r"""
        Construct an element of ``self`` from ``data``.

        INPUT:

        - ``data`` -- a multilist

        EXAMPLES::

            sage: GYW = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = GYW([[],[1,0],[2,1]]) # indirect doctest
            sage: y
            [[], [1, 0], [2, 1]]
        """
        return self.element_class(self, data, check=check)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(4)
            sage: Y
            Crystal of generalized Young walls of type ['A', 4, 1]
        """
        return "Crystal of key tableaux of type {} and shape {}".format(self._cartan_type, self._shape)

    def __contains__(self, data):
        """
        Check if ``data`` is an element of ``self``.
        """
        if isinstance(data, self.element_class):
            return data.parent() is self

        # Check the shape agrees
        if not isinstance(data, Sequence):
            return False
        if len(self._shape) != len(data):
            return False
        if any((not isinstance(row, Sequence)) or len(row) != ell
               for row, ell in zip(data, self._shape)):
            return False

        # Check entries in each column are distinct
        for i in range(self._width):
            col = [row[i] for row in data if len(row) > i]
            if len(col) != len(set(col)):
                return False

        # Check the other nontrivial content conditions
        return self._check_content(data)

    def _check_content(self, data):
        """
        Check the content of ``data`` is an element of ``self``, where
        we assume that ``data`` is a filling of ``self._shape`` with
        with distinct entries in the columns.
        """
        # Check rows are weakly decreasing
        for row in data:
            if any(row[i] < row[i+1] for i in range(len(row)-1)):
                return False

        # Check the other condition
        for ind, row in enumerate(data):
            for ri, k in enumerate(row):
                i = max((row[ri] for row in data[ind+1:] if len(row) > ri and row[ri] < k), default=None)
                if i is None:
                    continue
                if ri == len(row) - 1 or i >= row[ri+1]:
                    return False

        # Check the semistandard condition
        return all(not row or row[0] <= i for i, row in enumerate(data, start=1))

    def from_tableau(self, T):
        """
        Return the element of ``self`` corresponding to ``T`` if the result
        is in ``self`` and ``None`` otherwise.
        """
        if isinstance(T, self._ssyt.element_class):
            T = T.to_tableau()

        # Special case of the empty tableau
        if not T:
            if T not in self:
                return None
            return self.element_class(self, [], check=False)

        data = [[None] * ell for ell in self._shape]
        n = self._cartan_type.rank() + 2
        T = ([n - val for val in col] for col in reversed(T.conjugate()))

        for j in range(self._width-1,-1,-1):
            col = next(T)
            for ind, row in enumerate(data):
                if len(row) <= j:
                    continue

                if len(row) == j + 1:  # rightmost entry in the row
                    row[j] = col.pop()
                else:
                    # Necessarily there is an entry to the right
                    for i in range(len(col)-1,-1,-1):
                        if col[i] >= row[j+1]:
                            row[j] = col.pop(i)
                            break

                # Check the semistandard condition
                if row[j] > ind + 1:
                    return None
                if not col:
                    break

        assert data in self
        data = tuple(map(tuple, data))
        return self.element_class(self, data, check=False)

    def shape(self):
        """
        Return the shape of ``self``.
        """
        return self._shape