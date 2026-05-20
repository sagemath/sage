from sage.structure.parent import Parent
from sage.structure.element import Element, parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.rings.integer_ring import ZZ


class CylindricShape(Element):
    r"""
    An element representing a cylindric shape of period `(d, L)`.

    A cylindric shape `\lambda` is encoded by a length-`d` tuple
    ``values``, and is extended to all integer indices via
    `\lambda_{i+d} = \lambda_i - L`.

    EXAMPLES::

        sage: from sage.combinat.cylindric_shapes import CylindricShapes
        sage: CS = CylindricShapes(3, 2)
        sage: s = CS((2, 1, 0)); s
        [2, 1, 0]

        sage: [s[i] for i in range(-2, 5)]
        [3, 2, 2, 1, 0, 0, -1]

    """
    def __init__(self, parent, values):
        """
        Construct a cylindric shape.

        INPUT:

        - ``parent`` -- a :class:`CylindricShapes` instance of period `(d,L)`
        - ``values`` -- an iterable of length ``parent._d`` giving the profile
          on a fundamental domain

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShape, CylindricShapes
            sage: CS = CylindricShapes(2, 5)
            sage: CS((2, 1)).parent() is CS
            True
        """
        self._values = values
        Element.__init__(self, parent)

    def _repr_(self):
        """
        Return a string representation showing the fundamental-domain values.

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(4, 3)
            sage: CS((3,2,1,0))
            [3, 2, 1, 0]
        """
        return str(list(self._values))

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        We use bars to indicate negative values. If all letters are
        single-digit letters, we omit commas.

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(4, 11)
            sage: latex(CS((2,0,-1,-1)))
            20\bar{1}\bar{1}

            sage: latex(CS((10,0,-1,-1)))
            10,0,\bar{1},\bar{1}
        """
        if all(abs(v) < 10 for v in self._values):
            sep = ""
        else:
            sep = ","
        return sep.join(r"\bar{" + str(-v) + "}" if v < 0 else str(v)
                        for v in self._values)

    def __getitem__(self, i):
        r"""
        Return the value at index ``i``, accounting for periodicity.

        A cylindric shape `\lambda` in `C_{d, L}` satisfies
        `\lambda_{i+d} = \lambda_i - L`.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: s = CS((2, 1, 0))
            sage: s[0], s[1], s[2]
            (2, 1, 0)

        One period to the right subtracts `L`::

            sage: s[3]
            0

        One period to the right adds `L`::

            sage: s[-3]
            4
        """
        P = self.parent()
        q, r = divmod(i, P._d)
        return self._values[r] - q * P._L

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: a = CS((1, 0, 0))
            sage: a == CS((1, 0, 0))
            True
            sage: CT = CylindricShapes(3, 5)
            sage: a == CT((1, 0, 0))
            False
        """
        return self.parent() is parent(other) and self._values == other._values

    def __ne__(self, other):
        """
        Check inequality.

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: a = CS((1, 0, 0)); b = CS((0, 0, 0))
            sage: a != b
            True
            sage: a != CS((1, 0, 0))
            False
        """
        return not (self == other)

    def __le__(self, other):
        r"""
        Return whether ``self[i] <= other[i]`` for all ``i``.

        INPUT:

        - ``other`` -- a :class:`CylindricShape` in the same parent

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: CS((1, 1, 1)) <= CS((2, 1, 1))
            True
            sage: CS((1, 1, 1)) <= CS((2, 1, 0))
            False
        """
        P = self.parent()
        if P is not other.parent():
            raise TypeError("the parents of the cylindric shapes must coincide")

        return all(self._values[i] <= other._values[i] for i in range(self.parent()._d))

    def __lt__(self, other):
        r"""
        Return whether ``self[i] <= other[i]`` for all ``i``.

        INPUT:

        - ``other`` -- a :class:`CylindricShape` in the same parent

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: CS((1, 1, 1)) < CS((2, 1, 1))
            True
            sage: CS((1, 1, 1)) < CS((1, 1, 1))
            False
        """
        return self != other and self <= other

    def __ge__(self, other):
        r"""
        Return whether ``self[i] >= other[i]`` for all ``i``.

        INPUT:

        - ``other`` -- a :class:`CylindricShape` in the same parent

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: CS((2, 2, 1)) >= CS((2, 1, 1))
            True
            sage: CS((2, 2, 1)) >= CS((2, 2, 2))
            False
        """
        P = self.parent()
        if P is not other.parent():
            raise TypeError("the parents of the cylindric shapes must coincide")

        return all(self._values[i] >= other._values[i] for i in range(self.parent()._d))

    def __hash__(self):
        """
        Return a hash.

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: hash(CS((1, 0, 0))) == hash(CS((1, 1, 0)))
            False
        """
        return hash(self._values)

    def union(self, other):
        r"""
        Return the join (component-wise maximum) of two shapes.

        This is the shape with fundamental-domain values
        max(self.values, other.values) taken component-wise.

        INPUT:

        - ``other`` -- another :class:`CylindricShape`

        OUTPUT:

        - a :class:`CylindricShape`

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: a = CS((2, 2, 1)); b = CS((3, 1, 1))
            sage: a.union(b)
            [3, 2, 1]
        """
        P = self.parent()
        if P is not other.parent():
            raise TypeError("the parents of the cylindric shapes must coincide")

        new_vals = tuple(max(self._values[i], other._values[i])
                         for i in range(P._d))
        return P(new_vals)

    def intersection(self, other):
        r"""
        Return the component-wise minimum of ``self`` and
        ``other``.

        INPUT:

        - ``other`` -- another :class:`CylindricShape`

        OUTPUT:

        - a :class:`CylindricShape`

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: a = CS((3, 2, 1)); b = CS((1, 1, 1))
            sage: a.intersection(b)
            [1, 1, 1]
        """
        P = self.parent()
        if P is not other.parent():
            raise TypeError("the parents of the cylindric shapes must coincide")

        new_vals = tuple(min(self._values[i], other._values[i])
                         for i in range(P._d))
        return P(new_vals)

    def add_cell(self, i):
        r"""
        Return the shape obtained by adding one cell at residue index ``i``.

        Internally, we increment the component at ``i mod d`` by 1.

        INPUT:

        - ``i`` -- an integer; only the residue class modulo d is used

        OUTPUT:

        - a :class:`CylindricShape`

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: s = CS((1, 0, 0))
            sage: s.add_cell(1)
            [1, 1, 0]
            sage: s.add_cell(3)
            [2, 0, 0]
        """
        idx = i % self.parent()._d
        new_vals = list(self._values)
        new_vals[idx] += 1
        return self.parent()(tuple(new_vals))

    def remove_cell(self, i):
        r"""
        Return the shape obtained by removing one cell at residue index ``i``.

        Internally, we decrement the component at ``i mod d`` by 1.

        INPUT:

        - ``i`` -- an integer; only the residue class modulo d is used

        OUTPUT:

        - a :class:`CylindricShape`

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: s = CS((0, 0, -1))
            sage: s.remove_cell(1)
            [0, -1, -1]
            sage: s.remove_cell(2)
            [0, 0, -2]
        """
        idx = i % self.parent()._d
        new_vals = list(self._values)
        new_vals[idx] -= 1
        return self.parent()(tuple(new_vals))

    def complement(self):
        r"""
        Return the complement shape in the same parent.

        If ``self`` has fundamental-domain values `(\mu_1,\dots,
        \mu_d)`, the complement has values `(L - \mu_d,\dots, L -
        \mu_1)` in the same parent (d, L).

        OUTPUT:

        - a :class:`CylindricShape` in the same parent

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: s = CS((2, 1, 0))
            sage: s.complement()
            [2, 1, 0]
            sage: t = CS((2, 2, 1))
            sage: t.complement()
            [1, 0, 0]

        Complement is an involution::

            sage: (t.complement()).complement() == t
            True
        """
        P = self.parent()
        L = P._L
        return P(L - v for v in reversed(self._values))

    def conjugate(self):
        r"""Return the conjugate shape in the parent with `d` and `L`
        swapped.

        If ``self`` has fundamental-domain values `(\mu_1,\dots,
        \mu_d)`, the conjugate shape `\mu'` is in the parent with `L`
        and `d` swapped and has values defined by
        .. MATH::

            \mu'_j = \max\{ i \in \{1,\dots,d\} : \mu_i \ge j \}

        Equivalently, because `\mu` is weakly decreasing, `\mu'_j` is
        the number of parts of `\mu` that are at least `j`.

        OUTPUT:

        - a :class:`CylindricShape` in :class:`CylindricShapes(L, d)`

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: s = CS((2, 1, 0))
            sage: sc = s.conjugate()
            sage: sc, sc.parent()
            ([2, 1], Cylindric Shapes of period (2, 3))

        Conjugation is an involution::

            sage: s2 = sc.conjugate()
            sage: s2, s2.parent()
            ([2, 1, 0], Cylindric Shapes of period (3, 2))
            sage: s2 == s
            True
        """
        P = self.parent()
        d = P._d
        L = P._L
        new_P = CylindricShapes(L, d)
        vals = list(self._values)
        # Compute counts: for each j=1..L, how many vals are >= j
        new_vals = []
        for j in range(1, L + 1):
            cnt = 0
            for v in vals:
                if v >= j:
                    cnt += 1
                else:
                    # vals is weakly decreasing, so we can break
                    break
            new_vals.append(cnt)
        return new_P(new_vals)


class CylindricShapes(UniqueRepresentation, Parent):
    r"""
    The set of all cylindric shapes of period `(d, L)`.

    Elements are length-`d` integer tuples, interpreted on all of
    `\mathbb Z` via `\lambda_{i+d} = \lambda_i - L`.

    INPUT:

    - ``d`` -- positive integer, the fundamental-domain width
    - ``L`` -- integer, the decrement across one period

    EXAMPLES::

        sage: from sage.combinat.cylindric_shapes import CylindricShapes
        sage: CS = CylindricShapes(3, 2)
        sage: CS
        Cylindric Shapes of period (3, 2)
        sage: s = CS((2, 1, 0)); s
        [2, 1, 0]
        sage: s.parent() is CS
        True
    """
    def __init__(self, d, L):
        """
        Initialize the parent of cylindric shapes with period `(d, L)`.

        INPUT:

        - ``d`` -- positive integer
        - ``L`` -- integer

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CylindricShapes(1, 0)
            Cylindric Shapes of period (1, 0)
        """
        if d <= 0:
            raise ValueError(f"d must be a positive integer (got {d})")
        self._d = int(d)
        self._L = int(L)
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, values):
        r"""
        Construct an element in this parent from a tuple or list
        of length `d`.

        The entries are coerced to ZZ and must satisfy the
        cylindric-shape constraints `\lambda_d +
        L\ge\lambda_1\ge\lambda_2\ge\cdots\ge\lambda_d`.

        INPUT:

        - ``values`` -- iterable of length ``self._d``

        OUTPUT:

        - a :class:`CylindricShape` with these fundamental-domain values

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(3, 2)
            sage: CS((2, 1, 1))
            [2, 1, 1]

        The entries must be weakly decreasing::

            sage: CS((1, 2, 0))
            Traceback (most recent call last):
            ...
            ValueError: values (1, 2, 0) should be weakly decreasing

        The wrap-around constraint must hold::

            sage: CT = CylindricShapes(3, 1)
            sage: CT((3, 2, 0))   # 0 + L = 1 < 3
            Traceback (most recent call last):
            ...
            ValueError: wrap-around constraint violated: last + L < first for (3, 2, 0) with L=1
        """
        vals = tuple(ZZ(x) for x in values)
        if len(vals) != self._d:
            raise ValueError(f"{vals} should be an iterable of length-{self._d}")

        # Monotonicity: vals[0] >= ... >= vals[d-1]
        if any(vals[i] < vals[i + 1] for i in range(self._d - 1)):
            raise ValueError(f"values {vals} should be weakly decreasing")

        # Wrap-around: vals[-1] + L >= vals[0]
        if vals[-1] + self._L < vals[0]:
            raise ValueError(f"wrap-around constraint violated: last + L < first for {vals} with L={self._L}")

        return self.Element(self, vals)

    def zero(self):
        """
        Return the empty cylindric shape (all zero entries on the fundamental domain).

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CS = CylindricShapes(4, 3)
            sage: CS.zero()
            [0, 0, 0, 0]
        """
        return self((0,) * self._d)

    def _repr_(self):
        """
        String representation of this parent.

        EXAMPLES::

            sage: from sage.combinat.cylindric_shapes import CylindricShapes
            sage: CylindricShapes(5, 1)
            Cylindric Shapes of period (5, 1)
        """
        return f"Cylindric Shapes of period ({self._d}, {self._L})"

    Element = CylindricShape
