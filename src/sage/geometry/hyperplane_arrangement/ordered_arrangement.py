r"""
Ordered Hyperplane Arrangements

The :class:`HyperplaneArrangements` orders the hyperplanes in a arrangement
independently of the way the hyperplanes are introduced. The class
:class:`OrderedHyperplaneArrangements` fixes an order specified by
the user. This can be needed for certain properties, e.g., fundamental group with
information about meridians, braid monodromy with information about the strands;
in the future, it may be useful for combinatorial properties.
There are no other differences with usual hyperplane arrangements.

An ordered arrangement is an arrangement where the hyperplanes are sorted
by the user::

    sage: H0.<t0, t1, t2> = HyperplaneArrangements(QQ)
    sage: H0(t0 - t1, t1 - t2, t0 - t2)
    Arrangement <t1 - t2 | t0 - t1 | t0 - t2>
    sage: H.<t0, t1, t2> = OrderedHyperplaneArrangements(QQ)
    sage: H(t0 - t1, t1 - t2, t0 - t2)
    Arrangement <t0 - t1 | t1 - t2 | t0 - t2>

Some methods are adapted, e.g., :meth:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement.hyperplanes`,
and some new ones are created, regarding
hyperplane sections and fundamental groups::

    sage: H.<x,y> = HyperplaneArrangements(QQ)
    sage: H1.<x,y> = OrderedHyperplaneArrangements(QQ)
    sage: A1 = H1(x, y); A = H(A1)
    sage: A.hyperplanes()
    (Hyperplane 0*x + y + 0, Hyperplane x + 0*y + 0)
    sage: A1.hyperplanes()
    (Hyperplane x + 0*y + 0, Hyperplane 0*x + y + 0)

We see the differences in :meth:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement.union`::

    sage: H.<x,y> = HyperplaneArrangements(QQ)
    sage: H1.<x,y> = OrderedHyperplaneArrangements(QQ)
    sage: A = H([1,2,3], [0,1,1], [0,1,-1], [1,-1,0], [1,1,0])
    sage: B = H([1,1,1], [1,-1,1], [1,0,-1])
    sage: C = A.union(B)
    sage: A1 = H1(A); B1 = H1(B); C1 = A1.union(B1)
    sage: [C1.hyperplanes().index(h) for h in C.hyperplanes()]
    [0, 5, 6, 1, 2, 3, 7, 4]

Also in meth:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement.cone`::

    sage: # needs sage.combinat
    sage: a.<x,y,z> = hyperplane_arrangements.semiorder(3)
    sage: H.<x,y,z> = OrderedHyperplaneArrangements(QQ)
    sage: a1 = H(a)
    sage: b = a.cone(); b1 = a1.cone()
    sage: [b1.hyperplanes().index(h) for h in b.hyperplanes()]
    [0, 2, 4, 6, 1, 3, 5]

And in :meth:`~sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement.restriction`::

    sage: # needs sage.graphs
    sage: A.<u, x, y, z> = hyperplane_arrangements.braid(4)
    sage: L.<u, x, y, z> = OrderedHyperplaneArrangements(QQ)
    sage: A1 = L(A)
    sage: H = A[0]; H
    Hyperplane 0*u + 0*x + y - z + 0
    sage: A.restriction(H)
    Arrangement <x - z | u - x | u - z>
    sage: A1.restriction(H)
    Arrangement <x - z | u - x | u - z>
    sage: A1.restriction(H, repetitions=True)
    Arrangement of 5 hyperplanes of dimension 3 and rank 2

AUTHORS:

- Enrique Artal (2023-12): initial version

This module adds some features to the *unordered* one for some
properties which depend on the order.
"""

# *****************************************************************************
#       Copyright (C) 2023 Enrique Artal <artal@unizar.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangementElement
from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
from sage.geometry.hyperplane_arrangement.hyperplane import Hyperplane
from sage.matrix.constructor import matrix, vector
from sage.misc.misc_c import prod
from sage.groups.free_group import FreeGroup
from sage.rings.integer_ring import ZZ
from sage.rings.qqbar import QQbar
from sage.schemes.curves.plane_curve_arrangement import AffinePlaneCurveArrangements
from sage.schemes.curves.plane_curve_arrangement import ProjectivePlaneCurveArrangements


class OrderedHyperplaneArrangementElement(HyperplaneArrangementElement):
    """
    An ordered hyperplane arrangement.

    .. WARNING::

        You should never create
        :class:`OrderedHyperplaneArrangementElement` instances directly,
        always use the parent.
    """
    def __init__(self, parent, hyperplanes, check=True, backend=None):
        """
        Construct an ordered hyperplane arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`OrderedHyperplaneArrangements`

        - ``hyperplanes`` -- tuple of hyperplanes

        - ``check`` -- boolean (default: ``True``); whether
          to check input

        - ``backend`` -- string (default: ``None``); the backend to
          use for the related polyhedral objects

        EXAMPLES::

            sage: H.<x,y> = OrderedHyperplaneArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement <x | y>
            sage: TestSuite(elt).run()
        """
        super().__init__(parent, hyperplanes, check=check, backend=backend)
        self._affine_fundamental_group = None
        self._affine_meridians = None
        self._projective_fundamental_group = None
        self._projective_meridians = None

    def hyperplane_section(self, proj=True):
        r"""
        Compute a generic hyperplane section of ``self``.

        INPUT:

        - ``proj`` -- (default: ``True``) if the
          ambient space is affine or projective

        OUTPUT:

        An arrangement `\mathcal{A}` obtained by intersecting with a
        generic hyperplane

        EXAMPLES::

            sage: L.<x, y, z> = OrderedHyperplaneArrangements(QQ)
            sage: L(x, y - 1, z).hyperplane_section()
            Traceback (most recent call last):
            ...
            TypeError: the arrangement is not projective

            sage: # needs sage.graphs
            sage: A0.<u,x,y,z> = hyperplane_arrangements.braid(4); A0
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
            sage: L.<u,x,y,z> = OrderedHyperplaneArrangements(QQ)
            sage: A = L(A0)
            sage: M = A.matroid()
            sage: A1 = A.hyperplane_section()
            sage: A1
            Arrangement of 6 hyperplanes of dimension 3 and rank 3
            sage: M1 = A1.matroid()
            sage: A2 = A1.hyperplane_section(); A2
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
            sage: M2 = A2.matroid()
            sage: T1 = M1.truncation()
            sage: T1.is_isomorphic(M2)
            True
            sage: T1.isomorphism(M2)
            {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}

            sage: # needs sage.combinat
            sage: a0 = hyperplane_arrangements.semiorder(3); a0
            Arrangement of 6 hyperplanes of dimension 3 and rank 2
            sage: L.<t0, t1, t2> = OrderedHyperplaneArrangements(QQ)
            sage: a = L(a0)
            sage: ca = a.cone()
            sage: m = ca.matroid()
            sage: a1 = a.hyperplane_section(proj=False)
            sage: a1
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
            sage: ca1 = a1.cone()
            sage: m1 = ca1.matroid()
            sage: m.isomorphism(m1)
            {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
            sage: p0 = hyperplane_arrangements.Shi(4)
            sage: L.<t0, t1, t2, t3> = OrderedHyperplaneArrangements(QQ)
            sage: p = L(p0)
            sage: a = p.hyperplane_section(proj=False); a
            Arrangement of 12 hyperplanes of dimension 3 and rank 3
            sage: ca = a.cone()
            sage: m = ca.matroid().truncation()
            sage: a1 = a.hyperplane_section(proj=False); a1
            Arrangement of 12 hyperplanes of dimension 2 and rank 2
            sage: ca1 = a1.cone()
            sage: m1 = ca1.matroid()
            sage: m1.is_isomorphism(m, {j: j for j in range(13)})
            True
        """
        if proj and not self.is_linear():
            raise TypeError('the arrangement is not projective')
        n0 = self.dimension()
        if not proj:
            H = self.cone()
            H1 = H.hyperplane_section()
            mat = matrix(h.coefficients()[1:] for h in H1)
            m = mat.nrows()
            for j in range(mat.ncols()):
                if mat[m - 1, j] != 0:
                    mat.swap_columns(0, j)
                    break
            for j in range(1, mat.ncols()):
                mat.add_multiple_of_column(j, 0, -mat[m - 1, j] / mat[m - 1, 0])
            vrs = H1.parent().variable_names()[1:]
            A1 = OrderedHyperplaneArrangements(self.base_ring(), names=vrs)
            mat_rows = mat.rows()[:-1]
            H1b = A1(mat_rows)
            return H1b
        P = self.intersection_poset(element_label='subspace')
        center = P.maximal_elements()[0].linear_part()
        n1 = center.dimension()
        U = []
        for p in P:
            if p.dimension() == n1 + 1:
                B = [u for u in p.linear_part().basis() if u not in center]
                U.append(B[0])
        # U = [p.linear_part().basis()[0] for p in P if p.dimension() == n1 + 1]
        U0 = sum(U)
        # We look for a linear hyperplane with integer coefficients
        # defining a transversal hyperplane
        for v in ZZ**n0:
            v1 = v + U0
            if 0 not in [w * v1 for w in U]:
                break
        h0 = self.parent()((0,) + tuple(v1))
        H1 = self.add_hyperplane(h0)
        return H1.restriction(h0)

    def affine_fundamental_group(self):
        r"""
        Return the fundamental group of the complement of an affine
        hyperplane arrangement in `\CC^n` whose equations have
        coefficients in a subfield of `\QQbar`.

        OUTPUT: a finitely presented fundamental group

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: A.<x, y> = OrderedHyperplaneArrangements(QQ)
            sage: L = [y + x, y + x - 1]
            sage: H = A(L)
            sage: H.affine_fundamental_group()
            Finitely presented group < x0, x1 |  >
            sage: L = [x, y, x + 1, y + 1, x - y]
            sage: A(L).affine_fundamental_group()
            Finitely presented group
            < x0, x1, x2, x3, x4 | x4*x0*x4^-1*x0^-1,
                                   x0*x2*x3*x2^-1*x0^-1*x3^-1,
                                   x1*x2*x4*x2^-1*x1^-1*x4^-1,
                                   x2*x3*x0*x2^-1*x0^-1*x3^-1,
                                   x2*x4*x1*x2^-1*x1^-1*x4^-1,
                                   x4*x1*x4^-1*x3^-1*x2^-1*x1^-1*x2*x3 >
            sage: H = A(x, y, x + y)
            sage: H.affine_fundamental_group()
            Finitely presented group
            < x0, x1, x2 | x0*x1*x2*x1^-1*x0^-1*x2^-1, x1*x2*x0*x1^-1*x0^-1*x2^-1 >
            sage: H.affine_fundamental_group()  # repeat to use the attribute
            Finitely presented group
            < x0, x1, x2 | x0*x1*x2*x1^-1*x0^-1*x2^-1, x1*x2*x0*x1^-1*x0^-1*x2^-1 >
            sage: T.<t> = QQ[]
            sage: K.<a> = NumberField(t^3 + t + 1)
            sage: L.<x, y> = OrderedHyperplaneArrangements(K)
            sage: H = L(a*x + y -1, x + a*y + 1, x - 1, y - 1)
            sage: H.affine_fundamental_group()
            Traceback (most recent call last):
            ...
            TypeError: the base field is not in QQbar
            sage: L.<t> = OrderedHyperplaneArrangements(QQ)
            sage: L([t - j for j in range(4)]).affine_fundamental_group()
            Finitely presented group < x0, x1, x2, x3 |  >
            sage: L.<x, y, z> = OrderedHyperplaneArrangements(QQ)
            sage: L(L.gens() + (x + y + z + 1,)).affine_fundamental_group().sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3 | x3^-1*x2^-1*x3*x2, x3^-1*x1^-1*x3*x1,
                               x3^-1*x0^-1*x3*x0, x2^-1*x1^-1*x2*x1,
                               x2^-1*x0^-1*x2*x0, x1^-1*x0^-1*x1*x0 >
            sage: A = OrderedHyperplaneArrangements(QQ, names=())
            sage: H = A(); H
            Empty hyperplane arrangement of dimension 0
            sage: H.affine_fundamental_group()
            Finitely presented group <  |  >
        """
        K = self.base_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        if self._affine_fundamental_group:
            return self._affine_fundamental_group
        n = self.dimension()
        r = len(self)
        if n == 0:
            return FreeGroup(0) / []
        if n == 1:
            G = FreeGroup(r) / []
            dic = {j: G.gen(j) for j in range(r)}
            dic[r] = [prod(G.gens()) ** -1]
            self._affine_fundamental_group = G
            self._affine_meridians = dic
            return G
        if n == 2:
            S = self.parent().ambient_space().symmetric_space()
            coord = vector((1,) + S.gens())
            Af = AffinePlaneCurveArrangements(K, names=self.parent().variable_names())
            L = Af([vector(line.coefficients()) * coord for line in self])
            G = L.fundamental_group()
            self._affine_fundamental_group = G
            self._affine_meridians = L._meridians_simpl_vertical
            return G
        H1 = self.hyperplane_section(proj=False)
        G = H1.affine_fundamental_group()
        self._affine_fundamental_group = G
        self._affine_meridians = H1._affine_meridians
        return G

    def affine_meridians(self):
        r"""
        Return the meridians of each hyperplane (including the one at infinity).

        OUTPUT: a dictionary

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: A.<x, y> = OrderedHyperplaneArrangements(QQ)
            sage: L = [y + x, y + x - 1]
            sage: H = A(L)
            sage: g = H.affine_fundamental_group()
            sage: g
            Finitely presented group < x0, x1 |  >
            sage: H.affine_meridians()
            {0: [x0], 1: [x1], 2: [x1^-1*x0^-1]}
            sage: H1 = H.add_hyperplane(y - x)
            sage: H1.affine_meridians()
            {0: [x0], 1: [x1], 2: [x2], 3: [x2^-1*x1^-1*x0^-1]}
        """
        if self._affine_meridians is None:
            self.affine_fundamental_group()
        return dict(self._affine_meridians)

    def projective_fundamental_group(self):
        r"""
        Return the fundamental group of the complement of a projective
        hyperplane arrangement.

        OUTPUT:

        The finitely presented group of the complement
        in the projective space whose equations have
        coefficients in a subfield of `\QQbar`.

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: A.<x, y> = OrderedHyperplaneArrangements(QQ)
            sage: H = A(x, y, x + y)
            sage: H.projective_fundamental_group()
            Finitely presented group < x0, x1 |  >

            sage: # needs sirocco sage.graphs
            sage: A3.<x, y, z> = OrderedHyperplaneArrangements(QQ)
            sage: H = A3(hyperplane_arrangements.braid(4).essentialization())
            sage: G3 = H.projective_fundamental_group(); G3.sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3, x4 | x4^-1*x3^-1*x2^-1*x3*x4*x0*x2*x0^-1,
                                   x4^-1*x2^-1*x4*x2, x4^-1*x1^-1*x0^-1*x1*x4*x0,
                                   x4^-1*x1^-1*x0^-1*x4*x0*x1,
                                   x4^-1*x1^-1*x3*x0*x1*x3^-1*x2^-1*x4*x0^-1*x2,
                                   x3^-1*x2^-1*x1^-1*x0^-1*x3*x0*x1*x2,
                                   x3^-1*x1^-1*x3*x1 >
            sage: G3.abelian_invariants()
            (0, 0, 0, 0, 0)
            sage: A4.<t1, t2, t3, t4> = OrderedHyperplaneArrangements(QQ)
            sage: H = A4(hyperplane_arrangements.braid(4))
            sage: G4 = H.projective_fundamental_group(); G4.sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3, x4 | x4^-1*x3^-1*x2^-1*x3*x4*x0*x2*x0^-1,
                                   x4^-1*x2^-1*x4*x2, x4^-1*x1^-1*x0^-1*x1*x4*x0,
                                   x4^-1*x1^-1*x0^-1*x4*x0*x1,
                                   x4^-1*x1^-1*x3*x0*x1*x3^-1*x2^-1*x4*x0^-1*x2,
                                   x3^-1*x2^-1*x1^-1*x0^-1*x3*x0*x1*x2,
                                   x3^-1*x1^-1*x3*x1 >
            sage: G4.abelian_invariants()
            (0, 0, 0, 0, 0)

            sage: # needs sirocco
            sage: L.<t0, t1, t2, t3, t4> = OrderedHyperplaneArrangements(QQ)
            sage: H = hyperplane_arrangements.coordinate(5)
            sage: H = L(H)
            sage: g = H.projective_fundamental_group()
            sage: g.is_abelian(), g.abelian_invariants()
            (True, (0, 0, 0, 0))
            sage: L(t0, t1, t2, t3, t4, t0 - 1).projective_fundamental_group()
            Traceback (most recent call last):
            ...
            TypeError: the arrangement is not projective
            sage: T.<t> = QQ[]
            sage: K.<a> = NumberField(t^3 + t + 1)
            sage: L.<x, y, z> = OrderedHyperplaneArrangements(K)
            sage: H = L(a*x + y - z, x + a*y + z, x - z, y - z)
            sage: H.projective_fundamental_group()
            Traceback (most recent call last):
            ...
            TypeError: the base field is not in QQbar
            sage: A.<x> = OrderedHyperplaneArrangements(QQ)
            sage: H = A(); H
            Empty hyperplane arrangement of dimension 1
            sage: H.projective_fundamental_group()
            Finitely presented group <  |  >
        """
        K = self.base_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        if not self.is_linear():
            raise TypeError('the arrangement is not projective')
        if self._projective_fundamental_group:
            return self._projective_fundamental_group
        n = self.dimension()
        r = len(self)
        if n == 1:
            return FreeGroup(0) / []
        if n == 2:
            G = FreeGroup(r - 1) / []
            dic = {j: G.gen(j) for j in range(r - 1)}
            dic[r - 1] = [prod(G.gens()) ** -1]
            self._projective_fundamental_group = G
            self._projective_meridians = dic
            return G
        if n == 3:
            S = self.parent().ambient_space().symmetric_space()
            coord = vector(S.gens())
            Proj = ProjectivePlaneCurveArrangements(K,
                                                    names=self.parent().variable_names())
            L = Proj([vector(line.coefficients()[1:]) * coord for line in self])
            G = L.fundamental_group()
            self._projective_fundamental_group = G
            self._projective_meridians = L._meridians_simpl
            return G
        H1 = self.hyperplane_section()
        G = H1.projective_fundamental_group()
        self._projective_fundamental_group = G
        self._projective_meridians = H1._projective_meridians
        return G

    def projective_meridians(self):
        r"""
        Return the meridian of each hyperplane.

        OUTPUT: a dictionary

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: A.<x, y> = OrderedHyperplaneArrangements(QQ)
            sage: H = A(x, y, x + y)
            sage: H.projective_meridians()
            {0: x0, 1: x1, 2: [x1^-1*x0^-1]}

            sage: # needs sirocco sage.graphs
            sage: A3.<x, y, z> = OrderedHyperplaneArrangements(QQ)
            sage: H = A3(hyperplane_arrangements.braid(4).essentialization())
            sage: H.projective_meridians()
            {0: [x2^-1*x0^-1*x4^-1*x3^-1*x1^-1],
             1: [x3], 2: [x4], 3: [x1], 4: [x2], 5: [x0]}
            sage: A4.<t1, t2, t3, t4> = OrderedHyperplaneArrangements(QQ)
            sage: H = A4(hyperplane_arrangements.braid(4))
            sage: H.projective_meridians()
            {0: [x2^-1*x0^-1*x4^-1*x3^-1*x1^-1], 1: [x3],
             2: [x4], 3: [x0], 4: [x2], 5: [x1]}

            sage: # needs sirocco
            sage: L.<t0, t1, t2, t3, t4> = OrderedHyperplaneArrangements(QQ)
            sage: H = hyperplane_arrangements.coordinate(5)
            sage: H = L(H)
            sage: H.projective_meridians()
            {0: [x2], 1: [x3], 2: [x0], 3: [x3^-1*x2^-1*x1^-1*x0^-1], 4: [x1]}
        """
        if self._projective_meridians is None:
            self.projective_fundamental_group()
        return dict(self._projective_meridians)


class OrderedHyperplaneArrangements(HyperplaneArrangements):
    """
    Ordered Hyperplane arrangements.

    For more information on hyperplane arrangements, see
    :mod:`sage.geometry.hyperplane_arrangement.arrangement`.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::

        sage: H.<x,y> = HyperplaneArrangements(QQ)
        sage: x
        Hyperplane x + 0*y + 0
        sage: x + y
        Hyperplane x + y + 0
        sage: H(x, y, x-1, y-1)
        Arrangement <y - 1 | y | x - 1 | x>
    """
    Element = OrderedHyperplaneArrangementElement

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        INPUT:

        - ``*args`` -- positional arguments, each defining a
          hyperplane; alternatively, a single polytope or a single
          hyperplane arrangement

        - ``signed`` -- boolean (default: ``True``); whether to
          preserve signs of hyperplane equations

        - ``check`` -- boolean (default: ``True``); whether to
          perform argument checking

        EXAMPLES::

            sage: L.<x, y> = OrderedHyperplaneArrangements(QQ)
            sage: L(x)
            Arrangement <x>
            sage: L(x, y)
            Arrangement <x | y>
            sage: L([x, y])
            Arrangement <x | y>
            sage: L([0, 1, 0], [0, 0, 1])
            Arrangement <x | y>
            sage: L([[0, 0, 1], [0, 1, 0]])
            Arrangement <y | x>

            sage: L(polytopes.hypercube(2))
            Arrangement <-x + 1 | -y + 1 | x + 1 | y + 1>

            sage: L(-x, x + y - 1, signed=False)
            Arrangement <x | -x - y + 1>

        TESTS::

            sage: L()
            Empty hyperplane arrangement of dimension 2
            sage: L(0)        # zero is equivalent to no argument, Issue #8648
            Empty hyperplane arrangement of dimension 2
            sage: L(0*x)      # degenerate hyperplane is NOT allowed
            Traceback (most recent call last):
            ...
            ValueError: linear expression must be non-constant to define a hyperplane
            sage: L(0*x, y)   # ditto
            Traceback (most recent call last):
            ...
            ValueError: linear expression must be non-constant to define a hyperplane
        """
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, HyperplaneArrangementElement) and arg.parent() is self:
                # optimization if argument is already a hyperplane arrangement
                return arg
            if arg == 0 and not isinstance(arg, Hyperplane):
                # zero = neutral element under addition = the empty hyperplane arrangement
                args = []
        # process keyword arguments
        not_char2 = (self.base_ring().characteristic() != 2)
        signed = kwds.pop('signed', not_char2)
        check = kwds.pop('check', True)
        backend = kwds.pop('backend', None)
        if kwds:
            raise ValueError('unknown keyword argument')
        # process positional arguments
        AA = self.ambient_space()
        try:
            hyperplanes = [AA(a) for a in args]
        except (TypeError, ValueError, AttributeError):
            if len(args) > 1:
                raise
            arg = args[0]
            if hasattr(arg, 'Hrepresentation'):
                hyperplanes = [AA(h) for h in arg.Hrepresentation()]
            else:
                hyperplanes = [AA(a) for a in arg]
        hyperplanes = [h.primitive(signed) for h in hyperplanes]
        if check:
            if signed and not not_char2:
                raise ValueError('cannot be signed in characteristic 2')
            for h in hyperplanes:
                if h.A() == 0:
                    raise ValueError('linear expression must be non-constant to define a hyperplane')
                if not_char2 and -h in hyperplanes:
                    raise ValueError('arrangement cannot simultaneously have h and -h as hyperplane')
        return self.element_class(self, tuple(hyperplanes), backend=backend)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: L.<x, y> = OrderedHyperplaneArrangements(QQ);  L
            Ordered hyperplane arrangements in 2-dimensional linear space
            over Rational Field with coordinates x, y
        """
        return 'Ordered hyperplane arrangements in {0}'.format(self.ambient_space())
