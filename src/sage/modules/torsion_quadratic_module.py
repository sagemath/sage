r"""
Finite `\ZZ`-modules with with bilinear and quadratic forms.

AUTHORS:

- Simon Brandhorst (2017-09): First created
"""

# ****************************************************************************
#       Copyright (C) 2017 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.fg_pid.fgp_module import FGP_Module_class
from sage.modules.fg_pid.fgp_element import FGP_Element
from sage.modules.free_quadratic_module import FreeQuadraticModule
from sage.arith.misc import gcd
from sage.rings.all import ZZ, Zp, QQ, IntegerModRing
from sage.groups.additive_abelian.qmodnz import QmodnZ
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.integer_mod import mod
from sage.arith.misc import legendre_symbol
from sage.structure.unique_representation import CachedRepresentation

def TorsionQuadraticForm(q):
    r"""
    Create a torsion quadratic form module from a rational matrix.

    The resulting quadratic form takes values in `\QQ / \ZZ`
    or `\QQ / 2 \ZZ` (depending on ``q``).
    If it takes values modulo `2`, then it is non-degenerate.
    In any case the bilinear form is non-degenerate.

    INPUT:

    - ``q`` -- a symmetric rational matrix

    EXAMPLES::

        sage: q1 = Matrix(QQ,2,[1,1/2,1/2,1])
        sage: TorsionQuadraticForm(q1)
        Finite quadratic module over Integer Ring with invariants (2, 2)
        Gram matrix of the quadratic form with values in Q/2Z:
        [  1 1/2]
        [1/2   1]

    In the following example the quadratic form is degenerate.
    But the bilinear form is still non-degenerate::

        sage: q2 = diagonal_matrix(QQ,[1/4,1/3])
        sage: TorsionQuadraticForm(q2)
        Finite quadratic module over Integer Ring with invariants (12,)
        Gram matrix of the quadratic form with values in Q/Z:
        [7/12]

    TESTS::

        sage: TorsionQuadraticForm(matrix.diagonal([3/4,3/8,3/8]))
        Finite quadratic module over Integer Ring with invariants (4, 8, 8)
        Gram matrix of the quadratic form with values in Q/2Z:
        [3/4   0   0]
        [  0 3/8   0]
        [  0   0 3/8]
    """
    q = matrix(QQ, q)
    if q.nrows() != q.ncols():
        raise ValueError("the input must be a square matrix")
    if q != q.transpose():
        raise ValueError("the input must be a symmetric matrix")

    Q, d = q._clear_denom()
    S, U, V = Q.smith_form()
    D = U * q * V
    Q = FreeQuadraticModule(ZZ, q.ncols(), inner_product_matrix=d**2 * q)
    denoms = [D[i, i].denominator() for i in range(D.ncols())]
    rels = Q.span(diagonal_matrix(ZZ, denoms) * U)
    return TorsionQuadraticModule((1/d)*Q, (1/d)*rels, modulus=1)


class TorsionQuadraticModuleElement(FGP_Element):
    r"""
    An element of a torsion quadratic module.

    INPUT:

    - ``parent`` -- parent

    - ``x`` -- element of ``parent.V()``

    - ``check`` -- bool (default: ``True``)

    TESTS::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: T = TorsionQuadraticModule(ZZ^3, 6*ZZ^3)
        sage: loads(dumps(T)) == T
        True
        sage: t = T.gen(0)
        sage: loads(dumps(t)) == t
        True

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = span([[1/2,1,1], [3/2,2,1], [0,0,1]], ZZ)
        sage: b = V.basis()
        sage: W = V.span([2*b[0]+4*b[1], 9*b[0]+12*b[1], 4*b[2]])
        sage: Q = TorsionQuadraticModule(V, W)
        sage: x = Q(b[0] - b[1])
        sage: TestSuite(x).run()
        """

    def _mul_(self, other):
        r"""
        Compute the inner product of two elements.

        OUTPUT:

        - an element of `\QQ / m\ZZ` with `m\ZZ = (V, W)`

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = (1/2)*ZZ^2; W = ZZ^2
            sage: T = TorsionQuadraticModule(V, W)
            sage: g = T.gens()
            sage: x = g[0]
            sage: y = g[0] + g[1]
            sage: x
            (1, 0)
            sage: x*y
            1/4

        The inner product has further aliases::

            sage: x.inner_product(y)
            1/4
            sage: x.b(y)
            1/4
        """
        value_module = self.parent().value_module()
        return value_module( self.lift().inner_product(other.lift()) )

    inner_product = _mul_
    b = _mul_

    def quadratic_product(self):
        r"""
        Compute the quadratic_product of ``self``.

        OUTPUT:

        - an element of `\QQ / n\ZZ` where `n\ZZ = 2(V,W) +
          \ZZ \{ (w,w) | w \in W \}`

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: W = FreeQuadraticModule(ZZ, 2, 2*matrix.identity(2))
            sage: V = (1/2) * W
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.gen(0)
            sage: x
            (1, 0)
            sage: x.quadratic_product()
            1/2
            sage: x.quadratic_product().parent()
            Q/2Z
            sage: x*x
            1/2
            sage: (x*x).parent()
            Q/Z
        """
        value_module_qf = self.parent().value_module_qf()
        lift = self.lift()
        return value_module_qf(lift.inner_product(lift))

    q = quadratic_product


class TorsionQuadraticModule(FGP_Module_class, CachedRepresentation):
    r"""
    Finite quotients with a bilinear and a quadratic form.

    Let `V` be a symmetric FreeQuadraticModule and `W \subseteq V` a
    submodule of the same rank as `V`. The quotient `V / W` is a torsion
    quadratic module. It inherits a bilinear form `b` and a quadratic
    form `q`.

    `b: V \times V \to \QQ / m\ZZ`, where  `m\ZZ = (V,W)`
    and `b(x,y) = (x,y) + m\ZZ`

    `q: V \to \QQ / n\ZZ`, where `n\ZZ = 2(V,W) + \ZZ \{ (w,w) | w \in W \}`

    INPUT:

    - ``V`` -- a :class:`FreeModule` with a symmetric inner product matrix

    - ``W`` -- a submodule of ``V`` of the same rank as ``V``

    - ``check`` -- bool (default: ``True``)

    - ``modulus`` -- a rational number dividing `m` (default: `m`);
      the inner product `b` is defined in `\QQ /` ``modulus`` `\ZZ`

    - ``modulus_qf`` -- a rational number dividing `n` (default: `n`);
      the quadratic form `q` is defined in `\QQ /` ``modulus_qf`` `\ZZ`

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = FreeModule(ZZ, 3)
        sage: T = TorsionQuadraticModule(V, 5*V)
        sage: T
        Finite quadratic module over Integer Ring with invariants (5, 5, 5)
        Gram matrix of the quadratic form with values in Q/5Z:
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    Element = TorsionQuadraticModuleElement

    @staticmethod
    def __classcall__(cls, V, W, gens=None, modulus=None, modulus_qf=None, check=True):
        r"""
        Return a :class:`TorsionQuadraticModule`.

        This method does the preprocessing for :meth:`sage.structure.CachedRepresentation`.

        TESTS::

            sage: q = matrix([1/2])
            sage: D1 = TorsionQuadraticForm(q)
            sage: D2 = TorsionQuadraticForm(q)
            sage: D1 is D2
            True
        """
        if check:
            if V.rank() != W.rank():
                raise ValueError("modules must be of the same rank")
            if V.base_ring() is not ZZ:
                raise NotImplementedError("only currently implemented over ZZ")
            if V.inner_product_matrix() != V.inner_product_matrix().transpose():
                raise ValueError("the cover must have a symmetric inner product")

            if gens is not None and V.span(gens) + W != V:
                raise ValueError("provided gens do not generate the quotient")

        # compute the modulus - this may be expensive
        if modulus is None or check:
            # The inner product of two elements `b(v1+W,v2+W)`
            # is defined `mod (V,W)`
            num = V.basis_matrix() * V.inner_product_matrix() * W.basis_matrix().T
            max_modulus = gcd(num.list())

        if modulus is None:
            modulus = max_modulus
        elif check and max_modulus / modulus not in V.base_ring():
            raise ValueError("the modulus must divide (V, W)")

        if modulus_qf is None or check:
            # The quadratic_product of an element `q(v+W)` is defined
            # `\mod 2(V,W) + ZZ\{ (w,w) | w in w\}`
            norm = gcd(W.gram_matrix().diagonal())
            max_modulus_qf = gcd(norm, 2 * modulus)

        if modulus_qf is None:
            modulus_qf = max_modulus_qf
        elif check and max_modulus_qf / modulus_qf not in V.base_ring():
            raise ValueError("the modulus_qf must divide (V, W)")
        return super(TorsionQuadraticModule, cls).__classcall__(cls, V, W, gens, modulus, modulus_qf)

    def __init__(self, V, W, gens, modulus, modulus_qf):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule(ZZ^3, 6*ZZ^3)
            sage: TestSuite(T).run()
        """

        FGP_Module_class.__init__(self, V, W, check=True)
        if gens is not None:
            self._gens_user = tuple(self(g) for g in gens)
        else:
            # this is taken care of in the .gens method
            # we do not want this at initialization
            self._gens_user = None
        self._modulus = modulus
        self._modulus_qf = modulus_qf


    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ,3)
            sage: T = TorsionQuadraticModule(V, 5*V,modulus=1)
            sage: T
            Finite quadratic module over Integer Ring with invariants (5, 5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        return ( "Finite quadratic module over %s with invariants %s\n"
                 % (self.base_ring(), self.invariants()) +
                 "Gram matrix of the quadratic form with values in %r:\n%r"
                 % (self.value_module_qf(), self.gram_matrix_quadratic()))

    def _module_constructor(self, V, W, check=False):
        r"""
        Construct a torsion quadratic module ``V / W``.

        INPUT:

        - ``V`` -- an module

        - ``W`` -- an submodule of ``V`` over the same base ring

        - ``check`` -- bool (default: ``False``);

          * if ``False``, then the value modulus is inherited from ``self``
          * if ``True``, it figures it out on its own. But that is expensive

        OUTPUT:

        The quotient ``V / W`` as a :class:`TorsionQuadraticModule`.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = span([[1/2,1,1], [3/2,2,1], [0,0,1]], ZZ)
            sage: b = V.basis()
            sage: W = V.span([2*b[0]+4*b[1], 9*b[0]+12*b[1], 4*b[2]])
            sage: Q = TorsionQuadraticModule(V, W); Q
            Finite quadratic module over Integer Ring with invariants (4, 12)
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
            sage: Q._module_constructor(V,W)
            Finite quadratic module over Integer Ring with invariants (4, 12)
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
        """
        if check:
            # figuring out the modulus can be expensive
            return TorsionQuadraticModule(V, W, check=check)
        else:
            return TorsionQuadraticModule(V, W, check=check,
                                          modulus=self._modulus,
                                          modulus_qf=self._modulus_qf)

    def all_submodules(self):
        r"""
        Return a list of all submodules of ``self``.

        .. WARNING::

            This method creates all submodules in memory. The number of submodules
            grows rapidly with the number of generators. For example consider a
            vector space of dimension `n` over a finite field of prime order `p`.
            The number of subspaces is (very) roughly `p^{(n^2-n)/2}`.

        EXAMPLES::

            sage: D = IntegralLattice("D4").discriminant_group()
            sage: D.all_submodules()
            [Finite quadratic module over Integer Ring with invariants ()
              Gram matrix of the quadratic form with values in Q/2Z:
              [],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2,)
              Gram matrix of the quadratic form with values in Q/2Z:
              [1],
             Finite quadratic module over Integer Ring with invariants (2, 2)
              Gram matrix of the quadratic form with values in Q/2Z:
              [  1 1/2]
              [1/2   1]]
        """
        from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        invs = self.invariants()
        # knows how to compute all subgroups
        A = AbelianGroupGap(invs)
        S = A.all_subgroups()
        # over ZZ submodules and subgroups are the same thing.
        submodules = []
        for sub in S:
            gen = [A(g).exponents() for g in sub.gens()]
            gen = [self.linear_combination_of_smith_form_gens(g) for g in gen]
            submodules.append(self.submodule(gen))
        return submodules

    @cached_method
    def brown_invariant(self):
        r"""
        Return the Brown invariant of this torsion quadratic form.

        Let `(D,q)` be a torsion quadratic module with values in `\QQ / 2 \ZZ`.
        The Brown invariant `Br(D,q) \in \Zmod{8}` is defined by the equation

        .. MATH::

            \exp \left( \frac{2 \pi i }{8} Br(q)\right) =
            \frac{1}{\sqrt{D}} \sum_{x \in D} \exp(i \pi q(x)).

        The Brown invariant is additive with respect to direct sums of
        torsion quadratic modules.

        OUTPUT:

        - an element of `\Zmod{8}`

        EXAMPLES::

            sage: L = IntegralLattice("D4")
            sage: D = L.discriminant_group()
            sage: D.brown_invariant()
            4

        We require the quadratic form to be defined modulo `2 \ZZ`::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3))
            sage: T = TorsionQuadraticModule((1/10)*V, V)
            sage: T.brown_invariant()
            Traceback (most recent call last):
            ...
            ValueError: the torsion quadratic form must have values in QQ / 2 ZZ
        """
        if self._modulus_qf != 2:
            raise ValueError("the torsion quadratic form must have values in "
                             "QQ / 2 ZZ")
        from sage.quadratic_forms.genera.normal_form import collect_small_blocks
        brown = IntegerModRing(8).zero()
        for p in self.annihilator().gen().prime_divisors():
            q = self.primary_part(p).normal_form()
            q = q.gram_matrix_quadratic()
            L = collect_small_blocks(q)
            for qi in L:
                brown += _brown_indecomposable(qi, p)
        return brown

    @cached_method
    def gram_matrix_bilinear(self):
        r"""
        Return the Gram matrix with respect to the generators.

        OUTPUT:

        A rational matrix ``G`` with ``G[i,j]`` given by the inner product
        of the `i`-th and `j`-th generator. Its entries are only well
        defined `\mod (V, W)`.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ, 3, matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T.gram_matrix_bilinear()
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
        """
        gens = self.gens()
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i + 1):
                G[i, j] = G[j, i] = (gens[i] * gens[j]).lift()
        return G

    @cached_method
    def gram_matrix_quadratic(self):
        r"""
        The Gram matrix of the quadratic form with respect to the generators.

        OUTPUT:

        - a rational matrix ``Gq`` with ``Gq[i,j] = gens[i]*gens[j]``
          and ``G[i,i] = gens[i].q()``

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ, [[2,0,0,-1],[0,2,0,-1],[0,0,2,-1],[-1,-1,-1,2]])
            sage: D4 = FreeQuadraticModule(ZZ, 4, D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: discrForm = TorsionQuadraticModule(D4dual, D4)
            sage: discrForm.gram_matrix_quadratic()
            [  1 1/2]
            [1/2   1]
            sage: discrForm.gram_matrix_bilinear()
            [  0 1/2]
            [1/2   0]
        """
        gens = self.gens()
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i):
                G[i, j] = G[j, i] = (gens[i] * gens[j]).lift()
            G[i, i] = gens[i].q().lift()
        return G

    def gens(self):
        r"""
        Return generators of ``self``.

        There is no assumption on the generators except that they
        generate the module.

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ, 3)
            sage: T = TorsionQuadraticModule(V, 5*V)
            sage: T.gens()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        """
        if self._gens_user is None:
            return self.smith_form_gens()
        return self._gens_user

    def genus(self, signature_pair):
        r"""
        Return the genus defined by ``self`` and the ``signature_pair``.

        If no such genus exists, raise a ``ValueError``.

        REFERENCES:

        [Nik1977]_ Corollary 1.9.4 and 1.16.3.

        EXAMPLES::

            sage: L = IntegralLattice("D4").direct_sum(IntegralLattice("A2"))
            sage: D = L.discriminant_group()
            sage: genus = D.genus(L.signature_pair())
            sage: genus
            Genus of
            None
            Signature:  (6, 0)
            Genus symbol at 2:    1^4:2^-2
            Genus symbol at 3:     1^-5 3^-1
            sage: genus == L.genus()
            True

        Let `H` be an even unimodular lattice of signature `(9, 1)`.
        Then `L = D_4 + A_2` is primitively embedded in `H`. We compute the discriminant
        form of the orthogonal complement of `L` in `H`::

            sage: DK = D.twist(-1)
            sage: DK
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]

        We know that  `K` has signature `(5, 1)` and thus we can compute
        the genus of `K` as::

            sage: DK.genus((3,1))
            Genus of
            None
            Signature:  (3, 1)
            Genus symbol at 2:    1^2:2^-2
            Genus symbol at 3:     1^-3 3^1

        We can also compute the genus of an odd lattice
        from its discriminant form::

            sage: L = IntegralLattice(matrix.diagonal(range(1,5)))
            sage: D = L.discriminant_group()
            sage: D.genus((4,0))
            Genus of
            None
            Signature:  (4, 0)
            Genus symbol at 2:    [1^-2 2^1 4^1]_6
            Genus symbol at 3:     1^-3 3^1

        TESTS::

            sage: L.genus() == D.genus((4,0))
            True
            sage: D.genus((1,0))
            Traceback (most recent call last):
            ...
            ValueError: this discriminant form and signature do not define a genus

        A systematic test of lattices of small ranks and determinants::

            sage: from sage.quadratic_forms.genera.genus import genera
            sage: signatures = [(1,0),(1,1),(1,2),(3,0),(0,4)]
            sage: dets = range(1,33)
            sage: genera = flatten([genera(s, d, even=False) for d in dets for s in signatures])    # long time
            sage: all(g == g.discriminant_form().genus(g.signature_pair()) for g in genera)  # long time
            True
            """
        from sage.quadratic_forms.genera.genus import (Genus_Symbol_p_adic_ring,
                                                       GenusSymbol_global_ring,
                                                       p_adic_symbol,
                                                       is_GlobalGenus,
                                                       _blocks)
        from sage.misc.misc_c import prod
        s_plus = signature_pair[0]
        s_minus = signature_pair[1]
        rank = s_plus + s_minus
        if len(self.invariants()) > rank:
            raise ValueError("this discriminant form and " +
                             "signature do not define a genus")
        disc = self.cardinality()
        determinant = (-1)**s_minus * disc
        local_symbols = []
        for p in (2 * disc).prime_divisors():
            D = self.primary_part(p)
            if len(D.invariants()) != 0:
                G_p = D.gram_matrix_quadratic().inverse()
                # get rid of denominators without changing the local equivalence class
                G_p *= G_p.denominator()**2
                G_p = G_p.change_ring(ZZ)
                local_symbol = p_adic_symbol(G_p, p, D.invariants()[-1].valuation(p))
            else:
                local_symbol = []

            rk = rank - len(D.invariants())
            if rk > 0:
                if p == 2:
                    det = determinant.prime_to_m_part(2)
                    det *= prod([di[2] for di in local_symbol])
                    det = det % 8
                    local_symbol.append([ZZ(0), rk, det, ZZ(0), ZZ(0)])
                else:
                    det = legendre_symbol(determinant.prime_to_m_part(p), p)
                    det = (det * prod([di[2] for di in local_symbol]))
                    local_symbol.append([ZZ(0), rk, det])
            local_symbol.sort()
            local_symbol = Genus_Symbol_p_adic_ring(p, local_symbol)
            local_symbols.append(local_symbol)

        # This genus has the right discriminant group
        # but it may be empty
        genus = GenusSymbol_global_ring(signature_pair, local_symbols)
        sym2 = local_symbols[0].symbol_tuple_list()

        if sym2[0][0] != 0:
            sym2 = [[ZZ(0), ZZ(0), ZZ(1), ZZ(0), ZZ(0)]] + sym2
        if len(sym2) <= 1 or sym2[1][0] != 1:
            sym2 = sym2[:1] + [[ZZ(1), ZZ(0), ZZ(1), ZZ(0), ZZ(0)]] + sym2[1:]
        if len(sym2) <= 2 or sym2[2][0] != 2:
            sym2 = sym2[:2] + [[ZZ(2), ZZ(0), ZZ(1), ZZ(0), ZZ(0)]] + sym2[2:]

        if self.value_module_qf().n == 1:
            # in this case the blocks of scales 1, 2, 4 are under determined
            # make sure the first 3 symbols are of scales 1, 2, 4
            # i.e. their valuations are 0, 1, 2

            # the form is odd
            block0 = [b for b in _blocks(sym2[0]) if b[3] == 1]

            o = sym2[1][3]
            # no restrictions on determinant and
            # oddity beyond existence
            # but we know if even or odd
            block1 = [b for b in _blocks(sym2[1]) if b[3] == o]

            d = sym2[2][2]
            o = sym2[2][3]
            t = sym2[2][4]
            # if the jordan block of scale 2 is even we know it
            if o == 0:
                block2 = [sym2[2]]
            # if it is odd we know det and oddity mod 4 at least
            else:
                block2 = [b for b in _blocks(sym2[2]) if b[3] == o
                          and (b[2] - d) % 4 == 0
                          and (b[4] - t) % 4 == 0
                          and (b[2] - d) % 8 == (b[4] - t) % 8  # if the oddity is altered by 4 then so is the determinant
                         ]
        elif self.value_module_qf().n == 2:
            # the form is even
            block0 = [b for b in _blocks(sym2[0]) if b[3] == 0]

            # if the jordan block of scale 2 is even we know it
            d = sym2[1][2]
            o = sym2[1][3]
            t = sym2[1][4]
            if o == 0:
                block1 = [sym2[1]]
            else:
                # the block is odd and we know det and oddity mod 4
                block1 = [b for b in _blocks(sym2[1])
                          if b[3] == o
                          and (b[2] - d) % 4 == 0
                          and (b[4] - t) % 4 == 0
                          and (b[2] - d) % 8 == (b[4] - t) % 8 # if the oddity is altered by 4 then so is the determinant
                         ]
            # this is completely determined
            block2 = [sym2[2]]
        else:
            raise ValueError("this is not a discriminant form")

        # figure out which symbol defines a genus and return that
        for b0 in block0:
            for b1 in block1:
                for b2 in block2:
                    sym2[:3] = [b0, b1, b2]
                    local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
                    genus = GenusSymbol_global_ring(signature_pair, local_symbols)
                    if is_GlobalGenus(genus):
                        # make the symbol sparse again.
                        i = 0
                        k = 0
                        while i < 3:
                            if sym2[k][1] == 0:
                                sym2.pop(k)
                            else:
                                k = k + 1
                            i = i + 1
                        local_symbols[0] = Genus_Symbol_p_adic_ring(2, sym2)
                        genus = GenusSymbol_global_ring(signature_pair, local_symbols)
                        return genus
        raise ValueError("this discriminant form and signature do not define a genus")

    def is_genus(self, signature_pair, even=True):
        r"""
        Return ``True`` if there is a lattice with this signature and discriminant form.

        .. TODO::

            implement the same for odd lattices

        INPUT:

        - signature_pair -- a tuple of non negative integers ``(s_plus, s_minus)``
        - even -- bool (default: ``True``)

        EXAMPLES::

            sage: L = IntegralLattice("D4").direct_sum(IntegralLattice(3 * Matrix(ZZ,2,[2,1,1,2])))
            sage: D = L.discriminant_group()
            sage: D.is_genus((6,0))
            True

        Let us see if there is a lattice in the genus defined by the same discriminant form
        but with a different signature::

            sage: D.is_genus((4,2))
            False
            sage: D.is_genus((16,2))
            True
        """
        s_plus = ZZ(signature_pair[0])
        s_minus = ZZ(signature_pair[1])
        if s_plus < 0 or s_minus < 0:
            raise ValueError("signature invariants must be non negative")
        rank = s_plus + s_minus
        signature = s_plus - s_minus
        D = self.cardinality()
        det = (-1)**s_minus * D
        if rank < len(self.invariants()):
            return False
        if even and self._modulus_qf != 2:
            raise ValueError("the discriminant form of an even lattice has"
                             "values modulo 2.")
        if (not even) and not (self._modulus == self._modulus_qf == 1):
            raise ValueError("the discriminant form of an odd lattice has"
                             "values modulo 1.")
        if not even:
            raise NotImplementedError("at the moment sage knows how to do this only for even genera. " +
                                      " Help us to implement this for odd genera.")
        for p in D.prime_divisors():
            # check the determinant conditions
            Q_p = self.primary_part(p)
            gram_p = Q_p.gram_matrix_quadratic()
            length_p = len(Q_p.invariants())
            u = det.prime_to_m_part(p)
            up = gram_p.det().numerator().prime_to_m_part(p)
            if p != 2 and length_p == rank:
                if legendre_symbol(u, p) != legendre_symbol(up, p):
                    return False
            if p == 2:
                if rank % 2 != length_p % 2:
                    return False
                n = (rank - length_p) / 2
                if u % 4 != (-1)**(n % 2) * up % 4:
                    return False
                if rank == length_p:
                    a = QQ(1) / QQ(2)
                    b = QQ(3) / QQ(2)
                    diag = gram_p.diagonal()
                    if not (a in diag or b in diag):
                        if u % 8 != up % 8:
                            return False
        if self.brown_invariant() != signature:
            return False
        return True

    def orthogonal_group(self, gens=None, check=False):
        r"""
        Orthogonal group of the associated torsion quadratic form.

        .. WARNING::

            This is can be smaller than the orthogonal group of the bilinear form.

        INPUT:

        - ``gens`` --  a list of generators, for instance square matrices,
                       something that acts on ``self``, or an automorphism
                       of the underlying abelian group
        - ``check`` -- perform additional checks on the generators

        EXAMPLES:

        You can provide generators to obtain a subgroup of the full orthogonal group::

            sage: D = TorsionQuadraticForm(matrix.identity(2)/2)
            sage: f = matrix(2,[0,1,1,0])
            sage: D.orthogonal_group(gens=[f]).order()
            2

        If no generators are given a slow brute force approach is used to calculate the full orthogonal group::

            sage: D = TorsionQuadraticForm(matrix.identity(3)/2)
            sage: OD = D.orthogonal_group()
            sage: OD.order()
            6
            sage: fd = D.hom([D.1,D.0,D.2])
            sage: OD(fd)
            [0 1 0]
            [1 0 0]
            [0 0 1]

        We compute the kernel of the action of the orthogonal group of `L` on the discriminant group.

            sage: L = IntegralLattice('A4')
            sage: O = L.orthogonal_group()
            sage: D = L.discriminant_group()
            sage: Obar = D.orthogonal_group(O.gens())
            sage: O.order()
            240
            sage: Obar.order()
            2
            sage: phi = O.hom(Obar.gens())
            sage: phi.kernel().order()
            120
        """
        from sage.groups.fqf_orthogonal import FqfOrthogonalGroup,_isom_fqf
        from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap

        ambient = AbelianGroupGap(self.invariants()).aut()
        # slow brute force implementation
        if gens is None:
            try:
                gens = self._orthogonal_group_gens
            except AttributeError:
                gens = _isom_fqf(self)
                gens = tuple(ambient(g) for g in gens)
                self._orthogonal_group_gens = gens
        else:
            # see if there is an action
            try:
                gens = [matrix(x*g for x in self.smith_form_gens()) for g in gens]
            except TypeError:
                pass 
            # the ambient knows what to do with the generators
            gens = tuple(ambient(g) for g in gens)
        Oq =  FqfOrthogonalGroup(ambient, gens, self, check=check)
        return Oq

    def orthogonal_submodule_to(self, S):
        r"""
        Return the submodule orthogonal to ``S``.

        INPUT:

        - ``S`` -- a submodule, list, or tuple of generators

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ, 10)
            sage: T = TorsionQuadraticModule(V, 3*V)
            sage: S = T.submodule(T.gens()[:5])
            sage: O = T.orthogonal_submodule_to(S)
            sage: O
            Finite quadratic module over Integer Ring with invariants (3, 3, 3, 3, 3)
            Gram matrix of the quadratic form with values in Q/3Z:
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            sage: O.V() + S.V() == T.V()
            True
        """
        if not isinstance(S, TorsionQuadraticModule):
            S = self.submodule(S)
        else:
            if not S.is_submodule(self):
                raise ValueError("S must be a submodule of this module")

        G = self.V().inner_product_matrix()
        T = self.V().basis_matrix()
        S = S.V().basis_matrix()
        m = self._modulus

        Y = T * G * S.transpose()
        # Elements of the ambient module which pair integrally with self.V()
        integral = Y.inverse() * T
        # Element of the ambient module which pair in mZZ with self.V()
        orthogonal = m * integral
        orthogonal = self.V().span(orthogonal.rows())
        # We have to make sure we get a submodule
        orthogonal = orthogonal.intersection(self.V())
        orthogonal = self.submodule(orthogonal.gens())
        return orthogonal

    @cached_method
    def normal_form(self, partial=False):
        r"""
        Return the normal form of this torsion quadratic module.

        Two torsion quadratic modules are isomorphic if and only if they have
        the same value modules and the same normal form.

        A torsion quadratic module `(T,q)` with values in `\QQ/n\ZZ` is
        in normal form if the rescaled quadratic module `(T, q/n)`
        with values in `\QQ/\ZZ` is in normal form.

        For the definition of normal form see [MirMor2009]_ IV Definition 4.6.
        Below are some of its properties.
        Let `p` be odd and `u` be the smallest non-square modulo `p`.
        The normal form is a diagonal matrix with diagonal entries either `p^n`
        or `u p^n`.

        If `p = 2` is even, then the normal form consists of
        1 x 1 blocks of the form

        .. MATH::

            (0), \quad 2^n(1),\quad 2^n(3),\quad 2^n(5) ,\quad 2^n(7)

        or of `2 \times 2` blocks of the form

        .. MATH::

            2^n
            \left(\begin{matrix}
                2 & 1\\
                1 & 2
            \end{matrix}\right), \quad
            2^n
            \left(\begin{matrix}
                0 & 1\\
                1 & 0
            \end{matrix}\right).

       The blocks are ordered by their valuation.

        INPUT:

        - partial - bool (default: ``False``) return only a partial normal form
          it is not unique but still useful to extract invariants

        OUTPUT:

        - a torsion quadratic module

        EXAMPLES::

            sage: L1=IntegralLattice(matrix([[-2,0,0],[0,1,0],[0,0,4]]))
            sage: L1.discriminant_group().normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0]
            [  0 1/4]
            sage: L2=IntegralLattice(matrix([[-2,0,0],[0,1,0],[0,0,-4]]))
            sage: L2.discriminant_group().normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/2   0]
            [  0 1/4]

        We check that :trac:`24864` is fixed::

            sage: L1=IntegralLattice(matrix([[-4,0,0],[0,4,0],[0,0,-2]]))
            sage: AL1=L1.discriminant_group()
            sage: L2=IntegralLattice(matrix([[-4,0,0],[0,-4,0],[0,0,2]]))
            sage: AL2=L2.discriminant_group()
            sage: AL1.normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 4)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/4   0]
            [  0   0 5/4]
            sage: AL2.normal_form()
            Finite quadratic module over Integer Ring with invariants (2, 4, 4)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/2   0   0]
            [  0 1/4   0]
            [  0   0 5/4]

        Some exotic cases::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2])
            sage: D4 = FreeQuadraticModule(ZZ,4,D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: T = TorsionQuadraticModule((1/6)*D4dual,D4)
            sage: T
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [ 1/18  1/12  5/36  1/36]
            [ 1/12   1/6  1/36   1/9]
            [ 5/36  1/36  1/36 11/72]
            [ 1/36   1/9 11/72  1/36]
            sage: T.normal_form()
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [ 1/6 1/12    0    0    0    0    0    0]
            [1/12  1/6    0    0    0    0    0    0]
            [   0    0 1/12 1/24    0    0    0    0]
            [   0    0 1/24 1/12    0    0    0    0]
            [   0    0    0    0  1/9    0    0    0]
            [   0    0    0    0    0  1/9    0    0]
            [   0    0    0    0    0    0  1/9    0]
            [   0    0    0    0    0    0    0  1/9]

        TESTS:

        A degenerate case::

            sage: T = TorsionQuadraticModule((1/6)*D4dual, D4, modulus=1/36)
            sage: T.normal_form()
            Finite quadratic module over Integer Ring with invariants (6, 6, 12, 12)
            Gram matrix of the quadratic form with values in Q/(1/18)Z:
            [1/36 1/72    0    0    0    0    0    0]
            [1/72 1/36    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
            [   0    0    0    0    0    0    0    0]
        """
        gens = []
        from sage.quadratic_forms.genera.normal_form import p_adic_normal_form, _normalize
        for p in self.annihilator().gen().prime_divisors():
            D_p = self.primary_part(p)
            q_p = D_p.gram_matrix_quadratic()
            q_p = q_p / D_p._modulus_qf

            # continue with the non-degenerate part
            r = q_p.rank()
            if r != q_p.ncols():
                U = q_p._clear_denom()[0].hermite_form(transformation=True)[1]
            else:
                U = q_p.parent().identity_matrix()
            kernel = U[r:, :]
            nondeg = U[:r, :]
            q_p = nondeg * q_p * nondeg.T

            # the normal form is implemented for p-adic lattices
            # so we should work with the lattice q_p --> q_p^-1
            q_p1 = q_p.inverse()
            prec = self.annihilator().gen().valuation(p) + 5
            D, U = p_adic_normal_form(q_p1, p, precision=prec + 5, partial=partial)
            # if we compute the inverse in the p-adics everything explodes --> go to ZZ
            U = U.change_ring(ZZ).inverse().transpose()

            # the inverse is in normal form - so to get a normal form for the original one
            # it is enough to massage each 1x1 resp. 2x2 block.
            U = U.change_ring(Zp(p, type='fixed-mod', prec=prec)).change_ring(ZZ)
            D = U * q_p * U.T * p**q_p.denominator().valuation(p)
            D = D.change_ring(Zp(p, type='fixed-mod', prec=prec))
            _, U1 = _normalize(D, normal_odd=False)
            U = U1.change_ring(ZZ) * U

            # reattach the degenerate part
            nondeg = U * nondeg
            U = nondeg.stack(kernel)

            # apply U to the generators
            n = U.ncols()
            gens_p = []
            for i in range(n):
                g = self.V().zero()
                for j in range(n):
                    g += D_p.gens()[j].lift() * U[i, j]
                gens_p.append(g)
            gens += gens_p
        return self.submodule_with_gens(gens)

    def primary_part(self, m):
        r"""
        Return the ``m``-primary part of this torsion quadratic module
        as a submodule.

        INPUT:

        - ``m`` -- an integer

        OUTPUT:

        - a submodule

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule((1/6)*ZZ^3,ZZ^3)
            sage: T
            Finite quadratic module over Integer Ring with invariants (6, 6, 6)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/36    0    0]
            [   0 1/36    0]
            [   0    0 1/36]
            sage: T.primary_part(2)
            Finite quadratic module over Integer Ring with invariants (2, 2, 2)
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/4   0   0]
            [  0 1/4   0]
            [  0   0 1/4]

        TESTS::

            sage: T == T.primary_part(T.annihilator().gen())
            True
        """
        annihilator = self.annihilator().gen()
        a = annihilator.prime_to_m_part(m)
        return self.submodule((a * self.V()).gens())

    def submodule_with_gens(self, gens):
        r"""
        Return a submodule with generators given by ``gens``.

        INPUT:

        - ``gens`` -- a list of generators that convert into ``self``

        OUTPUT:

        - a submodule with the specified generators

        EXAMPLES::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*10)
            sage: T = TorsionQuadraticModule((1/10)*V, V)
            sage: g = T.gens()
            sage: new_gens = [2*g[0], 5*g[0]]
            sage: T.submodule_with_gens(new_gens)
            Finite quadratic module over Integer Ring with invariants (10,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [2/5   0]
            [  0 1/2]

        The generators do not need to be independent::

            sage: new_gens = [g[0], 2*g[1], g[0], g[1]]
            sage: T.submodule_with_gens(new_gens)
            Finite quadratic module over Integer Ring with invariants (10, 10)
            Gram matrix of the quadratic form with values in Q/2Z:
            [1/10    0 1/10    0]
            [   0  2/5    0  1/5]
            [1/10    0 1/10    0]
            [   0  1/5    0 1/10]

        TESTS:

        Test that things work without specified gens too::

            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T
            Finite quadratic module over Integer Ring with invariants (5, 5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
            sage: T.submodule(T.gens()[:2])
            Finite quadratic module over Integer Ring with invariants (5, 5)
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0]
            [  0 1/5]
        """
        gens = tuple(self(v) for v in gens)
        V = self.V().submodule([v.lift() for v in gens]) + self._W
        W = self.W()
        return TorsionQuadraticModule(V, W, gens=gens, modulus=self._modulus,
                                      modulus_qf=self._modulus_qf, check=False)

    def twist(self, s):
        r"""
        Return the torsion quadratic module with quadratic form scaled by ``s``.

        If the old form was defined modulo `n`, then the new form is defined
        modulo `n s`.

        INPUT:

        - ``s`` - a rational number

        EXAMPLES::

            sage: q = TorsionQuadraticForm(matrix.diagonal([3/9, 1/9]))
            sage: q.twist(-1)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/Z:
            [2/3   0]
            [  0 8/9]

        This form is defined modulo `3`::

            sage: q.twist(3)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/3Z:
            [  1   0]
            [  0 1/3]

        The next form is defined modulo `4`::

            sage: q.twist(4)
            Finite quadratic module over Integer Ring with invariants (3, 9)
            Gram matrix of the quadratic form with values in Q/4Z:
            [4/3   0]
            [  0 4/9]
        """
        s = self.base_ring().fraction_field()(s)
        n = self.V().degree()
        inner_product_matrix = s * self.V().inner_product_matrix()
        ambient = FreeQuadraticModule(self.base_ring(), n, inner_product_matrix)
        V = ambient.span(self.V().basis())
        W = ambient.span(self.W().basis())
        return TorsionQuadraticModule(V, W)

    def value_module(self):
        r"""
        Return `\QQ / m\ZZ` with `m = (V, W)`.

        This is where the inner product takes values.

        EXAMPLES::

            sage: A2 = Matrix(ZZ, 2, 2, [2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module()
            Q/Z
        """
        return QmodnZ(self._modulus)

    def value_module_qf(self):
        r"""
        Return `\QQ / n\ZZ` with `n\ZZ = (V,W) + \ZZ \{ (w,w) | w \in W \}`.

        This is where the torsion quadratic form takes values.

        EXAMPLES::

            sage: A2 = Matrix(ZZ, 2, 2, [2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module over Integer Ring with invariants (2, 6)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module_qf()
            Q/2Z
        """
        return QmodnZ(self._modulus_qf)


def _brown_indecomposable(q, p):
    r"""
    Return the Brown invariant of the indecomposable form ``q``.

    The values are taken from Table 2.1 in [Shim2016]_.

    INPUT:

    - ``q`` - an indecomposable quadratic form represented by a
      rational `1 \times 1` or `2 \times 2` matrix
    - ``p`` - a prime number

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import _brown_indecomposable
        sage: q = Matrix(QQ, [1/3])
        sage: _brown_indecomposable(q,3)
        6
        sage: q = Matrix(QQ, [2/3])
        sage: _brown_indecomposable(q,3)
        2
        sage: q = Matrix(QQ, [5/4])
        sage: _brown_indecomposable(q,2)
        5
        sage: q = Matrix(QQ, [7/4])
        sage: _brown_indecomposable(q,2)
        7
        sage: q = Matrix(QQ, 2, [0,1,1,0])/2
        sage: _brown_indecomposable(q,2)
        0
        sage: q = Matrix(QQ, 2, [2,1,1,2])/2
        sage: _brown_indecomposable(q,2)
        4
    """
    v = q.denominator().valuation(p)
    if p == 2:
        # brown(U) = 0
        if q.ncols() == 2:
            if q[0, 0].valuation(2) > v + 1 and q[1, 1].valuation(2) > v + 1:
                # type U
                return mod(0, 8)
            else:
                # type V
                return mod(4 * v, 8)
        u = q[0, 0].numerator()
        return mod(u + v * (u**2 - 1) / 2, 8)
    if p % 4 == 1:
        e = -1
    if p % 4 == 3:
        e = 1
    if v % 2 == 1:
        u = q[0, 0].numerator() // 2
        if legendre_symbol(u, p) == 1:
            return mod(1 + e, 8)
        else:
            return mod(-3 + e, 8)
    return mod(0, 8)
