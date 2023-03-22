"""
The Fusion Ring of the Drinfeld Double of a Finite Group
"""
# ****************************************************************************
#  Copyright (C) 2023 Wenqi Li
#                     Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.all import Algebras, AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.misc.misc import inject_variable
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ideal import Ideal

class FusionDouble(CombinatorialFreeModule):
    r"""
    This constructs the Fusion Ring of the modular
    tensor category of modules over the Drinfeld
    Double of a finite group. Usage is similar
    to :class:`FusionRing`::.

    INPUT:

    - ``G`` -- a finite group
    - ``prefix`` (optional: defaults to 's') a prefix for the names of simple objects.
    - ``inject_varables`` (optional): set TRUE to create variables for the simple objects.

    REFERENCES:

    - [BaKi2001]_ Chapter 3
    - [Mas1995]_

    EXAMPLES ::

        sage: G = DihedralGroup(5)
        sage: H = FusionDouble(G, inject_variables=True)
        sage: H.basis()
        Finite family {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8, 9: s9, 10: s10, 11: s11, 12: s12, 13: s13, 14: s14, 15: s15}
        sage: for x in H.basis():
        ....:     print ("%s : %s"%(x,x^2))
        ....:
        s0 : s0
        s1 : s0
        s2 : s0 + s1 + s3
        s3 : s0 + s1 + s2
        s4 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s5 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s6 : s0 + s1 + s11
        s7 : s0 + s1 + s13
        s8 : s0 + s1 + s15
        s9 : s0 + s1 + s12
        s10 : s0 + s1 + s14
        s11 : s0 + s1 + s6
        s12 : s0 + s1 + s9
        s13 : s0 + s1 + s7
        s14 : s0 + s1 + s10
        s15 : s0 + s1 + s8
        sage: s4*s5
        s1 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        sage: s4.ribbon()
        1
        sage: s5.ribbon()
        -1
        sage: s8.ribbon()
        zeta5^3
    """
    def __init__(self, G, prefix="s",inject_variables=False):
        self._G = G
        self._prefix = prefix
        self._names = {}
        self._elt = {}
        self._chi = {}
        count = 0
        for g in G.conjugacy_classes_representatives():
            for chi in G.centralizer(g).irreducible_characters():
                self._names[count] = "%s%s"%(prefix, count)
                self._elt[count] = g
                self._chi[count] = chi
                count += 1
        self._rank = count
        self._cyclotomic_order = G.exponent()
        cat = AlgebrasWithBasis(ZZ).Subobjects()
        CombinatorialFreeModule.__init__(self, ZZ, [k for k in self._names], category=cat)
        if inject_variables:
            self.inject_variables()

    def _repr_(self):
        return "The Fusion Ring of the Drinfeld Double of %s"%self._G

    def __call__(self, *args):
        if len(args) > 1:
            args = (args,)
        return super(GAlg, self).__call__(*args)

    def _element_constructor(self, k):
        return self.monomial(k)

    def inject_variables(self):
        for i in range(self._rank):
            inject_variable(self._names[i],self.monomial(i))

    @cached_method
    def s_ij(self, i, j, unitary=False):
        sum = 0
        G = self._G
        [i, j] = [x.support_of_term() for x in [i,j]]
        [a, chi_1] = [self._elt[i], self._chi[i]]
        [b, chi_2] = [self._elt[j], self._chi[j]]
        for g in G:
            if a*g*b*g.inverse() == g*b*g.inverse()*a:
                sum += chi_1(g*b*g.inverse()) * chi_2(g.inverse()*a*g)
        if unitary:
            coef = 1 / (G.centralizer(a).order() * G.centralizer(b).order())
        else:
            coef = G.order() / (G.centralizer(a).order() * G.centralizer(b).order())
        return coef * sum

    def s_ijconj(self, i, j, unitary=False):
        return self.s_ij(i, j, unitary=unitary).conjugate()

    @cached_method
    def N_ijk(self, i, j, k):
        """
        The symmetric invariant of three simple objects,
        this returns the dimension of

        .. MATH::
           Hom(i \\otimes j\\otimes k, s_0)

        where `s_0` is the unit element (assuming prefix='s').
        Method of computation is through the Verlinde formula.
        """
        sz = self.one()
        return ZZ(sum(self.s_ij(i, r, unitary=True) * self.s_ij(j, r, unitary=True) * self.s_ij(k, r, unitary=True)/self.s_ij(sz, r, unitary=True) for r in self.basis()))

    def Nk_ij(self, i, j, k):
        r"""
        Returns the fusion coefficient `N^k_{ij}`, computed using the Verlinde formula:

        .. MATH::
            N^k_{ij} = \sum_l \frac{s(i, \ell)\, s(j, \ell)\, \overline{s(k, \ell)}}{s(I, \ell)},

        """
        return self.N_ijk(i, j, self.dual(k))

    def char_Nk_ij(self, i, j, k):
        r"""
        Use character theoretic method to compute the fusion coefficient `N_{ij}^k`.
        Each simple object, for example `i` corresponds to a conjugacy class `\mathcal{C}_i`
        of the underlying group `G`, and an irreducible character `\chi_i` of the
        centralizer `C(g_i)` of a representative `g_i` of `\mathcal{C}_i`. In addition
        to a fixed representative `g_i` and `g_j` of the classes `\mathcal{C}_i`
        and `\mathcal{C}_j`, the formula will make use of variable elements `h_i` and
        `h_j` that are subject to the condition `h_ih_j=g_k`.

        .. MATH::

            \frac{|\mathcal{C}_i|}{|G|}\sum_{\substack{h_i\in\mathcal{C}_i\\ h_j\in\mathcal{C}_j\\ h_ih_j=g_k}}|C(h_i)\cap C(h_j)|\,
            \langle\chi_i^{(h_i)}\chi_j^{(h_j)},\chi_k\rangle_{C(h_i)\cap C(h_j)},

        where `\chi_i^{(h_i)}` is the character `\chi_i` of `C(g_i)` conjugated to a
        character of `C(h_i)`, when `h_i` is a conjugate of the fixed representative `g_i`.
        More exactly, there exists `r_i` such that `r_i g_i r_i^{-1}=h_i`, and then
        `\chi_i^{(h_i)}(x)=\chi_i(r_i^{-1}xr_i)`, and this definition does not
        depend on the choice of `r_i`.
        This formula is due to Christopher Goff when the centralizers are normal, and
        to Wenqi Li in the general case.
        """
        G = self._G
        I = G.conjugacy_class(i.g())
        J = G.conjugacy_class(j.g())
        IJ = Set(I_elem * J_elem for I_elem in I for J_elem in J)
        if k.g() not in IJ:
            return 0

        K = G.conjugacy_class(k.g())
        CI = G.centralizer(i.g())
        CJ = G.centralizer(j.g())
        CK = G.centralizer(k.g())

        c = K.cardinality() / G.order()
        summands = [(I_elem, J_elem) for I_elem in I for J_elem in J if I_elem * J_elem == k.g()]
        res = 0
        for p in summands:
            I_elem, J_elem = p
            for g in G:
                if g.inverse() * i.g() * g == I_elem:
                    i_twist = g
                if g.inverse() * j.g() * g == J_elem:
                    j_twist = g
            A = Set(i_twist.inverse() * zi * i_twist for zi in CI)
            B = Set(j_twist.inverse() * zj * j_twist for zj in CJ)
            inner_summands = A.intersection(B).intersection(Set(CK))
            for x in inner_summands:
                res += i.chi()(i_twist * x * i_twist.inverse()) * j.chi()(j_twist * x * j_twist.inverse()) * k.chi()(x).conjugate()
        return c * res

    def field(self):
        return CyclotomicField(4 * self._cyclotomic_order)

    def root_of_unity(self, r, base_coercion=True):
        r"""
        Return `e^{i\pi r}` as an element of ``self.field()`` if possible.

        INPUT:

        - ``r`` -- a rational number

        EXAMPLES::

            sage: H = FusionDouble(DihedralGroup(6))
            sage: H.field()
            Cyclotomic Field of order 24 and degree 8
            sage: for n in [1..7]:
            ....:     try:
            ....:         print (n,H.root_of_unity(2/n))
            ....:     except ValueError as err:
            ....:         print (n,err)
            ....: 
            1 1
            2 -1
            3 zeta24^4 - 1
            4 zeta24^6
            5 not a root of unity in the field
            6 zeta24^4
            7 not a root of unity in the field
        """
        n = 2 * r * self._cyclotomic_order
        if n not in ZZ:
            raise ValueError("not a root of unity in the field")
        return  self.field().gen() ** n

    @cached_method
    def r_matrix(self, i, j, k):
        r"""
        Return the R-matrix entry corresponding to the subobject ``k``
        in the tensor product of ``i`` with ``j``. This method is only
        correct if the fusion coefficient ``N_{ij}^k\leq 1``. See the
        :class:`FusionRing` method for more information, including
        the reason for this caveat, and the algorithm.
        """
        if self.Nk_ij(i, j, k) == 0:
            return self.field().zero()
        if i != j:
            ret = self.root_of_unity((k.twist() - i.twist() - j.twist()) / 2)
        else:
            i0 = self.one()
            B = self.basis()
            ret = sum(y.ribbon()**2 / (i.ribbon() * x.ribbon()**2)
                   * self.s_ij(i0, y) * self.s_ij(i, z) * self.s_ijconj(x, z)
                   * self.s_ijconj(k, x) * self.s_ijconj(y, z) / self.s_ij(i0, z)
                   for x in B for y in B for z in B) / (self.total_q_order(base_coercion=False)**4)

        return ret

    def total_q_order(self):
        r"""
        Return the positive square root of :meth:`self.global_q_dimension()
        <global_q_dimension>` as an element of :meth:`self.field() <field>`.

        For the Drinfeld double of a finite group `G`, this equals the
        cardinality of `G`.
        """
        return self._G.order()

    def is_multiplicity_free(self, verbose=False):
        """
        Returns True if all fusion coefficients are at most 1.
        """
        for i in self.basis():
            for j in self.basis():
                for k in self.basis():
                    if self.N_ijk(i,j,k) > 1:
                        if verbose:
                            print ("N(%s,%s,%s)=%s"%(i,j,k,self.N_ijk(i,j,k)))
                        return False
        return True

    def one(self):
        """
        The unit element of the ring.
        """
        return self.basis()[0]

    def dual(self,i):
        sz = self.one()
        for j in self.basis():
            if self.N_ijk(i,j,sz) > 0:
                return j

    def product_on_basis(self, a, b):
        d = {k.support_of_term() : self.N_ijk(self.monomial(a),self.monomial(b),self.dual(k)) for k in self.basis()}
        return self._from_dict(d)

    def _repr_term(self, t):
        return self._names[t]

    def group(self):
        """
        Returns the name of the underlying group.

        EXAMPLE::
            sage: FusionDouble(DiCyclicGroup(4)).group()
            Diyclic group of order 16 as a permutation group
        """
        return self._G

    def IdGroup(self):
        """
        returns the Gap Group ID

        EXAMPLE::
            sage: FusionDouble(DiCyclicGroup(4)).IdGroup()
            [ 16, 9 ]

        """
        return self._G._libgap_().IdGroup()

    def get_fmatrix(self, *args, **kwargs):
        """
        Construct an :class:`FMatrix` factory to solve the pentagon relations
        and organize the resulting F-symbols.
        """
        if not hasattr(self, 'fmats') or kwargs.get('new', False):
            kwargs.pop('new', None)
            from sage.algebras.fusion_rings.f_matrix import FMatrix
            self.fmats = FMatrix(self, *args, **kwargs)
        return self.fmats

    class Element(CombinatorialFreeModule.Element):
        def is_simple_object(self):
            r"""
            Determine whether ``self`` is a simple object of the fusion ring.
            """
            return self in self.parent().basis()

        def g(self):
           """
           Returns the conjugacy class representative of the underlying
           group corresponding to a simple object.
           """
           return self.parent()._elt[self.support_of_term()]

        def chi(self):
           """
           Returns the character of the centralizer of a conjugacy class 
           representative of the underlying group corresponding to a simple object.
           """
           return self.parent()._chi[self.support_of_term()]

        def ribbon(self):
            """
            The twist or ribbon of the simple object.

            EXAMPLE::
                sage: H = FusionDouble(CyclicPermutationGroup(3))
                sage: [i.ribbon() for i in H.basis()]
                [1, 1, 1, 1, zeta3, -zeta3 - 1, 1, -zeta3 - 1, zeta3]
            """
            return self.chi()(self.g()) / self.chi()(self.parent()._G.one())

        def twist(self, reduced=True):
            r"""
            Return a rational number `h` such that `\theta = e^{i \pi h}`
            is the twist of ``self``. The quantity `e^{i \pi h}` is
            also available using :meth:`ribbon`.

            This method is only available for simple objects.
            """
            if not self.is_simple_object():
                raise ValueError("quantum twist is only available for simple objects of a FusionRing")
            zeta = self.parent().field().gen()
            rib = self.ribbon()
            for k in range(4*self.parent()._cyclotomic_order):
                if zeta**k == rib:
                    return k/(2*self.parent()._cyclotomic_order)

        def dual(self):
            """
            Return the dual of self.

            EXAMPLE::
                sage: G = CyclicPermutationGroup(4)
                sage: H = FusionDouble(G, prefix="b")
                sage: [x for x in H.basis() if x==x.dual()]
                [b0, b1, b8, b9]
            """
            if not self.is_simple_object():
                raise ValueError("dual is only available for simple objects of a FusionRing")
            return self.parent().dual(self)
            return

        @cached_method
        def q_dimension(self):
            """
            Return the q-dimension of self.

            EXAMPLE::
                sage: G = AlternatingGroup(4)
                sage: H = FusionDouble(G)
                sage: [x.q_dimension() for x in H.basis()]
                [1, 1, 1, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4]
                sage: sum(x.q_dimension()^2 for x in H.basis()) == G.order()^2
                True
            """
            if not self.is_simple_object():
                raise ValueError("quantum dimension is only available for simple objects of a FusionRing")
            return self.parent().s_ij(self,self.parent().one())
