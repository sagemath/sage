r"""
LLT symmetric functions

REFERENCES:

.. [LLT1997] Alain Lascoux, Bernard Leclerc, Jean-Yves Thibon,
   Ribbon tableaux, Hall-Littlewood functions, quantum affine algebras, and unipotent varieties,
   J. Math. Phys. 38 (1997), no. 2, 1041-1068,
   arXiv:q-alg/9512-31v1 [math.q.alg]

.. [LT2000] Bernard Leclerc and Jean-Yves Thibon,
   Littlewood-Richardson coefficients and Kazhdan-Lusztig polynomials,
   in: Combinatorial methods in representation theory (Kyoto)
   Adv. Stud. Pure Math., vol. 28, Kinokuniya, Tokyo, 2000, pp 155-220
   arXiv:math/9809122v3 [math.q-alg]
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.calculus.var import var
import sfa
import sage.combinat.ribbon_tableau as ribbon_tableau
import sage.combinat.skew_partition
from sage.rings.all import ZZ
import sage.combinat.partition
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.rings.rational_field import QQ

# cache for H spin basis
hsp_to_m_cache={}
m_to_hsp_cache={}

# cache for H cospin basis
hcosp_to_m_cache={}
m_to_hcosp_cache={}

QQt = QQ['t'].fraction_field()
# This is to become the "abstract algebra" for llt polynomials

class LLT_class(UniqueRepresentation):
    r"""
    A class for working with LLT symmetric functions.

    EXAMPLES::

        sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
        sage: L3 = Sym.llt(3); L3
        level 3 LLT polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: L3.cospin([3,2,1])
        (t+1)*m[1, 1] + m[2]
        sage: HC3 = L3.hcospin(); HC3
        Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT cospin basis
        sage: m = Sym.monomial()
        sage: m( HC3[1,1] )
        (t+1)*m[1, 1] + m[2]

    We require that the parameter `t` must be in the base ring::

        sage: Symxt = SymmetricFunctions(QQ['x','t'].fraction_field())
        sage: (x,t) = Symxt.base_ring().gens()
        sage: LLT3x = Symxt.llt(3,t=x)
        sage: LLT3 = Symxt.llt(3)
        sage: HS3x = LLT3x.hspin()
        sage: HS3t = LLT3.hspin()
        sage: s = Symxt.schur()
        sage: s(HS3x[2,1])
        s[2, 1] + x*s[3]
        sage: s(HS3t[2,1])
        s[2, 1] + t*s[3]
        sage: HS3x(HS3t[2,1])
        HSp3[2, 1] + (-x+t)*HSp3[3]
        sage: s(HS3x(HS3t[2,1]))
        s[2, 1] + t*s[3]
        sage: LLT3t2 = Symxt.llt(3,t=2)
        sage: HC3t2 = LLT3t2.hcospin()
        sage: HS3x(HC3t2[3,1])
        2*HSp3[3, 1] + (-2*x+1)*HSp3[4]
    """

    def __init__(self, Sym, k, t='t'):
        r"""
        Class of LLT symmetric function bases

        INPUT:

        - ``self`` -- a family of LLT symmetric function bases
        - ``k`` -- a positive integer (the level)
        - ``t`` -- a parameter (default: `t`)

        EXAMPLES::

            sage: L3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3)
            sage: L3 == loads(dumps(L3))
            True
            sage: TestSuite(L3).run(skip=["_test_associativity","_test_distributivity","_test_prod"])

        TESTS::

            sage: L3 != SymmetricFunctions(FractionField(QQ['t'])).llt(2)
            True
            sage: L3p = SymmetricFunctions(FractionField(QQ['t'])).llt(3,t=1)
            sage: L3 != L3p
            True
            sage: L3p != SymmetricFunctions(QQ).llt(3,t=1)
            True
        """
        self._k = k
        self._sym = Sym
        self._name = "level %s LLT polynomials"%self._k
        if not (t in Sym.base_ring() or var(t) in Sym.base_ring()):
            raise ValueError, "parameter t must be in the base ring"
        self.t = Sym.base_ring()(t)
        self._name_suffix = ""
        if str(t) !='t':
            self._name_suffix += " with t=%s"%self.t
        self._name += self._name_suffix+" over %s"%self._sym.base_ring()
        self._m = sage.combinat.sf.sf.SymmetricFunctions(QQt).monomial()

    def __repr__(self):
        r"""
        Representation of the LLT symmetric functions

        INPUT:

        - ``self`` -- a family of LLT symmetric function bases

        OUTPUT:

        - returns a string representing the LLT symmetric functions

        EXAMPLES::

            sage: SymmetricFunctions(FractionField(QQ['t'])).llt(3)
            level 3 LLT polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: SymmetricFunctions(QQ).llt(3,t=2)
            level 3 LLT polynomials with t=2 over Rational Field
        """
        return self._name

    def symmetric_function_ring( self ):
        r"""
        The symmetric function algebra associated to the family of LLT
        symmetric function bases

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases

        OUTPUT:

        - returns the symmetric function ring associated to ``self``.

        EXAMPLES ::

            sage: L3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3)
            sage: L3.symmetric_function_ring()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._sym

    def base_ring(self):
        r"""
        Returns the base ring of ``self``.

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases

        OUTPUT:

        - returns the base ring of the symmetric function ring associated to ``self``

        EXAMPLES::

            sage: SymmetricFunctions(FractionField(QQ['t'])).llt(3).base_ring()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._sym.base_ring()

    def level(self):
        r"""
        Returns the level of ``self``.

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases

        OUTPUT:

        - the level is the parameter of `k` in the basis

        EXAMPLES::

            sage: SymmetricFunctions(FractionField(QQ['t'])).llt(3).level()
            3
        """
        return self._k

    def _llt_generic(self, skp, stat):
        r"""
        Takes in partition, list of partitions, or a list of skew
        partitions as well as a function which takes in two partitions and
        a level and returns a coefficient.

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases
        - ``skp`` -- a partition or a list of partitions or a list of skew partitions
        - ``stat`` -- a function which accepts two partitions and a value
          for the level and returns a coefficient which is a polynomial
          in a parameter `t`.  The first partition is the index of the
          LLT function, the second partition is the index of the monomial
          basis element.

        OUTPUT:

        - returns the monomial expansion of the LLT symmetric function
          indexed by ``skp``

        EXAMPLES::

            sage: L3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3)
            sage: f = lambda skp,mu,level: QQ(1)
            sage: L3._llt_generic([3,2,1],f)
            m[1, 1] + m[2]
            sage: L3._llt_generic([[2,1],[1],[2]],f)
            m[1, 1, 1, 1, 1, 1] + m[2, 1, 1, 1, 1] + m[2, 2, 1, 1] + m[2, 2, 2] + m[3, 1, 1, 1] + m[3, 2, 1] + m[3, 3] + m[4, 1, 1] + m[4, 2] + m[5, 1] + m[6]
            sage: L3._llt_generic([[[2,2],[1]],[[2,1],[]]],f)
            m[1, 1, 1, 1] + m[2, 1, 1] + m[2, 2] + m[3, 1] + m[4]
        """
        if skp in sage.combinat.partition.Partitions():
            m = (sum(skp) / self.level()).floor()
            if m == 0:
                raise ValueError, "level (%=) must divide %s "%(sum(skp), self.level())
            mu = sage.combinat.partition.Partitions( ZZ(sum(skp) / self.level()) )

        elif isinstance(skp, list) and skp[0] in sage.combinat.skew_partition.SkewPartitions():
            #skp is a list of skew partitions
            skp2 =  [sage.combinat.partition.Partition(core=[], quotient=[skp[i][0] for i in range(len(skp))])]
            skp2 += [sage.combinat.partition.Partition(core=[], quotient=[skp[i][1] for i in range(len(skp))])]
            mu = sage.combinat.partition.Partitions(ZZ((skp2[0].size()-skp2[1].size()) / self.level()))
            skp = skp2
        elif isinstance(skp, list) and skp[0] in sage.combinat.partition.Partitions():
            #skp is a list of partitions
            skp = sage.combinat.partition.Partition(core=[], quotient=skp)
            mu = sage.combinat.partition.Partitions( ZZ(sum(skp) / self.level() ))
        else:
            raise ValueError, "LLT polynomials not defined for %s"%skp

        BR = self.base_ring()
        return sum([ BR(stat(skp,nu,self.level()).subs(t=self.t))*self._m(nu) for nu in mu])

    def spin_square(self, skp):
         r"""
         Calculates a single instance of a spin squared LLT symmetric function
         associated with a partition, list of partitions, or a list of skew partitions.
         This family of symmetric functions is defined in [LT2000]_ equation (43).

         INPUT:

         - ``self`` -- a family of LLT symmetric functions bases
         - ``skp`` -- a partition of a list of partitions or a list of skew
           partitions

         OUTPUT:

         - returns the monomial expansion of the LLT symmetric function spin-square
           functions indexed by ``skp``

         EXAMPLES::

             sage: L3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3)
             sage: L3.spin_square([2,1])
             t*m[1]
             sage: L3.spin_square([3,2,1])
             (t^3+t)*m[1, 1] + t^3*m[2]
             sage: L3.spin_square([[1],[1],[1]])
             (t^6+2*t^4+2*t^2+1)*m[1, 1, 1] + (t^6+t^4+t^2)*m[2, 1] + t^6*m[3]
             sage: L3.spin_square([[[2,2],[1]],[[2,1],[]]])
             (2*t^4+3*t^2+1)*m[1, 1, 1, 1] + (t^4+t^2)*m[2, 1, 1] + t^4*m[2, 2]
         """
         return self._llt_generic(skp, ribbon_tableau.spin_polynomial_square)

    def cospin(self, skp):
        r"""
        Calculates a single instance of the cospin symmetric functions.
        These are the functions defined in [LLT1997]_ equation (26).

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases
        - ``skp`` -- a partition or a list of partitions or a list of skew
          partitions

        OUTPUT:

        - returns the monomial expansion of the LLT symmetric function cospin
          functions indexed by ``skp``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: L3 = Sym.llt(3)
            sage: L3.cospin([2,1])
            m[1]
            sage: L3.cospin([3,2,1])
            (t+1)*m[1, 1] + m[2]
            sage: s = Sym.schur()
            sage: s(L3.cospin([[2],[1],[2]]))
            t^4*s[2, 2, 1] + t^3*s[3, 1, 1] + (t^3+t^2)*s[3, 2] + (t^2+t)*s[4, 1] + s[5]
        """
        return self._llt_generic(skp, ribbon_tableau.cospin_polynomial)

#### Is it safe to delete this function?
##     def llt_inv(self, skp):
##         """
##         """
##         l = sage.combinat.partitions( sum( [ p.size() for p in skp ] ) ).list()
##         res = m(0)
##         for p in l:
##             inv_p = [ ktuple.inversions() for ktuple in kTupleTableaux(skp, p) ]
##             res += sum([t**x for x in inv_p])*m(p)
##         return res

    def hcospin(self):
        r"""
        Returns the HCospin basis.
        This basis is defined [LLT1997]_ equation (27).

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases

        OUPUT:

        - returns the h-cospin basis of the LLT symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HCosp3 = Sym.llt(3).hcospin(); HCosp3
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT cospin basis
            sage: HCosp3([1])^2
            1/t*HCosp3[1, 1] + ((t-1)/t)*HCosp3[2]

            sage: s = Sym.schur()
            sage: HCosp3(s([2]))
            HCosp3[2]
            sage: HCosp3(s([1,1]))
            1/t*HCosp3[1, 1] - 1/t*HCosp3[2]
            sage: s(HCosp3([2,1]))
            t*s[2, 1] + s[3]
        """
        return LLT_cospin(self)

    def hspin(self):
        r"""
        Returns the HSpin basis.
        This basis is defined [LLT1997]_ equation (28).

        INPUT:

        - ``self`` -- a family of LLT symmetric functions bases

        OUPUT:

        - returns the h-spin basis of the LLT symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HSp3 = Sym.llt(3).hspin(); HSp3
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT spin basis
            sage: HSp3([1])^2
            HSp3[1, 1] + (-t+1)*HSp3[2]

            sage: s = Sym.schur()
            sage: HSp3(s([2]))
            HSp3[2]
            sage: HSp3(s([1,1]))
            HSp3[1, 1] - t*HSp3[2]
            sage: s(HSp3([2,1]))
            s[2, 1] + t*s[3]
        """
        return LLT_spin(self)



class LLT_generic(sfa.SymmetricFunctionAlgebra_generic):

    def __init__(self, llt, prefix):
        r"""
        A class of methods which are common to both the hspin and hcospin
        of the LLT symmetric functions.

        INPUT:

        - ``self`` -- an instance of the LLT hspin or hcospin basis
        - ``llt`` -- a family of LLT symmetric functions

        EXAMPLES::

            sage: SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT spin basis
            sage: SymmetricFunctions(QQ).llt(3,t=2).hspin()
            Symmetric Functions over Rational Field in the level 3 LLT spin with t=2 basis
            sage: QQz = FractionField(QQ['z']); z = QQz.gen()
            sage: SymmetricFunctions(QQz).llt(3,t=z).hspin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in z over Rational Field in the level 3 LLT spin with t=z basis
        """
        s = self.__class__.__name__[4:]
        sfa.SymmetricFunctionAlgebra_generic.__init__(
            self, llt._sym,
            basis_name = "level %s LLT "%llt.level() + s + llt._name_suffix,
            prefix = prefix)

        self.t = llt.t
        self._sym = llt._sym
        self._llt = llt
        self._k = llt._k

        sfa.SymmetricFunctionAlgebra_generic.__init__(self, self._sym)

        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        category = sage.categories.all.ModulesWithBasis(self._sym.base_ring())
        self._m = llt._sym.m()
        self   .register_coercion(SetMorphism(Hom(self._m, self, category), self._m_to_self))
        self._m.register_coercion(SetMorphism(Hom(self, self._m, category), self._self_to_m))

    def _m_to_self(self, x):
        r"""
        Isomorphism from the monomial basis into ``self``

        INPUT:

        - ``self`` - an instance of the LLT hspin or hcospin basis
        - ``x`` - an element of the monomial basis

        OUTPUT:

        - returns ``x`` expanded in the basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HSp3 = Sym.llt(3).hspin()
            sage: m = Sym.monomial()
            sage: HSp3._m_to_self(m[2,1])
            -2*HSp3[1, 1, 1] + (2*t^2+2*t+1)*HSp3[2, 1] + (-2*t^2-t)*HSp3[3]

        This is for internal use only. Please use instead::

            sage: HSp3(m[2,1])
            -2*HSp3[1, 1, 1] + (2*t^2+2*t+1)*HSp3[2, 1] + (-2*t^2-t)*HSp3[3]
        """
        return self._from_cache(x, self._m_cache, self._m_to_self_cache, t = self.t)

    def _self_to_m(self, x):
        r"""
        Isomorphism from self to the monomial basis

        INPUT:

        - ``self`` -- an instance of the LLT hspin or hcospin basis
        - ``x`` -- an element of ``self``

        OUTPUT:

        - returns ``x`` expanded in the monomial basis.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HSp3 = Sym.llt(3).hspin()
            sage: m = Sym.monomial()
            sage: HSp3._self_to_m(HSp3[2,1])
            (t+2)*m[1, 1, 1] + (t+1)*m[2, 1] + t*m[3]

        This is for internal use only. Please use instead::

            sage: m(HSp3[2,1])
            (t+2)*m[1, 1, 1] + (t+1)*m[2, 1] + t*m[3]
        """
        return self._m._from_cache(x, self._m_cache, self._self_to_m_cache, t = self.t)


    def level(self):
        r"""
        Returns the level of ``self``.

        INPUT:

        - ``self`` -- an instance of the LLT hspin or hcospin basis

        OUTPUT:

        - returns the level associated to the basis ``self``.

        EXAMPLES::

            sage: HSp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            sage: HSp3.level()
            3
        """
        return self._k

    def llt_family( self ):
        r"""
        The family of the llt bases of the symmetric functions.

        INPUT:

        - ``self`` -- an instance of the LLT hspin or hcospin basis

        OUTPUT:

        - returns an instance of the family of LLT bases associated to ``self``.

        EXAMPLES::

            sage: HSp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            sage: HSp3.llt_family()
            level 3 LLT polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._llt

    def _multiply(self, left, right):
        r"""
        Convert to the monomial basis, do the multiplication there, and
        convert back to the basis ``self``.

        INPUT:

        - ``self`` -- an instance of the LLT hspin or hcospin basis
        - ``left``, ``right`` -- elements of the symmetric functions

        OUTPUT:

        - returns the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES::

            sage: HSp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            sage: HSp3._multiply(HSp3([1]), HSp3([2]))
            HSp3[2, 1] + (-t+1)*HSp3[3]
            sage: HCosp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hcospin()
            sage: HCosp3._multiply(HCosp3([1]), HSp3([2]))
            1/t*HCosp3[2, 1] + ((t-1)/t)*HCosp3[3]
        """
        return self( self._m(left) * self._m(right) )

    def _m_cache(self, n):
        r"""
        Compute the change of basis from the monomial symmetric functions
        to ``self``.

        INPUT:

        - ``self`` -- an instance of the LLT hspin or hcospin basis
        - ``n`` -- a positive integer representing the degree

        EXAMPLES::

            sage: HSp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            sage: HSp3._m_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( HSp3._self_to_m_cache[2] )
            [([1, 1], [([1, 1], t + 1), ([2], t)]), ([2], [([1, 1], 1), ([2], 1)])]
            sage: l( HSp3._m_to_self_cache[2] )
            [([1, 1], [([1, 1], 1), ([2], -t)]), ([2], [([1, 1], -1), ([2], t + 1)])]
            sage: HCosp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hcospin()
            sage: HCosp3._m_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( HCosp3._self_to_m_cache[2] )
            [([1, 1], [([1, 1], t + 1), ([2], 1)]), ([2], [([1, 1], 1), ([2], 1)])]
            sage: l( HCosp3._m_to_self_cache[2] )
            [([1, 1], [([1, 1], 1/t), ([2], -1/t)]),
             ([2], [([1, 1], -1/t), ([2], (t + 1)/t)])]
        """
        self._invert_morphism(n, QQt, self._self_to_m_cache, \
                              self._m_to_self_cache, to_other_function = self._to_m)

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        pass

# the H-spin basis
class LLT_spin(LLT_generic):

    def __init__(self, llt):
        r"""
        A class of methods for the h-spin LLT basis of the symmetric functions.

        INPUT:

        - ``self`` -- an instance of the LLT hcospin basis
        - ``llt`` -- a family of LLT symmetric function bases

        TESTS::

            sage: HSp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            sage: TestSuite(HSp3).run(skip = ["_test_associativity", "_test_distributivity", "_test_prod"]) # products are too expensive, long time (10s on sage.math, 2012)
            sage: TestSuite(HSp3).run(elements = [HSp3.t*HSp3[1,1]+HSp3.t*HSp3[2], HSp3[1]+(1+HSp3.t)*HSp3[1,1]])  # long time (depends on previous)

        ::

            sage: HS3t2 = SymmetricFunctions(QQ).llt(3,t=2).hspin()
            sage: TestSuite(HS3t2).run() # products are too expensive, long time (7s on sage.math, 2012)

        ::

            sage: HS3x = SymmetricFunctions(FractionField(QQ['x'])).llt(3,t=x).hspin()
            sage: TestSuite(HS3x).run(skip = ["_test_associativity", "_test_distributivity", "_test_prod"]) # products are too expensive, long time (4s on sage.math, 2012)
            sage: TestSuite(HS3x).run(elements = [HS3x.t*HS3x[1,1]+HS3x.t*HS3x[2], HS3x[1]+(1+HS3x.t)*HS3x[1,1]])  # long time (depends on previous)
        """
        level = llt._k
        if level not in hsp_to_m_cache:
            hsp_to_m_cache[level] = {}
            m_to_hsp_cache[level] = {}
        self._self_to_m_cache = hsp_to_m_cache[level]
        self._m_to_self_cache = m_to_hsp_cache[level]

        LLT_generic.__init__(self, llt, prefix="HSp%s"%level)


    def _to_m(self, part):
        r"""
        Returns a function which gives the coefficient of a partition
        in the monomial expansion of self(part).

        INPUT:

        - ``self`` -- an instance of the LLT hspin basis
        - ``part`` -- a partition

        OUTPUT:

        - returns a function which accepts a partition and returns the coefficient
          in the expansion of the monomial basis

        EXAMPLES::

            sage: HSp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hspin()
            sage: f21 = HSp3._to_m(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [t, t + 1, t + 2]
            sage: HSp3.symmetric_function_ring().m()( HSp3[2,1] )
            (t+2)*m[1, 1, 1] + (t+1)*m[2, 1] + t*m[3]
        """
        level = self.level()
        f = lambda part2: QQt(ribbon_tableau.spin_polynomial([level*i for i in part], part2, level))
        return f

    class Element(LLT_generic.Element):
        pass


# the h-cospin basis
class LLT_cospin(LLT_generic):
    def __init__(self, llt):
        r"""
        A class of methods for the h-cospin LLT basis of the symmetric functions.

        INPUT:

        - ``self`` -- an instance of the LLT hcospin basis
        - ``llt`` -- a family of LLT symmetric function bases

        TESTS::

            sage: HCosp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hcospin()
            sage: TestSuite(HCosp3).run(skip = ["_test_associativity", "_test_distributivity", "_test_prod"]) # products are too expensive, long time (11s on sage.math, 2012)
            sage: TestSuite(HCosp3).run(elements = [HCosp3.t*HCosp3[1,1]+HCosp3.t*HCosp3[2], HCosp3[1]+(1+HCosp3.t)*HCosp3[1,1]])  # long time (depends on previous)

        ::

            sage: HC3t2 = SymmetricFunctions(QQ).llt(3,t=2).hcospin()
            sage: TestSuite(HC3t2).run() # products are too expensive, long time (6s on sage.math, 2012)

        ::

            sage: HC3x = SymmetricFunctions(FractionField(QQ['x'])).llt(3,t=x).hcospin()
            sage: TestSuite(HC3x).run(skip = ["_test_associativity", "_test_distributivity", "_test_prod"]) # products are too expensive, long time (5s on sage.math, 2012)
            sage: TestSuite(HC3x).run(elements = [HC3x.t*HC3x[1,1]+HC3x.t*HC3x[2], HC3x[1]+(1+HC3x.t)*HC3x[1,1]])  # long time (depends on previous)
        """
        level = llt._k
        if level not in hcosp_to_m_cache:
            hcosp_to_m_cache[level] = {}
            m_to_hcosp_cache[level] = {}
        self._self_to_m_cache = hcosp_to_m_cache[level]
        self._m_to_self_cache = m_to_hcosp_cache[level]
        LLT_generic.__init__(self, llt, prefix= "HCosp%s"%level)

    def _to_m(self, part):
        r"""
        Returns a function which gives the coefficient of part2 in the
        monomial expansion of self(part).

        INPUT:

        - ``self`` -- an instance of the LLT hcospin basis
        - ``part`` -- a partition

        OUTPUT:

        - returns a function which accepts a partition and returns the coefficient
          in the expansion of the monomial basis

        EXAMPLES::

            sage: HCosp3 = SymmetricFunctions(FractionField(QQ['t'])).llt(3).hcospin()
            sage: f21 = HCosp3._to_m(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [1, t + 1, 2*t + 1]
            sage: HCosp3.symmetric_function_ring().m()( HCosp3[2,1] )
            (2*t+1)*m[1, 1, 1] + (t+1)*m[2, 1] + m[3]
        """
        level = self.level()
        f = lambda part2: QQt(ribbon_tableau.cospin_polynomial([level*i for i in part], part2, level))
        return f

    class Element(LLT_generic.Element):
        pass

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.llt', 'LLTElement_spin',  LLT_spin.Element)
register_unpickle_override('sage.combinat.sf.llt', 'LLTElement_cospin',  LLT_cospin.Element)
