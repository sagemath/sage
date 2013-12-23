"""
Classical symmetric functions.
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

from sage.rings.integer import Integer

import sage.combinat.skew_partition
import sage.libs.symmetrica.all as symmetrica # used in eval()

from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField



import hall_littlewood
import sfa
import llt
import macdonald
import jack
import orthotriang
import kschur

ZZ = IntegerRing()
QQ = RationalField()

translate = {'monomial':'MONOMIAL', 'homogeneous':'HOMSYM', 'powersum':'POWSYM', 'elementary':'ELMSYM', 'Schur':'SCHUR'}

conversion_functions = {}

def init():
    """
    Set up the conversion functions between the classical bases.

    EXAMPLES::

        sage: from sage.combinat.sf.classical import init
        sage: sage.combinat.sf.classical.conversion_functions = {}
        sage: init()
        sage: sage.combinat.sf.classical.conversion_functions[('Schur', 'powersum')]
        <built-in function t_SCHUR_POWSYM_symmetrica>

    The following checks if the bug described in :trac:`15312` is fixed.::

        sage: change = sage.combinat.sf.classical.conversion_functions[('powersum', 'Schur')]
        sage: hideme = change({Partition([1]*47):ZZ(1)}) # long time
        sage: change({Partition([2,2]):QQ(1)})
        s[1, 1, 1, 1] - s[2, 1, 1] + 2*s[2, 2] - s[3, 1] + s[4]
    """
    for other_basis in translate:
        for basis in translate:
            try:
                conversion_functions[(other_basis, basis)] = eval('symmetrica.t_' + translate[other_basis] + '_' +  translate[basis])
            except AttributeError:
                pass


init()


###################################
#                                 #
#  Classical Symmetric Functions  #
#                                 #
###################################
class SymmetricFunctionAlgebra_classical(sfa.SymmetricFunctionAlgebra_generic):
    """
    The class of classical symmetric functions.

    .. TODO:: delete this class once all coercions will be handled by Sage's coercion model

    TESTS::

        sage: TestSuite(SymmetricFunctions(QQ).s()).run()
        sage: TestSuite(SymmetricFunctions(QQ).h()).run()
        sage: TestSuite(SymmetricFunctions(QQ).m()).run()
        sage: TestSuite(SymmetricFunctions(QQ).e()).run()
        sage: TestSuite(SymmetricFunctions(QQ).p()).run()
    """

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``, if coercion failed

        INPUT:

        - ``self`` -- a basis of the symmetric functions
        - ``x`` -- an element of the symmetric functions

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s(2)
            2*s[]
            sage: s([2,1]) # indirect doctest
            s[2, 1]

            sage: s([[2,1],[1]])
            s[1, 1] + s[2]
            sage: s([[],[]])
            s[]

            sage: McdJ = SymmetricFunctions(QQ['q','t'].base_ring())
            sage: s = SymmetricFunctions(McdJ.base_ring()).s()
            sage: s._element_constructor_(McdJ(s[2,1]))
            s[2, 1]
        """
        R = self.base_ring()

        eclass = self.element_class
        if isinstance(x, int):
            x = Integer(x)


        ##############
        # Partitions #
        ##############
        if x in sage.combinat.partition.Partitions():
            return eclass(self, {sage.combinat.partition.Partition(filter(lambda x: x!=0, x)): R.one()})

        # Todo: discard all of this which is taken care by Sage's coercion
        # (up to changes of base ring)

        ##############
        # Dual bases #
        ##############
        elif sfa.is_SymmetricFunction(x) and hasattr(x, 'dual'):
            #Check to see if it is the dual of some other basis
            #If it is, try to coerce its corresponding element
            #in the other basis
            return self(x.dual())

        ##################################################################
        # Symmetric Functions, same basis, possibly different coeff ring #
        ##################################################################

        # self.Element is used below to test if another symmetric
        # function is expressed in the same basis but in another
        # ground ring.  This idiom is fragile and depends on the
        # internal (unstable) specifications of parents and categories
        #
        # TODO: find the right idiom
        #
        # One cannot use anymore self.element_class: it is build by
        # the category mechanism, and depends on the coeff ring.

        elif isinstance(x, self.Element):
            P = x.parent()
            #same base ring
            if P is self:
                return x
            #different base ring
            else:
                return eclass(self, dict([ (e1,R(e2)) for e1,e2 in x._monomial_coefficients.items()]))

        ##################################################
        # Classical Symmetric Functions, different basis #
        ##################################################
        elif isinstance(x, SymmetricFunctionAlgebra_classical.Element):


            R = self.base_ring()
            xP = x.parent()
            xm = x.monomial_coefficients()

            #determine the conversion function.
            try:
                t = conversion_functions[(xP.basis_name(),self.basis_name())]
            except AttributeError:
                raise TypeError, "do not know how to convert from %s to %s"%(xP.basis_name(), self.basis_name())

            if R == QQ and xP.base_ring() == QQ:
                if xm:
                    return self._from_dict(t(xm)._monomial_coefficients, coerce=True)
                else:
                    return self.zero()
            else:
                f = lambda part: self._from_dict(t( {part: ZZ.one()} )._monomial_coefficients)
                return self._apply_module_endomorphism(x, f)


        ###############################
        # Hall-Littlewood Polynomials #
        ###############################
        elif isinstance(x, hall_littlewood.HallLittlewood_generic.Element):
            #
            #Qp: Convert to Schur basis and then convert to self
            #
            if isinstance(x, hall_littlewood.HallLittlewood_qp.Element):
                Qp = x.parent()
                sx = Qp._s._from_cache(x, Qp._s_cache, Qp._self_to_s_cache, t=Qp.t)
                return self(sx)
            #
            #P: Convert to Schur basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewood_p.Element):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            #
            #Q: Convert to P basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewood_q.Element):
                return self( x.parent()._P( x ) )

        #######
        # LLT #
        #######
        #Convert to m and then to self.
        elif isinstance(x, llt.LLT_generic.Element):
            P = x.parent()
            BR = self.base_ring()
            zero = BR.zero()
            PBR = P.base_ring()
            if not BR.has_coerce_map_from(PBR):
                raise TypeError, "no coerce map from x's parent's base ring (= %s) to self's base ring (= %s)"%(PBR, self.base_ring())

            z_elt = {}
            for m, c in x._monomial_coefficients.iteritems():
                n = sum(m)
                P._m_cache(n)
                for part in P._self_to_m_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*P._self_to_m_cache[n][m][part].subs(t=P.t))

            m = P._sym.monomial()
            return self( m._from_dict(z_elt) )

        #########################
        # Macdonald Polynomials #
        #########################
        elif isinstance(x, macdonald.MacdonaldPolynomials_generic.Element):
            if isinstance(x, macdonald.MacdonaldPolynomials_j.Element):
                J = x.parent()
                sx = J._s._from_cache(x, J._s_cache, J._self_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomials_q.Element, macdonald.MacdonaldPolynomials_p.Element)):
                J = x.parent()._J
                jx = J(x)
                sx = J._s._from_cache(jx, J._s_cache, J._self_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomials_h.Element,macdonald.MacdonaldPolynomials_ht.Element)):
                H = x.parent()
                sx = H._s._from_cache(x, H._s_cache, H._self_to_s_cache, q=H.q, t=H.t)
                return self(sx)
            elif isinstance(x, macdonald.MacdonaldPolynomials_s.Element):
                S = x.parent()
                sx = S._s._from_cache(x, S._s_cache, S._self_to_s_cache, q=S.q, t=S.t)
                return self(sx)
            else:
                raise TypeError

        ####################
        # Jack Polynomials #
        ####################
        elif isinstance(x, jack.JackPolynomials_generic.Element):
            if isinstance(x, jack.JackPolynomials_p.Element):
                P = x.parent()
                mx = P._m._from_cache(x, P._m_cache, P._self_to_m_cache, t=P.t)
                return self(mx)
            if isinstance(x, (jack.JackPolynomials_j.Element, jack.JackPolynomials_q.Element)):
                return self( x.parent()._P(x) )
            else:
                raise TypeError

        #####################
        # k-Schur Functions #
        #####################
        if isinstance(x, kschur.kSchurFunctions_generic.Element):
            if isinstance(x, kschur.kSchurFunctions_t.Element):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            else:
                raise TypeError

        ####################################################
        # Bases defined by orthogonality and triangularity #
        ####################################################
        elif isinstance(x, orthotriang.SymmetricFunctionAlgebra_orthotriang.Element):
            #Convert to its base and then to self
            xp = x.parent()
            if self is xp._sf_base:
                return xp._sf_base._from_cache(x, xp._base_cache, xp._self_to_base_cache)
            else:
                return self( xp._sf_base(x) )

        ###################
        # Skew Partitions #
        ###################
        elif x in sage.combinat.skew_partition.SkewPartitions():
            import sage.libs.lrcalc.lrcalc as lrcalc
            skewschur = lrcalc.skew(x[0], x[1])
            return self._from_dict(skewschur)


        #############################
        # Elements of the base ring #
        #############################
        elif x.parent() is R:
            return eclass(self, {sage.combinat.partition.Partition([]):x})

        ###########################################
        # Elements that coerce into the base ring #
        ###########################################
        elif R.has_coerce_map_from(x.parent()):
            return eclass(self, {sage.combinat.partition.Partition([]):R(x)})

        #################################
        # Last shot -- try calling R(x) #
        #################################
        else:
            try:
                return eclass(self, {sage.combinat.partition.Partition([]):R(x)})
            except StandardError:
                raise TypeError, "do not know how to make x (= %s) an element of self"%(x)

    # This subclass is currently needed for the test above:
    #    isinstance(x, SymmetricFunctionAlgebra_classical.Element):
    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        """
        A symmetric function.
        """
        pass
