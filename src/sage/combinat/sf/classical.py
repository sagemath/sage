# sage.doctest: needs sage.combinat sage.modules
"""
Classical symmetric functions
"""
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.combinat.partition import _Partitions


from . import hall_littlewood
from . import sfa
from . import llt
from . import macdonald
from . import jack
from . import orthotriang

translate = {'monomial': 'MONOMIAL',
             'homogeneous': 'HOMSYM',
             'powersum': 'POWSYM',
             'elementary': 'ELMSYM',
             'Schur': 'SCHUR'}

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

    The following checks if the bug described in :issue:`15312` is fixed. ::

        sage: change = sage.combinat.sf.classical.conversion_functions[('powersum', 'Schur')]
        sage: hideme = change({Partition([1]*47):ZZ(1)}) # long time
        sage: change({Partition([2,2]):QQ(1)})
        s[1, 1, 1, 1] - s[2, 1, 1] + 2*s[2, 2] - s[3, 1] + s[4]
    """
    import sage.libs.symmetrica.all as symmetrica
    for other_basis, other_name in translate.items():
        for basis, name in translate.items():
            try:
                conversion_functions[(other_basis, basis)] = getattr(symmetrica,
                        f't_{other_name}_{name}')
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
        Convert ``x`` into ``self``, if coercion failed.

        INPUT:

        - ``x`` -- an element of the symmetric functions

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s(2)
            2*s[]
            sage: s([2,1]) # indirect doctest
            s[2, 1]

            sage: McdJ = SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().J()
            sage: s = SymmetricFunctions(McdJ.base_ring()).s()
            sage: s._element_constructor_(McdJ(s[2,1]))
            s[2, 1]

        TESTS:

        Check that non-Schur bases raise an error when given skew partitions
        (:issue:`19218`)::

            sage: e = SymmetricFunctions(QQ).e()
            sage: e([[2,1],[1]])
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= [[2, 1], [1]]) an element of self

        Check that :issue:`34576` is fixed::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: f = s(0/2); f
            0
            sage: f == 0
            True
            sage: f._monomial_coefficients
            {}

            sage: s2 = SymmetricFunctions(GF(2)).s()
            sage: f = s2(2*s[2,1]); f
            0
            sage: f == 0
            True
            sage: f._monomial_coefficients
            {}
        """
        R = self.base_ring()

        eclass = self.element_class
        if isinstance(x, int):
            x = Integer(x)

        ##############
        # Partitions #
        ##############
        if x in _Partitions:
            return eclass(self, {_Partitions(x): R.one()})

        # Todo: discard all of this which is taken care by Sage's coercion
        # (up to changes of base ring)

        ##############
        # Dual bases #
        ##############
        elif isinstance(x, sfa.SymmetricFunctionAlgebra_generic.Element) and hasattr(x, 'dual'):
            # Check to see if it is the dual of some other basis
            # If it is, try to coerce its corresponding element
            # in the other basis
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
            # same base ring
            if P is self:
                return x
            # different base ring
            else:
                return eclass(self, {la: rc for la, c in x._monomial_coefficients.items()
                                     if (rc := R(c))})

        ##################################################
        # Classical Symmetric Functions, different basis #
        ##################################################
        elif isinstance(x, SymmetricFunctionAlgebra_classical.Element):

            P = x.parent()
            m = x.monomial_coefficients()

            # determine the conversion function.
            try:
                t = conversion_functions[(P.basis_name(), self.basis_name())]
            except AttributeError:
                raise TypeError("do not know how to convert from %s to %s"
                                % (P.basis_name(), self.basis_name()))

            if R == QQ and P.base_ring() == QQ:
                if m:
                    return self._from_dict(t(m)._monomial_coefficients,
                                           coerce=True)
                return self.zero()
            else:
                f = lambda part: self._from_dict(t({part: ZZ.one()})._monomial_coefficients)
                return self._apply_module_endomorphism(x, f)

        ###############################
        # Hall-Littlewood Polynomials #
        ###############################
        elif isinstance(x, hall_littlewood.HallLittlewood_generic.Element):
            #
            # Qp: Convert to Schur basis and then convert to self
            #
            if isinstance(x, (hall_littlewood.HallLittlewood_qp.Element,
                              hall_littlewood.HallLittlewood_p.Element)):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            #
            # Q: Convert to P basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewood_q.Element):
                return self(x.parent()._P(x))

        #######
        # LLT #
        #######
        # Convert to m and then to self.
        elif isinstance(x, llt.LLT_generic.Element):
            P = x.parent()
            Rx = P.base_ring()
            zero = R.zero()
            if not R.has_coerce_map_from(Rx):
                raise TypeError("no coerce map from x's parent's base ring (= %s) to self's base ring (= %s)"
                                % (Rx, R))

            z_elt = {}
            for m, c in x._monomial_coefficients.items():
                n = sum(m)
                P._m_cache(n)
                for part in P._self_to_m_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + R(c*P._self_to_m_cache[n][m][part].subs(t=P.t))

            m = P._sym.monomial()
            return self(m._from_dict(z_elt))

        #########################
        # Macdonald Polynomials #
        #########################
        elif isinstance(x, macdonald.MacdonaldPolynomials_generic.Element):
            if isinstance(x, (macdonald.MacdonaldPolynomials_j.Element,
                              macdonald.MacdonaldPolynomials_s.Element)):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, q=P.q, t=P.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomials_q.Element,
                                macdonald.MacdonaldPolynomials_p.Element)):
                J = x.parent()._J
                jx = J(x)
                sx = J._s._from_cache(jx, J._s_cache, J._self_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomials_h.Element,
                                macdonald.MacdonaldPolynomials_ht.Element)):
                P = x.parent()
                sx = P._self_to_s(x)
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
            if isinstance(x, (jack.JackPolynomials_j.Element,
                              jack.JackPolynomials_q.Element)):
                return self(x.parent()._P(x))
            else:
                raise TypeError

        ####################################################
        # Bases defined by orthogonality and triangularity #
        ####################################################
        elif isinstance(x, orthotriang.SymmetricFunctionAlgebra_orthotriang.Element):
            # Convert to its base and then to self
            P = x.parent()
            if self is P._sf_base:
                return P._sf_base._from_cache(x, P._base_cache, P._self_to_base_cache)
            else:
                return self( P._sf_base(x) )

        #################################
        # Last shot -- try calling R(x) #
        #################################
        else:
            try:
                c = R(x)
            except (TypeError, ValueError):
                raise TypeError("do not know how to make x (= {}) an element of self".format(x))
            else:
                if not c:
                    return self.zero()
                return eclass(self, {_Partitions([]): c})

    # This subclass is currently needed for the test above:
    #    isinstance(x, SymmetricFunctionAlgebra_classical.Element):
    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        """
        A symmetric function.
        """
        pass
