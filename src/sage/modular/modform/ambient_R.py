"""
Modular Forms over a Non-minimal Base Ring
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from . import ambient
from .cuspidal_submodule import CuspidalSubmodule_R
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

class ModularFormsAmbient_R(ambient.ModularFormsAmbient):
    def __init__(self, M, base_ring):
        """
        Ambient space of modular forms over a ring other than QQ.

        EXAMPLES::

            sage: M = ModularForms(23,2,base_ring=GF(7)) ## indirect doctest
            sage: M
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(23) of weight 2 over Finite Field of size 7
            sage: M == loads(dumps(M))
            True
        """
        self.__M = M
        if M.character() is not None:
            self.__R_character = M.character().change_ring(base_ring)
        else:
            self.__R_character = None
        ambient.ModularFormsAmbient.__init__(self, M.group(), M.weight(), base_ring, M.character(), M._eis_only)

    @cached_method(key=lambda self,sign: ZZ(sign)) # convert sign to an Integer before looking this up in the cache
    def modular_symbols(self,sign=0):
        r"""
        Return the space of modular symbols attached to this space, with the given sign (default 0).

        TESTS::

            sage: K.<i> = QuadraticField(-1)
            sage: chi = DirichletGroup(5, base_ring = K).0
            sage: L.<c> = K.extension(x^2 - 402*i)
            sage: M = ModularForms(chi, 7, base_ring = L)
            sage: symbs = M.modular_symbols()
            sage: symbs.character() == chi
            True
            sage: symbs.base_ring() == L
            True
        """
        sign = ZZ(sign)
        return self.__M.modular_symbols(sign).change_ring(self.base_ring())

    def _repr_(self):
        """
        String representation for self.

        EXAMPLES::

            sage: M = ModularForms(23,2,base_ring=GF(7)) ## indirect doctest
            sage: M._repr_()
            'Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(23) of weight 2 over Finite Field of size 7'

            sage: chi = DirichletGroup(109).0 ** 36
            sage: ModularForms(chi, 2, base_ring = chi.base_ring())
            Modular Forms space of dimension 9, character [zeta3] and weight 2 over Cyclotomic Field of order 108 and degree 36
        """
        s = str(self.__M)
        i = s.find('over')
        if i != -1:
            s = s[:i]
        return s + 'over %s'%self.base_ring()

    def _compute_q_expansion_basis(self, prec=None):
        """
        Compute q-expansions for a basis of self to precision prec.

        EXAMPLES::

            sage: M = ModularForms(23,2,base_ring=GF(7))
            sage: M._compute_q_expansion_basis(10)
            [q + 6*q^3 + 6*q^4 + 5*q^6 + 2*q^7 + 6*q^8 + 2*q^9 + O(q^10),
            q^2 + 5*q^3 + 6*q^4 + 2*q^5 + q^6 + 2*q^7 + 5*q^8 + O(q^10),
            1 + 5*q^3 + 5*q^4 + 5*q^6 + 3*q^8 + 5*q^9 + O(q^10)]

        TESTS:

        This checks that :trac:`13445` is fixed::

            sage: M = ModularForms(Gamma1(29), base_ring=GF(29))
            sage: S = M.cuspidal_subspace()
            sage: 0 in [f.valuation() for f in S.basis()]
            False
            sage: from sage.modular.dims import dimension_cusp_forms
            sage: len(S.basis()) == dimension_cusp_forms(Gamma1(29), 2)
            True
        """
        if prec is None:
            prec = self.prec()
        R = self._q_expansion_ring()
        c = self.base_ring().characteristic()
        if c == 0:
            B = self.__M.q_expansion_basis(prec)
            return [R(f) for f in B]
        elif c.is_prime_power():
            K = self.base_ring()
            p = K.characteristic().prime_factors()[0]
            from sage.rings.finite_rings.finite_field_constructor import GF
            Kp = GF(p)
            newB = [f.change_ring(K) for f in list(self.__M.cuspidal_subspace().q_integral_basis(prec))]
            A = Kp**prec
            gens = [f.padded_list(prec) for f in newB]
            V = A.span(gens)
            B = [f.change_ring(K) for f in self.__M.q_integral_basis(prec)]
            for f in B:
                fc = f.padded_list(prec)
                gens.append(fc)
                if not A.span(gens) == V:
                    newB.append(f)
                    V = A.span(gens)
            if len(newB) != self.dimension():
                raise RuntimeError("The dimension of the space is %s but the basis we computed has %s elements"%(self.dimension(), len(newB)))
            lst = [R(f) for f in newB]
            return [f/f[f.valuation()] for f in lst]
        else:
            # this returns a basis of q-expansions, without guaranteeing that
            # the first vectors form a basis of the cuspidal subspace
            # TODO: bring this in line with the other cases
            # simply using the above code fails because free modules over
            # general rings do not have a .span() method
            B = self.__M.q_integral_basis(prec)
            return [R(f) for f in B]

    def cuspidal_submodule(self):
        r"""
        Return the cuspidal subspace of this space.

        EXAMPLES::

            sage: C = CuspForms(7, 4, base_ring=CyclotomicField(5)) # indirect doctest
            sage: type(C)
            <class 'sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_R_with_category'>
        """
        return CuspidalSubmodule_R(self)


    def change_ring(self, R):
        r"""
        Return this modular forms space with the base ring changed to the ring R.

        EXAMPLES::

            sage: chi = DirichletGroup(109, CyclotomicField(3)).0
            sage: M9 = ModularForms(chi, 2, base_ring = CyclotomicField(9))
            sage: M9.change_ring(CyclotomicField(15))
            Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2 over Cyclotomic Field of order 15 and degree 8
            sage: M9.change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: Space cannot be defined over Rational Field
        """
        if not R.has_coerce_map_from(self.__M.base_ring()):
            raise ValueError("Space cannot be defined over %s" % R)
        return ModularFormsAmbient_R(self.__M, R)
