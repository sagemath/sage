# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
"""
Abelian varieties attached to newforms

TESTS::

    sage: A = AbelianVariety('23a')
    sage: loads(dumps(A)) == A
    True
"""
###########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  https://www.gnu.org/licenses/                           #
###########################################################################
from sage.misc.lazy_import import lazy_import

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.modular.modform.element import Newform
from sage.modular.arithgroup.all import Gamma0_class, Gamma1_class, GammaH_class

from .abvar import ModularAbelianVariety_modsym_abstract
from sage.modular.abvar import homspace
lazy_import('sage.databases.cremona', 'cremona_letter_code')


class ModularAbelianVariety_newform(ModularAbelianVariety_modsym_abstract):
    """
    A modular abelian variety attached to a specific newform.
    """
    def __init__(self, f, internal_name=False):
        """
        Create the modular abelian variety `A_f` attached to the
        newform `f`.

        INPUT:

        - ``f`` -- a newform

        EXAMPLES::

            sage: f = CuspForms(37).newforms('a')[0]
            sage: f.abelian_variety()
            Newform abelian subvariety 37a of dimension 1 of J0(37)

            sage: AbelianVariety(Newforms(1, 12)[0])
            Traceback (most recent call last):
            ...
            TypeError: f must have weight 2
        """
        if not isinstance(f, Newform):
            raise TypeError("f must be a newform")
        if f.weight() != 2:
            raise TypeError("f must have weight 2")
        self.__f = f
        self._is_hecke_stable = True
        K = f.qexp().base_ring()
        if K == QQ:
            variable_name = None
        else:
            variable_name = K.variable_name()
        self.__named_newforms = {variable_name: self.__f}
        if not internal_name:
            self.__named_newforms[None] = self.__f
        ModularAbelianVariety_modsym_abstract.__init__(self, (f.group(),), QQ,
            is_simple=True, newform_level=(f.level(), f.group()),
            isogeny_number=f.number(), number=0)

    def _modular_symbols(self, sign=0):
        """
        EXAMPLES::

            sage: f = CuspForms(52).newforms('a')[0]
            sage: A = f.abelian_variety()
            sage: A._modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 15 for Gamma_0(52) of weight 2 with sign 0 over Rational Field
            sage: A._modular_symbols(1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 10 for Gamma_0(52) of weight 2 with sign 1 over Rational Field
            sage: A._modular_symbols(-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(52) of weight 2 with sign -1 over Rational Field
        """
        return self.__f.modular_symbols(sign=sign)

    def newform(self, names=None):
        r"""
        Return the newform that this modular abelian variety is attached to.

        EXAMPLES::

            sage: f = Newform('37a')
            sage: A = f.abelian_variety()
            sage: A.newform()
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: A.newform() is f
            True

        If the a variable name has not been specified, we must specify one::

            sage: A = AbelianVariety('67b')
            sage: A.newform()
            Traceback (most recent call last):
            ...
            TypeError: You must specify the name of the generator.
            sage: A.newform('alpha')
            q + alpha*q^2 + (-alpha - 3)*q^3 + (-3*alpha - 3)*q^4 - 3*q^5 + O(q^6)

        If the eigenform is actually over `\QQ` then we don't have to specify
        the name::

            sage: A = AbelianVariety('67a')
            sage: A.newform()
            q + 2*q^2 - 2*q^3 + 2*q^4 + 2*q^5 + O(q^6)
        """
        try:
            return self.__named_newforms[names]
        except KeyError:
            self.__named_newforms[names] = Newform(self.__f.parent().change_ring(QQ), self.__f.modular_symbols(1), names=names, check=False)
            return self.__named_newforms[names]

    def label(self) -> str:
        """
        Return canonical label that defines this newform modular
        abelian variety.

        OUTPUT: string

        EXAMPLES::

            sage: A = AbelianVariety('43b')
            sage: A.label()
            '43b'
        """
        G = self.__f.group()
        if isinstance(G, Gamma0_class):
            group = ''
        elif isinstance(G, Gamma1_class):
            group = 'G1'
        elif isinstance(G, GammaH_class):
            group = 'GH[' + ','.join(str(z) for z in G._generators_for_H()) + ']'
        return '%s%s%s' % (self.level(), cremona_letter_code(self.factor_number()), group)

    def factor_number(self):
        """
        Return factor number.

        OUTPUT: int

        EXAMPLES::

            sage: A = AbelianVariety('43b')
            sage: A.factor_number()
            1
        """
        try:
            return self.__factor_number
        except AttributeError:
            self.__factor_number = self.__f.number()
            return self.__factor_number

    def _repr_(self) -> str:
        """
        String representation of this modular abelian variety.

        EXAMPLES::

            sage: AbelianVariety('37a')._repr_()
            'Newform abelian subvariety 37a of dimension 1 of J0(37)'
        """
        return "Newform abelian subvariety %s of dimension %s of %s" % (
            self.newform_label(), self.dimension(), self._ambient_repr())

    def endomorphism_ring(self):
        """
        Return the endomorphism ring of this newform abelian variety.

        EXAMPLES::

            sage: A = AbelianVariety('23a')
            sage: E = A.endomorphism_ring(); E
            Endomorphism ring of Newform abelian subvariety 23a of dimension 2 of J0(23)

        We display the matrices of these two basis matrices::

            sage: E.0.matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: E.1.matrix()
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]

        The result is cached::

            sage: E is A.endomorphism_ring()
            True
        """
        try:
            return self.__endomorphism_ring
        except AttributeError:
            pass

        E = homspace.EndomorphismSubring(self)
        self.__endomorphism_ring = E
        return self.__endomorphism_ring

    def _calculate_endomorphism_generators(self):
        """
        EXAMPLES::

            sage: A = AbelianVariety('43b')
            sage: B = A.endomorphism_ring(); B   # indirect doctest
            Endomorphism ring of Newform abelian subvariety 43b of dimension 2 of J0(43)
            sage: [b.matrix() for b in B.gens()]
            [
            [1 0 0 0]  [ 0  1  0  0]
            [0 1 0 0]  [ 1 -2  0  0]
            [0 0 1 0]  [ 1  0 -2 -1]
            [0 0 0 1], [ 0  1 -1  0]
            ]
        """
        M = self.modular_symbols()
        bound = M.sturm_bound()

        d = self.dimension()
        T1list = self.hecke_operator(1).matrix().list()
        EndVecZ = ZZ**(len(T1list))
        V = EndVecZ.submodule([T1list])
        n = 2

        while V.dimension() < d:
            W = EndVecZ.submodule([((self.hecke_operator(n).matrix())**i).list()
                                   for i in range(1, d + 1)])
            V = V + W
            n += 1
            if n > bound:
                raise ArithmeticError("Error computing endomorphism generators")

        return V.saturation().basis()
