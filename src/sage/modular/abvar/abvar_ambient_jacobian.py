"""
Ambient Jacobian Abelian Variety

TESTS::

    sage: loads(dumps(J0(37))) == J0(37)
    True
    sage: loads(dumps(J1(13))) == J1(13)
    True
"""

import weakref
from sage.structure.sequence import Sequence

from .abvar import (ModularAbelianVariety_modsym_abstract,
                    simple_factorization_of_modsym_space, modsym_lattices,
                    ModularAbelianVariety_modsym)
from sage.rings.rational_field import QQ

from sage.modular.modsym.modsym import ModularSymbols
from sage.modular.modform.constructor import Newforms
from sage.modular.arithgroup.all import is_Gamma0, is_Gamma1
from . import morphism


_cache = {}


def ModAbVar_ambient_jacobian(group):
    """
    Return the ambient Jacobian attached to a given congruence
    subgroup.

    The result is cached using a weakref. This function is called
    internally by modular abelian variety constructors.

    INPUT:


    -  ``group`` - a congruence subgroup.


    OUTPUT: a modular abelian variety attached

    EXAMPLES::

        sage: import sage.modular.abvar.abvar_ambient_jacobian as abvar_ambient_jacobian
        sage: A = abvar_ambient_jacobian.ModAbVar_ambient_jacobian(Gamma0(11))
        sage: A
        Abelian variety J0(11) of dimension 1
        sage: B = abvar_ambient_jacobian.ModAbVar_ambient_jacobian(Gamma0(11))
        sage: A is B
        True

    You can get access to and/or clear the cache as follows::

        sage: abvar_ambient_jacobian._cache = {}
        sage: B = abvar_ambient_jacobian.ModAbVar_ambient_jacobian(Gamma0(11))
        sage: A is B
        False
    """
    try:
        X = _cache[group]()
        if X is not None:
            return X
    except KeyError:
        pass
    X = ModAbVar_ambient_jacobian_class(group)
    _cache[group] = weakref.ref(X)
    return X


class ModAbVar_ambient_jacobian_class(ModularAbelianVariety_modsym_abstract):
    """
    An ambient Jacobian modular abelian variety attached to a
    congruence subgroup.
    """
    def __init__(self, group):
        """
        Create an ambient Jacobian modular abelian variety.

        EXAMPLES::

            sage: A = J0(37); A
            Abelian variety J0(37) of dimension 2
            sage: type(A)
            <class 'sage.modular.abvar.abvar_ambient_jacobian.ModAbVar_ambient_jacobian_class_with_category'>
            sage: A.group()
            Congruence Subgroup Gamma0(37)
        """
        ModularAbelianVariety_modsym_abstract.__init__(self, (group,), QQ)
        self.__group = group
        self._is_hecke_stable = True

    def _modular_symbols(self):
        """
        Return the modular symbols space associated to this ambient
        Jacobian.

        OUTPUT: modular symbols space

        EXAMPLES::

            sage: M = J0(33)._modular_symbols(); M
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: J0(33)._modular_symbols() is M
            True
        """
        try:
            return self.__modsym
        except AttributeError:
            self.__modsym = ModularSymbols(self.__group, weight=2).cuspidal_submodule()
            return self.__modsym

    def _repr_(self):
        """
        Return string representation of this Jacobian modular abelian
        variety.

        EXAMPLES::

            sage: A = J0(11); A
            Abelian variety J0(11) of dimension 1
            sage: A._repr_()
            'Abelian variety J0(11) of dimension 1'
            sage: A.rename("J_0(11)")
            sage: A
            J_0(11)

        We now clear the cache to get rid of our renamed
        `J_0(11)`.

        ::

            sage: import sage.modular.abvar.abvar_ambient_jacobian as abvar_ambient_jacobian
            sage: abvar_ambient_jacobian._cache = {}
        """
        return 'Abelian variety %s of dimension %s%s' % (self._ambient_repr(),
                                                         self.dimension(),
                                    '' if self.base_field() == QQ else ' over %s' % self.base_field())

    def _latex_(self):
        """
        Return Latex representation of ``self``.

        EXAMPLES::

            sage: latex(J0(37))
            J_0(37)
            sage: J1(13)._latex_()
            'J_1(13)'
            sage: latex(JH(389,[16]))
            J_H(389,[16])
        """
        return self._ambient_latex_repr()

    def ambient_variety(self):
        """
        Return the ambient modular abelian variety that contains self.
        Since self is a Jacobian modular abelian variety, this is just
        self.

        OUTPUT: abelian variety

        EXAMPLES::

            sage: A = J0(17)
            sage: A.ambient_variety()
            Abelian variety J0(17) of dimension 1
            sage: A is A.ambient_variety()
            True
        """
        return self

    def group(self):
        """
        Return the group that this Jacobian modular abelian variety is
        attached to.

        EXAMPLES::

            sage: J1(37).group()
            Congruence Subgroup Gamma1(37)
            sage: J0(5077).group()
            Congruence Subgroup Gamma0(5077)
            sage: J = GammaH(11,[3]).modular_abelian_variety(); J
            Abelian variety JH(11,[3]) of dimension 1
            sage: J.group()
            Congruence Subgroup Gamma_H(11) with H generated by [3]
        """
        return self.__group

    def groups(self):
        """
        Return the tuple of congruence subgroups attached to this ambient
        Jacobian. This is always a tuple of length 1.

        OUTPUT: tuple

        EXAMPLES::

            sage: J0(37).groups()
            (Congruence Subgroup Gamma0(37),)
        """
        return (self.__group,)

    def _calculate_endomorphism_generators(self):
        """
        Calculate generators for the endomorphism ring of self.

        EXAMPLES::

            sage: J0(11)._calculate_endomorphism_generators()
            [Abelian variety endomorphism of Abelian variety J0(11) of dimension 1]
            sage: ls = J0(46)._calculate_endomorphism_generators() ; ls
            [Abelian variety endomorphism of Abelian variety J0(46) of dimension 5,
             Abelian variety endomorphism of Abelian variety J0(46) of dimension 5,
             Abelian variety endomorphism of Abelian variety J0(46) of dimension 5,
             Abelian variety endomorphism of Abelian variety J0(46) of dimension 5,
             Abelian variety endomorphism of Abelian variety J0(46) of dimension 5]
            sage: len(ls) == J0(46).dimension()
            True
        """
        D = self.decomposition()
        phi = self._isogeny_to_product_of_simples()
        psi = phi.complementary_isogeny()

        m1 = phi.matrix()
        m2 = psi.matrix()

        H = self.Hom(self)
        M = H.matrix_space()

        ls = []
        ind = 0
        for d in D:
            to_newform = d._isogeny_to_newform_abelian_variety()
            n1 = to_newform.matrix()
            n2 = to_newform.complementary_isogeny().matrix()
            f_gens = to_newform.codomain()._calculate_endomorphism_generators()
            small_space = to_newform.parent().matrix_space()
            f_gens = [small_space(x.list()) for x in f_gens]
            for m in f_gens:
                mat = H.matrix_space()(0)
                mat.set_block(ind, ind, n1 * m * n2)
                ls.append((m1 * mat * m2).list())
            ind += 2 * d.dimension()

        return [H(morphism.Morphism(H, M(x))) for x in ls]

    def degeneracy_map(self, level, t=1, check=True):
        """
        Return the t-th degeneracy map from self to J(level). Here t must
        be a divisor of either level/self.level() or self.level()/level.

        INPUT:


        -  ``level`` - integer (multiple or divisor of level of
           self)

        -  ``t`` - divisor of quotient of level of self and
           level

        -  ``check`` - bool (default: True); if True do some
           checks on the input


        OUTPUT: a morphism

        EXAMPLES::

            sage: J0(11).degeneracy_map(33)
            Degeneracy map from Abelian variety J0(11) of dimension 1 to Abelian variety J0(33) of dimension 3 defined by [1]
            sage: J0(11).degeneracy_map(33).matrix()
            [ 0 -3  2  1 -2  0]
            [ 1 -2  0  1  0 -1]
            sage: J0(11).degeneracy_map(33,3).matrix()
            [-1  0  0  0  1 -2]
            [-1 -1  1 -1  1  0]
            sage: J0(33).degeneracy_map(11,1).matrix()
            [ 0  1]
            [ 0 -1]
            [ 1 -1]
            [ 0  1]
            [-1  1]
            [ 0  0]
            sage: J0(11).degeneracy_map(33,1).matrix() * J0(33).degeneracy_map(11,1).matrix()
            [4 0]
            [0 4]
        """
        if check:
            if (level % self.level()) and (self.level() % level):
                raise ValueError("level must be divisible by level of self")
            if (max(level, self.level()) // min(self.level(), level)) % t:
                raise ValueError("t must divide the quotient of the two levels")

        Mself = self.modular_symbols()
        Jdest = (type(Mself.group()))(level).modular_abelian_variety()
        Mdest = Jdest.modular_symbols()

        symbol_map = Mself.degeneracy_map(level, t).restrict_codomain(Mdest)
        H = self.Hom(Jdest)

        return H(morphism.DegeneracyMap(H, symbol_map.matrix(), [t]))

    def dimension(self):
        """
        Return the dimension of this modular abelian variety.

        EXAMPLES::

            sage: J0(2007).dimension()
            221
            sage: J1(13).dimension()
            2
            sage: J1(997).dimension()
            40920
            sage: J0(389).dimension()
            32
            sage: JH(389,[4]).dimension()
            64
            sage: J1(389).dimension()
            6112
        """
        try:
            return self._dimension
        except AttributeError:
            d = self.group().genus()
            self._dimension = d
            return d

    def decomposition(self, simple=True, bound=None):
        """
        Decompose this ambient Jacobian as a product of abelian
        subvarieties, up to isogeny.

        EXAMPLES::

            sage: J0(33).decomposition(simple=False)
            [
            Abelian subvariety of dimension 2 of J0(33),
            Abelian subvariety of dimension 1 of J0(33)
            ]
            sage: J0(33).decomposition(simple=False)[1].is_simple()
            True
            sage: J0(33).decomposition(simple=False)[0].is_simple()
            False
            sage: J0(33).decomposition(simple=False)
            [
            Abelian subvariety of dimension 2 of J0(33),
            Simple abelian subvariety 33a(None,33) of dimension 1 of J0(33)
            ]
            sage: J0(33).decomposition(simple=True)
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ]
        """
        try:
            return self.__decomposition[simple]
        except KeyError:
            pass
        except AttributeError:
            self.__decomposition = {}

        M = self.modular_symbols().ambient_module()
        level = M.level()
        group = M.group()
        factors = simple_factorization_of_modsym_space(M, simple=simple)
        factors = modsym_lattices(M, factors)

        D = []
        is_simple = True if simple else None
        for newform_level, isogeny_number, number, modsym, lattice in factors:
            A = ModularAbelianVariety_modsym(modsym, lattice=lattice,
                               newform_level=(newform_level, group),
                                             is_simple=is_simple,
                                             isogeny_number=isogeny_number,
                                             number=(number, level),
                                             check=False)
            D.append(A)

            # This line below could be safely deleted.  It basically creates a circular
            # reference so that say J0(389)[0] + J0(389)[1] doesn't do two separate
            # decompositions.  Memory will be freed though, at least if you do
            # import gc; gc.collect().
            A._ambient = self

        D.sort()
        D = Sequence(D, immutable=True, cr=True, universe=self.category())
        self.__decomposition[simple] = D
        return D

    def newform_decomposition(self, names=None):
        """
        Return the newforms of the simple subvarieties in the decomposition of
        self as a product of simple subvarieties, up to isogeny.

        OUTPUT:

        - an array of newforms

        EXAMPLES::

            sage: J0(81).newform_decomposition('a')
            [q - 2*q^4 + O(q^6), q - 2*q^4 + O(q^6), q + a0*q^2 + q^4 - a0*q^5 + O(q^6)]

            sage: J1(19).newform_decomposition('a')
            [q - 2*q^3 - 2*q^4 + 3*q^5 + O(q^6),
             q + a1*q^2 + (-1/9*a1^5 - 1/3*a1^4 - 1/3*a1^3 + 1/3*a1^2 - a1 - 1)*q^3 + (4/9*a1^5 + 2*a1^4 + 14/3*a1^3 + 17/3*a1^2 + 6*a1 + 2)*q^4 + (-2/3*a1^5 - 11/3*a1^4 - 10*a1^3 - 14*a1^2 - 15*a1 - 9)*q^5 + O(q^6)]
        """
        if self.dimension() == 0:
            return []
        G = self.group()
        if not (is_Gamma0(G) or is_Gamma1(G)):
            return [S.newform(names=names) for S in self.decomposition()]
        Gtype = G.parent()
        N = G.level()
        preans = [Newforms(Gtype(d), names=names) * len((N // d).divisors())
                  for d in N.divisors()]
        return [newform for l in preans for newform in l]
