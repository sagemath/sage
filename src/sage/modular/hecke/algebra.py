"""
Hecke algebras

In Sage a "Hecke algebra" always refers to an algebra of endomorphisms of some
explicit module, rather than the abstract Hecke algebra of double cosets
attached to a subgroup of the modular group.

We distinguish between "anemic Hecke algebras", which are algebras of Hecke
operators whose indices do not divide some integer N (the level), and "full
Hecke algebras", which include Hecke operators coprime to the level. Morphisms
in the category of Hecke modules are not required to commute with the action of
the full Hecke algebra, only with the anemic algebra.
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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


import weakref

import sage.arith.all as arith
import sage.rings.infinity
import sage.misc.latex as latex
import sage.rings.commutative_algebra
from sage.matrix.constructor import matrix
from sage.arith.all import lcm
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ, QQ
from sage.structure.element import Element

def is_HeckeAlgebra(x):
    r"""
    Return True if x is of type HeckeAlgebra.

    EXAMPLES::

        sage: from sage.modular.hecke.algebra import is_HeckeAlgebra
        sage: is_HeckeAlgebra(CuspForms(1, 12).anemic_hecke_algebra())
        True
        sage: is_HeckeAlgebra(ZZ)
        False
    """
    return isinstance(x, HeckeAlgebra_base)

# The basis_matrix stuff here is a workaround for a subtle bug discovered by
# me (David Loeffler) 2009-04-13. The problem is that if one creates two
# subspaces of a Hecke module which are equal as subspaces but have different
# bases, then the caching machinery needs to distinguish between them.

# The AttributeError occurs in two distinct ways: if M is not a free module
# over its base ring, it might not have a basis_matrix method; and for
# SupersingularModule objects, the basis_matrix method exists but raises an
# error -- this is a known bug (#4306).

# See the doctest for the __call__ method below, which tests that this caching
# is working as it should.

_anemic_cache = {}
def AnemicHeckeAlgebra(M):
    r"""
    Return the anemic Hecke algebra associated to the Hecke module
    M. This checks whether or not the object already exists in memory,
    and if so, returns the existing object rather than a new one.

    EXAMPLES::

        sage: CuspForms(1, 12).anemic_hecke_algebra() # indirect doctest
        Anemic Hecke algebra acting on Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field

    We test uniqueness::

        sage: CuspForms(1, 12).anemic_hecke_algebra() is CuspForms(1, 12).anemic_hecke_algebra()
        True

    We can't ensure uniqueness when loading and saving objects from files, but we can ensure equality::

        sage: CuspForms(1, 12).anemic_hecke_algebra() is loads(dumps(CuspForms(1, 12).anemic_hecke_algebra()))
        False
        sage: CuspForms(1, 12).anemic_hecke_algebra() == loads(dumps(CuspForms(1, 12).anemic_hecke_algebra()))
        True
    """

    try:
        k = (M, M.basis_matrix())
    except AttributeError:
        k = M
    if k in _anemic_cache:
        T = _anemic_cache[k]()
        if not (T is None):
            return T
    T = HeckeAlgebra_anemic(M)
    _anemic_cache[k] = weakref.ref(T)
    return T

_cache = {}
def HeckeAlgebra(M):
    """
    Return the full Hecke algebra associated to the Hecke module
    M. This checks whether or not the object already exists in memory,
    and if so, returns the existing object rather than a new one.

    EXAMPLES::

        sage: CuspForms(1, 12).hecke_algebra() # indirect doctest
        Full Hecke algebra acting on Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field

    We test uniqueness::

        sage: CuspForms(1, 12).hecke_algebra() is CuspForms(1, 12).hecke_algebra()
        True

    We can't ensure uniqueness when loading and saving objects from files, but we can ensure equality::

        sage: CuspForms(1, 12).hecke_algebra() is loads(dumps(CuspForms(1, 12).hecke_algebra()))
        False
        sage: CuspForms(1, 12).hecke_algebra() == loads(dumps(CuspForms(1, 12).hecke_algebra()))
        True
    """
    try:
        k = (M, M.basis_matrix())
    except AttributeError:
        k = M
    if k in _cache:
        T = _cache[k]()
        if not (T is None):
            return T
    T = HeckeAlgebra_full(M)
    _cache[k] = weakref.ref(T)
    return T


def _heckebasis(M):
    r"""
    Gives a basis of the hecke algebra of M as a ZZ-module

    INPUT:

    - ``M`` - a hecke module

    OUTPUT:

    - a list of hecke algebra elements represented as matrices

    EXAMPLES::

        sage: M = ModularSymbols(11,2,1)
        sage: sage.modular.hecke.algebra._heckebasis(M)
        [Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
        [1 0]
        [0 1],
        Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
        [0 1]
        [0 5]]
    """
    d = M.rank()
    VV = QQ**(d**2)
    WW = ZZ**(d**2)
    MM = MatrixSpace(QQ,d)
    MMZ = MatrixSpace(ZZ,d)
    S = []; Denom = []; B = []; B1 = []
    for i in xrange(1, M.hecke_bound() + 1):
        v = M.hecke_operator(i).matrix()
        den = v.denominator()
        Denom.append(den)
        S.append(v)
    den = lcm(Denom)
    for m in S:
        B.append(WW((den*m).list()))
    UU = WW.submodule(B)
    B = UU.basis()
    for u in B:
        u1 = u.list()
        m1 = M.hecke_algebra()(MM(u1), check=False)
        #m1 = MM(u1)
        B1.append((1/den)*m1)
    return B1


class HeckeAlgebra_base(sage.rings.commutative_algebra.CommutativeAlgebra):
    """
    Base class for algebras of Hecke operators on a fixed Hecke module.
    """
    def __init__(self, M):
        """
        INPUT:

        -  ``M`` - a Hecke module

        EXAMPLE::

            sage: CuspForms(1, 12).hecke_algebra() # indirect doctest
            Full Hecke algebra acting on Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
        """
        from . import module
        if not module.is_HeckeModule(M):
            raise TypeError("M (=%s) must be a HeckeModule"%M)
        self.__M = M
        sage.rings.commutative_algebra.CommutativeAlgebra.__init__(self, M.base_ring())

    def _an_element_impl(self):
        r"""
        Return an element of this algebra. Used by the coercion machinery.

        EXAMPLE::

            sage: CuspForms(1, 12).hecke_algebra().an_element() # indirect doctest
            Hecke operator T_2 on Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
        """
        return self.hecke_operator(self.level() + 1)

    def __call__(self, x, check=True):
        r"""
        Convert x into an element of this Hecke algebra. Here x is either:

        - an element of a Hecke algebra equal to this one

        - an element of the corresponding anemic Hecke algebra, if x is a full
          Hecke algebra

        - an element of the corresponding full Hecke algebra of the
          form `T_i` where i is coprime to ``self.level()``, if self
          is an anemic Hecke algebra

        - something that can be converted into an element of the
          underlying matrix space.

        In the last case, the parameter ``check'' controls whether or
        not to check that this element really does lie in the
        appropriate algebra. At present, setting ``check=True'' raises
        a NotImplementedError unless x is a scalar (or a diagonal
        matrix).

        EXAMPLES::

            sage: T = ModularSymbols(11).hecke_algebra()
            sage: T.gen(2) in T
            True
            sage: 5 in T
            True
            sage: T.gen(2).matrix() in T
            Traceback (most recent call last):
            ...
            NotImplementedError: Membership testing for '...' not implemented
            sage: T(T.gen(2).matrix(), check=False)
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [ 3  0 -1]
            [ 0 -2  0]
            [ 0  0 -2]
            sage: A = ModularSymbols(11).anemic_hecke_algebra()
            sage: A(T.gen(3))
            Hecke operator T_3 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: A(T.gen(11))
            Traceback (most recent call last):
            ...
            TypeError: Don't know how to construct an element of Anemic Hecke algebra acting on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field from Hecke operator T_11 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field

        TESTS:

        We test that coercion is OK between the Hecke algebras associated to two submodules which are equal but have different bases::

            sage: M = CuspForms(Gamma0(57))
            sage: f1,f2,f3 = M.newforms()
            sage: N1 = M.submodule(M.free_module().submodule_with_basis([f1.element().element(), f2.element().element()]))
            sage: N2 = M.submodule(M.free_module().submodule_with_basis([f1.element().element(), (f1.element() + f2.element()).element()]))
            sage: N1.hecke_operator(5).matrix_form()
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [ 0  1]
            sage: N2.hecke_operator(5).matrix_form()
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [-4  1]
            sage: N1.hecke_algebra()(N2.hecke_operator(5)).matrix_form()
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [ 0  1]
            sage: N1.hecke_algebra()(N2.hecke_operator(5).matrix_form())
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [ 0  1]
        """
        from . import hecke_operator
        try:
            if not isinstance(x, Element):
                x = self.base_ring()(x)
            if x.parent() is self:
                return x
            elif hecke_operator.is_HeckeOperator(x):
                if x.parent() == self \
                        or (self.is_anemic() == False and x.parent() == self.anemic_subalgebra()) \
                        or (self.is_anemic() == True and x.parent().anemic_subalgebra() == self and arith.gcd(x.index(), self.level()) == 1):
                    return hecke_operator.HeckeOperator(self, x.index())
                else:
                    raise TypeError
            elif hecke_operator.is_HeckeAlgebraElement(x):
                if x.parent() == self or (self.is_anemic() == False and x.parent() == self.anemic_subalgebra()):
                    if x.parent().module().basis_matrix() == self.module().basis_matrix():
                        return hecke_operator.HeckeAlgebraElement_matrix(self, x.matrix())
                    else:
                        A = matrix([self.module().coordinate_vector(x.parent().module().gen(i)) \
                            for i in xrange(x.parent().module().rank())])
                        return hecke_operator.HeckeAlgebraElement_matrix(self, ~A * x.matrix() * A)
                elif x.parent() == self.anemic_subalgebra():
                    pass

                else:
                    raise TypeError
            else:
                A = self.matrix_space()(x)
                if check:
                    if not A.is_scalar():
                        raise NotImplementedError("Membership testing for '%s' not implemented" % self)
                return hecke_operator.HeckeAlgebraElement_matrix(self, A)

        except TypeError:
            raise TypeError("Don't know how to construct an element of %s from %s" % (self, x))

    def _coerce_impl(self, x):
        r"""
        Implicit coercion of x into this Hecke algebra. The only things that
        coerce implicitly into self are: elements of Hecke algebras which are
        equal to self, or to the anemic subalgebra of self if self is not
        anemic; and elements that coerce into the base ring of self.  Bare
        matrices do *not* coerce implicitly into self.

        EXAMPLE::

            sage: C = CuspForms(3, 12)
            sage: A = C.anemic_hecke_algebra()
            sage: F = C.hecke_algebra()
            sage: F.coerce(A.2) # indirect doctest
            Hecke operator T_2 on Cuspidal subspace of dimension 3 of Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(3) of weight 12 over Rational Field
        """
        if x.parent() == self or (self.is_anemic() == False and x.parent() == self.anemic_subalgebra()):
            return self(x)
        else:
            return self(self.matrix_space()(1) * self.base_ring().coerce(x))
        #return self._coerce_try(x, self.matrix_space())

    def gen(self, n):
        """
        Return the `n`-th Hecke operator.

        EXAMPLES::

            sage: T = ModularSymbols(11).hecke_algebra()
            sage: T.gen(2)
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
        """
        return self.hecke_operator(n)

    def ngens(self):
        r"""
        The size of the set of generators returned by gens(), which is clearly
        infinity. (This is not necessarily a minimal set of generators.)

        EXAMPLES::

            sage: CuspForms(1, 12).anemic_hecke_algebra().ngens()
            +Infinity
        """
        return sage.rings.infinity.infinity

    def is_noetherian(self):
        """
        Return True if this Hecke algebra is Noetherian as a ring. This is true
        if and only if the base ring is Noetherian.

        EXAMPLES::

            sage: CuspForms(1, 12).anemic_hecke_algebra().is_noetherian()
            True
        """
        return self.base_ring().is_noetherian()

    def matrix_space(self):
        r"""
        Return the underlying matrix space of this module.

        EXAMPLES::

            sage: CuspForms(3, 24, base_ring=Qp(5)).anemic_hecke_algebra().matrix_space()
            Full MatrixSpace of 7 by 7 dense matrices over 5-adic Field with capped relative precision 20
        """
        try:
            return self.__matrix_space_cache
        except AttributeError:
            self.__matrix_space_cache = sage.matrix.matrix_space.MatrixSpace(
                self.base_ring(), self.module().rank())
            return self.__matrix_space_cache

    def _latex_(self):
        r"""
        LaTeX representation of self.

        EXAMPLES::

            sage: latex(CuspForms(3, 24).hecke_algebra()) # indirect doctest
            \mathbf{T}_{\text{\texttt{Cuspidal...Gamma0(3)...24...}
        """
        from sage.misc.latex import latex
        return "\\mathbf{T}_{%s}" % latex(self.__M)

    def level(self):
        r"""
        Return the level of this Hecke algebra, which is (by definition) the
        level of the Hecke module on which it acts.

        EXAMPLE::

            sage: ModularSymbols(37).hecke_algebra().level()
            37
        """
        return self.module().level()

    def module(self):
        """
        The Hecke module on which this algebra is acting.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.module()
            Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        return self.__M

    def rank(self):
        r"""
        The rank of this Hecke algebra as a module over its base
        ring. Not implemented at present.

        EXAMPLE::

            sage: ModularSymbols(Gamma1(3), 3).hecke_algebra().rank()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def basis(self):
        r"""
        Return a basis for this Hecke algebra as a free module over
        its base ring.

        EXAMPLE::

            sage: ModularSymbols(Gamma1(3), 3).hecke_algebra().basis()
            [Hecke operator on Modular Symbols space of dimension 2 for Gamma_1(3) of weight 3 with sign 0 and over Rational Field defined by:
            [1 0]
            [0 1],
            Hecke operator on Modular Symbols space of dimension 2 for Gamma_1(3) of weight 3 with sign 0 and over Rational Field defined by:
            [0 0]
            [0 2]]
        """
        try:
            return self.__basis_cache
        except AttributeError:
            pass
        level = self.level()
        bound = self.__M.hecke_bound()
        dim = self.__M.rank()
        if dim == 0:
            basis = []
        elif dim == 1:
            basis = [self.hecke_operator(1)]
        else:
            span = [self.hecke_operator(n) for n in range(1, bound+1) if not self.is_anemic() or gcd(n, level) == 1]
            rand_max = 5
            while True:
                # Project the full Hecke module to a random submodule to ease the HNF reduction.
                v = (ZZ**dim).random_element(x=rand_max)
                proj_span = matrix([T.matrix()*v for T in span])._clear_denom()[0]
                proj_basis = proj_span.hermite_form()
                if proj_basis[dim-1] == 0:
                    # We got unlucky, choose another projection.
                    rand_max *= 2
                    continue
                # Lift the projected basis to a basis in the Hecke algebra.
                trans = proj_span.solve_left(proj_basis)
                basis = [sum(c*T for c,T in zip(row,span) if c != 0) for row in trans[:dim]]
                break

        self.__basis_cache = tuple(basis)
        return basis

    def discriminant(self):
        r"""
        Return the discriminant of this Hecke algebra, i.e. the
        determinant of the matrix `{\rm Tr}(x_i x_j)` where `x_1,
        \dots,x_d` is a basis for self, and `{\rm Tr}(x)` signifies
        the trace (in the sense of linear algebra) of left
        multiplication by `x` on the algebra (*not* the trace of the
        operator `x` acting on the underlying Hecke module!). For
        further discussion and conjectures see Calegari + Stein,
        *Conjectures about discriminants of Hecke algebras of prime
        level*, Springer LNCS 3076.

        EXAMPLE::

            sage: BrandtModule(3, 4).hecke_algebra().discriminant()
            1
            sage: ModularSymbols(65, sign=1).cuspidal_submodule().hecke_algebra().discriminant()
            6144
            sage: ModularSymbols(1,4,sign=1).cuspidal_submodule().hecke_algebra().discriminant()
            1
            sage: H = CuspForms(1, 24).hecke_algebra()
            sage: H.discriminant()
            83041344
        """
        try:
            return self.__disc
        except AttributeError:
            pass
        basis = self.basis()
        d = len(basis)
        if d <= 1:
            self.__disc = ZZ(1)
            return self.__disc
        trace_matrix = matrix(ZZ, d)
        for i in range(d):
            for j in range(i+1):
                trace_matrix[i,j] = trace_matrix[j,i] = basis[i].matrix().trace_of_product(basis[j].matrix())
        self.__disc = trace_matrix.det()
        return self.__disc

    def gens(self):
        r"""
        Return a generator over all Hecke operator `T_n` for
        `n = 1, 2, 3, \ldots`. This is infinite.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: g = T.gens()
            sage: next(g)
            Hecke operator T_1 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
            sage: next(g)
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        n = 1
        while True:
            yield self.hecke_operator(n)
            n += 1

    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator `T_n`.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.hecke_operator(2)
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        try:
            return self.__hecke_operator[n]
        except AttributeError:
            self.__hecke_operator = {}
        except KeyError:
            pass
        n = int(n)
        T = self.__M._hecke_operator_class()(self, n)
        self.__hecke_operator[n] = T
        return T

    def hecke_matrix(self, n, *args, **kwds):
        """
        Return the matrix of the n-th Hecke operator `T_n`.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.hecke_matrix(2)
            [ -24    0    0]
            [   0  -24    0]
            [4860    0 2049]
        """
        return self.__M.hecke_matrix(n, *args, **kwds)

    def diamond_bracket_matrix(self, d):
        r"""
        Return the matrix of the diamond bracket operator `\langle d \rangle`.

        EXAMPLE::

            sage: T = ModularSymbols(Gamma1(7), 4).hecke_algebra()
            sage: T.diamond_bracket_matrix(3)
            [    0     0     1     0     0     0     0     0     0     0     0     0]
            [    1     0     0     0     0     0     0     0     0     0     0     0]
            [    0     1     0     0     0     0     0     0     0     0     0     0]
            [    0     0     0 -11/9  -4/9     1   2/3   7/9   2/9   7/9  -5/9  -2/9]
            [    0     0     0  58/9  17/9    -5 -10/3   4/9   5/9 -50/9  37/9  13/9]
            [    0     0     0 -22/9  -8/9     2   4/3   5/9   4/9  14/9 -10/9  -4/9]
            [    0     0     0  44/9  16/9    -4  -8/3   8/9   1/9 -28/9  20/9   8/9]
            [    0     0     0     0     0     0     0     0     0     0     1     0]
            [    0     0     0     0     0     0     0     0     0     0     0     1]
            [    0     0     0     1     0     0     0     0     0     0     0     0]
            [    0     0     0     2     0    -1     0     0     0     0     0     0]
            [    0     0     0    -4     0     4     1     0     0     0     0     0]

        """
        return self.__M.diamond_bracket_matrix(d)

    def diamond_bracket_operator(self, d):
        r"""
        Return the diamond bracket operator `\langle d \rangle`.

        EXAMPLE::

            sage: T = ModularSymbols(Gamma1(7), 4).hecke_algebra()
            sage: T.diamond_bracket_operator(3)
            Diamond bracket operator <3> on Modular Symbols space of dimension 12 for Gamma_1(7) of weight 4 with sign 0 and over Rational Field
        """
        d = int(d) % self.__M.level()
        try:
            return self.__diamond_operator[d]
        except AttributeError:
            self.__diamond_operator = {}
        except KeyError:
            pass
        D = self.__M._diamond_operator_class()(self, d)
        self.__diamond_operator[d] = D
        return D

class HeckeAlgebra_full(HeckeAlgebra_base):
    r"""
    A full Hecke algebra (including the operators `T_n` where `n` is not
    assumed to be coprime to the level).
    """
    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: ModularForms(37).hecke_algebra()._repr_()
            'Full Hecke algebra acting on Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(37) of weight 2 over Rational Field'
        """
        return "Full Hecke algebra acting on %s"%self.module()

    def __cmp__(self, other):
        r"""
        Compare self to other.

        EXAMPLES::

            sage: A = ModularForms(37).hecke_algebra()
            sage: A == QQ
            False
            sage: A == ModularForms(37).anemic_hecke_algebra()
            False
            sage: A == A
            True
        """
        if not isinstance(other, HeckeAlgebra_full):
            return -1
        return cmp(self.module(), other.module())

    def is_anemic(self):
        """
        Return False, since this the full Hecke algebra.

        EXAMPLES::

            sage: H = CuspForms(3, 12).hecke_algebra()
            sage: H.is_anemic()
            False
        """
        return False

    def anemic_subalgebra(self):
        r"""
        The subalgebra of self generated by the Hecke operators of
        index coprime to the level.

        EXAMPLE::

            sage: H = CuspForms(3, 12).hecke_algebra()
            sage: H.anemic_subalgebra()
            Anemic Hecke algebra acting on Cuspidal subspace of dimension 3 of Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(3) of weight 12 over Rational Field
        """
        return self.module().anemic_hecke_algebra()


class HeckeAlgebra_anemic(HeckeAlgebra_base):
    r"""
    An anemic Hecke algebra, generated by Hecke operators with index coprime to the level.
    """
    def _repr_(self):
        r"""
        EXAMPLE::

            sage: H = CuspForms(3, 12).anemic_hecke_algebra()._repr_()
        """
        return "Anemic Hecke algebra acting on %s"%self.module()

    def __cmp__(self, other):
        r"""
        Compare self to other.

        EXAMPLES::

            sage: A = ModularForms(23).anemic_hecke_algebra()
            sage: A == QQ
            False
            sage: A == ModularForms(23).hecke_algebra()
            False
            sage: A == A
            True

        """
        if not isinstance(other, HeckeAlgebra_anemic):
            return -1
        return cmp(self.module(), other.module())

    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator, for `n` any
        positive integer coprime to the level.

        EXAMPLES::

            sage: T = ModularSymbols(Gamma1(5),3).anemic_hecke_algebra()
            sage: T.hecke_operator(2)
            Hecke operator T_2 on Modular Symbols space of dimension 4 for Gamma_1(5) of weight 3 with sign 0 and over Rational Field
            sage: T.hecke_operator(5)
            Traceback (most recent call last):
            ...
            IndexError: Hecke operator T_5 not defined in the anemic Hecke algebra
        """
        n = int(n)
        if arith.gcd(self.module().level(), n) != 1:
            raise IndexError("Hecke operator T_%s not defined in the anemic Hecke algebra"%n)
        return self.module()._hecke_operator_class()(self, n)

    def is_anemic(self):
        """
        Return True, since this the anemic Hecke algebra.

        EXAMPLES::

            sage: H = CuspForms(3, 12).anemic_hecke_algebra()
            sage: H.is_anemic()
            True
        """
        return True

    def gens(self):
        """
        Return a generator over all Hecke operator `T_n` for
        `n = 1, 2, 3, \ldots`, with `n` coprime to the
        level. This is an infinite sequence.

        EXAMPLES::

            sage: T = ModularSymbols(12,2).anemic_hecke_algebra()
            sage: g = T.gens()
            sage: next(g)
            Hecke operator T_1 on Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Rational Field
            sage: next(g)
            Hecke operator T_5 on Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Rational Field
        """
        level = self.level()
        n = 1
        while True:
            if arith.gcd(n, level) == 1:
                yield self.hecke_operator(n)
            n += 1




