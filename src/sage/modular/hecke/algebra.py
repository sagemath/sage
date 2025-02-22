# sage.doctest: needs sage.libs.flint sage.libs.pari
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
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from collections.abc import Iterator

from sage.rings.infinity import infinity
from sage.categories.algebras import Algebras
from sage.matrix.constructor import matrix
from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.misc.latex import latex
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.element import Element
from sage.structure.unique_representation import CachedRepresentation
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_method, richcmp


def is_HeckeAlgebra(x) -> bool:
    r"""
    Return ``True`` if x is of type HeckeAlgebra.

    EXAMPLES::

        sage: from sage.modular.hecke.algebra import is_HeckeAlgebra
        sage: is_HeckeAlgebra(CuspForms(1, 12).anemic_hecke_algebra())
        doctest:warning...
        DeprecationWarning: the function is_HeckeAlgebra is deprecated;
        use 'isinstance(..., HeckeAlgebra_base)' instead
        See https://github.com/sagemath/sage/issues/37895 for details.
        True
        sage: is_HeckeAlgebra(ZZ)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37895, "the function is_HeckeAlgebra is deprecated; use 'isinstance(..., HeckeAlgebra_base)' instead")
    return isinstance(x, HeckeAlgebra_base)


def _heckebasis(M):
    r"""
    Return a basis of the Hecke algebra of M as a ZZ-module.

    INPUT:

    - ``M`` -- a Hecke module

    OUTPUT: list of Hecke algebra elements represented as matrices

    EXAMPLES::

        sage: M = ModularSymbols(11,2,1)
        sage: sage.modular.hecke.algebra._heckebasis(M)
        [Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
        [1 0]
        [0 1],
        Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
        [0 2]
        [0 5]]
    """
    d = M.rank()
    WW = ZZ**(d**2)
    MM = MatrixSpace(QQ, d)
    S = []
    Denom = []
    B1 = []
    for i in range(1, M.hecke_bound() + 1):
        v = M.hecke_operator(i).matrix()
        den = v.denominator()
        Denom.append(den)
        S.append(v)
    den = lcm(Denom)
    B = [WW((den * m).list()) for m in S]
    UU = WW.submodule(B)
    for u in UU.basis():
        u1 = u.list()
        m1 = M.hecke_algebra()(MM(u1), check=False)
        B1.append((1 / den) * m1)
    return B1


@richcmp_method
class HeckeAlgebra_base(CachedRepresentation, Parent):
    """
    Base class for algebras of Hecke operators on a fixed Hecke module.

    INPUT:

    - ``M`` -- a Hecke module

    EXAMPLES::

        sage: CuspForms(1, 12).hecke_algebra() # indirect doctest
        Full Hecke algebra acting on Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
    """
    @staticmethod
    def __classcall__(cls, M):
        r"""
        We need to work around a problem originally discovered by David
        Loeffler on 2009-04-13: The problem is that if one creates two
        subspaces of a Hecke module which are equal as subspaces but have
        different bases, then the caching machinery needs to distinguish
        between them. So we need to add ``basis_matrix`` to the cache key even
        though it is not looked at by the constructor.

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
        if isinstance(M, tuple):
            M = M[0]
        try:
            M = (M, M.basis_matrix())
        except AttributeError:
            # The AttributeError occurs if M is not a free module; then it might not have a basis_matrix method
            pass
        return super().__classcall__(cls, M)

    def __init__(self, M) -> None:
        """
        Initialization.

        EXAMPLES::

            sage: from sage.modular.hecke.algebra import HeckeAlgebra_base
            sage: type(HeckeAlgebra_base(CuspForms(1, 12)))
            <class 'sage.modular.hecke.algebra.HeckeAlgebra_base_with_category'>
        """
        if isinstance(M, tuple):
            M = M[0]
        from .module import HeckeModule_generic
        if not isinstance(M, HeckeModule_generic):
            msg = f"M (={M}) must be a HeckeModule"
            raise TypeError(msg)
        self.__M = M
        cat = Algebras(M.base_ring()).Commutative()
        Parent.__init__(self, base=M.base_ring(), category=cat)

    def _an_element_(self):
        r"""
        Return an element of this algebra.

        Used by the coercion machinery.

        EXAMPLES::

            sage: CuspForms(1, 12).hecke_algebra().an_element() # indirect doctest
            Hecke operator T_2 on Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
        """
        return self.hecke_operator(self.level() + 1)

    def _element_constructor_(self, x, check=True):
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

        In the last case, the parameter ``check`` controls whether or
        not to check that this element really does lie in the
        appropriate algebra. At present, setting ``check=True`` raises
        a :exc:`NotImplementedError` unless x is a scalar (or a diagonal
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
            NotImplementedError: membership testing for '...' not implemented
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
        """
        from .hecke_operator import HeckeAlgebraElement_matrix, HeckeOperator, HeckeAlgebraElement

        if not isinstance(x, Element):
            x = self.base_ring()(x)

        parent = x.parent()

        if parent is self.base_ring():
            return HeckeAlgebraElement_matrix(self, x * self.matrix_space().one())

        if parent is self:
            return x

        if isinstance(x, HeckeOperator):
            if x.parent() == self \
                    or (not self.is_anemic() and x.parent() == self.anemic_subalgebra()) \
                    or (self.is_anemic() and x.parent().anemic_subalgebra() == self and gcd(x.index(), self.level()) == 1):
                return HeckeOperator(self, x.index())

        if isinstance(x, HeckeAlgebraElement):
            if x.parent() == self or (not self.is_anemic() and x.parent() == self.anemic_subalgebra()):
                if x.parent().module().basis_matrix() == self.module().basis_matrix():
                    return HeckeAlgebraElement_matrix(self, x.matrix())
                else:
                    A = matrix([self.module().coordinate_vector(x.parent().module().gen(i))
                                for i in range(x.parent().module().rank())])
                    return HeckeAlgebraElement_matrix(self, ~A * x.matrix() * A)

        try:
            A = self.matrix_space()(x)
            if check:
                if not A.is_scalar():
                    msg = f"membership testing for '{self}' not implemented"
                    raise NotImplementedError(msg)
            return HeckeAlgebraElement_matrix(self, A)
        except TypeError:
            raise TypeError("Don't know how to construct an element of %s from %s" % (self, x))

    def _coerce_map_from_(self, R):
        """
        Coercion of a parent ``R`` into this Hecke algebra.

        The parents that coerce into ``self`` are: Hecke
        algebras which are equal to ``self``, or to the anemic subalgebra
        of ``self`` if ``self`` is not anemic; and parents that coerce into
        the base ring of ``self``.

        EXAMPLES::

            sage: C = CuspForms(3, 12)
            sage: A = C.anemic_hecke_algebra()
            sage: F = C.hecke_algebra()
            sage: F.coerce(A.2) # indirect doctest
            Hecke operator T_2 on Cuspidal subspace of dimension 3 of Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(3) of weight 12 over Rational Field
        """
        if R == self:
            return True
        if not self.is_anemic() and R == self.anemic_subalgebra():
            return True
        return self.base_ring().has_coerce_map_from(R)

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
        infinity.

        (This is not necessarily a minimal set of generators.)

        EXAMPLES::

            sage: CuspForms(1, 12).anemic_hecke_algebra().ngens()
            +Infinity
        """
        return infinity

    def one(self):
        """
        Return the unit of the Hecke algebra.

        EXAMPLES::

            sage: M = ModularSymbols(11,2,1)
            sage: H = M.hecke_algebra()
            sage: H.one()
            Hecke operator on Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field defined by:
            [1 0]
            [0 1]
        """
        from .hecke_operator import HeckeAlgebraElement_matrix
        A = self.matrix_space()
        return HeckeAlgebraElement_matrix(self, A.one())

    def is_noetherian(self) -> bool:
        """
        Return ``True`` if this Hecke algebra is Noetherian as a ring.

        This is true if and only if the base ring is Noetherian.

        EXAMPLES::

            sage: CuspForms(1, 12).anemic_hecke_algebra().is_noetherian()
            True
        """
        return self.base_ring().is_noetherian()

    @cached_method
    def matrix_space(self):
        r"""
        Return the underlying matrix space of this module.

        EXAMPLES::

            sage: CuspForms(3, 24, base_ring=Qp(5)).anemic_hecke_algebra().matrix_space()
            Full MatrixSpace of 7 by 7 dense matrices over 5-adic Field with capped relative precision 20
        """
        return MatrixSpace(self.base_ring(), self.module().rank())

    def _latex_(self) -> str:
        r"""
        LaTeX representation of ``self``.

        EXAMPLES::

            sage: latex(CuspForms(3, 24).hecke_algebra()) # indirect doctest
            \mathbf{T}_{\text{\texttt{Cuspidal...Gamma0(3)...24...}
        """
        return "\\mathbf{T}_{%s}" % latex(self.__M)

    def level(self):
        r"""
        Return the level of this Hecke algebra, which is (by definition) the
        level of the Hecke module on which it acts.

        EXAMPLES::

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

        EXAMPLES::

            sage: ModularSymbols(Gamma1(3), 3).hecke_algebra().rank()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    @cached_method
    def basis(self):
        r"""
        Return a basis for this Hecke algebra as a free module over
        its base ring.

        EXAMPLES::

            sage: ModularSymbols(Gamma1(3), 3).hecke_algebra().basis()
            (Hecke operator on Modular Symbols space of dimension 2 for Gamma_1(3) of weight 3 with sign 0 over Rational Field defined by:
            [1 0]
            [0 1],
            Hecke operator on Modular Symbols space of dimension 2 for Gamma_1(3) of weight 3 with sign 0 over Rational Field defined by:
            [0 0]
            [0 2])

            sage: M = ModularSymbols(Gamma0(22), sign=1)
            sage: H = M.hecke_algebra()
            sage: B = H.basis()
            sage: len(B)
            5
            sage: all(b in H for b in B)
            True
            sage: [B[0, 0] for B in M.anemic_hecke_algebra().basis()]
            Traceback (most recent call last):
            ...
            NotImplementedError: basis not implemented for anemic Hecke algebra
        """
        bound = self.__M.hecke_bound()
        dim = self.__M.rank()

        # current implementation gets stuck in an infinite loop when dimension of Hecke alg != dimension of space
        if self.is_anemic():
            raise NotImplementedError("basis not implemented for anemic Hecke algebra")

        if dim == 0:
            basis = []
        elif dim == 1:
            basis = [self.hecke_operator(1)]
        else:
            span = [self.hecke_operator(n) for n in range(1, bound + 1)]
            rand_max = 5
            while True:
                # Project the full Hecke module to a random submodule to ease the HNF reduction.
                v = (ZZ**dim).random_element(x=rand_max)
                proj_span = matrix([T.matrix() * v for T in span])._clear_denom()[0]
                proj_basis = proj_span.hermite_form()
                if proj_basis[dim - 1] == 0:
                    # We got unlucky, choose another projection.
                    rand_max *= 2
                    continue
                # Lift the projected basis to a basis in the Hecke algebra.
                trans = proj_span.solve_left(proj_basis)
                basis = [sum(c * T for c, T in zip(row, span) if c != 0)
                         for row in trans[:dim]]
                break

        return tuple(basis)

    @cached_method
    def discriminant(self):
        r"""
        Return the discriminant of this Hecke algebra.

        This is the
        determinant of the matrix `{\rm Tr}(x_i x_j)` where `x_1,
        \dots,x_d` is a basis for self, and `{\rm Tr}(x)` signifies
        the trace (in the sense of linear algebra) of left
        multiplication by `x` on the algebra (*not* the trace of the
        operator `x` acting on the underlying Hecke module!). For
        further discussion and conjectures see Calegari + Stein,
        *Conjectures about discriminants of Hecke algebras of prime
        level*, Springer LNCS 3076.

        EXAMPLES::

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
        basis = self.basis()
        d = len(basis)
        if d <= 1:
            return ZZ.one()
        trace_matrix = matrix(ZZ, d)
        for i in range(d):
            for j in range(i + 1):
                trace_matrix[i, j] = trace_matrix[j, i] = basis[i].matrix().trace_of_product(basis[j].matrix())
        return trace_matrix.det()

    def gens(self) -> Iterator:
        r"""
        Return a generator over all Hecke operator `T_n`
        for `n = 1, 2, 3, \ldots`.

        This is infinite.

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

    @cached_method(key=lambda self, n: int(n))
    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator `T_n`.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.hecke_operator(2)
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        return self.__M._hecke_operator_class()(self, n)

    def hecke_matrix(self, n, *args, **kwds):
        """
        Return the matrix of the `n`-th Hecke operator `T_n`.

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

        EXAMPLES::

            sage: T = ModularSymbols(Gamma1(7), 4).hecke_algebra()
            sage: d3 = T.diamond_bracket_matrix(3)
            sage: x = d3.charpoly().variables()[0]
            sage: d3.charpoly() == (x^3-1)^4
            True
        """
        return self.__M.diamond_bracket_matrix(d)

    @cached_method(key=lambda self, d: int(d) % self.__M.level())
    def diamond_bracket_operator(self, d):
        r"""
        Return the diamond bracket operator `\langle d \rangle`.

        EXAMPLES::

            sage: T = ModularSymbols(Gamma1(7), 4).hecke_algebra()
            sage: T.diamond_bracket_operator(3)
            Diamond bracket operator <3> on Modular Symbols space of dimension 12 for Gamma_1(7) of weight 4 with sign 0 over Rational Field
        """
        return self.__M._diamond_operator_class()(self, d)


class HeckeAlgebra_full(HeckeAlgebra_base):
    r"""
    A full Hecke algebra (including the operators `T_n` where `n` is not
    assumed to be coprime to the level).
    """
    def _repr_(self) -> str:
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: ModularForms(37).hecke_algebra()._repr_()
            'Full Hecke algebra acting on Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(37) of weight 2 over Rational Field'
        """
        return "Full Hecke algebra acting on %s" % self.module()

    def __richcmp__(self, other, op) -> bool:
        r"""
        Compare ``self`` to ``other``.

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
            return NotImplemented
        return richcmp(self.module(), other.module(), op)

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
        The subalgebra of ``self`` generated by the Hecke operators of
        index coprime to the level.

        EXAMPLES::

            sage: H = CuspForms(3, 12).hecke_algebra()
            sage: H.anemic_subalgebra()
            Anemic Hecke algebra acting on Cuspidal subspace of dimension 3 of Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(3) of weight 12 over Rational Field
        """
        return self.module().anemic_hecke_algebra()


HeckeAlgebra = HeckeAlgebra_full


class HeckeAlgebra_anemic(HeckeAlgebra_base):
    r"""
    An anemic Hecke algebra, generated by Hecke operators with index coprime to the level.
    """
    def _repr_(self) -> str:
        r"""
        EXAMPLES::

            sage: H = CuspForms(3, 12).anemic_hecke_algebra()._repr_()
        """
        return "Anemic Hecke algebra acting on %s" % self.module()

    def __richcmp__(self, other, op) -> bool:
        r"""
        Compare ``self`` to ``other``.

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
            return NotImplemented
        return richcmp(self.module(), other.module(), op)

    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator, for `n` any
        positive integer coprime to the level.

        EXAMPLES::

            sage: T = ModularSymbols(Gamma1(5),3).anemic_hecke_algebra()
            sage: T.hecke_operator(2)
            Hecke operator T_2 on Modular Symbols space of dimension 4 for Gamma_1(5) of weight 3 with sign 0 over Rational Field
            sage: T.hecke_operator(5)
            Traceback (most recent call last):
            ...
            IndexError: Hecke operator T_5 not defined in the anemic Hecke algebra
        """
        n = int(n)
        if gcd(self.module().level(), n) != 1:
            raise IndexError("Hecke operator T_%s not defined in the anemic Hecke algebra" % n)
        return self.module()._hecke_operator_class()(self, n)

    def is_anemic(self) -> bool:
        """
        Return ``True``, since this is the anemic Hecke algebra.

        EXAMPLES::

            sage: H = CuspForms(3, 12).anemic_hecke_algebra()
            sage: H.is_anemic()
            True
        """
        return True

    def gens(self) -> Iterator:
        r"""
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
            if gcd(n, level) == 1:
                yield self.hecke_operator(n)
            n += 1


AnemicHeckeAlgebra = HeckeAlgebra_anemic
