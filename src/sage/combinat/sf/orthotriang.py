# sage.doctest: needs sage.combinat sage.modules
"""
Symmetric functions defined by orthogonality and triangularity

One characterization of Schur functions is that they are upper
triangularly related to the monomial symmetric functions and
orthogonal with respect to the Hall scalar product. We can use the
class SymmetricFunctionAlgebra_orthotriang to obtain the Schur
functions from this definition.

::

    sage: from sage.combinat.sf.sfa import zee
    sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
    sage: Sym = SymmetricFunctions(QQ)
    sage: m = Sym.m()
    sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
    sage: s([2,1])^2
    s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

::

    sage: s2 = SymmetricFunctions(QQ).s()
    sage: s2([2,1])^2
    s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]
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

from . import sfa
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom


class SymmetricFunctionAlgebra_orthotriang(sfa.SymmetricFunctionAlgebra_generic):

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        pass

    @staticmethod
    def __classcall__(cls, Sym, base, scalar, prefix, basis_name, leading_coeff=None):
        """
        Normalize the arguments.

        TESTS::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: B1 = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur')
            sage: B2 = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur', None)
            sage: B1 is B2
            True
        """
        return super().__classcall__(cls, Sym, base, scalar, prefix, basis_name, leading_coeff)

    def __init__(self, Sym, base, scalar, prefix, basis_name, leading_coeff):
        r"""
        Initialization of the symmetric function algebra defined via orthotriangular rules.

        INPUT:

        - ``self`` -- a basis determined by an orthotriangular definition
        - ``Sym`` -- ring of symmetric functions
        - ``base`` -- an instance of a basis of the ring of symmetric functions
          (e.g. the Schur functions)
        - ``scalar`` -- a function ``zee`` on partitions. The function
          ``zee`` determines the scalar product on the power sum basis
          with normalization `\langle p_{\mu}, p_{\mu} \rangle =
          \mathrm{zee}(\mu)`.
        - ``prefix`` -- the prefix used to display the basis
        - ``basis_name`` -- the name used for the basis

        .. NOTE::

            The base ring is required to be a `\QQ`-algebra for this
            method to be usable, since the scalar product is defined by
            its values on the power sum basis.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur'); s
            Symmetric Functions over Rational Field in the Schur basis

        TESTS::

            sage: TestSuite(s).run(elements=[s[1,1]+2*s[2], s[1]+3*s[1,1]])
            sage: TestSuite(s).run(skip=["_test_associativity", "_test_prod", '_test_construction'])  # long time (7s on sage.math, 2011)

        Note: ``s.an_element()`` is of degree 4; so we skip
        ``_test_associativity`` and ``_test_prod`` which involve
        (currently?) expensive calculations up to degree 12.
        """
        self._sym = Sym
        self._sf_base = base
        self._scalar = scalar
        self._leading_coeff = leading_coeff
        sfa.SymmetricFunctionAlgebra_generic.__init__(self, Sym, prefix=prefix, basis_name=basis_name)

        self._self_to_base_cache = {}
        self._base_to_self_cache = {}
        self.register_coercion(SetMorphism(Hom(base, self), self._base_to_self))
        base.register_coercion(SetMorphism(Hom(self, base), self._self_to_base))

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a
        :class:`SymmetricFunctionsFunctor` and `R` is a ring, such
        that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur')
            sage: s.construction()
            (SymmetricFunctionsFunctor[Schur], Rational Field)
        """
        return OrthotriangBasisFunctor(self), self.base_ring()

    def _base_to_self(self, x):
        """
        Coerce a symmetric function in base ``x`` into ``self``.

        INPUT:

        - ``self`` -- a basis determined by an orthotriangular definition
        - ``x`` -- an element of the basis `base` to which ``self`` is triangularly related

        OUTPUT: an element of ``self`` equivalent to ``x``

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
            sage: s._base_to_self(m([2,1]))
            -2*s[1, 1, 1] + s[2, 1]
        """
        return self._from_cache(x, self._base_cache, self._base_to_self_cache)

    def _self_to_base(self, x):
        """
        Coerce a symmetric function in ``self`` into base ``x``.

        INPUT:

        - ``self`` -- a basis determined by an orthotriangular definition
        - ``x`` -- an element of ``self`` as a basis of the ring of symmetric functions

        OUTPUT: the element ``x`` expressed in the basis to which ``self`` is triangularly related

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
            sage: s._self_to_base(s([2,1]))
            2*m[1, 1, 1] + m[2, 1]
        """
        return self._sf_base._from_cache(x, self._base_cache, self._self_to_base_cache)

    def _base_cache(self, n):
        """
        Compute the change of basis between ``self`` and base for the
        homogeneous component of size ``n``

        INPUT:

        - ``self`` -- a basis determined by an orthotriangular definition
        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
            sage: s._base_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(s._base_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -1), ([2], 1)])]
            sage: l(s._self_to_base_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], 1), ([2], 1)])]
        """
        if n in self._self_to_base_cache:
            return
        else:
            self._self_to_base_cache[n] = {}

        self._gram_schmidt(n, self._sf_base, self._scalar,
                           self._self_to_base_cache,
                           leading_coeff=self._leading_coeff,
                           upper_triangular=True)
        self._invert_morphism(n, self.base_ring(), self._self_to_base_cache,
                              self._base_to_self_cache,
                              to_other_function=self._to_base)

    def _to_base(self, part):
        r"""
        Return a function which takes in a partition `\mu` and returns the
        coefficient of a partition in the expansion of ``self`` `(part)` in base.

        INPUT:

        - ``self`` -- a basis determined by an orthotriangular definition
        - ``part`` -- a partition

        .. NOTE::

           We assume that self._gram_schmidt has been called before
           self._to_base is called.

        OUTPUT:

        - a function which accepts a partition ``mu`` and returns the coefficients
          in the expansion of ``self(part)`` in the triangularly related basis.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
            sage: m(s([2,1]))
            2*m[1, 1, 1] + m[2, 1]
            sage: f = s._to_base(Partition([2,1]))
            sage: [f(p) for p in Partitions(3)]
            [0, 1, 2]
        """
        f = lambda mu: self._self_to_base_cache[part].get(mu, 0)
        return f

    def product(self, left, right):
        """
        Return ``left`` * ``right`` by converting both to the base and then
        converting back to ``self``.

        INPUT:

        - ``self`` -- a basis determined by an orthotriangular definition
        - ``left``, ``right`` -- elements in ``self``

        OUTPUT: the expansion of the product of ``left`` and ``right`` in the basis ``self``

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
            sage: s([1])*s([2,1]) #indirect doctest
            s[2, 1, 1] + s[2, 2] + s[3, 1]
        """
        return self(self._sf_base(left) * self._sf_base(right))


from sage.combinat.sf.sfa import SymmetricFunctionsFunctor


class OrthotriangBasisFunctor(SymmetricFunctionsFunctor):
    """
    A constructor for algebras of symmetric functions constructed by
    orthogonality and triangularity.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import zee
        sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
        sage: Sym = SymmetricFunctions(QQ)
        sage: m = Sym.m()
        sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur')
        sage: s.construction()
        (SymmetricFunctionsFunctor[Schur], Rational Field)
    """
    def __init__(self, basis):
        r"""
        Initialize the functor.

        INPUT:

        - ``basis`` -- the basis of the symmetric function algebra

        TESTS::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur')
            sage: F, R = s.construction()
            sage: TestSuite(F).run()
        """
        self._basis_name = basis.basis_name()
        self._sf_base = basis._sf_base
        self._scalar = basis._scalar
        self._leading_coeff = basis._leading_coeff
        self._prefix = basis._prefix
        super().__init__(basis, self._basis_name)

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import zee
            sage: from sage.combinat.sf.orthotriang import SymmetricFunctionAlgebra_orthotriang
            sage: Sym = SymmetricFunctions(ZZ)
            sage: m = Sym.m()
            sage: s = SymmetricFunctionAlgebra_orthotriang(Sym, m, zee, 's', 'Schur functions')
            sage: F, R = s.construction()  # indirect doctest
            sage: F(QQbar)
            Symmetric Functions over Algebraic Field in the Schur functions basis
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        return self._basis(SymmetricFunctions(R), self._sf_base.change_ring(R),
                           self._scalar, self._prefix, self._basis_name,
                           self._leading_coeff)


# Backward compatibility for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.sf.orthotriang', 'SymmetricFunctionAlgebraElement_orthotriang', SymmetricFunctionAlgebra_orthotriang.Element)
