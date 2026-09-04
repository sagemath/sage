r"""
Bases of Supersymmetric functions

AUTHORS:

- Shriya M
"""
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.combinat.partition import Partition, _Partitions
from sage.rings.integer_ring import ZZ
from sage.categories.tensor import TensorProductsCategory
from sage.categories.tensor import tensor
from sage.combinat.sf.sf import SymmetricFunctions

class SuperSymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of bases of supersymmetric functions.

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: from sage.combinat.super_sf.super_sfa import SuperSymmetricFunctionsBases
        sage: s1 = SuperSymmetricFunctionsBases(s)
        sage: s1
        Category of bases of Supersymmetric functions over Rational Field
    """
    def _repr_(self):
        r"""
        Return the string representation of this category.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import SuperSymmetricFunctionsBases
            sage: s1 = SuperSymmetricFunctionsBases(s)
            sage: s1._repr_()
            'Category of bases of Supersymmetric functions over Rational Field'
        """
        return "Category of bases of %s" % self.base()

    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import SuperSymmetricFunctionsBases
            sage: s1 = SuperSymmetricFunctionsBases(s)
            sage: s1
            Category of bases of Supersymmetric functions over Rational Field
            sage: s1.super_categories()
            [Category of realizations of Supersymmetric functions over Rational Field,
             Category of commutative Hopf algebras with basis over Rational Field,
             Join of Category of realizations of Hopf algebras over Rational Field
                 and Category of graded algebras over Rational Field
                 and Category of graded coalgebras over Rational Field,
             Category of unique factorization domains]
        """
        R = self.base().base_ring()
        cat = HopfAlgebras(R)
        categories = [self.base().Realizations(),
                      cat.Commutative().WithBasis(),
                      cat.Graded().Realizations()]
        if R in PrincipalIdealDomains:
            categories.append(UniqueFactorizationDomains())
        return categories

    class ParentMethods:
        def one_basis(self):
            r"""
            Return the empty partition, as per ``AlgebrasWithBasis.ParentMethods.one_basis``.

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ)
                sage: h = s.h()
                sage: h.one_basis()
                []
            """
            return _Partitions([])

        def is_integral_domain(self, proof=True):
            r"""
            Return whether ``self`` is an integral domain.

            INPUT:

            - ``proof`` -- (default: ``True``) when ``False`` is returned, if ``True`` then
              this is not an integral domain;  otherwise the test might not be able
              to determine if it an integral domain or not

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ).h()
                sage: s.is_integral_domain()
                True

                sage: s = SuperSymmetricFunctions(Zmod(14)).p()
                sage: s.is_integral_domain()
                False
            """
            return self.base_ring().is_integral_domain(proof=proof)

        def fraction_field(self):
            r"""
            Return the fraction field of ``self``.

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ).p()
                sage: s.fraction_field()
                Fraction Field of Supersymmetric functions over Rational Field in the powersum basis
            """
            if not self.is_integral_domain():
                raise TypeError("self must be an integral domain")
            from sage.rings.fraction_field import FractionField_generic
            return FractionField_generic(self)

        def lift(self, x):
            r"""
            Return value of ``x`` under plethysm from
            ``SupersymFunctionAlgebra_powersum`` to
            ``SymmetricFunctionAlgebra_power # SymmetricFunctionAlgebra_power``.

            INPUT:

            - ``x`` -- element of ``self``

            OUTPUT:

            The value of ``x`` under the plethysm to the tensor product of Sym.

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ)
                sage: h = s.h()
                sage: h.lift(h[5])
                -1/120*p[] # p[1, 1, 1, 1, 1] - 1/12*p[] # p[2, 1, 1, 1] -
                1/8*p[] # p[2, 2, 1] - 1/6*p[] # p[3, 1, 1] - 1/6*p[] # p[3, 2] -
                1/4*p[] # p[4, 1] - 1/5*p[] # p[5] +
                1/24*p[1] # p[1, 1, 1, 1] + 1/4*p[1] # p[2, 1, 1] +
                1/8*p[1] # p[2, 2] + 1/3*p[1] # p[3, 1] + 1/4*p[1] # p[4] -
                1/12*p[1, 1] # p[1, 1, 1] - 1/4*p[1, 1] # p[2, 1] -
                1/6*p[1, 1] # p[3] + 1/12*p[1, 1, 1] # p[1, 1] +
                1/12*p[1, 1, 1] # p[2] - 1/24*p[1, 1, 1, 1] # p[1] +
                1/120*p[1, 1, 1, 1, 1] # p[] - 1/12*p[2] # p[1, 1, 1] -
                1/4*p[2] # p[2, 1] - 1/6*p[2] # p[3] + 1/4*p[2, 1] # p[1, 1] +
                1/4*p[2, 1] # p[2] - 1/4*p[2, 1, 1] # p[1] +
                1/12*p[2, 1, 1, 1] # p[] - 1/8*p[2, 2] # p[1] +
                1/8*p[2, 2, 1] # p[] + 1/6*p[3] # p[1, 1] + 1/6*p[3] # p[2] -
                1/3*p[3, 1] # p[1] + 1/6*p[3, 1, 1] # p[] + 1/6*p[3, 2] # p[] -
                1/4*p[4] # p[1] + 1/4*p[4, 1] # p[] + 1/5*p[5] # p[]
            """
            p = self.realization_of().a_realization()
            x = p(x)
            return p.lift(x)

        def retract(self, x):
            r"""
            Return the value under retraction of the plethysm from
            ``SupersymFunctionAlgebra_powersum`` to
            ``SymmetricFunctionAlgebra_power # SymmetricFunctionAlgebra_power``.

            INPUT:

            - ``x`` -- element of tensor product of powersum symmetric functions

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ)
                sage: e = s.e()
                sage: e.retract(e.lift(e[5]))
                e[5]
            """
            p = self.realization_of().a_realization()
            return self(p.retract(x))

    class TensorProducts(TensorProductsCategory):
        class ParentMethods:
            def antipode_on_basis(self, x):
                r"""
                Return antipode of ``x``.

                EXAMPLES::

                    sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                    sage: s = SuperSymmetricFunctions(QQ)
                    sage: p = s.p()
                    sage: T = p.tensor_square()
                    sage: from sage.combinat.partition import _Partitions
                    sage: f = T.sum_of_monomials((_Partitions([4,3]), _Partitions([7,5,4])))
                    sage: T.antipode_on_basis(f)
                    -p[4, 3] # p[7, 5, 4]
                """
                TF = self.tensor_factors()
                return tensor([F.antipode_on_basis(c[0]) for F, c in zip(TF, x)])

class GradedSuperSymBases(Category_realization_of_parent):
    r"""
    The category of graded bases of supersymmetric functions.

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: from sage.combinat.super_sf.super_sfa import GradedSuperSymBases
        sage: f = GradedSuperSymBases(s)
        sage: f
        Category of graded bases of Supersymmetric functions over Rational Field
    """
    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import GradedSuperSymBases
            sage: f = GradedSuperSymBases(s)
            sage: f._repr_()
            'Category of graded bases of Supersymmetric functions over Rational Field'
        """
        return "Category of graded bases of %s" % self.base()

    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import GradedSuperSymBases
            sage: g = GradedSuperSymBases(s)
            sage: g.super_categories()
            [Category of bases of Supersymmetric functions over Rational Field,
             Category of commutative graded Hopf algebras with basis over Rational Field]
        """
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Graded()
        return [SuperSymmetricFunctionsBases(self.base()), cat]

    class ParentMethods:
        def counit(self, element):
            r"""
            Return the counit of ``element``.

            The counit is the constant term of ``element``.

            INPUT:

            - ``element`` -- element in a basis of the ring of supersymmetric functions

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: Sym = SuperSymmetricFunctions(QQ)
                sage: p = Sym.powersum()
                sage: f = 2*p[2,1] + 3*p[[]]
                sage: f.counit()
                3
            """
            return element.coefficient(_Partitions([]))

class SuperSymAlgebra_generic(CombinatorialFreeModule):
    r"""
    Abstract class for Supersymmetric Function Algebras.
    """
    def __init__(self, SuperSym=None, graded=True, prefix=None, basis_name=None):
        r"""
        Initialize a supersymmetric function algebra.

        INPUT:

        - ``SuperSym`` -- ring of supersymmetric functions (default: ``None``)
        - ``graded`` -- boolean (default: ``True``); if ``True``, then the basis is
          considered to be graded, otherwise the basis is filtered

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: from sage.combinat.super_sf.super_sfa import SuperSymAlgebra_generic
            sage: p = s.p()
            sage: isinstance(p, SuperSymAlgebra_generic)
            True
        """
        self.basis_name = basis_name
        if prefix is not None:
            self._prefix = prefix
        R = SuperSym.base_ring()
        from sage.categories.commutative_rings import CommutativeRings
        if R not in CommutativeRings():
            raise TypeError("argument R must be a commutative ring")
        if graded:
            cat = GradedSuperSymBases(SuperSym)
        else:
            cat = SuperSymmetricFunctionsBases(SuperSym)
        CombinatorialFreeModule.__init__(self, R, basis_keys=_Partitions, category=cat, bracket='', prefix=prefix)

    def _repr_(self):
            r"""
            Return a string representation of ``self``.

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ)
                sage: e = s.e()
                sage: e._repr_()
                'Supersymmetric functions over Rational Field in the elementary basis'
            """
            return "%s in the %s basis" % (self.realization_of(), self.basis_name)

    def __getitem__(self, c):
        r"""
        Return the monomial corresponding to the given list/partition/integer.

        INPUT:

        - ``c`` -- list, integer or partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p[5,4,3]
            p[5, 4, 3]
            sage: from sage.combinat.partition import _Partitions
            sage: part =_Partitions([3,2])
            sage: p[part]
            p[3, 2]
        """
        if not isinstance(c, Partition):
            if c not in ZZ:
                list(c).sort(reverse=True)
            else:
                c = [c]
            c = _Partitions(c)
        return self.monomial(c)

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` to ``self``.

        INPUT:

        - ``x`` -- an element of supersymmetric functions

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: h = s.h()
            sage: p(h[5])
            1/120*p[1, 1, 1, 1, 1] + 1/12*p[2, 1, 1, 1] + 1/8*p[2, 2, 1] +
            1/6*p[3, 1, 1] + 1/6*p[3, 2] + 1/4*p[4, 1] + 1/5*p[5]
            sage: e = s.e()
            sage: h(e[3,1])
            h[1, 1, 1, 1] - 2*h[2, 1, 1] + h[3, 1]
        """
        old_basis = x.parent().basis_name
        new_basis = self.basis_name
        R = self.base_ring()
        Sym = SymmetricFunctions(R)
        old_basis_sym = getattr(Sym, old_basis)()
        new_basis_sym = getattr(Sym, new_basis)()
        new_x = old_basis_sym._from_dict(x.monomial_coefficients())
        res = new_basis_sym(new_x)
        return self._from_dict(d=res.monomial_coefficients())

class SuperSymAlgebra_multiplicative(SuperSymAlgebra_generic):
    r"""
    Abstract class of multiplicative supersymmetric function algebras.

    A realization `h` of the ring of supersymmetric functions is multiplicative
    if for a partition `\lambda = (\lambda_1,\lambda_2,\ldots)` we have
    `h_\lambda = h_{\lambda_1} h_{\lambda_2} \cdots`.
    """
    def product_on_basis(self, left, right):
        r"""
        Return the product of ``left`` and ``right``.

        INPUT:

        - ``left``, ``right`` -- partitions

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: e = s.e()
            sage: e.product_on_basis([3,2], [6,5])
            e[6, 5, 3, 2]
        """
        m = list(left) + list(right)
        m.sort(reverse=True)
        return self.monomial(_Partitions(m))

    def coproduct_on_basis(self, mu):
        r"""
        Return the coproduct on a basis element for multiplicative bases.

        INPUT:

        - ``mu`` -- a partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: from sage.combinat.partition import Partition, _Partitions
            sage: part = _Partitions([4,4,2])
            sage: h.coproduct_on_basis(part)
            h[] # h[4, 4, 2] + 2*h[1] # h[4, 3, 2] + h[1] # h[4, 4, 1] +
            h[1, 1] # h[3, 3, 2] + 2*h[1, 1] # h[4, 3, 1] +
            h[1, 1, 1] # h[3, 3, 1] + 2*h[2] # h[4, 2, 2] + h[2] # h[4, 4] +
            2*h[2, 1] # h[3, 2, 2] + 2*h[2, 1] # h[4, 2, 1] +
            2*h[2, 1] # h[4, 3] + 2*h[2, 1, 1] # h[3, 2, 1] +
            h[2, 1, 1] # h[3, 3] + h[2, 2] # h[2, 2, 2] + 2*h[2, 2] # h[4, 2] +
            h[2, 2, 1] # h[2, 2, 1] + 2*h[2, 2, 1] # h[3, 2] +
            h[2, 2, 2] # h[2, 2] + 2*h[3] # h[4, 2, 1] +
            2*h[3, 1] # h[3, 2, 1] + 2*h[3, 1] # h[4, 1, 1] +
            2*h[3, 1, 1] # h[3, 1, 1] + 2*h[3, 2] # h[2, 2, 1] +
            2*h[3, 2] # h[4, 1] + 2*h[3, 2, 1] # h[2, 1, 1] +
            2*h[3, 2, 1] # h[3, 1] + 2*h[3, 2, 2] # h[2, 1] +
            h[3, 3] # h[2, 1, 1] + h[3, 3, 1] # h[1, 1, 1] +
            h[3, 3, 2] # h[1, 1] + 2*h[4] # h[4, 2] + 2*h[4, 1] # h[3, 2] +
            2*h[4, 1] # h[4, 1] + 2*h[4, 1, 1] # h[3, 1] +
            2*h[4, 2] # h[2, 2] + 2*h[4, 2] # h[4] + 2*h[4, 2, 1] # h[2, 1] +
            2*h[4, 2, 1] # h[3] + 2*h[4, 2, 2] # h[2] + 2*h[4, 3] # h[2, 1] +
            2*h[4, 3, 1] # h[1, 1] + 2*h[4, 3, 2] # h[1] + h[4, 4] # h[2] +
            h[4, 4, 1] # h[1] + h[4, 4, 2] # h[]
        """
        T = self.tensor_square()
        return T.prod(self.coproduct_on_generators(p) for p in mu)
