r"""
Ore modules.

Let `R` be a commutative ring, `\theta : K \to K` by a ring
endomorphism and `\partial : K \to K` be a `\theta`-derivation,
that is an additive map satisfying the following axiom

.. MATH::

    \partial(x y) = \theta(x) \partial(y) + \partial(x) y

A Ore module over `(R, \theta, \partial)` is a `R`-module `M`
equipped with a additive `f : M \to M` such that

.. MATH::

    f(a x) = \theta(a) f(x) + \partial(a) x

Such a map `f` is called a pseudomorphism.

Equivalently, a Ore module is a module over the (noncommutative)
Ore polynomial ring `\mathcal S = R[X; \theta, \partial]`.

.. RUBRIC:: Defining Ore modules

SageMath provides support for creating and manipulating Ore
modules that are finite free over the base ring `R`.

To start with, the method :meth:`quotient_module` creates the
quotient `\mathcal S/ \mathcal S P`, endowed with its structure
of `\mathcal S`-module, that is its structure of Ore module::

    sage: K.<z> = GF(5^3)
    sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
    sage: M = S.quotient_module(X^2 + z)
    sage: M
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

Classical methods are available and we can work with elements in
`M` as we do usually for vectors in finite free modules::

    sage: M.basis()
    [(1, 0), (0, 1)]

    sage: v = M((z, z^2)); v
    (z, z^2)
    sage: z*v
    (z^2, 2*z + 2)

The Ore action (or equivalently the structure of `mathcal A`-module)
is also easily accessible)::

    sage: X*v
    (3*z^2 + 2*z, 2*z^2 + 4*z + 4)

The method :meth:`pseudohom` returns the map `f` defining the action
of `X`::

    sage: M.pseudohom()
    Free module pseudomorphism (twisted by z |--> z^5) defined by the matrix
    [  0   1]
    [4*z   0]
    Domain: Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    Codomain: Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

A useful feature is the possibility to give chosen names to the vectors
of the canonical basis. This is easily done as follows::

    sage: N.<u,v,w> = S.quotient_module(X^3 + z*X + 1)
    sage: N
    Ore module <u, v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: N.basis()
    [u, v, w]

Alternatively, one can pass in the argument ``names``; this could be
useful in particular when we want to name them `e_0, e_1, \ldots`::

    sage: A = S.quotient_module(X^11 + z, names='e')
    sage: A
    Ore module <e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: A.basis()
    [e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10]

Do not forget to use the method :meth:`inject_variables` to get the
`e_i` in your namespace::

    sage: e0
    Traceback (most recent call last):
    ...
    NameError: name 'e0' is not defined
    sage: A.inject_variables()
    Defining e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10
    sage: e0
    e0

.. RUBRIC:: Submodules and quotients

SageMath provides facilities for creating submodules and quotient
modules of Ore modules.
First of all, we define the Ore module `\mathcal S/\mathcal S P^2`
(for some Ore polynomials `P`), which is obviously not simple::

    sage: P = X^2 + z*X + 1
    sage: U = S.quotient_module(P^2, names='u')
    sage: U.inject_variables()
    Defining u0, u1, u2, u3

We now build the submodule `\mathcal S Q / \mathcal S PQ` using
the method :meth:`span`::

    sage: V = U.span(P*u0)
    sage: V
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: V.basis()
    [u0 + (z^2 + 2*z + 2)*u2 + 4*z*u3,
     u1 + (2*z^2 + 4*z + 4)*u2 + u3]

We underline that the span is really the `\mathcal S`-span and
not the `R`-span (as otherwise, it will not be a Ore module).

As before, one can use the attributes ``names`` to give explicit
names to the basis vectors::

    sage: V = U.span(P*u0, names='v')
    sage: V
    Ore module <v0, v1> over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: V.inject_variables()
    Defining v0, v1
    sage: v0
    v0
    sage: U(v0)
    u0 + (z^2 + 2*z + 2)*u2 + 4*z*u3

A coercion map from `V` to `U` is automatically created and set up.
Hence, we can safely combine vectors in `V` and vectors in `U` in a
single expression::

    sage: v0 - u0
    (z^2 + 2*z + 2)*u2 + 4*z*u3

We can create the quotient `U/V` using a similar syntax::

    sage: W = U.quo(P*u0)
    sage: W
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: W.basis()
    [u2, u3]

We see that SageMath reuses by default the names of the representatives
to denote the vectors in the quotient `U/V`. This behaviour can be overrided
by providing explicit names using the attributes ``names``.

Shortcuts for creating quotients are also available::

    sage: U / (P*u0)
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
    sage: U/V
    Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

.. RUBRIC:: Morphisms of Ore modules

For a tutorial on morphisms of Ore modules, we refer to
:module:`sage.modules.ore_modules_morphism`.

AUTHOR:

- Xavier Caruso (2024-10)
"""

# ***************************************************************************
#    Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import operator
from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name
from sage.structure.sequence import Sequence
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.action import Action
from sage.categories.fields import Fields
from sage.categories.ore_modules import OreModules

from sage.matrix.matrix0 import Matrix
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.modules.free_module import FreeModule_ambient
from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.modules.ore_module_element import OreModuleElement

# Action by left multiplication on Ore modules
##############################################

class ScalarAction(Action):
    r"""
    Action by scalar multiplication on Ore modules.
    """
    def _act_(self, a, x):
        r"""
        Return the result of the action of `a` on `x`.

        INPUT:

        - ``a`` -- a scalar in the base ring

        - ``x`` -- a vector in a Ore module

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<e0,e1> = S.quotient_module(X^2 + z)  # indirect doctest
            sage: z*e0  # indirect doctest
            z*e0
        """
        return x._rmul_(a)

class OreAction(Action):
    r"""
    Action by left multiplication of Ore polynomial rings
    over Ore modules.
    """
    def _act_(self, P, x):
        r"""
        Return the result of the action of `P` on `x`.

        INPUT:

        - ``P`` -- a Ore polynomial

        - ``x`` -- a vector in a Ore module

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<e0,e1> = S.quotient_module(X^2 + z)  # indirect doctest
            sage: X*e0  # indirect doctest
            e1
        """
        M = x.parent()
        ans = P[0]*x
        y = x
        for i in range(1, P.degree() + 1):
            y = y.image()
            ans += y._rmul_(P[i])
        return ans

# Generic class for Ore modules
###############################

def normalize_names(names, rank):
    r"""
    Return a normalized form of ``names``.

    INPUT:

    - ``names`` -- a string, a list of strings or ``None``

    - ``rank`` -- the number of names to normalize

    EXAMPLES::

        sage: from sage.modules.ore_module import normalize_names

    When ``names`` is a string, indices are added::

        sage: normalize_names('e', 3)
        ('e0', 'e1', 'e2')

    When ``names`` is a list or a tuple, it remains untouched
    except that it is always casted to a tuple (in order to be
    hashable and serve as a key)::

        sage: normalize_names(['u', 'v', 'w'], 3)
        ('u', 'v', 'w')

    Similarly, when ``names`` is ``None``, nothing is returned::

        sage: normalize_names(None, 3)

    If the number of names is not equal to ``rank``, an error
    is raised::

        sage: normalize_names(['u', 'v', 'w'], 2)
        Traceback (most recent call last):
        ...
        ValueError: the number of given names does not match the rank of the Ore module
    """
    if names is None:
        pass
    elif isinstance(names, (list, tuple)):
        if rank != len(names):
            raise ValueError("the number of given names does not match the rank of the Ore module")
        names = tuple([str(name) for name in names])
    elif isinstance(names, str):
        names = tuple([names + str(i) for i in range(rank)])
    else:
        raise ValueError("names must be a string or a list/tuple of strings")
    return names

class OreModule(UniqueRepresentation, FreeModule_ambient):
    r"""
    Generic class for Ore modules.
    """
    Element = OreModuleElement

    def __classcall_private__(cls, mat, twist, names=None, category=None):
        r"""
        Normalize the input before passing it to the init function
        (useful to ensure the uniqueness assuption).

        INPUT:

        - ``mat`` -- the matrix defining the action of the Ore variable

        - ``twist`` -- the twisting morphism/derivation

        - ``names`` (default: ``None``) -- a string of a list of stings,
          the names of the vector of the canonical basis; if ``None``,
          elements are represented as vectors in `K^d`

        - ``category`` (default: ``None``) -- the category of this
          Ore module

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^2 + z

            sage: M1 = S.quotient_module(P)
            sage: M2 = S.quotient_module(P, names='e')
            sage: M3.<e0,e1> = S.quotient_module(P)

            sage: M1 is M2
            False
            sage: M2 is M3
            True
        """
        base = mat.base_ring()
        if category is None:
            category = OreModules(base, twist)
        rank = mat.nrows()
        if mat.ncols() != rank:
            raise ValueError("matrix must be square")
        names = normalize_names(names, rank)
        return cls.__classcall__(cls, mat, category._ore, names, category)

    def __init__(self, mat, ore, names, category):
        r"""
        Initialize this Ore module.

        INPUT:

        - ``mat`` -- the matrix defining the action of the Ore variable

        - ``ore`` -- the underlying Ore polynomial ring

        - ``names`` -- a string of a list of stings,
          the names of the vector of the canonical basis; if ``None``,
          elements are represented as vectors in `K^d`

        - ``category`` -- the category of this Ore module

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)  # indirect doctest
            sage: type(M)
            <class 'sage.modules.ore_module.OreModule_with_category'>

            sage: TestSuite(M).run()
        """
        base = mat.base_ring()
        rank = mat.nrows()
        FreeModule_ambient.__init__(self, base, rank, category=category)
        self.register_action(ScalarAction(base, self, True, operator.mul))
        self._ore = ore
        self._ore_category = category
        self._names = names
        if names is not None:
            self._latex_names = [latex_variable_name(name) for name in names]
        self._submodule_class = OreSubmodule
        self._quotientModule_class = OreQuotientModule
        self._pseudohom = FreeModule_ambient.pseudohom(self, mat, ore, codomain=self)

    def _repr_(self):
        r"""
        Return a string representation of this Ore module.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: M._repr_()
            'Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5'
            sage: N = S.quotient_module(X^2 + z, names='e')
            sage: N._repr_()
            'Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5'

        ::

            sage: K.<z> = Frac(GF(17)['z'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M = S.quotient_module(X^2 + z)
            sage: M._repr_()
            'Ore module of rank 2 over Fraction Field of Univariate Polynomial Ring in z over Finite Field of size 17 twisted by d/dz'
        """
        s = "Ore module "
        if self._names is None:
            s += "of rank %s " % self.rank()
        else:
            s += "<" + ", ".join(self._names) + "> "
        s += "over %s %s" % (self.base_ring(), self._ore._repr_twist())
        return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this Ore module.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: S.<X> = OrePolynomialRing(K, Frob)
            sage: M = S.quotient_module(X^2 + z)
            sage: latex(M)
            \texttt{Ore module of rank } 2\texttt{ over } \Bold{F}_{5^{3}} \texttt{ twisted by } z \mapsto z^{5}
            sage: N = S.quotient_module(X^2 + z, names='e')
            sage: latex(N)
            \left<e_{0}, e_{1}\right>_{\Bold{F}_{5^{3}} , z \mapsto z^{5} }

        ::

            sage: T.<Y> = OrePolynomialRing(K, Frob^3, polcast=False)
            sage: M = T.quotient_module(Y^2 + z^2)
            sage: latex(M)
            \texttt{Ore module of rank } 2\texttt{ over } \Bold{F}_{5^{3}}\texttt{ untwisted}
        """
        if self._names is None:
            s = "\\texttt{Ore module of rank } %s" % self.rank()
            s += "\\texttt{ over } %s" % latex(self.base_ring())
            twist = self._ore._latex_twist()
            if twist == "":
                s += "\\texttt{ untwisted}"
            else:
                s += "\\texttt{ twisted by }" + twist
        else:
            s = "\\left<" + ", ".join(self._latex_names) + "\\right>"
            s += "_{%s" % latex(self.base_ring())
            twist = self._ore._latex_twist()
            if twist != "":
                s += "," + twist
            s += "}"
        return s

    def _repr_element(self, x):
        r"""
        Return a string representation of the element `x` in
        this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: M((z, z^2))  # indirect doctest
            (z, z^2)
        """
        return FreeModuleElement_generic_dense._repr_(x)

    def _latex_element(self, x):
        r"""
        Return a LaTeX representation of the element `x` in
        this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: v = M((z, z^2))
            sage: latex(v)  # indirect doctest
            \left(z,\,z^{2}\right)
        """
        return FreeModuleElement_generic_dense._latex_(x)

    def _coerce_map_from_(self, S):
        r"""
        Return a coercion map from `M` to ``self``, or ``None``.

        This method always returns ``None``; all coercions between
        Ore modules are currently handled using the method
        :meth:`register_coercion`.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)

            sage: M._coerce_map_from_(N)
            sage: M.coerce_map_from(N)
            Ore module morphism:
              From: Ore module of rank 1 over Finite Field in z of size 5^3 twisted by z |--> z^5
              To:   Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
        """
        pass

    def is_zero(self):
        r"""
        Return ``True`` if this Ore module is reduced to zero.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: M
            Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: M.is_zero()
            False

            sage: Q = M.quo(M)
            sage: Q.is_zero()
            True
        """
        return self.rank() == 0

    def rename_basis(self, names, coerce=False):
        r"""
        Return the same Ore module with the given naming
        for the vectors in its distinguished basis.

        INPUT:

        - ``names`` -- a string or a list of strings, the
          new names

        - ``coerce`` (default: ``False``) -- a boolean; if
          ``True``, a coercion map from this Ore module to
          renamed version is set

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z)
            sage: M
            Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

            sage: Me = M.rename_basis('e')
            sage: Me
            Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5

        Now compare how elements are displayed::

            sage: M.random_element()   # random
            (3*z^2 + 4*z + 2, 3*z^2 + z)
            sage: Me.random_element()  # random
            (2*z + 4)*e0 + (z^2 + 4*z + 4)*e1

        At this point, there is no coercion map between ``M``
        and ``Me``. Therefore, adding elements in both parents
        results in an error:

            sage: M.random_element() + Me.random_element()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +:
                'Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5' and
                'Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5'

        In order to set this coercion, one should define ``Me``
        by passing the extra argument ``coerce=True``::

            sage: Me = M.rename_basis('e', coerce=True)
            sage: M.random_element() + Me.random_element()  # random
            2*z^2*e0 + (z^2 + z + 4)*e1

        .. WARNING::

            Use ``coerce=True`` with extreme caution. Indeed,
            setting inappropriate coercion maps may result in a
            circular path in the coercion graph which, in turn,
            could eventually break the coercion system.

        Note that the bracket construction also works::

            sage: M.<v,w> = M.rename_basis()
            sage: M
            Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5

        In this case, `v` and `w` are automatically defined::

            sage: v + w
            v + w
        """
        rank = self.rank()
        names = normalize_names(names, rank)
        cls = self.__class__
        M = cls.__classcall__(cls, self._pseudohom._matrix,
                              self._ore, names, self._ore_category)
        if coerce:
            mat = identity_matrix(self.base_ring(), rank)
            id = self.hom(mat, codomain=M)
            M._unset_coercions_used()
            M.register_coercion(id)
        return M

    def pseudohom(self):
        r"""
        Return the pseudomorphism giving the action of the Ore
        variable on this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 - z^2*X + (z + 2)
            sage: M = S.quotient_module(P)
            sage: M.pseudohom()
            Free module pseudomorphism (twisted by z |--> z^5) defined by the matrix
            [      0       1       0]
            [      0       0       1]
            [4*z + 3     z^2     4*z]
            Domain: Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5
            Codomain: Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5

        .. SEEALSO:

            :meth:`matrix`
        """
        return self._pseudohom

    def ore_ring(self, names='x', action=True):
        r"""
        Return the underlying Ore polynomial ring.

        INPUT:

        - ``names`` (default: ``x``) -- a string, the name
          of the variable

        - ``action`` (default: ``True``) -- a boolean; if
          ``True``, an action of the Ore polynomial ring on
          the Ore module is set

        EXAMPLES::

            sage: K.<a> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<e1,e2> = S.quotient_module(X^2 - a)
            sage: M.ore_ring()
            Ore Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5

        We can use a different variable name::

            sage: M.ore_ring('Y')
            Ore Polynomial Ring in Y over Finite Field in a of size 5^3 twisted by a |--> a^5

        Alternatively, one can use the following shortcut::

            sage: T.<Z> = M.ore_ring()
            sage: T
            Ore Polynomial Ring in Z over Finite Field in a of size 5^3 twisted by a |--> a^5

        In all the above cases, an action of the returned Ore polynomial
        ring on `M` is registered::

            sage: Z*e1
            e2
            sage: Z*e2
            a*e1

        Specifying ``action=False`` prevents this to happen::

            sage: T.<U> = M.ore_ring(action=False)
            sage: U*e1
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
                'Ore Polynomial Ring in U over Finite Field in a of size 5^3 twisted by a |--> a^5' and
                'Ore module <e1, e2> over Finite Field in a of size 5^3 twisted by a |--> a^5'
        """
        S = self._ore_category.ore_ring(names)
        if action:
            self._unset_coercions_used()
            self.register_action(OreAction(S, self, True, operator.mul))
        return S

    def twisting_morphism(self):
        r"""
        Return the twisting morphism corresponding to this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X + z)
            sage: M.twisting_morphism()
            Frobenius endomorphism z |--> z^5 on Finite Field in z of size 5^3

        When the twisting morphism is trivial (that is, the identity),
        nothing is returned::

            sage: R.<t> = QQ[]
            sage: T.<Y> = OrePolynomialRing(R, R.derivation())
            sage: M = T.quotient_module(Y + t^2)
            sage: M.twisting_morphism()

        .. SEEALSO:

            :meth:`twisting_derivation`
        """
        return self._ore.twisting_morphism()

    def twisting_derivation(self):
        r"""
        Return the twisting derivation corresponding to this Ore module.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: T.<Y> = OrePolynomialRing(R, R.derivation())
            sage: M = T.quotient_module(Y + t^2)
            sage: M.twisting_derivation()
            d/dt

        When the twisting derivation in zero, nothing is returned::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X + z)
            sage: M.twisting_derivation()

        .. SEEALSO:

            :meth:`twisting_morphism`
        """
        return self._ore.twisting_derivation()

    def matrix(self):
        r"""
        Return the matrix giving the action of the Ore variable
        on this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + z*X^2 - z^2*X + (z + 2)
            sage: M = S.quotient_module(P)
            sage: M.matrix()
            [      0       1       0]
            [      0       0       1]
            [4*z + 3     z^2     4*z]

        We recognize the companion matrix attached to the Ore
        polynomial `P`. This is of course not a coincidence given
        that the pseudomorphism corresponds to the left multiplication

        .. SEEALSO:

            :meth:`pseudohom`
        """
        return self._pseudohom.matrix()

    def basis(self):
        r"""
        Return the canonical basis of this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^3 - z)
            sage: M.basis()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        rank = self.rank()
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        coeffs = [zero] * rank
        B = [ ]
        for i in range(rank):
            coeffs[i] = one
            B.append(self(coeffs))
            coeffs[i] = zero
        return B

    def gens(self):
        r"""
        Return the canonical basis of this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^3 - z)
            sage: M.gens()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return self.basis()

    def gen(self, i):
        r"""
        Return the `i`th vector in the canonical basis
        of this Ore module.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^3 - z)
            sage: M.gen(0)
            (1, 0, 0)
            sage: M.gen(1)
            (0, 1, 0)
            sage: M.gen(2)
            (0, 0, 1)
            sage: M.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: generator is not defined
        """
        rank = self.rank()
        if i < 0 or i >= rank:
            raise IndexError("generator is not defined")
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        coeffs = [zero] * rank
        coeffs[i] = one
        return self(coeffs)

    def an_element(self):
        r"""
        Return an element of this Ore module.

        EXAMPLES:

        When the Ore module is not zero, the returned element
        is the first vector of the distinguished basis::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<u,v> = S.quotient_module(X^2 - t)
            sage: M.an_element()
            u

       On the contrary, when the Ore module vanishes, the
       returned element is of course zero::

            sage: N = M / u
            sage: N
            Ore module of rank 0 over Fraction Field of Univariate Polynomial Ring in t over Rational Field twisted by d/dt
            sage: N.an_element()
            0
        """
        if self.rank() > 0:
            return self.gen(0)
        else:
            return self.zero()

    def random_element(self, *args, **kwds):
        r"""
        Return a random element in this Ore module.

        Extra arguments are passed to the random generator
        of the base ring.

        EXAMPLES::

            sage: A.<t> = QQ['t']
            sage: S.<X> = OrePolynomialRing(A, A.derivation())
            sage: M = S.quotient_module(X^3 - t, names='e')
            sage: M.random_element()   # random
            (-1/2*t^2 - 3/4*t + 3/2)*e0 + (-3/2*t^2 - 3*t + 4)*e1 + (-6*t + 2)*e2

            sage: M.random_element(degree=5)   # random
            (4*t^5 - 1/2*t^4 + 3/2*t^3 + 6*t^2 - t - 1/10)*e0 + (19/3*t^5 - t^3 - t^2 + 1)*e1 + (t^5 + 4*t^4 + 4*t^2 + 1/3*t - 33)*e2
        """
        K = self.base_ring()
        r = self.rank()
        vs = [K.random_element(*args, **kwds) for _ in range(r)]
        return self(vs)

    def module(self):
        r"""
        Return the underlying free module of this Ore module.

        EXAMPLES::

            sage: A.<t> = QQ['t']
            sage: S.<X> = OrePolynomialRing(A, A.derivation())
            sage: M = S.quotient_module(X^3 - t)
            sage: M
            Ore module of rank 3 over Univariate Polynomial Ring in t over Rational Field twisted by d/dt

            sage: M.module()
            Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in t over Rational Field
        """
        return self.base_ring() ** self.rank()

    def _Hom_(self, codomain, category):
        r"""
        Return the space of Ore morphisms from this Ore module
        to ``codomain``.

        INPUT:

        - ``codomain`` -- a Ore module

        - ``category`` -- the category in which the morphisms are

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 - z)
            sage: N = S.quotient_module(X^3 - z)

            sage: Hom(M, N)  # indirect doctest
            Set of Morphisms
            from Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
            to Ore module of rank 3 over Finite Field in z of size 5^3 twisted by z |--> z^5
            in Category of enumerated finite dimensional Ore modules with basis over Finite Field in z of size 5^3 twisted by z |--> z^5

        ::

            sage: End(M)  # indirect doctest
            Set of Morphisms
            from Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
            to Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5
            in Category of enumerated finite dimensional Ore modules with basis over Finite Field in z of size 5^3 twisted by z |--> z^5
        """
        from sage.modules.ore_module_homspace import OreModule_homspace
        return OreModule_homspace(self, codomain)

    def hom(self, im_gens, codomain=None):
        r"""
        Return the morphism from this Ore module to ``codomain``
        defined by ``im_gens``.

        INPUT:

        - ``im_gens`` -- a datum defining the morphism to build;
          it could either a list, a tuple, a dictionary or a morphism
          of Ore modules

        - ``codomain`` (default: ``None``) -- a Ore module, the
          codomain of the morphism; if ``None``, it is inferred from
          ``im_gens``

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: P = X^3 + 2*t*X^2 + (t^2 + 2)*X + t
            sage: Q = t*X^2 - X + 1

            sage: U = S.quotient_module(P, names='u')
            sage: U.inject_variables()
            Defining u0, u1, u2
            sage: V = S.quotient_module(P*Q, names='v')
            sage: V.inject_variables()
            Defining v0, v1, v2, v3, v4

        The first method for creating a morphism from `U` to `V` is
        to explicitely write down its matrix in the canonical bases::

            sage: mat = matrix(3, 5, [1, 4, t, 0, 0,
            ....:                     0, 1, 0, t, 0,
            ....:                     0, 0, 1, 1, t])
            sage: f = U.hom(mat, codomain=V)
            sage: f
            Ore module morphism:
              From: Ore module <u0, u1, u2> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module <v0, v1, v2, v3, v4> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt

        This method is however not really convenient because it
        requires to compute beforehand all the entries of the
        defining matrix.
        Instead, we can pass the list of images of the generators::

            sage: g = U.hom([Q*v0, X*Q*v0, X^2*Q*v0])
            sage: g
            Ore module morphism:
              From: Ore module <u0, u1, u2> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module <v0, v1, v2, v3, v4> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: g.matrix()
            [1 4 t 0 0]
            [0 1 0 t 0]
            [0 0 1 1 t]

        One can even give the values of the morphism on a smaller
        set as soon as the latter generates the domain as Ore module.
        The syntax uses dictionaries as follows::

            sage: h = U.hom({u0: Q*v0})
            sage: h
            Ore module morphism:
              From: Ore module <u0, u1, u2> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module <v0, v1, v2, v3, v4> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: g == h
            True

        Finally ``im_gens`` can also be itself a Ore morphism, in which
        case SageMath tries to cast it into a morphism with the requested
        domains and codomains.
        As an example below, we restrict `g` to a subspace::

            sage: C.<c0,c1> = U.span((X + t)*u0)
            sage: gC = C.hom(g)
            sage: gC
            Ore module morphism:
              From: Ore module <c0, c1> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module <v0, v1, v2, v3, v4> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt

            sage: g(c0) == gC(c0)
            True
            sage: g(c1) == gC(c1)
            True

        TESTS::

            sage: U.hom(0)
            Traceback (most recent call last):
            ...
            ValueError: im_gens must be a list, a tuple, a dictionary, a matrix or a Ore module morphism

            sage: U.hom([Q*v0])
            Traceback (most recent call last):
            ...
            ValueError: wrong number of generators

            sage: U.hom({u0: Q*v0, u1: Q*v0})
            Traceback (most recent call last):
            ...
            ValueError: does not define a morphism of Ore modules

            sage: U.hom({(X+t)*u0: (X+t)*Q*v0})
            Traceback (most recent call last):
            ...
            ValueError: does not define a morphism of Ore modules
        """
        from sage.modules.ore_module_morphism import OreModuleMorphism
        if codomain is None:
            if isinstance(im_gens, Matrix):
                codomain = self
            elif isinstance(im_gens, OreModuleMorphism):
                codomain = im_gens.codomain()
            elif isinstance(im_gens, (list, tuple)):
                codomain = Sequence(im_gens).universe()
            elif isinstance(im_gens, dict):
                codomain = Sequence(im_gens.values()).universe()
            else:
                raise ValueError("im_gens must be a list, a tuple, a dictionary, a matrix or a Ore module morphism")
        H = self.Hom(codomain)
        return H(im_gens)

    def multiplication_map(self, P):
        r"""
        Return the multiplication by `P` acting on this Ore module.

        INPUT:

        - ``P`` -- a scalar in the base ring, or a Ore polynomial

        EXAMPLES::

            sage: K.<a> = GF(7^5)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^3 + a*X^2 + X - a^2
            sage: M = S.quotient_module(P)

        We define the scalar multiplication by an element in the base ring::

            sage: f = M.multiplication_map(3)
            sage: f
            Ore module endomorphism of Ore module of rank 3 over Finite Field in a of size 7^5 twisted by a |--> a^7
            sage: f.matrix()
            [3 0 0]
            [0 3 0]
            [0 0 3]

        Be careful that an element in the base ring defines a Ore morphism
        if and only if it is fixed by the twisting morphisms and killed by
        the derivation (otherwise the multiplication by this element does
        not commute with the Ore action).
        In SageMath, attempting to create the multiplication by an element
        which does not fulfill these requirements leads to an error::

            sage: M.multiplication_map(a)
            Traceback (most recent call last):
            ...
            ValueError: does not define a morphism of Ore modules

        As soon as it defines a Ore morphism, one can also build the left
        multiplication by an Ore polynomial::

            sage: g = M.multiplication_map(X^5)
            sage: g
            Ore module endomorphism of Ore module of rank 3 over Finite Field in a of size 7^5 twisted by a |--> a^7
            sage: g.matrix()
            [    3*a^4 + 3*a^3 + 6*a^2 + 5*a       4*a^4 + 5*a^3 + 2*a^2 + 6         6*a^4 + 6*a^3 + a^2 + 4]
            [                        a^2 + 3 5*a^4 + 5*a^3 + 6*a^2 + 4*a + 1                 a^3 + 5*a^2 + 4]
            [6*a^4 + 6*a^3 + 3*a^2 + 3*a + 1         4*a^4 + 2*a^3 + 3*a + 5 6*a^4 + 6*a^3 + 2*a^2 + 5*a + 2]

        We check that the characteristic polynomial of `g` is the reduced
        norm of the Ore polynomial `P` we started with (this is a classical
        property)::

            sage: g.charpoly()
            x^3 + 4*x^2 + 2*x + 5
            sage: P.reduced_norm(var='x')
            x^3 + 4*x^2 + 2*x + 5
        """
        if isinstance(P, OrePolynomial):
            S = P.parent()
            ore = self._ore
            if S._morphism != ore._morphism or S._derivation != ore._derivation:
                raise ValueError("twist does not match")
            action = OreAction(S, self, True, operator.mul)
            M = matrix([action._act_(P, x).list() for x in self.basis()])
        else:
            P = self.base_ring()(P)
            r = self.rank()
            M = matrix(r, r, P)
        H = self.Hom(self)
        return H(M)

    def identity_morphism(self):
        r"""
        Return the identity morphism of this Ore module.

        EXAMPLES::

            sage: K.<a> = GF(7^5)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<u,v> = S.quotient_module(X^2 + a*X + a^2)
            sage: id = M.identity_morphism()
            sage: id
            Ore module endomorphism of Ore module <u, v> over Finite Field in a of size 7^5 twisted by a |--> a^7

            sage: id(u)
            u
            sage: id(v)
            v
        """
        H = self.Hom(self)
        M = identity_matrix(self.base_ring(), self.rank())
        return H(M)

    def _span(self, gens):
        r"""
        Return a matrix defining a basis of the submodule of this
        Ore modules generated by ``gens``.

        INPUT:

        - ``gens`` -- a list of vectors or submodules of this Ore module

        TESTS::

            sage: K.<a> = GF(7^5)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X^2 + a
            sage: Q = X^3 + a^2*X + 1
            sage: M = S.quotient_module(P*Q, names='e')
            sage: M.inject_variables()
            Defining e0, e1, e2, e3, e4

            sage: M._span([Q*e0])
            [                      1                       0 a^4 + 5*a^3 + 6*a^2 + a                       1                   6*a^2]
            [                      0                       1 2*a^4 + 6*a^2 + 2*a + 3                       0                       1]

            sage: M._span([Q*e0, e1])
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

            sage: N = M.span(Q*e0)
            sage: M._span([N, e2])
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
        """
        base = self.base_ring()
        rank = self.rank()
        f = self._pseudohom
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        rows = []
        for gen in gens:
            if isinstance(gen, OreModule):
                incl = self.coerce_map_from(gen)
                if incl is None:
                    raise ValueError("not canonically a submodule")
                rows += incl._matrix.rows()
            elif isinstance(gen, OreModuleElement):
                rows.append(self(gen).list())
        if len(rows) < 2*rank:
            zero = rank * [base.zero()]
            rows += (2*rank - len(rows)) * [rank*[0]]
        M = matrix(base, rows)
        M.echelonize()
        oldr = 0
        r = M.rank()
        iter = 1
        while r > oldr:
            for i in range(r):
                v = M.row(i)
                for _ in range(iter):
                    v = f(v)
                v = v.list()
                for j in range(rank):
                    M[i+r,j] = v[j]
            M.echelonize()
            oldr = r
            r = M.rank()
            iter *= 2
        return M.matrix_from_rows(range(r))

    def span(self, gens, names=None):
        r"""
        Return the submodule of this Ore module generated by ``gens``.

        INPUT:

        - ``gens`` -- a list of vectors or submodules of this Ore module

        - ``names`` (default: ``None``) -- the name of the vectors in a
          basis of this submodule

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: P = X^2 + t*X + 1
            sage: M = S.quotient_module(P^3, names='e')
            sage: M.inject_variables()
            Defining e0, e1, e2, e3, e4, e5

        We create the submodule `M P`::

            sage: MP = M.span([P*e0])
            sage: MP
            Ore module of rank 4 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: MP.basis()
            [e0 + (t^4 + t^2 + 3)*e4 + t^3*e5,
             e1 + (4*t^3 + 2*t)*e4 + (4*t^2 + 3)*e5,
             e2 + (2*t^2 + 2)*e4 + 2*t*e5,
             e3 + 4*t*e4 + 4*e5]

        When there is only one generated, encapsulating it in a list in
        not necessary; one can equally write::

            sage: MP = M.span(P*e0)

        If we want, one can give names to the basis of the submodule using
        the attribute ``names``::

            sage: MP2 = M.span(P^2*e0, names='u')
            sage: MP2.inject_variables()
            Defining u0, u1
            sage: MP2.basis()
            [u0, u1]

            sage: M(u0)
            e0 + (t^2 + 4)*e2 + 3*t^3*e3 + (t^2 + 1)*e4 + 3*t*e5

        Note that a coercion map from the submodule to the ambient module
        is automatically set::

            sage: M.has_coerce_map_from(MP2)
            True

        Therefore, combining elements of ``M`` and ``MP2`` in the same
        expression perfectly works::

            sage: t*u0 + e1
            t*e0 + e1 + (t^3 + 4*t)*e2 + 3*t^4*e3 + (t^3 + t)*e4 + 3*t^2*e5

        Here is an example with multiple generators::

            sage: MM = M.span([MP2, P*e1])
            sage: MM.basis()
            [e0, e1, e2, e3, e4, e5]

        In this case, we obtain the whole space.

        Creating submodules of submodules is also allowed::

            sage: N = MP.span(P^2*e0)
            sage: N
            Ore module of rank 2 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: N.basis()
            [e0 + (t^2 + 4)*e2 + 3*t^3*e3 + (t^2 + 1)*e4 + 3*t*e5,
             e1 + (4*t^2 + 4)*e3 + 3*t*e4 + 4*e5]

        .. SEEALSO::

            :meth:`quotient`
        """
        gens = self._span(gens)
        return self._submodule_class(self, gens, names=names)

    def quotient(self, sub, names=None, check=True):
        r"""
        Return the quotient of this Ore module by the submodule
        generated by ``gens``.

        INPUT:

        - ``gens`` -- a list of vectors or submodules of this Ore module

        - ``names`` (default: ``None``) -- the name of the vectors in a
          basis of the quotient

        - ``check`` (default: ``True``) -- a boolean, ignored

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: P = X^2 + t*X + 1
            sage: M = S.quotient_module(P^3, names='e')
            sage: M.inject_variables()
            Defining e0, e1, e2, e3, e4, e5

        We create the quotient `M/MP`::

            sage: modP = M.quotient(P*e0)
            sage: modP
            Ore module of rank 2 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt

        As a shortcut, we can write ``quo`` instead of ``quotient`` or even
        use the `/` operator::

            sage: modP = M / (P*e0)
            sage: modP
            Ore module of rank 2 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt

        By default, the vectors in the quotient have the same names as their
        representatives in `M`.

            sage: modP.basis()
            [e4, e5]

        One can override this by setting the attributes ``names``::

            sage: modP = M.quo(P*e0, names='u')
            sage: modP.inject_variables()
            Defining u0, u1
            sage: modP.basis()
            [u0, u1]

        Note that a coercion map from the initial Ore module to its quotient
        is automatically set. As a consequence, combining elements of ``M``
        and ``modP`` in the same formula works::

            sage: t*u0 + e1
            (t^3 + 4*t)*u0 + (t^2 + 2)*u1

        One can combine the construction of quotients and submodules without
        trouble. For instance, here we build the space `M P / M P^2`::

            sage: modP2 = M / (P^2*e0)
            sage: N = modP2.span(P*e0)
            sage: N
            Ore module of rank 2 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: N.basis()
            [e2 + (2*t^2 + 2)*e4 + 2*t*e5,
             e3 + 4*t*e4 + 4*e5]

        .. SEEALSO::

            :meth:`quo`, :meth:`span`
        """
        gens = self._span(sub)
        return self._quotientModule_class(self, gens, names=names)

    quo = quotient

    def __eq__(self, other):
        r"""
        Return ``True`` if this Ore module is the same than ``other``.

        TESTS:

        Different names lead to different parents::

            sage: K.<a> = GF(7^5)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<u,v> = S.quotient_module(X^2 + a*X + a^2)
            sage: N.<e0,e1> = S.quotient_module(X^2 + a*X + a^2)
            sage: M == N
            False

        However, different syntaxes resulting in the same names lead
        to the same parent::

            sage: N2 = S.quotient_module(X^2 + a*X + a^2, names='e')
            sage: N == N2
            True
        """
        return self is other

    def __hash__(self):
        r"""
        Return a hash of this Ore module.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^3 - z)
            sage: hash(M)  # random
            128873304640624
        """
        return id(self)

# Submodules
############

class OreSubmodule(OreModule):
    r"""
    Class for submodules of Ore modules.
    """
    def __classcall_private__(cls, ambient, gens, names):
        r"""
        Normalize the input before passing it to the init function
        (useful to ensure the uniqueness assuption).

        INPUT:

        - ``ambient`` -- a Ore module, the ambient module where
          this submodule sits

        - ``gens`` -- a list of generators (formatted as coordinates
          vectors) of this submodule

        - ``names`` -- the name of the vectors of the basis of
          the submodule, or ``None``

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N1 = M.span((X + z)*v)
            sage: N2 = M.span((X + z^5)*w)
            sage: N1 is N2
            True

        ::

            sage: R.<t> = QQ[]
            sage: S.<X> = OrePolynomialRing(R, R.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: M.span((X + t)*v)
            Traceback (most recent call last):
            ...
            NotImplementedError: Ore submodules are currently only implemented over fields
        """
        base = ambient.base_ring()
        if base not in Fields():
            raise NotImplementedError("Ore submodules are currently only implemented over fields")
        if isinstance(gens, Matrix):
            basis = gens
        else:
            basis = matrix(base, gens)
        basis = basis.echelon_form()
        basis.set_immutable()
        rank = basis.rank()
        if basis.nrows() != rank:
            basis = basis.matrix_from_rows(range(rank))
        names = normalize_names(names, rank)
        return cls.__classcall__(cls, ambient, basis, names)

    def __init__(self, ambient, basis, names):
        r"""
        Initialize this Ore submodule.

        INPUT:

        - ``ambient`` -- a Ore module, the ambient module where
          this submodule sits

        - ``basis`` -- the echelon basis of this submodule

        - ``names`` -- the name of the vectors of the basis of
          the submodule, or ``None``

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)  # indirect doctest
            sage: type(N)
            <class 'sage.modules.ore_module.OreSubmodule_with_category'>

            sage: TestSuite(N).run()
        """
        from sage.modules.ore_module_morphism import OreModuleRetraction
        base = ambient.base_ring()
        self._ambient = ambient
        self._basis = basis
        rows = [basis.solve_left(ambient(x).image()) for x in basis.rows()]
        OreModule.__init__(self, matrix(base, rows),
                           ambient.ore_ring(action=False),
                           names, ambient._ore_category)
        coerce = self.hom(basis, codomain=ambient)
        ambient.register_coercion(coerce)
        self._inject = coerce.__copy__()
        self.register_conversion(OreModuleRetraction(ambient, self))

    def _repr_element(self, x):
        r"""
        Return a string representation of ``x``.

        By default, elements in a Ore submodule are printed as their
        images in the ambient module.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)
            sage: N.an_element()  # indirect doctest
            v + (3*z^2 + 4)*w
        """
        return self._ambient(x)._repr_()

    def _latex_element(self, x):
        r"""
        Return a LaTeX representation of ``x``.

        By default, elements in a Ore submodule are rendered as their
        images in the ambient module.

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)
            sage: latex(N.an_element())  # indirect doctest
            v + \left(3 z^{2} + 4\right) w
        """
        return self._ambient(x)._latex_()

    def ambient(self):
        r"""
        Return the ambient Ore module in which this submodule lives.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)
            sage: N.ambient()
            Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: N.ambient() is M
            True
        """
        return self._ambient

    def rename_basis(self, names, coerce=False):
        r"""
        Return the same Ore module with the given naming
        for the vectors in its distinguished basis.

        INPUT:

        - ``names`` -- a string or a list of strings, the
          new names

        - ``coerce`` (default: ``False``) -- a boolean; if
          ``True``, a coercion map from this Ore module to
          renamed version is set

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z^2)
            sage: M
            Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

            sage: Me = M.rename_basis('e')
            sage: Me
            Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5

        Now compare how elements are displayed::

            sage: M.random_element()   # random
            (3*z^2 + 4*z + 2, 3*z^2 + z)
            sage: Me.random_element()  # random
            (2*z + 4)*e0 + (z^2 + 4*z + 4)*e1

        At this point, there is no coercion map between ``M``
        and ``Me``. Therefore, adding elements in both parents
        results in an error:

            sage: M.random_element() + Me.random_element()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +:
                'Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5' and
                'Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5'

        In order to set this coercion, one should define ``Me``
        by passing the extra argument ``coerce=True``::

            sage: Me = M.rename_basis('e', coerce=True)
            sage: M.random_element() + Me.random_element()  # random
            2*z^2*e0 + (z^2 + z + 4)*e1

        .. WARNING::

            Use ``coerce=True`` with extreme caution. Indeed,
            setting inappropriate coercion maps may result in a
            circular path in the coercion graph which, in turn,
            could eventually break the coercion system.

        Note that the bracket construction also works::

            sage: M.<v,w> = M.rename_basis()
            sage: M
            Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5

        In this case, `v` and `w` are automatically defined::

            sage: v + w
            v + w

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X + z
            sage: A.<v,w> = S.quotient_module(P^2)
            sage: M = A.span(P*v)
            sage: Me = M.rename_basis('e', coerce=True)
            sage: M.an_element() + Me.an_element()
            2*e0

        """
        rank = self.rank()
        names = normalize_names(names, rank)
        cls = self.__class__
        M = cls.__classcall__(cls, self._ambient, self._basis, names)
        if coerce:
            mat = identity_matrix(self.base_ring(), rank)
            id = self.hom(mat, codomain=M)
            M._unset_coercions_used()
            M.register_coercion(id)
        return M

    def injection_morphism(self):
        r"""
        Return the inclusion of this submodule in the ambient space.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)
            sage: N.injection_morphism()
            Ore module morphism:
              From: Ore module of rank 1 over Finite Field in z of size 5^3 twisted by z |--> z^5
              To:   Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
        """
        return self._inject

    def morphism_restriction(self, f):
        r"""
        Return the restriction of `f` to this submodule.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M.<v,w> = S.quotient_module((X + z)^2)
            sage: N = M.span((X + z)*v)

            sage: f = M.multiplication_map(X^3)
            sage: f
            Ore module endomorphism of Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5

            sage: g = N.morphism_restriction(f)
            sage: g
            Ore module morphism:
              From: Ore module of rank 1 over Finite Field in z of size 5^3 twisted by z |--> z^5
              To:   Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: g.matrix()
            [        3 4*z^2 + 2]
        """
        if f.domain() is not self._ambient:
            raise ValueError
        return f * self._inject

    def morphism_corestriction(self, f):
        r"""
        If the image of `f` is contained in this submodule,
        return the corresponding corestriction of `f`.

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X + z
            sage: M.<v,w> = S.quotient_module(P^2)
            sage: N = M.span(P*v)

            sage: f = M.hom({v: P*v})
            sage: f
            Ore module endomorphism of Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5

            sage: g = N.morphism_corestriction(f)
            sage: g
            Ore module morphism:
              From: Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5
              To:   Ore module of rank 1 over Finite Field in z of size 5^3 twisted by z |--> z^5
            sage: g.matrix()
            [    z]
            [4*z^2]

        When the image of the morphism is not contained in this submodule,
        an error is raised::

            sage: h = M.multiplication_map(X^3)
            sage: N.morphism_corestriction(h)
            Traceback (most recent call last):
            ...
            ValueError: the image of the morphism is not contained in this submodule
        """
        if f.codomain() is not self._ambient:
            raise ValueError
        rows = []
        basis = self._basis
        try:
            rows = [basis.solve_left(y) for y in f._matrix.rows()]
        except ValueError:
            raise ValueError("the image of the morphism is not contained in this submodule")
        return f.domain().hom(rows, codomain=self)

    _hom_change_domain = morphism_restriction
    _hom_change_codomain = morphism_corestriction

# Quotients
###########

class OreQuotientModule(OreModule):
    r"""
    Class for quotients of Ore modules.
    """
    def __classcall_private__(cls, cover, gens, names):
        r"""
        Normalize the input before passing it to the init function
        (useful to ensure the uniqueness assuption).

        INPUT:

        - ``cover`` -- a Ore module, the cover module of this
          quotient

        - ``gens`` -- a list of generators (formatted as coordinates
          vectors) of the submodule by which we quotient out

        - ``names`` -- the name of the vectors of the basis of
          the quotient, or ``None``

        TESTS::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: Q1 = M.quo((X + t)*v)
            sage: Q2 = M.quo(v + (X + t)*w)
            sage: Q1 is Q2
            True

        ::

            sage: R.<t> = QQ[]
            sage: S.<X> = OrePolynomialRing(R, R.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: M.quo((X + t)*v)
            Traceback (most recent call last):
            ...
            NotImplementedError: quotient of Ore modules are currently only implemented over fields
        """
        base = cover.base_ring()
        if base not in Fields():
            raise NotImplementedError("quotient of Ore modules are currently only implemented over fields")
        if isinstance(gens, Matrix):
            basis = gens
        else:
            basis = matrix(base, gens)
        basis = basis.echelon_form()
        basis.set_immutable()
        rank = basis.rank()
        if basis.nrows() != rank:
            basis = basis.matrix_from_rows(range(rank))
        names = normalize_names(names, cover.rank() - rank)
        return cls.__classcall__(cls, cover, basis, names)

    def __init__(self, cover, basis, names):
        r"""
        Initialize this Ore quotient.

        INPUT:

        - ``cover`` -- a Ore module, the cover module of this
          quotient

        - ``basis`` -- the echelon basis of the submodule
          defining the quotient

        - ``names`` -- the name of the vectors of the basis of
          the submodule, or ``None``

        TESTS::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: Q = M.quo((X + t)*v)  # indirect doctest
            sage: type(Q)
            <class 'sage.modules.ore_module.OreQuotientModule_with_category'>

            sage: TestSuite(N).run()
        """
        from sage.modules.ore_module_morphism import OreModuleSection
        self._cover = cover
        d = cover.rank()
        base = cover.base_ring()
        self._relations = basis
        pivots = basis.pivots()
        r = basis.rank()
        coerce = matrix(base, d, d-r)
        indices = []
        i = 0
        for j in range(d):
            if i < r and pivots[i] == j:
                i += 1
            else:
                indices.append(j)
                coerce[j,j-i] = base.one()
        for i in range(r):
            for j in range(d-r):
                coerce[pivots[i],j] = -basis[i,indices[j]]
        rows = [cover.gen(i).image() * coerce for i in indices]
        OreModule.__init__(self, matrix(base, rows),
                           cover.ore_ring(action=False),
                           names, cover._ore_category)
        self._indices = indices
        self._project = coerce = cover.hom(coerce, codomain=self)
        self.register_coercion(coerce)
        cover.register_conversion(OreModuleSection(self, cover))

    def _repr_element(self, x):
        r"""
        Return a string representation of `x`.

        By default, elements in a Ore quotient are printed as
        their (canonical) representatives.

        TESTS::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: Q = M.quo((X + t)*v)
            sage: Q.an_element()  # indirect doctest
            w
        """
        M = self._cover
        indices = self._indices
        base = self.base_ring()
        coords = M.rank() * [base.zero()]
        for i in range(self.rank()):
            coords[indices[i]] = x[i]
        return M(coords)._repr_()

    def _latex_element(self, x):
        r"""
        Return a LaTeX representation of `x`.

        By default, elements in a Ore quotient are rendered as
        their (canonical) representatives with a bar.

        TESTS::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: Q = M.quo((X + t)*v)
            sage: latex(Q.an_element())  # indirect doctest
            \overline{w}
        """
        M = self._cover
        indices = self._indices
        base = self.base_ring()
        coords = M.rank() * [base.zero()]
        for i in range(self.rank()):
            coords[indices[i]] = x[i]
        return "\\overline{%s}" % M(coords)._latex_()

    def cover(self):
        r"""
        If this quotient in `M/N`, return `M`.

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: N = M.quo((X + t)*v)

            sage: N.cover()
            Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: N.cover() is M
            True

        .. SEEALSO::

            :meth:`relations`
        """
        return self._cover

    def relations(self, names=None):
        r"""
        If this quotient in `M/N`, return `N`.

        INPUT:

        - ``names`` -- the names of the vectors of the basis
          of `N`, or ``None``

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: Q = M.quo((X + t)*v)

            sage: N = Q.relations()
            sage: N
            Ore module of rank 1 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: (X + t)*v in N
            True
            sage: Q == M/N
            True

        It is also possible to define names for the basis elements
        of `N`::

            sage: N.<u> = Q.relations()
            sage: N
            Ore module <u> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: M(u)
            v + 1/t*w

        .. SEEALSO::

            :meth:`relations`
        """
        return self._submodule_class(self._cover, self._relations, names=names)

    def rename_basis(self, names, coerce=False):
        r"""
        Return the same Ore module with the given naming
        for the vectors in its distinguished basis.

        INPUT:

        - ``names`` -- a string or a list of strings, the
          new names

        - ``coerce`` (default: ``False``) -- a boolean; if
          ``True``, a coercion map from this Ore module to
          renamed version is set

        EXAMPLES::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: M = S.quotient_module(X^2 + z*X + 1)
            sage: M
            Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5

            sage: Me = M.rename_basis('e')
            sage: Me
            Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5

        Now compare how elements are displayed::

            sage: M.random_element()   # random
            (3*z^2 + 4*z + 2, 3*z^2 + z)
            sage: Me.random_element()  # random
            (2*z + 4)*e0 + (z^2 + 4*z + 4)*e1

        At this point, there is no coercion map between ``M``
        and ``Me``. Therefore, adding elements in both parents
        results in an error:

            sage: M.random_element() + Me.random_element()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +:
                'Ore module of rank 2 over Finite Field in z of size 5^3 twisted by z |--> z^5' and
                'Ore module <e0, e1> over Finite Field in z of size 5^3 twisted by z |--> z^5'

        In order to set this coercion, one should define ``Me``
        by passing the extra argument ``coerce=True``::

            sage: Me = M.rename_basis('e', coerce=True)
            sage: M.random_element() + Me.random_element()  # random
            2*z^2*e0 + (z^2 + z + 4)*e1

        .. WARNING::

            Use ``coerce=True`` with extreme caution. Indeed,
            setting inappropriate coercion maps may result in a
            circular path in the coercion graph which, in turn,
            could eventually break the coercion system.

        Note that the bracket construction also works::

            sage: M.<v,w> = M.rename_basis()
            sage: M
            Ore module <v, w> over Finite Field in z of size 5^3 twisted by z |--> z^5

        In this case, `v` and `w` are automatically defined::

            sage: v + w
            v + w

        TESTS::

            sage: K.<z> = GF(5^3)
            sage: S.<X> = OrePolynomialRing(K, K.frobenius_endomorphism())
            sage: P = X + z
            sage: A.<v,w> = S.quotient_module(P^2)
            sage: M = A.quo(P*v)
            sage: Me = M.rename_basis('e', coerce=True)
            sage: M.an_element() + Me.an_element()
            2*e0

        """
        rank = self.rank()
        names = normalize_names(names, rank)
        cls = self.__class__
        M = cls.__classcall__(cls, self._cover, self._relations, names)
        if coerce:
            mat = identity_matrix(self.base_ring(), rank)
            id = self.hom(mat, codomain=M)
            M._unset_coercions_used()
            M.register_coercion(id)
        return M

    def projection_morphism(self):
        r"""
        Return the projection from the cover module to this quotient.

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((X + t)^2)
            sage: Q = M.quo((X + t)*v)
            sage: Q.projection_morphism()
            Ore module morphism:
              From: Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module of rank 1 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
        """
        return self._project

    def morphism_quotient(self, f):
        r"""
        If this quotient in `M/N` and `f : M \to X` is a morphism
        vanishing on `N`, return the induced map `M/N \to X`.

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: P = X + t
            sage: M.<v,w> = S.quotient_module(P^2)
            sage: Q.<wbar> = M.quo(P*v)

            sage: f = M.hom({v: P*v})
            sage: f
            Ore module endomorphism of Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: g = Q.morphism_quotient(f)
            sage: g
            Ore module morphism:
              From: Ore module <wbar> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt

        When the given morphism does not vanish on `N`, an error is raised::

            sage: h = M.multiplication_map(X^5)
            sage: h
            Ore module endomorphism of Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: Q.morphism_quotient(h)
            Traceback (most recent call last):
            ...
            ValueError: the morphism does not factor through this quotient
        """
        if f.domain() is not self._cover:
            raise ValueError
        Z = self._relations * f._matrix
        if not Z.is_zero():
            raise ValueError("the morphism does not factor through this quotient")
        mat = f._matrix.matrix_from_rows(self._indices)
        return self.hom(mat, codomain=f.codomain())

    def morphism_modulo(self, f):
        r"""
        If this quotient in `M/N` and `f : X \to M` is a morphism,
        return the induced map `X \to M/N`.

        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: P = X + t
            sage: M.<v,w> = S.quotient_module(P^2)
            sage: Q.<wbar> = M.quo(P*v)

            sage: f = M.multiplication_map(X^5)
            sage: f
            Ore module endomorphism of Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
            sage: g = Q.morphism_modulo(f)
            sage: g
            Ore module morphism:
              From: Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
              To:   Ore module <wbar> over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 twisted by d/dt
        """
        if f.codomain() is not self._cover:
            raise ValueError
        return self._project * f

    _hom_change_domain = morphism_quotient
    _hom_change_codomain = morphism_modulo
