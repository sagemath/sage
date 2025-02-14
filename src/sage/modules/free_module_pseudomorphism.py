"""
Pseudomorphisms of free modules

AUTHORS:

- Xavier Caruso, Yossef Musleh (2024-09): initial version
"""
# ****************************************************************************
#  Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#                     Yossef Musleh <specialholonomy@gmail.com>
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
# ****************************************************************************

from sage.categories.morphism import Morphism
from sage.structure.richcmp import rich_to_bool, richcmp
from sage.modules.free_module_morphism import FreeModuleMorphism


class FreeModulePseudoMorphism(Morphism):
    r"""
    Let `M, M'` be modules over a ring `R`, `\theta: R \to R` a
    ring homomorphism, and `\delta: R \to R` a `\theta`-derivation,
    which is a map such that:

    .. MATH::

        \delta(xy) = \theta(x)\delta(y) + \delta(x)y.

    A pseudomorphism `f : M \to M` is an additive map such that

    .. MATH::

        f(\lambda x) = \theta(\lambda)f(x) + \delta(\lambda) x

    The map `\theta` (resp. `\delta`) is referred to as the
    twisting endomorphism (resp. the twisting derivation) of `f`.

    .. NOTE::

        The implementation currently requires that `M` and `M'`
        are free modules.

    We represent pseudomorphisms by matrices with coefficient in the
    base ring `R`. The matrix `\mathcal M_f` representing `f` is such
    that its lines (resp. columns if ``side`` is ``"right"``) are the
    coordinates of the images of the distinguished basis of the domain
    (see also method :meth:`matrix`).
    More concretely, let `n` (resp. `n'`) be the dimension of `M`
    (resp. `M'`), let `(e_1, \dots, e_n)` be a basis of `M`.
    For any `x = \sum_{i=1}^n x_i e_i \in M`, we have

    .. MATH::

        f(x) = \begin{pmatrix}
                 \theta(x_1) & \cdots & \theta(x_n)
               \end{pmatrix}
               \mathcal M_f
               +
               \begin{pmatrix}
                 \delta(x_1) & \cdots & \theta(x_n)
               \end{pmatrix}
               .

    When ``side`` is ``"right"``, the formula is

    .. MATH::

        f(x) = \mathcal M_f
               \begin{pmatrix}
                 \theta(x_1) \\ \vdots \\ \theta(x_n)
               \end{pmatrix}
               +
               \begin{pmatrix}
                 \delta(x_1) \\ \vdots \\ \theta(x_n)
               \end{pmatrix}
               .

    This class is not supposed to be instantiated directly; the user
    should use instead the method
    :meth:`sage.rings.module.free_module.FreeModule_generic.pseudohom`
    to create a pseudomorphism.

    TESTS::

        sage: P.<x> = ZZ[]
        sage: d = P.derivation()
        sage: M = P^2
        sage: f = M.pseudohom([[1, 2*x], [x, 1]], d); f
        Free module pseudomorphism (twisted by d/dx) defined by the matrix
        [  1 2*x]
        [  x   1]
        Domain: Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring
        Codomain: Ambient free module of rank 2 over the integral domain Univariate Polynomial Ring in x over Integer Ring
        sage: e = M((2*x^2 + 3*x + 1, x^3 + 7*x + 4))
        sage: f(e)
        (x^4 + 9*x^2 + 11*x + 4, 5*x^3 + 9*x^2 + 9*x + 11)
        sage: f = M.pseudohom([[1, 2], [1, 1]], d)
        sage: f(e)
        (x^3 + 2*x^2 + 14*x + 8, x^3 + 7*x^2 + 13*x + 13)

    ::

        sage: Fq.<z> = GF(7^3)
        sage: Frob = Fq.frobenius_endomorphism()
        sage: M = Fq^3
        sage: N = Fq^2
        sage: phi = M.pseudohom([[2, 3, 1], [1, 4, 6]], Frob, codomain=N, side="right")
        sage: phi
        Free module pseudomorphism (twisted by z |--> z^7) defined as left-multiplication by the matrix
        [2 3 1]
        [1 4 6]
        Domain: Vector space of dimension 3 over Finite Field in z of size 7^3
        Codomain: Vector space of dimension 2 over Finite Field in z of size 7^3
        sage: v = (4*z^2 + 4*z + 3, 2, z + 5)
        sage: phi(v)
        (2*z + 1, 6*z^2 + 4*z + 5)
    """
    def __init__(self, parent, f, side):
        """
        Constructs a pseudomorphism of free modules.

        INPUT:

        - ``parent`` -- the parent space of pseudomorphisms

        - ``f`` -- a pseudomorphism or a matrix defining this
          pseudomorphism

        - ``side`` -- side of the vectors acted on by the matrix

        TESTS::

            sage: F.<z> = GF(5^3)
            sage: Frob = F.frobenius_endomorphism()
            sage: M = F^2
            sage: H = M.pseudoHom(Frob)
            sage: H
            Set of Pseudoendomorphisms (twisted by z |--> z^5) of Vector space of dimension 2 over Finite Field in z of size 5^3

        The attribute ``f`` can be a matrix::

            sage: mat = matrix(F, 2, [1, z, z^2, z^3])
            sage: f = H(mat)
            sage: f
            Free module pseudomorphism (twisted by z |--> z^5) defined by the matrix
            [      1       z]
            [    z^2 2*z + 2]
            Domain: Vector space of dimension 2 over Finite Field in z of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 5^3

            sage: type(f)
            <class 'sage.modules.free_module_pseudohomspace.FreeModulePseudoHomspace_with_category.element_class'>

        or a pseudomorphism with the same parent::

            sage: H(f)
            Free module pseudomorphism (twisted by z |--> z^5) defined by the matrix
            [      1       z]
            [    z^2 2*z + 2]
            Domain: Vector space of dimension 2 over Finite Field in z of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 5^3

        When the twisting morphism and the twisting derivation are both trivial,
        pseudomorphisms are just linear applications and coercion between those
        works::

            sage: id = End(F).identity()
            sage: g = M.hom(mat)
            sage: g2 = M.pseudoHom(id)(g)
            sage: g2
            Free module pseudomorphism (untwisted) defined by the matrix
            [      1       z]
            [    z^2 2*z + 2]
            Domain: Vector space of dimension 2 over Finite Field in z of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 5^3

        An example with ``side=right``::

            sage: h = M.pseudohom(mat, Frob, side="right")
            sage: h
            Free module pseudomorphism (twisted by z |--> z^5) defined as left-multiplication by the matrix
            [      1       z]
            [    z^2 2*z + 2]
            Domain: Vector space of dimension 2 over Finite Field in z of size 5^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 5^3

        ::

            sage: M.pseudohom(mat, Frob, side="middle")
            Traceback (most recent call last):
            ...
            ValueError: the side must be either 'left' or 'right'

        ::

            sage: TestSuite(f).run()
            sage: TestSuite(g2).run()
            sage: TestSuite(h).run()
        """
        Morphism.__init__(self, parent)
        dom = parent.domain()
        codom = parent.codomain()
        if side != "left" and side != "right":
            raise ValueError("the side must be either 'left' or 'right'")
        matrix_space = parent.matrix_space()
        if ((isinstance(f, FreeModulePseudoMorphism) and f.parent() is parent)
         or (isinstance(f, FreeModuleMorphism)
          and f.domain() is dom and f.codomain() is codom
          and parent._morphism is None and parent._derivation is None)):
            if f.side() == 'right':
                self._matrix = f.matrix().transpose()
            else:
                self._matrix = f.matrix()
        else:
            if side == "right":
                self._matrix = matrix_space.transposed(f).transpose()
            else:
                self._matrix = matrix_space(f)
        self._morphism = parent._morphism
        self._derivation = parent._derivation
        self._side = side

    def _call_(self, x):
        r"""
        Return the result of applying this pseudomorphism to ``x``.

        TESTS::

            sage: Fq.<z> = GF(7^3)
            sage: M = Fq^3
            sage: Frob = Fq.frobenius_endomorphism()
            sage: f = M.pseudohom([[1, z, 3], [0, 1, 1], [2, 1, 1]], Frob)
            sage: e = M((3*z^2 + 5*z + 2, 6*z^2 + 2*z + 2, z + 4))
            sage: f(e)
            (3*z^2 + 4*z + 4, 6*z^2 + 5*z + 6, 6*z^2 + 5*z + 3)

        ::

            sage: g = M.pseudohom([[1, z, 3], [0, 1, 1], [2, 1, 1]], Frob, side="right")
            sage: g(e)
            (z^2 + 6*z + 2, z^2 + 2*z + 1, 2*z^2 + 4*z)
        """
        D = self.domain()
        C = self.codomain()
        if D.is_ambient():
            x = x.element()
        else:
            x = D.coordinate_vector(x)
        if self._morphism is None:
            x_twist = x
        else:
            x_twist = D(list(map(self._morphism, x)))
        v = x_twist * self._matrix
        if C.is_ambient():
            v = C(v.list())
        else:
            v = C.linear_combination_of_basis(v)
        if self._derivation is not None:
            v += D(list(map(self._derivation, x)))
        return v

    def _repr_(self):
        r"""
        Return the string representation of a pseudomorphism.

        TESTS::

            sage: Fq.<z> = GF(7^3)
            sage: M = Fq^3
            sage: Frob = Fq.frobenius_endomorphism()

            sage: f = M.pseudohom([[1,1,1], [2,2,2], [3,3,3]], Frob)
            sage: f  # indirect doctest
            Free module pseudomorphism (twisted by z |--> z^7) defined by the matrix
            [1 1 1]
            [2 2 2]
            [3 3 3]
            Domain: Vector space of dimension 3 over Finite Field in z of size 7^3
            Codomain: Vector space of dimension 3 over Finite Field in z of size 7^3

            sage: g = M.pseudohom([[1,1,1], [2,2,2], [3,3,3]], Frob, side="right")
            sage: g  # indirect doctest
            Free module pseudomorphism (twisted by z |--> z^7) defined as left-multiplication by the matrix
            [1 1 1]
            [2 2 2]
            [3 3 3]
            Domain: Vector space of dimension 3 over Finite Field in z of size 7^3
            Codomain: Vector space of dimension 3 over Finite Field in z of size 7^3
        """
        twist = self.parent()._ore._repr_twist()
        s = "Free module pseudomorphism (%s) defined " % twist
        if self._side == "right":
            s += "as left-multiplication "
        s += "by the matrix\n%s\n" % self.matrix()
        s += "Domain: %s\n" % self.domain()
        s += "Codomain: %s" % self.codomain()
        return s

    def matrix(self):
        r"""
        Return the underlying matrix of this pseudomorphism.

        It is defined as the matrix `M` whose lines (resp. columns if
        ``side`` is ``"right"``) are the coordinates of the images of
        the distinguished basis of the domain.

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: M = Fq^3
            sage: f = M.pseudohom([[1, z, 3], [0, 1, z^2], [z+1, 1, 1]], Frob)
            sage: f.matrix()
            [    1     z     3]
            [    0     1   z^2]
            [z + 1     1     1]

        ::

            sage: e1, e2, e3 = M.basis()
            sage: f(e1)
            (1, z, 3)
            sage: f(e2)
            (0, 1, z^2)
            sage: f(e3)
            (z + 1, 1, 1)

        TESTS::

            sage: v = M.random_element()
            sage: f(v) == vector([Frob(c) for c in v]) * f.matrix()
            True
        """
        if self._side == "left":
            return self._matrix.__copy__()
        else:
            return self._matrix.transpose()

    def twisting_derivation(self):
        r"""
        Return the twisting derivation of the pseudomorphism
        (or ``None`` if the twisting derivation is zero).

        EXAMPLES::

            sage: P.<x> = ZZ[]
            sage: d = P.derivation()
            sage: M = P^2
            sage: f = M.pseudohom([[1, 2*x], [x, 1]], d)
            sage: f.twisting_derivation()
            d/dx

        ::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: f = V.pseudohom([[1, z], [0, z^2]], Frob)
            sage: f.twisting_derivation()
        """
        return self._derivation

    def twisting_morphism(self):
        r"""
        Return the twisting morphism of the pseudomorphism
        (or ``None`` if the twisting morphism is the identity).

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: f = V.pseudohom([[1, z], [0, z^2]], Frob)
            sage: f.twisting_morphism()
            Frobenius endomorphism z |--> z^7 on Finite Field in z of size 7^3

        ::

            sage: P.<x> = ZZ[]
            sage: d = P.derivation()
            sage: M = P^2
            sage: f = M.pseudohom([[1, 2*x], [x, 1]], d)
            sage: f.twisting_morphism()
        """
        return self._morphism

    def side(self):
        """
        Return the side of vectors acted on, relative to the matrix.

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2

            sage: m = matrix(2, [1, z, z^2, z^3])
            sage: h1 = V.pseudohom(m, Frob)
            sage: h1.side()
            'left'
            sage: h1([1, 0])
            (1, z)

            sage: h2 = V.pseudohom(m, Frob, side="right")
            sage: h2.side()
            'right'
            sage: h2([1, 0])
            (1, z^2)
        """
        return self._side

    def side_switch(self):
        """
        Return the same morphism, acting on vectors on the opposite side.

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2

            sage: m = matrix(2, [1, z, z^2, z^3])
            sage: h1 = V.pseudohom(m, Frob)
            sage: h1
            Free module pseudomorphism (twisted by z |--> z^7) defined by the matrix
            [      1       z]
            [    z^2 z^2 + 3]
            Domain: Vector space of dimension 2 over Finite Field in z of size 7^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 7^3

            sage: h2 = h1.side_switch()
            sage: h2
            Free module pseudomorphism (twisted by z |--> z^7) defined as left-multiplication by the matrix
            [      1     z^2]
            [      z z^2 + 3]
            Domain: Vector space of dimension 2 over Finite Field in z of size 7^3
            Codomain: Vector space of dimension 2 over Finite Field in z of size 7^3

        We check that ``h1`` and ``h2`` are the same::

            sage: v = V.random_element()
            sage: h1(v) == h2(v)
            True
        """
        if self._side == "left":
            side = "right"
            mat = self._matrix.transpose()
        else:
            side = "left"
            mat = self._matrix
        return self.parent()(mat, side)

    def __nonzero__(self):
        return not (self._derivation is None and self._matrix)

    def __eq__(self, other):
        r"""
        Compare this morphism with ``other``.

        TESTS::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: m = random_matrix(Fq, 2)

            sage: f = V.pseudohom(m, Frob)
            sage: g = V.pseudohom(m.transpose(), Frob, side="right")
            sage: f == g
            True

            sage: g = V.pseudohom(m.transpose(), Frob)
            sage: f == g
            False

            sage: g = V.pseudohom(m, Frob^2)
            sage: f == g
            False

            sage: g = V.pseudohom(m, Frob^3)
            sage: h = V.hom(m)
            sage: g == h
            True
        """
        if isinstance(other, FreeModuleMorphism):
            try:
                other = self.parent()(other)
            except ValueError:
                return False
        if isinstance(other, FreeModulePseudoMorphism):
            return self.parent() is other.parent() and self._matrix == other._matrix

    def ore_module(self, names=None):
        r"""
        Return the Ore module over which the Ore variable acts
        through this pseudomorphism.

        INPUT:

        - ``names`` -- a string, a list of strings or ``None``,
          the names of the vector of the canonical basis of the
          Ore module; if ``None``, elements are represented as
          vectors in `K^d` (where `K` is the base ring)

        EXAMPLES::

            sage: Fq.<z> = GF(7^3)
            sage: Frob = Fq.frobenius_endomorphism()
            sage: V = Fq^2
            sage: mat = matrix(2, [1, z, z^2, z^3])
            sage: f = V.pseudohom(mat, Frob)

            sage: M = f.ore_module()
            sage: M
            Ore module of rank 2 over Finite Field in z of size 7^3 twisted by z |--> z^7

        Here `M` is a module over the Ore ring `\FF_q[X; \\text{Frob}]`
        and the variable `X` acts on `M` through `f`::

            sage: S.<X> = M.ore_ring()
            sage: S
            Ore Polynomial Ring in X over Finite Field in z of size 7^3 twisted by z |--> z^7
            sage: v = M((1,0))
            sage: X*v
            (1, z)

        The argument ``names`` can be used to give chosen names
        to the vectors in the canonical basis::

            sage: M = f.ore_module(names=('v', 'w'))
            sage: M.basis()
            [v, w]

        or even::

            sage: M = f.ore_module(names='e')
            sage: M.basis()
            [e0, e1]

        Note that the bracket construction also works::

            sage: M.<v,w> = f.ore_module()
            sage: M.basis()
            [v, w]
            sage: v + w
            v + w

        We refer to :mod:`sage.modules.ore_module` for a
        tutorial on Ore modules in SageMath.

        .. SEEALSO::

            :mod:`sage.modules.ore_module`
        """
        from sage.modules.ore_module import OreModule
        return OreModule(self._matrix, self.parent()._ore, names)

    def _test_nonzero_equal(self, tester):
        pass
