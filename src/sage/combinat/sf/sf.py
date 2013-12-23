"""
Symmetric functions, with their multiple realizations
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2009-2012 Jason Bandlow <jbandlow@gmail.com>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
#                     2009-2012 Nicolas M. Thiery <nthiery at users.sf.net>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.all import Rings, GradedHopfAlgebras
from sage.combinat.partition import Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.rational_field import QQ

import schur
import monomial
import powersum
import elementary
import homogeneous
import hall_littlewood
import jack
import macdonald
import llt

class SymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The abstract algebra of commutative symmetric functions

    .. rubric:: Symmetric Functions in Sage

    .. MODULEAUTHOR:: Jason Bandlow, Anne Schilling, Nicolas M. Thiery, Mike Zabrocki

    This document is an introduction to working with symmetric function
    theory in Sage.
    It is not intended to be an introduction to the theory
    of symmetric functions ([MAC]_ and [STA]_, Chapter 7, are two excellent
    references.)  The reader is also expected to be familiar with Sage.

    .. rubric:: The algebra of symmetric functions

    The algebra of symmetric functions is the unique free commutative graded
    connected algebra over the given ring, with one generator in each degree.  It
    can also be thought of as the inverse limit (in the category of graded
    algebras) of the algebra of symmetric polynomials in `n` variables as `n \rightarrow \infty`.
    Sage allows us to construct the algebra of symmetric functions over
    any ring.  We will use a base ring of rational numbers in these first
    examples::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym
        Symmetric Functions over Rational Field

    Sage knows certain categorical information about this algebra::

        sage: Sym.category()
        Join of Category of graded hopf algebras over Rational Field and Category of monoids with realizations and Category of coalgebras over Rational Field with realizations

    Notice that ``Sym`` is an *abstract* algebra.  This reflects the fact that
    there are multiple natural bases.  To work with specific
    elements, we need a *realization* of this algebra.  In practice, this
    means we need to specify a basis.

    .. rubric:: An example basis - power sums

    Here is an example of how one might use the power sum realization::

        sage: p = Sym.powersum()
        sage: p
        Symmetric Functions over Rational Field in the powersum basis

    ``p`` now represents the realization of the symmetric function algebra on
    the power sum basis.  The basis itself is accessible through::

        sage: p.basis()
        Lazy family (Term map from Partitions to Symmetric Functions over Rational Field in the powersum basis(i))_{i in Partitions}
        sage: p.basis().keys()
        Partitions

    This last line means that ``p.basis()`` is an association between the set
    of Partitions and the basis elements of the algebra ``p``. To construct a
    specific element one can therefore do::

        sage: p.basis()[Partition([2,1,1])]
        p[2, 1, 1]

    As this is rather cumbersome, realizations of the symmetric function
    algebra allow for the following abuses of notation::

        sage: p[Partition([2, 1, 1])]
        p[2, 1, 1]
        sage: p[[2, 1, 1]]
        p[2, 1, 1]
        sage: p[2, 1, 1]
        p[2, 1, 1]

    or even::

        sage: p[(i for i in [2, 1, 1])]
        p[2, 1, 1]

    In the special case of the empty partition, due to a limitation in
    Python syntax, one cannot use::

        sage: p[]       # todo: not implemented

    Please use instead::

        sage: p[[]]
        p[]

    .. note:: When elements are constructed using the ``p[something ]`` syntax ,
       an error will be raised if the input cannot be interpreted as a partition.
       This is *not* the case when ``p.basis()`` is used::

        sage: p['something']
        Traceback (most recent call last):
        ...
        ValueError: ['s', 'o', 'm', 'e', 't', 'h', 'i', 'n', 'g'] is not an element of Partitions
        sage: p.basis()['something']
        p'something'

    Elements of ``p`` are linear combinations of such compositions::

      sage: p.an_element()
      2*p[] + 2*p[1] + 3*p[2]

    .. rubric:: Algebra structure

    Algebraic combinations of basis elements can be entered in a natural way::

        sage: p[2,1,1] + 2 * p[1] * (p[4] + p[2,1])
        3*p[2, 1, 1] + 2*p[4, 1]

    Let us explore the other operations of ``p``. We can ask for
    the mathematical properties of ``p``::

        sage: p.categories()
        [Category of bases of Symmetric Functions over Rational Field, Category of graded hopf algebras with basis over Rational Field, ...]

    To start with, ``p`` is a graded algebra, the grading being induced
    by the size of the partitions. Due to this, the one is the basis
    element indexed by the empty partition::

        sage: p.one()
        p[]

    The ``p`` basis is multiplicative; that is, multiplication is induced by
    linearity from the (nonincreasingly sorted) concatenation of partitions::

        sage: p[3,1] * p[2,1]
        p[3, 2, 1, 1]

        sage: (p.one() + 2 * p[3,1]) * p[4, 2]
        p[4, 2] + 2*p[4, 3, 2, 1]

    .. rubric:: The classical bases

    In addition to the power sum basis, the other classical bases of the
    symmetric function algebra are the elementary, complete homogeneous,
    monomial, and Schur bases.  These can be defined as follows::

        sage: e = Sym.elementary()
        sage: h = Sym.homogeneous()
        sage: m = Sym.monomial()
        sage: s = Sym.schur()

    These can be defined all at once with the single command
    ::

        sage: Sym.inject_shorthands()
        doctest:...: RuntimeWarning: redefining global value `h`
        doctest:...: RuntimeWarning: redefining global value `s`
        doctest:...: RuntimeWarning: redefining global value `e`
        doctest:...: RuntimeWarning: redefining global value `m`
        doctest:...: RuntimeWarning: redefining global value `p`

    We can then do conversions from one basis to another::

        sage: s(p[2,1])
        -s[1, 1, 1] + s[3]

        sage: m(p[3])
        m[3]
        sage: m(p[3,2])
        m[3, 2] + m[5]

    For computations which mix bases, Sage will return a result with respect
    to a single (not necessarily predictable) basis::

        sage: p[2] * s[2] - m[4]
        1/2*p[2, 1, 1] + 1/2*p[2, 2] - p[4]

        sage: p( m[1] * ( e[3]*s[2] + 1 ))
        p[1] + 1/12*p[1, 1, 1, 1, 1, 1] - 1/6*p[2, 1, 1, 1, 1] - 1/4*p[2, 2, 1, 1] + 1/6*p[3, 1, 1, 1] + 1/6*p[3, 2, 1]


    The one for different bases such as the power sum and Schur function is the same::

        sage: s.one() == p.one()
        True

    .. rubric:: Basic computations

    In this section, we explore some of the many methods that can be applied
    to an arbitrary symmetric function::

        sage: f = s[2]^2; f
        s[2, 2] + s[3, 1] + s[4]

    For more methods than discussed here, create a symmetric function as
    above, and use ``f.<tab>``.

    .. _`Representation theory of the symmetric group`:

    .. rubric:: Representation theory of the symmetric group

    The Schur functions `s_\lambda` can also be interpreted as irreducible characters of the symmetric
    group `S_n`, where `n` is the size of the partition `\lambda`. Since the Schur functions of
    degree `n` form a basis of the symmetric functions of degree `n`, it
    follows that an arbitrary symmetric function (homogeneous of degree
    `n`) may be interpreted as a function on the symmetric group. In this
    interpretation the power sum symmetric function `p_\lambda` is the characteristic
    function of the conjugacy class with shape `\lambda`, multiplied by the order of
    the centralizer of an element. Hence the irreducible characters can be computed
    as follows::

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.schur()
        sage: p = Sym.power()
        sage: P = Partitions(5).list()
        sage: P = [P[i] for i in range(len(P)-1,-1,-1)]
        sage: M = matrix([[s[P[i]].scalar(p[P[j]]) for j in range(len(P))] for i in range(len(P))])
        sage: M
        [ 1 -1  1  1 -1 -1  1]
        [ 4 -2  0  1  1  0 -1]
        [ 5 -1  1 -1 -1  1  0]
        [ 6  0 -2  0  0  0  1]
        [ 5  1  1 -1  1 -1  0]
        [ 4  2  0  1 -1  0 -1]
        [ 1  1  1  1  1  1  1]

    We can indeed check that this agrees with the character table of `S_5`::

        sage: SymmetricGroup(5).character_table() == M
        True

    In this interpretation of symmetric functions as characters on the
    symmetric group, the multiplication and comultiplication are
    interpreted as induction (from `S_n\times S_m` to `S_{n+m}`)
    and restriction, respectively. The Schur functions can also be interpreted
    as characters of `GL_n`, see `Partitions and Schur functions`__.

    __ ../../../../thematic_tutorials/lie/lie_basics.html#partitions-and-schur-polynomials

    .. rubric:: The omega involution

    The `\omega` involution is the linear extension of the map which sends
    `e_\lambda` to `h_{\lambda}`::

        sage: h(f)
        h[2, 2]
        sage: e(f.omega())
        e[2, 2]

    .. rubric:: The Hall scalar product

    The Hall scalar product on the algebra of symmetric functions makes the
    Schur functions into an orthonormal basis::

        sage: f.scalar(f)
        3

    .. rubric:: Skewing

    *Skewing* is the adjoint operation to multiplication with respect to
    this scalar product::

        sage: f.skew_by(s[1])
        2*s[2, 1] + 2*s[3]

    In general, ``s[la].skew_by(s[mu])`` is the symmetric function typically
    denoted `s_{\lambda \setminus \mu}` or `s_{\lambda / \mu}`.

    .. rubric:: Expanding into variables

    We can expand a symmetric function into a symmetric polynomial in a
    specified number of variables::

        sage: f.expand(2)
        x0^4 + 2*x0^3*x1 + 3*x0^2*x1^2 + 2*x0*x1^3 + x1^4

    See the documentation for ``expand`` for more examples.

    .. rubric:: The Kronecker product

    As in the section on the `Representation theory of
    the symmetric group`_, a symmetric function may be considered as a
    class function on the symmetric group where the elements
    `p_\mu/z_\mu` are the indicators of a permutation having
    cycle structure `\mu`.  The Kronecker product of two
    symmetric functions corresponds to the pointwise product
    of these class functions.

    Since the Schur functions are the irreducible characters
    of the symmetric group under this identification, the Kronecker
    product of two Schur functions corresponds to the internal
    tensor product of two irreducible symmetric group representations.

    Under this identification, the Kronecker
    product of `p_\mu/z_\mu` and `p_\nu/z_\nu` is `p_\mu/z_\mu`
    if `\mu=\nu`, and the result is equal to `0` otherwise.

    ``internal_product``, ``kronecker_product``, ``inner_tensor`` and
    ``itensor`` are different names for the same function.

    ::

        sage: f.kronecker_product(f)
        s[1, 1, 1, 1] + 3*s[2, 1, 1] + 4*s[2, 2] + 5*s[3, 1] + 3*s[4]

    .. rubric:: Plethysm

    The *plethysm* of symmetric functions is the operation corresponding to
    composition of representations of the general linear group.  See [STA]_
    Chapter 7, Appendix 2 for details.
    ::

        sage: s[2].plethysm(s[2])
        s[2, 2] + s[4]

    Plethysm can also be written as a composition of functions::

        sage: s[2]( s[2] )
        s[2, 2] + s[4]

    If the coefficient ring contains degree 1 elements, these are handled
    properly by plethysm::

        sage: R.<t> = QQ[]; s = SymmetricFunctions(R).schur()
        sage: s[2]( (1-t)*s[1] )
        (t^2-t)*s[1, 1] + (-t+1)*s[2]

    See the documentation for ``plethysm`` for more information.

    .. rubric:: Inner plethysm

    The operation of inner plethysm ``f.inner_plethysm(g)`` models the
    composition of the `S_n` representation represented by `g` with the
    `GL_m` representation whose character is `f`.  See the documentation of
    ``inner_plethysm``, [ST94]_ or [STA]_, exercise 7.74 solutions for more
    information::

        sage: s = SymmetricFunctions(QQ).schur()
        sage: f = s[2]^2
        sage: f.inner_plethysm(s[2])
        s[2]

    .. rubric:: Hopf algebra structure

    The ring of symmetric functions is further endowed with a coalgebra
    structure. The coproduct is an algebra morphism, and therefore
    determined by its values on the generators; the power sum generators
    are primitive::

        sage: p[1].coproduct()
        p[] # p[1] + p[1] # p[]
        sage: p[2].coproduct()
        p[] # p[2] + p[2] # p[]

    The coproduct, being cocommutative on the generators, is cocommutative everywhere::

        sage: p[2, 1].coproduct()
        p[] # p[2, 1] + p[1] # p[2] + p[2] # p[1] + p[2, 1] # p[]

    This coproduct, along with the counit which sends every symmetric function
    to its `0`-th homogeneous component, makes the ring of symmetric functions
    into a graded connected bialgebra. It is known that every graded connected
    bialgebra has an antipode. For the ring of symmetric functions, the antipode
    can be characterized explicitly: The antipode is an anti-algebra morphism
    (thus an algebra morphism, since our algebra is commutative) which sends
    `p_{\lambda}` to `(-1)^{\mathrm{length}(\lambda)} p_{\lambda}` for every
    partition `\lambda`. Thus, in particular, it sends the generators on the
    ``p`` basis to their opposites::

        sage: p[3].antipode()
        -p[3]
        sage: p[3,2,1].antipode()
        -p[3, 2, 1]

    The graded connected bialgebra of symmetric functions over a `\QQ`-algebra
    has a rather simply-understood structure: It is (isomorphic to) the
    symmetric algebra of its space of primitives (which is spanned by the
    power-sum symmetric functions).

    Here are further examples::

        sage: f = s[2]^2
        sage: f.antipode()
        s[1, 1, 1, 1] + s[2, 1, 1] + s[2, 2]
        sage: f.coproduct()
        s[] # s[2, 2] + s[] # s[3, 1] + s[] # s[4] + 2*s[1] # s[2, 1] + 2*s[1] # s[3] + s[1, 1] # s[1, 1]
        + s[1, 1] # s[2] + s[2] # s[1, 1] + 3*s[2] # s[2] + 2*s[2, 1] # s[1] + s[2, 2] # s[] + 2*s[3] # s[1]
        + s[3, 1] # s[] + s[4] # s[]
        sage: f.coproduct().apply_multilinear_morphism( lambda x,y: x*y.antipode() )
        0

    .. rubric:: Transformations of symmetric functions

    There are many methods in Sage which make it easy to manipulate symmetric
    functions.  For example, if we have some function which acts on partitions
    (say, conjugation), it is a simple matter to apply it to the support of a
    symmetric function.  Here is an example::

        sage: conj = lambda mu: mu.conjugate()
        sage: f = h[4] + 2*h[3,1]
        sage: f.map_support(conj)
        h[1, 1, 1, 1] + 2*h[2, 1, 1]

    We can also easily modify the coefficients::

        sage: def foo(mu, coeff): return mu.conjugate(), -coeff
        sage: f.map_item(foo)
        -h[1, 1, 1, 1] - 2*h[2, 1, 1]

    See also ``map_coefficients``.

    There are also methods for building functions directly::

        sage: s.sum_of_monomials(mu for mu in Partitions(3))
        s[1, 1, 1] + s[2, 1] + s[3]
        sage: s.sum_of_monomials(Partitions(3))
        s[1, 1, 1] + s[2, 1] + s[3]
        sage: s.sum_of_terms( (mu, mu[0]) for mu in Partitions(3))
        s[1, 1, 1] + 2*s[2, 1] + 3*s[3]

    These are the preferred way to build elements within a program;
    the result will usually be faster than using :func:`sum`. It also
    guarantees that empty sums yields the zero of ``s`` (see also
    ``s.sum``).

    Note also that it is a good idea to use::

        sage: s.one()
        s[]
        sage: s.zero()
        0

    instead of ``s(1)`` and ``s(0)`` within programs where speed is important,
    in order to prevent unnecessary coercions.

    .. rubric:: Different base rings

    Depending on the base ring, the different realizations of the symmetric
    function algebra may not span the same space::

        sage: SZ = SymmetricFunctions(ZZ)
        sage: p = SZ.power(); s = SZ.schur()
        sage: p(s[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Because of this, some functions may not behave as expected when working over
    the integers, even though they make mathematical sense::

        sage: s[1,1,1].plethysm(s[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    It is possible to work over different base rings simultaneously::

        sage: s = SymmetricFunctions(QQ).schur()
        sage: p = SymmetricFunctions(QQ).power()
        sage: sz = SymmetricFunctions(ZZ).schur(); sz._prefix = 'sz'
        sage: pz = SymmetricFunctions(ZZ).power(); pz._prefix = 'pz'
        sage: p(sz[1,1,1])
        1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3]
        sage: sz( 1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3] )
        sz[1, 1, 1]

    As shown in this example, if you are working over multiple base rings
    simultaneously, it is a good idea to change the prefix in some cases, so that
    you can tell from the output which realization your result is in.

    Let us change the notation back for the remainder of this tutorial::

        sage: sz._prefix = 's'
        sage: pz._prefix = 'p'

    One can also use the Sage standard renaming idiom to get shorter outputs::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym.rename("Sym")
        sage: Sym
        Sym
        sage: Sym.rename()

    And we name it back::

        sage: Sym.rename("Symmetric Functions over Rational Field"); Sym
        Symmetric Functions over Rational Field

    .. rubric:: Other bases

    There are two additional basis of the symmetric functions which are not
    considered as classical bases:

    * forgotten basis
    * Witt basis

    The forgotten basis is the dual basis of the elementary symmetric
    functions basis with respect to the Hall scalar product. The Witt basis
    can be constructed by

    .. MATH::

        \prod_{d=1}^{\infty} (1 - w_d t^d)^{-1} = \sum_{n=0}^{\infty} h_n t^n

    where `t` is a formal variable.

    There are further bases of the ring of symmetric functions, in general over
    fields with parameters such as `q` and `t`:

    * Hall-Littlewood bases
    * Jack bases
    * Macdonald bases
    * `k`-Schur functions

    We briefly demonstrate how to access these bases. For more information, see
    the documentation of the individual bases.

    The *Jack polynomials* can be obtained as::

        sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
        sage: Jack = Sym.jack()
        sage: P = Jack.P(); J = Jack.J(); Q = Jack.Q()
        sage: J(P[2,1])
        (1/(t+2))*JackJ[2, 1]

    The parameter `t` can be specialized as follows::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Jack = Sym.jack(t = 1)
        sage: P = Jack.P(); J = Jack.J(); Q = Jack.Q()
        sage: J(P[2,1])
        1/3*JackJ[2, 1]

    Similarly one can access the Hall-Littlewood and Macdonald polynomials, etc::

        sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
        sage: Mcd = Sym.macdonald()
        sage: P = Mcd.P(); J = Mcd.J(); Q = Mcd.Q()
        sage: J(P[2,1])
        (1/(-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1))*McdJ[2, 1]

    .. rubric:: `k`-Schur functions

    The `k`-Schur functions live in the `k`-bounded subspace of the ring of
    symmetric functions. It is possible to compute in the `k`-bounded subspace
    directly::

        sage: Sym = SymmetricFunctions(QQ)
        sage: ks = Sym.kschur(3,1)
        sage: f = ks[2,1]*ks[2,1]; f
        ks3[2, 2, 1, 1] + ks3[2, 2, 2] + ks3[3, 1, 1, 1]

    or to lift to the ring of symmetric functions::

        sage: f.lift()
        s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

    However, it is not always possible to convert a symmetric function to the `k`-bounded subspace::

        sage: s = Sym.schur()
        sage: ks(s[2,1,1])
        Traceback (most recent call last):
        ...
        ValueError: s[2, 1, 1] is not in the image of Generic morphism:
          From: 3-bounded Symmetric Functions over Rational Field with t=1 in the 3-Schur basis also with t=1
          To:   Symmetric Functions over Rational Field in the Schur basis

    The `k`-Schur functions are more generally defined with a parameter `t` and they are
    a basis of the subspace spanned by the Hall-Littlewood ``Qp`` symmetric functions
    indexed by partitions whose first part is less than or equal to `k`::

        sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
        sage: SymS3 = Sym.kBoundedSubspace(3) # default t='t'
        sage: ks = SymS3.kschur()
        sage: Qp = Sym.hall_littlewood().Qp()
        sage: ks(Qp[2,1,1,1])
        ks3[2, 1, 1, 1] + (t^2+t)*ks3[2, 2, 1] + (t^3+t^2)*ks3[3, 1, 1] + t^4*ks3[3, 2]

    The subspace spanned by the `k`-Schur functions with a parameter `t` are not known
    to form a natural algebra.  However it is known that the product of a `k`-Schur
    function and an `\ell`-Schur function is in the linear span of the `k+\ell`-Schur
    functions::

        sage: ks(ks[2,1]*ks[1,1])
        Traceback (most recent call last):
        ...
        ValueError: s[2, 1, 1, 1] + s[2, 2, 1] + s[3, 1, 1] + s[3, 2] is not in the image of Generic morphism:
          From: 3-bounded Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the 3-Schur basis
          To:   Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Schur basis
        sage: ks[2,1]*ks[1,1]
        s[2, 1, 1, 1] + s[2, 2, 1] + s[3, 1, 1] + s[3, 2]
        sage: ks6 = Sym.kBoundedSubspace(6).kschur()
        sage: ks6(ks[3,1,1]*ks[3])
        ks6[3, 3, 1, 1] + ks6[4, 2, 1, 1] + (t+1)*ks6[4, 3, 1] + t*ks6[4, 4]
        + ks6[5, 1, 1, 1] + ks6[5, 2, 1] + t*ks6[5, 3] + ks6[6, 1, 1]

    .. rubric:: dual `k`-Schur functions

    The dual space to the subspace spanned by the `k`-Schur functions is most naturally
    realized as a quotient of the ring of symmetric functions by an ideal.  When `t=1`
    the ideal is generated by the monomial symmetric functions indexed by partitions
    whose first part is greater than `k`.::

        sage: Sym = SymmetricFunctions(QQ)
        sage: SymQ3 = Sym.kBoundedQuotient(3,t=1)
        sage: km = SymQ3.kmonomial()
        sage: km[2,1]*km[2,1]
        4*m3[2, 2, 1, 1] + 6*m3[2, 2, 2] + 2*m3[3, 2, 1] + 2*m3[3, 3]
        sage: F = SymQ3.affineSchur()
        sage: F[2,1]*F[2,1]
        2*F3[1, 1, 1, 1, 1, 1] + 4*F3[2, 1, 1, 1, 1] + 4*F3[2, 2, 1, 1] + 4*F3[2, 2, 2]
        + 2*F3[3, 1, 1, 1] + 4*F3[3, 2, 1] + 2*F3[3, 3]

    When `t` is not equal to `1`, the subspace spanned by the `k`-Schur functions is
    realized as a quotient of the ring of symmetric functions by the ideal generated by
    the Hall-Littlewood symmetric functions in the P basis indexed by partitions with
    first part greater than `k`.::

        sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
        sage: SymQ3 = Sym.kBoundedQuotient(3)
        sage: kHLP = SymQ3.kHallLittlewoodP()
        sage: kHLP[2,1]*kHLP[2,1]
        (t^2+2*t+1)*HLP3[2, 2, 1, 1] + (t^3+2*t^2+2*t+1)*HLP3[2, 2, 2]
        + (-t^4-t^3+t+1)*HLP3[3, 1, 1, 1] + (-t^2+t+2)*HLP3[3, 2, 1] + (t+1)*HLP3[3, 3]
        sage: HLP = Sym.hall_littlewood().P()
        sage: kHLP(HLP[3,1])
        HLP3[3, 1]
        sage: kHLP(HLP[4])
        0

    In this space, the basis which is dual to the `k`-Schur functions conjecturally
    expands positively in the `k`-bounded Hall-Littlewood functions and has positive
    structure coefficients.::

        sage: dks = SymQ3.dual_k_Schur()
        sage: kHLP(dks[2,2])
        (t^4+t^2)*HLP3[1, 1, 1, 1] + t*HLP3[2, 1, 1] + HLP3[2, 2]
        sage: dks[2,1]*dks[1,1]
        (t^2+t)*dks3[1, 1, 1, 1, 1] + (t+1)*dks3[2, 1, 1, 1] + (t+1)*dks3[2, 2, 1]
        + dks3[3, 1, 1] + dks3[3, 2]

    At `t=1` the `k`-bounded Hall-Littlewood basis is equal to the `k`-bounded monomial
    basis and the dual `k`-Schur elements are equal to the affine Schur basis.  The
    `k`-bounded monomial basis and affine Schur functions are faster and should be used
    instead of the `k`-bounded Hall-Littlewood P basis and dual `k`-Schur functions when
    `t=1`.::

        sage: SymQ3 = Sym.kBoundedQuotient(3,t=1)
        sage: dks = SymQ3.dual_k_Schur()
        sage: F = SymQ3.affineSchur()
        sage: F[3,1]==dks[3,1]
        True

    .. rubric:: Implementing new bases

    .. todo:: to be described

    .. rubric:: Acknowledgements

    The design is heavily inspired from the implementation of
    symmetric functions in MuPAD-Combinat (see [HT04]_ and [FD06]_).

    REFERENCES:

        .. [FD06] Francois Descouens, Making research on symmetric functions using MuPAD-Combinat.
                 In Andres Iglesias and Nobuki Takayama, editors, 2nd International Congress on Mathematical Software (ICMS'06),
                 volume 4151 of LNCS, pages 407-418, Castro Urdiales, Spain, September 2006. Springer-Verlag.
                 :arXiv:`0806.1873`

        .. [HT04] Florent Hivert and Nicolas M. Thiery,
                 MuPAD-Combinat, an open-source package for research in algebraic combinatorics.
                 Sem. Lothar. Combin., 51 :Art. B51z, 70 pp. (electronic), 2004.
                 http://mupad-combinat.sf.net/.

        .. [MAC] Ian Macdonald, Symmetric Functions and Orthogonal Polynomials,
                 Second edition. With contributions by A. Zelevinsky. Oxford Mathematical Monographs.
                 Oxford Science Publications. The Clarendon Press, Oxford University Press, New York, 1995. x+475 pp.
                 ISBN: 0-19-853489-2

        .. [STA] Richard Stanley, Enumerative combinatorics. Vol. 2.
                 With a foreword by Gian-Carlo Rota and appendix 1 by Sergey Fomin.
                 Cambridge Studies in Advanced Mathematics, 62. Cambridge University Press, Cambridge, 1999. xii+581 pp.
                 ISBN: 0-521-56069-1; 0-521-78987-7

        .. [ST94]  Scharf, Thomas, Thibon, Jean-Yves,
                 A Hopf-algebra approach to inner plethysm.
                 Adv. Math.  104  (1994),  no. 1, 30-58.
                 :doi:`10.1006/aima.1994.1019`

    .. rubric:: Further tests

    TESTS::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym
        Symmetric Functions over Rational Field
        sage: h = Sym.h(); e = Sym.e(); s = Sym.s(); m = Sym.m(); p = Sym.p()
        sage: ( ( h[2,1] * ( 1 + 3 * h[2,1]) ) + s[2]. antipode()) . coproduct()
        h[] # h[1, 1] - h[] # h[2] + h[] # h[2, 1] + 3*h[] # h[2, 2, 1, 1] + h[1] # h[1] + h[1] # h[1, 1]
        + h[1] # h[2] + 6*h[1] # h[2, 1, 1, 1] + 6*h[1] # h[2, 2, 1] + h[1, 1] # h[] + h[1, 1] # h[1]
        + 3*h[1, 1] # h[1, 1, 1, 1] + 12*h[1, 1] # h[2, 1, 1] + 3*h[1, 1] # h[2, 2] + 6*h[1, 1, 1] # h[1, 1, 1]
        + 6*h[1, 1, 1] # h[2, 1] + 3*h[1, 1, 1, 1] # h[1, 1] - h[2] # h[] + h[2] # h[1] + 6*h[2] # h[2, 1, 1]
        + h[2, 1] # h[] + 6*h[2, 1] # h[1, 1, 1] + 12*h[2, 1] # h[2, 1] + 12*h[2, 1, 1] # h[1, 1]
        + 6*h[2, 1, 1] # h[2] + 6*h[2, 1, 1, 1] # h[1] + 3*h[2, 2] # h[1, 1] + 6*h[2, 2, 1] # h[1] + 3*h[2, 2, 1, 1] # h[]

    .. TODO::

        - Introduce fields with degree 1 elements as in
          MuPAD-Combinat, to get proper plethysm.
        - Use UniqueRepresentation to get rid of all the manual cache
          handling for the bases
        - Devise a mechanism so that pickling bases of symmetric
          functions pickles the coercions which have a cache.
    """

    def __init__(self, R):
        r"""
        Initialization of ``self``.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)

        TESTS::

            sage: TestSuite(Sym).run()

        """
        assert(R in Rings())
        self._base = R # Won't be needed when CategoryObject won't override anymore base_ring
        Parent.__init__(self, category = GradedHopfAlgebras(R).WithRealizations())

    def a_realization(self):
        r"""
        Returns a particular realization of ``self`` (the Schur basis).

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: Sym.a_realization()
            Symmetric Functions over Rational Field in the Schur basis
        """
        return self.schur()

    def _repr_(self): # could be taken care of by the category
        r"""
        Representation of ``self``

        TESTS::

            sage: SymmetricFunctions(RR) # indirect doctest
            Symmetric Functions over Real Field with 53 bits of precision
        """
        return "Symmetric Functions over %s"%self.base_ring()

    def schur(self):
        r"""
        The Schur basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).schur()
            Symmetric Functions over Rational Field in the Schur basis
        """
        return schur.SymmetricFunctionAlgebra_schur(self)
    s = schur
    Schur = schur # Currently needed by SymmetricFunctions.__init_extra__
                  # and sfa.SymmetricFunctionsBases.corresponding_basis_over

    def powersum(self):
        r"""
        The power sum basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).powersum()
            Symmetric Functions over Rational Field in the powersum basis
        """
        return powersum.SymmetricFunctionAlgebra_power(self)
    p = powersum
    power = powersum # Todo: get rid of this one when it won't be needed anymore

    def complete(self):
        r"""
        The complete basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).complete()
            Symmetric Functions over Rational Field in the homogeneous basis
        """
        return homogeneous.SymmetricFunctionAlgebra_homogeneous(self)
    h = complete
    homogeneous = complete

    def elementary(self):
        r"""
        The elementary basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).elementary()
            Symmetric Functions over Rational Field in the elementary basis
        """
        return elementary.SymmetricFunctionAlgebra_elementary(self)
    e = elementary

    def monomial(self):
        r"""
        The monomial basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).monomial()
            Symmetric Functions over Rational Field in the monomial basis
        """
        return monomial.SymmetricFunctionAlgebra_monomial(self)
    m = monomial

    def witt(self, coerce_h=True, coerce_e=False, coerce_p=False):
        r"""
        The Witt basis of the symmetric functions.

        EXAMPLES::

            sage: SymmetricFunctions(QQ).witt()
            Symmetric Functions over Rational Field in the Witt basis
            sage: SymmetricFunctions(QQ).witt(coerce_p=True)
            Symmetric Functions over Rational Field in the Witt basis
            sage: SymmetricFunctions(QQ).witt(coerce_h=False, coerce_e=True, coerce_p=True)
            Symmetric Functions over Rational Field in the Witt basis
        """
        import witt
        return witt.SymmetricFunctionAlgebra_witt(self, coerce_h=coerce_h, coerce_e=coerce_e, coerce_p=coerce_p)
    w = witt
    # Currently needed by sfa.SymmetricFunctionsBases.corresponding_basis_over
    Witt = witt

    def forgotten(self):
        r"""
        The forgotten basis of the Symmetric Functions (or the basis dual to
        the elementary basis with respect to the Hall scalar product).

        EXAMPLES::

            sage: SymmetricFunctions(QQ).forgotten()
            Symmetric Functions over Rational Field in the forgotten basis

        TESTS:

        Over the rationals::

            sage: Sym = SymmetricFunctions(QQ)
            sage: e = Sym.e()
            sage: f = Sym.f()
            sage: h = Sym.h()
            sage: p = Sym.p()
            sage: s = Sym.s()
            sage: m = Sym.m()
            sage: e(f([2,1]))
            -2*e[1, 1, 1] + 5*e[2, 1] - 3*e[3]
            sage: f(e([2,1]))
            3*f[1, 1, 1] + 2*f[2, 1] + f[3]
            sage: h(f([2,1]))
            h[2, 1] - 3*h[3]
            sage: f(h([2,1]))
            3*f[1, 1, 1] + f[2, 1]
            sage: p(f([2,1]))
            -p[2, 1] - p[3]
            sage: f(p([2,1]))
            -f[2, 1] - f[3]
            sage: s(f([2,1]))
            s[2, 1] - 2*s[3]
            sage: f(s([2,1]))
            2*f[1, 1, 1] + f[2, 1]
            sage: m(f([2,1]))
            -m[2, 1] - 2*m[3]
            sage: f(m([2,1]))
            -f[2, 1] - 2*f[3]

        Over the integers::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: e = Sym.e()
            sage: f = Sym.f()
            sage: h = Sym.h()
            sage: p = Sym.p()
            sage: s = Sym.s()
            sage: m = Sym.m()
            sage: e(f([2,1]))
            -2*e[1, 1, 1] + 5*e[2, 1] - 3*e[3]
            sage: f(e([2,1]))
            3*f[1, 1, 1] + 2*f[2, 1] + f[3]
            sage: h(f([2,1]))
            h[2, 1] - 3*h[3]
            sage: f(h([2,1]))
            3*f[1, 1, 1] + f[2, 1]
            sage: f(p([2,1]))
            -f[2, 1] - f[3]
            sage: s(f([2,1]))
            s[2, 1] - 2*s[3]
            sage: f(s([2,1]))
            2*f[1, 1, 1] + f[2, 1]
            sage: m(f([2,1]))
            -m[2, 1] - 2*m[3]
            sage: f(m([2,1]))
            -f[2, 1] - 2*f[3]

        Conversion from the forgotten basis to the power-sum basis over the
        integers is not well-defined in general, even if the result happens
        to have integral coefficients::

            sage: p(f([2,1]))
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        Fun exercise: prove that `p(f_{\lambda})` and `p(m_{\lambda})` have
        integral coefficients whenever `\lambda` is a strict partition.
        """
        return self.elementary().dual_basis()
    f = forgotten

    def macdonald(self, q='q', t='t'):
        r"""
        Returns the entry point for the various Macdonald bases.

        INPUT:

        - ``q``, ``t`` -- parameters

        Macdonald symmetric functions including bases `P`, `Q`, `J`, `H`, `Ht`.
        This also contains the `S` basis which is dual to the Schur basis with
        respect to the `q,t` scalar product.

        The parameters `q` and `t` must be in the base_ring of parent.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P(); P
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
            sage: P[2]
            McdP[2]
            sage: Q = Sym.macdonald().Q(); Q
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald Q basis
            sage: S = Sym.macdonald().S()
            sage: s = Sym.schur()
            sage: matrix([[S(la).scalar_qt(s(mu)) for la in Partitions(3)] for mu in Partitions(3)])
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: H = Sym.macdonald().H()
            sage: s(H[2,2])
            q^2*s[1, 1, 1, 1] + (q^2*t+q*t+q)*s[2, 1, 1] + (q^2*t^2+1)*s[2, 2] + (q*t^2+q*t+t)*s[3, 1] + t^2*s[4]

            sage: Sym = SymmetricFunctions(QQ['z','q'].fraction_field())
            sage: (z,q) = Sym.base_ring().gens()
            sage: Hzq = Sym.macdonald(q=z,t=q).H()
            sage: H1z = Sym.macdonald(q=1,t=z).H()
            sage: s = Sym.schur()
            sage: s(H1z([2,2]))
            s[1, 1, 1, 1] + (2*z+1)*s[2, 1, 1] + (z^2+1)*s[2, 2] + (z^2+2*z)*s[3, 1] + z^2*s[4]
            sage: s(Hzq[2,2])
            z^2*s[1, 1, 1, 1] + (z^2*q+z*q+z)*s[2, 1, 1] + (z^2*q^2+1)*s[2, 2] + (z*q^2+z*q+q)*s[3, 1] + q^2*s[4]
            sage: s(H1z(Hzq[2,2]))
            z^2*s[1, 1, 1, 1] + (z^2*q+z*q+z)*s[2, 1, 1] + (z^2*q^2+1)*s[2, 2] + (z*q^2+z*q+q)*s[3, 1] + q^2*s[4]
        """
        return macdonald.Macdonald(self, q=q, t=t)

    def hall_littlewood(self, t='t'):
        """
        Returns the entry point for the various Hall-Littlewood bases.

        INPUT:

        - ``t`` -- parameter

        Hall-Littlewood symmetric functions including bases `P`, `Q`, `Qp`.
        The Hall-Littlewood `P` and `Q` functions at `t=-1` are the
        Schur-P and Schur-Q functions when indexed by strict partitions.

        The parameter `t` must be in the base ring of parent.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: P = Sym.hall_littlewood().P(); P
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood P basis
            sage: P[2]
            HLP[2]
            sage: Q = Sym.hall_littlewood().Q(); Q
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood Q basis
            sage: Q[2]
            HLQ[2]
            sage: Qp = Sym.hall_littlewood().Qp(); Qp
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood Qp basis
            sage: Qp[2]
            HLQp[2]
        """
        return hall_littlewood.HallLittlewood(self, t=t)

    def jack(self, t='t'):
        """
        Returns the entry point for the various Jack bases.

        INPUT:

        - ``t`` -- parameter

        Jack symmetric functions including bases `P`, `Q`, `Qp`.

        The parameter `t` must be in the base ring of parent.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JP = Sym.jack().P(); JP
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack P basis
            sage: JQ = Sym.jack().Q(); JQ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack Q basis
            sage: JJ = Sym.jack().J(); JJ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack J basis
            sage: JQp = Sym.jack().Qp(); JQp
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack Qp basis
        """
        return jack.Jack( self, t=t )

    def zonal(self):
        """
        The zonal basis of the Symmetric Functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).zonal()
            Symmetric Functions over Rational Field in the zonal basis
        """
        return jack.SymmetricFunctionAlgebra_zonal( self )

    def llt(self, k, t='t'):
        """
        The LLT symmetric functions.

        INPUT:

        - ``k`` -- a positive integer indicating the level
        - ``t`` -- a parameter (default: `t`)

        LLT polynomials in `hspin` and `hcospin` bases.

        EXAMPLES::

            sage: llt3 = SymmetricFunctions(QQ['t'].fraction_field()).llt(3); llt3
            level 3 LLT polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: llt3.hspin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT spin basis
            sage: llt3.hcospin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT cospin basis
            sage: llt3.hcospin()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the level 3 LLT cospin basis
        """
        return llt.LLT_class( self, k, t=t )

    def from_polynomial(self, f):
        """
        Converts a symmetric polynomial ``f`` to a symmetric function.

        INPUT:

        - ``f`` -- a symmetric polynomial

        This function converts a symmetric polynomial `f` in a polynomial ring in finitely
        many variables to a symmetric function in the monomial
        basis of the ring of symmetric functions over the same base ring.

        EXAMPLES::

            sage: P = PolynomialRing(QQ, 'x', 3)
            sage: x = P.gens()
            sage: f = x[0] + x[1] + x[2]
            sage: S = SymmetricFunctions(QQ)
            sage: S.from_polynomial(f)
            m[1]

            sage: f = x[0] + 2*x[1] + x[2]
            sage: S.from_polynomial(f)
            Traceback (most recent call last):
            ...
            ValueError: x0 + 2*x1 + x2 is not a symmetric polynomial
        """
        return self.m().from_polynomial(f)

    def register_isomorphism(self, morphism, only_conversion=False):
        """
        Register an isomorphism between two bases of ``self``, as a canonical coercion
        (unless the optional keyword ``only_conversion`` is set to ``True``,
        in which case the isomorphism is registered as conversion only).

        EXAMPLES:

        We override the canonical coercion from the Schur basis to the
        powersum basis by a (stupid!) map `s_\lambda\mapsto 2p_\lambda`.
        ::

            sage: Sym = SymmetricFunctions(QQ['zorglub']) # make sure we are not going to screw up later tests
            sage: s = Sym.s(); p = Sym.p().dual_basis()
            sage: phi = s.module_morphism(diagonal = lambda t: 2, codomain = p)
            sage: phi(s[2, 1])
            2*d_p[2, 1]
            sage: Sym.register_isomorphism(phi)
            sage: p(s[2,1])
            2*d_p[2, 1]

        The map is supposed to implement the canonical isomorphism
        between the two bases. Otherwise, the results will be
        mathematically wrong, as above. Use with care!
        """
        if only_conversion:
            morphism.codomain().register_conversion(morphism)
        else:
            morphism.codomain().register_coercion(morphism)

    _shorthands = set(['e', 'h', 'm', 'p', 's'])

    def inject_shorthands(self, shorthands = _shorthands):
        """
        Imports standard shorthands into the global namespace

        INPUT:

        - ``shorthands`` -- a list (or iterable) of strings (default: ['e', 'h', 'm', 'p', 's'])

        EXAMPLES::

            sage: S = SymmetricFunctions(ZZ)
            sage: S.inject_shorthands()
            sage: s[1] + e[2] * p[1,1] + 2*h[3] + m[2,1]
            s[1] - 2*s[1, 1, 1] + s[1, 1, 1, 1] + s[2, 1] + 2*s[2, 1, 1] + s[2, 2] + 2*s[3] + s[3, 1]
            sage: e
            Symmetric Functions over Integer Ring in the elementary basis
            sage: p
            Symmetric Functions over Integer Ring in the powersum basis
            sage: s
            Symmetric Functions over Integer Ring in the Schur basis

            sage: e == S.e(), h == S.h(), m == S.m(), p == S.p(), s == S.s()
            (True, True, True, True, True)

        One can also just import a subset of the shorthands::

            sage: S = SymmetricFunctions(QQ)
            sage: S.inject_shorthands(['p', 's'])
            sage: p
            Symmetric Functions over Rational Field in the powersum basis
            sage: s
            Symmetric Functions over Rational Field in the Schur basis

        Note that ``e`` is left unchanged::

            sage: e
            Symmetric Functions over Integer Ring in the elementary basis
        """
        from sage.misc.misc import inject_variable
        for shorthand in shorthands:
            assert shorthand in self._shorthands
            inject_variable(shorthand, getattr(self, shorthand)())

    def __init_extra__(self):
        """
        Sets up the coercions between the different bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ) # indirect doctest
            sage: s = Sym.s(); p = Sym.p()
            sage: s.coerce_map_from(p)
            Generic morphism:
              From: Symmetric Functions over Rational Field in the powersum basis
              To:   Symmetric Functions over Rational Field in the Schur basis
        """
        #powersum   = self.powersum  ()
        #complete   = self.complete  ()
        #elementary = self.elementary()
        #schur      = self.schur     ()
        #monomial   = self.monomial  ()

        iso = self.register_isomorphism

        from sage.combinat.sf.classical import conversion_functions

        for (basis1_name, basis2_name) in conversion_functions.keys():
            basis1 = getattr(self, basis1_name)()
            basis2 = getattr(self, basis2_name)()
            on_basis = SymmetricaConversionOnBasis(t = conversion_functions[basis1_name,basis2_name], domain = basis1, codomain = basis2)
            from sage.rings.rational_field import RationalField
            if basis2_name != "powersum" or self._base.has_coerce_map_from(RationalField()):
                iso(basis1._module_morphism(on_basis, codomain = basis2))
            else:
                # Don't register conversions to powersums as coercions,
                # unless the base ring is a `\QQ`-algebra
                # (otherwise the coercion graph loses commutativity).
                iso(basis1._module_morphism(on_basis, codomain = basis2), only_conversion = True)

        # Todo: fill in with other conversion functions on the classical bases

    def kBoundedSubspace(self, k, t='t'):
        r"""
        Return the `k`-bounded subspace of the ring of symmetric functions.

        INPUT:

        - ``k`` - a positive integer
        - ``t`` a formal parameter; `t=1` yields a subring

        The subspace of the ring of symmetric functions spanned by
        `\{ s_{\lambda}[X/(1-t)] \}_{\lambda_1\le k} = \{ s_{\lambda}^{(k)}[X,t]\}_{\lambda_1 \le k}`
        over the base ring `\mathbb{Q}[t]`. When `t=1`, this space is in fact a subalgebra of
        the ring of symmetric functions generated by the complete homogeneous symmetric functions
        `h_i` for `1\le i \le k`.

        .. seealso:: :meth:`sage.combinat.sf.new_kschur.KBoundedSubspace`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: KB = Sym.kBoundedSubspace(3,1); KB
            3-bounded Symmetric Functions over Rational Field with t=1

            sage: Sym = SymmetricFunctions(QQ['t'])
            sage: Sym.kBoundedSubspace(3)
            3-bounded Symmetric Functions over Univariate Polynomial Ring in t over Rational Field

            sage: Sym = SymmetricFunctions(QQ['z'])
            sage: z = Sym.base_ring().gens()[0]
            sage: Sym.kBoundedSubspace(3,t=z)
            3-bounded Symmetric Functions over Univariate Polynomial Ring in z over Rational Field with t=z
        """
        from sage.combinat.sf.new_kschur import KBoundedSubspace
        return KBoundedSubspace(self, k, t=t)

    def kschur(self, k, t ='t'):
        r"""
        Returns the `k`-Schur functions.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ks = Sym.kschur(3,1)
            sage: ks[2]*ks[2]
            ks3[2, 2] + ks3[3, 1]
            sage: ks[2,1,1].lift()
            s[2, 1, 1] + s[3, 1]

            sage: Sym = SymmetricFunctions(QQ['t'])
            sage: ks = Sym.kschur(3)
            sage: ks[2,2,1].lift()
            s[2, 2, 1] + t*s[3, 2]
        """
        return self.kBoundedSubspace(k, t=t).kschur()

    def khomogeneous(self, k):
        r"""
        Returns the homogeneous symmetric functions in the `k`-bounded subspace.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: kh = Sym.khomogeneous(4)
            sage: kh[3]*kh[4]
            h4[4, 3]
            sage: kh[4].lift()
            h[4]
        """
        return self.kBoundedSubspace(k, t=1).khomogeneous()

    def kBoundedQuotient(self, k, t='t'):
        r"""
        Returns the `k`-bounded quotient space of the ring of symmetric functions.

        INPUT:

        - ``k`` - a positive integer

        The quotient of the ring of symmetric functions ...

        .. seealso:: :meth:`sage.combinat.sf.k_dual.KBoundedQuotient`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: KQ = Sym.kBoundedQuotient(3); KQ
            Traceback (most recent call last):
            ...
            TypeError: unable to convert t to a rational
            sage: KQ = Sym.kBoundedQuotient(3,t=1); KQ
            3-Bounded Quotient of Symmetric Functions over Rational Field with t=1
            sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
            sage: KQ = Sym.kBoundedQuotient(3); KQ
            3-Bounded Quotient of Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        from sage.combinat.sf.k_dual import KBoundedQuotient
        return KBoundedQuotient(self, k, t)

class SymmetricaConversionOnBasis:
    def __init__(self, t, domain, codomain):
        """
        Initialization of ``self``.

        INPUT:

        - ``t`` -- a function taking a monomial in CombinatorialFreeModule(QQ, Partitions()),
           and returning a (partition, coefficient) list.

        - ``domain``, ``codomain`` -- parents

        Construct a function mapping a partition to an element of ``codomain``.

        This is a temporary quick hack to wrap around the existing
        symmetrica conversions, without changing their specs.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ[x])
            sage: p = Sym.p(); s = Sym.s()
            sage: def t(x) : [(p,c)] = x; return [ (p,2*c), (p.conjugate(), c) ]
            sage: f = sage.combinat.sf.sf.SymmetricaConversionOnBasis(t, p, s)
            sage: f(Partition([3,1]))
            s[2, 1, 1] + 2*s[3, 1]
        """
        self._domain = domain
        self.fake_sym = CombinatorialFreeModule(QQ, Partitions())
        self._codomain = codomain
        self._t = t

    def __call__(self, partition):
        """
            sage: Sym = SymmetricFunctions(QQ[x])
            sage: p = Sym.p(); s = Sym.s()
            sage: p[1] + s[1]                           # indirect doctest
            2*p[1]
        """
        # TODO: use self._codomain.sum_of_monomials, when the later
        # will have an optional optimization for the case when there
        # is no repetition in the support
        return self._codomain._from_dict(dict(self._t(self.fake_sym.monomial(partition))), coerce = True)

