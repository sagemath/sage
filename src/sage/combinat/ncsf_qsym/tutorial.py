# -*- coding: utf-8 -*-
r"""
Introduction to Quasisymmetric Functions

In this document we briefly explain the quasisymmetric function bases and
related functionality in Sage.   We assume the reader is familiar with the
package :class:`SymmetricFunctions`.

Quasisymmetric functions, denoted `QSym`, form a subring of the power
series ring in countably many variables. `QSym` contains the symmetric
functions.  These functions first arose in the theory of
`P`-partitions.  The initial ideas in this field are attributed to
MacMahon, Knuth, Kreweras, Glânffrwd Thomas, Stanley. In 1984, Gessel
formalized the study of quasisymmetric functions and introduced the
basis of fundamental quasisymmetric functions [Ges]_. In 1995, Gelfand,
Krob, Lascoux, Leclerc, Retakh, and Thibon showed that the ring of
quasisymmetric functions is Hopf dual to the noncommutative symmetric
functions [NCSF]_.  Many results have built on these.

One advantage of working in `QSym` is that many interesting families of
symmetric functions have explicit expansions in fundamental quasisymmetric
functions such as Schur functions [Ges]_, Macdonald polynomials
[HHL05]_, and plethysm of Schur functions [LW12]_.

For more background see :wikipedia:`Quasisymmetric_function`.

To begin, initialize the ring.  Below we chose to use the rational
numbers `\QQ`.  Other options include the integers `\ZZ` and `\CC`::

    sage: QSym = QuasiSymmetricFunctions(QQ)
    sage: QSym
    Quasisymmetric functions over the Rational Field

    sage: QSym = QuasiSymmetricFunctions(CC); QSym
    Quasisymmetric functions over the Complex Field with 53 bits of precision

    sage: QSym = QuasiSymmetricFunctions(ZZ); QSym
    Quasisymmetric functions over the Integer Ring

All bases of `QSym` are indexed by compositions e.g. `[3,1,1,4]`. The
convention is to use capital letters for bases of `QSym` and lowercase
letters for bases of the symmetric functions `Sym`.  Next set up names for the
known bases by running ``inject_shorthands()``. As with symmetric functions,
you do not need to run this command and you could assign these bases other
names. ::

    sage: QSym = QuasiSymmetricFunctions(QQ)
    sage: QSym.inject_shorthands()
    Defining M as shorthand for Quasisymmetric functions over the Rational Field in the Monomial basis
    Defining F as shorthand for Quasisymmetric functions over the Rational Field in the Fundamental basis
    Defining E as shorthand for Quasisymmetric functions over the Rational Field in the Essential basis
    Defining dI as shorthand for Quasisymmetric functions over the Rational Field in the dualImmaculate basis
    Defining QS as shorthand for Quasisymmetric functions over the Rational Field in the Quasisymmetric Schur basis
    Defining YQS as shorthand for Quasisymmetric functions over the Rational Field in the Young Quasisymmetric Schur basis
    Defining phi as shorthand for Quasisymmetric functions over the Rational Field in the phi basis
    Defining psi as shorthand for Quasisymmetric functions over the Rational Field in the psi basis

Now one can start constructing quasisymmetric functions.

.. NOTE::

    It is best to use variables other than ``M`` and ``F``.

::

    sage: x = M[2,1] + M[1,2]
    sage: x
    M[1, 2] + M[2, 1]

    sage: y = 3*M[1,2] + M[3]^2; y
    3*M[1, 2] + 2*M[3, 3] + M[6]

    sage: F[3,1,3] + 7*F[2,1]
    7*F[2, 1] + F[3, 1, 3]

    sage: 3*F[2,1,2] + F[3]^2
    F[1, 2, 2, 1] + F[1, 2, 3] + 2*F[1, 3, 2] + F[1, 4, 1] + F[1, 5] + 3*F[2, 1, 2]
     + 2*F[2, 2, 2] + 2*F[2, 3, 1] + 2*F[2, 4] + F[3, 2, 1] + 3*F[3, 3] + 2*F[4, 2] + F[5, 1] + F[6]

To convert from one basis to another is easy::

    sage: z = M[1,2,1]
    sage: z
    M[1, 2, 1]

    sage: F(z)
    -F[1, 1, 1, 1] + F[1, 2, 1]

    sage: M(F(z))
    M[1, 2, 1]

To expand in variables, one can specify a finite size alphabet `x_1, x_2,
\ldots, x_m`::

    sage: y = M[1,2,1]
    sage: y.expand(4)
    x0*x1^2*x2 + x0*x1^2*x3 + x0*x2^2*x3 + x1*x2^2*x3

The usual methods on free modules are available such as coefficients,
degrees, and the support::

    sage: z=3*M[1,2]+M[3]^2; z
    3*M[1, 2] + 2*M[3, 3] + M[6]

    sage: z.coefficient([1,2])
    3

    sage: z.degree()
    6

    sage: sorted(z.coefficients())
    [1, 2, 3]

    sage: sorted(z.monomials(), key=lambda x: x.support())
    [M[1, 2], M[3, 3], M[6]]

    sage: z.monomial_coefficients()
    {[1, 2]: 3, [3, 3]: 2, [6]: 1}

As with the symmetric functions package, the quasisymmetric function ``1``
has several instantiations. However, the most obvious way to write ``1``
leads to an error (this is due to the semantics of python)::

    sage: M[[]]
    M[]
    sage: M.one()
    M[]
    sage: M(1)
    M[]
    sage: M[[]] == 1
    True
    sage: M[]
    Traceback (most recent call last):
    ...
    SyntaxError: invalid ...


Working with symmetric functions
--------------------------------

The quasisymmetric functions are a ring which contains the symmetric
functions as a subring.  The Monomial quasisymmetric functions are
related to the monomial symmetric functions by `m_\lambda =
\sum_{\mathrm{sort}(c) = \lambda} M_c`, where `\mathrm{sort}(c)`
means the partition obtained by sorting the composition `c`::

    sage: SymmetricFunctions(QQ).inject_shorthands()
    Defining e as shorthand for Symmetric Functions over Rational Field in the elementary basis
    Defining f as shorthand for Symmetric Functions over Rational Field in the forgotten basis
    Defining h as shorthand for Symmetric Functions over Rational Field in the homogeneous basis
    Defining m as shorthand for Symmetric Functions over Rational Field in the monomial basis
    Defining p as shorthand for Symmetric Functions over Rational Field in the powersum basis
    Defining s as shorthand for Symmetric Functions over Rational Field in the Schur basis

    sage: m[2,1]
    m[2, 1]
    sage: M(m[2,1])
    M[1, 2] + M[2, 1]
    sage: M(s[2,1])
    2*M[1, 1, 1] + M[1, 2] + M[2, 1]

There are methods to test if an expression `f` in the quasisymmetric functions
is a symmetric function::

    sage: f = M[1,1,2] + M[1,2,1]
    sage: f.is_symmetric()
    False
    sage: f = M[3,1] + M[1,3]
    sage: f.is_symmetric()
    True

If `f` is symmetric, there are methods to convert `f` to an expression in the
symmetric functions::

    sage: f.to_symmetric_function()
    m[3, 1]

The expansion of the Schur function in terms of the Fundamental quasisymmetric
functions is due to [Ges]_. There is one term in the expansion for each
standard tableau of shape equal to the partition indexing the Schur function.
::

    sage: f = F[3,2] + F[2,2,1] + F[2,3] + F[1,3,1] + F[1,2,2]
    sage: f.is_symmetric()
    True
    sage: f.to_symmetric_function()
    5*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + 2*m[2, 2, 1] + m[3, 1, 1] + m[3, 2]
    sage: s(f.to_symmetric_function())
    s[3, 2]

It is also possible to convert any symmetric function to the quasisymmetric
function expansion in any known basis. The converse is not true::

    sage: M( m[3,1,1] )
    M[1, 1, 3] + M[1, 3, 1] + M[3, 1, 1]
    sage: F( s[2,2,1] )
    F[1, 1, 2, 1] + F[1, 2, 1, 1] + F[1, 2, 2] + F[2, 1, 2] + F[2, 2, 1]

    sage: s(M[2,1])
    Traceback (most recent call last):
    ...
    TypeError: do not know how to make x (= M[2, 1]) an element of self

It is possible to experiment with the quasisymmetric function expansion of other
bases, but it is important that the base ring be the same for both algebras.
::

    sage: R = QQ['t']
    sage: Qp = SymmetricFunctions(R).hall_littlewood().Qp()
    sage: QSymt = QuasiSymmetricFunctions(R)
    sage: Ft = QSymt.F()
    sage: Ft( Qp[2,2] )
    F[1, 2, 1] + t*F[1, 3] + (t+1)*F[2, 2] + t*F[3, 1] + t^2*F[4]

::

    sage: K = QQ['q','t'].fraction_field()
    sage: Ht = SymmetricFunctions(K).macdonald().Ht()
    sage: Fqt = QuasiSymmetricFunctions(Ht.base_ring()).F()
    sage: Fqt(Ht[2,1])
    q*t*F[1, 1, 1] + (q+t)*F[1, 2] + (q+t)*F[2, 1] + F[3]

The following will raise an error because the base ring of ``F`` is not
equal to the base ring of ``Ht``::

    sage: F(Ht[2,1])
    Traceback (most recent call last):
    ...
    TypeError: do not know how to make x (= McdHt[2, 1]) an element of self (=Quasisymmetric functions over the Rational Field in the Fundamental basis)

QSym is a Hopf algebra
----------------------

The product on `QSym` is commutative and is inherited from the
product by the realization within the polynomial ring::

    sage: M[3]*M[1,1] == M[1,1]*M[3]
    True
    sage: M[3]*M[1,1]
    M[1, 1, 3] + M[1, 3, 1] + M[1, 4] + M[3, 1, 1] + M[4, 1]
    sage: F[3]*F[1,1]
    F[1, 1, 3] + F[1, 2, 2] + F[1, 3, 1] + F[1, 4] + F[2, 1, 2] + F[2, 2, 1] + F[2, 3] + F[3, 1, 1] + F[3, 2] + F[4, 1]
    sage: M[3]*F[2]
    M[1, 1, 3] + M[1, 3, 1] + M[1, 4] + M[2, 3] + M[3, 1, 1] + M[3, 2] + M[4, 1] + M[5]
    sage: F[2]*M[3]
    F[1, 1, 1, 2] - F[1, 2, 2] + F[2, 1, 1, 1] - F[2, 1, 2] - F[2, 2, 1] + F[5]

There is a coproduct on this ring as well, which in the Monomial basis acts by
cutting the composition into a left half and a right half. The co-product is
non-co-commutative::

    sage: M[1,3,1].coproduct()
    M[] # M[1, 3, 1] + M[1] # M[3, 1] + M[1, 3] # M[1] + M[1, 3, 1] # M[]
    sage: F[1,3,1].coproduct()
    F[] # F[1, 3, 1] + F[1] # F[3, 1] + F[1, 1] # F[2, 1] + F[1, 2] # F[1, 1] + F[1, 3] # F[1] + F[1, 3, 1] # F[]

.. rubric:: The Duality Pairing with Non-Commutative Symmetric Functions

These two operations endow `QSym` with the structure of a Hopf algebra. It is
the dual Hopf algebra of the non-commutative symmetric functions `NCSF`. Under
this duality, the Monomial basis of `QSym` is dual to the Complete basis of
`NCSF`, and the Fundamental basis of `QSym` is dual to the Ribbon basis of
`NCSF` (see [MR]_)::

    sage: S = M.dual(); S
    Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
    sage: M[1,3,1].duality_pairing( S[1,3,1] )
    1
    sage: M.duality_pairing_matrix( S, degree=4 )
    [1 0 0 0 0 0 0 0]
    [0 1 0 0 0 0 0 0]
    [0 0 1 0 0 0 0 0]
    [0 0 0 1 0 0 0 0]
    [0 0 0 0 1 0 0 0]
    [0 0 0 0 0 1 0 0]
    [0 0 0 0 0 0 1 0]
    [0 0 0 0 0 0 0 1]
    sage: F.duality_pairing_matrix( S, degree=4 )
    [1 0 0 0 0 0 0 0]
    [1 1 0 0 0 0 0 0]
    [1 0 1 0 0 0 0 0]
    [1 1 1 1 0 0 0 0]
    [1 0 0 0 1 0 0 0]
    [1 1 0 0 1 1 0 0]
    [1 0 1 0 1 0 1 0]
    [1 1 1 1 1 1 1 1]
    sage: NCSF = M.realization_of().dual()
    sage: R = NCSF.Ribbon()
    sage: F.duality_pairing_matrix( R, degree=4 )
    [1 0 0 0 0 0 0 0]
    [0 1 0 0 0 0 0 0]
    [0 0 1 0 0 0 0 0]
    [0 0 0 1 0 0 0 0]
    [0 0 0 0 1 0 0 0]
    [0 0 0 0 0 1 0 0]
    [0 0 0 0 0 0 1 0]
    [0 0 0 0 0 0 0 1]
    sage: M.duality_pairing_matrix( R, degree=4 )
    [ 1  0  0  0  0  0  0  0]
    [-1  1  0  0  0  0  0  0]
    [-1  0  1  0  0  0  0  0]
    [ 1 -1 -1  1  0  0  0  0]
    [-1  0  0  0  1  0  0  0]
    [ 1 -1  0  0 -1  1  0  0]
    [ 1  0 -1  0 -1  0  1  0]
    [-1  1  1 -1  1 -1 -1  1]

Let `H` and `G` be elements of `QSym` and `h` an element of `NCSF`. Then if
we represent the duality pairing with the mathematical notation `[ \cdot,
\cdot ]`, we have:

.. MATH::

    [H \cdot G, h] = [H \otimes G, \Delta(h)].

For example, the coefficient of ``M[2,1,4,1]`` in ``M[1,3]*M[2,1,1]`` may be
computed with the duality pairing::

    sage: I, J = Composition([1,3]), Composition([2,1,1])
    sage: (M[I]*M[J]).duality_pairing(S[2,1,4,1])
    1

And the coefficient of ``S[1,3] # S[2,1,1]`` in ``S[2,1,4,1].coproduct()`` is
equal to this result::

    sage: S[2,1,4,1].coproduct()
    S[] # S[2, 1, 4, 1] + ... + S[1, 3] # S[2, 1, 1] + ... + S[4, 1] # S[2, 1]

The duality pairing on the tensor space is another way of getting this
coefficient, but currently the method
:meth:`~sage.combinat.ncsf_qsym.generic_basis_code.BasesOfQSymOrNCSF.ParentMethods.duality_pairing()`
is not defined on the tensor squared space. However, we can extend this
functionality by applying a linear morphism to the terms in the coproduct,
as follows::

    sage: X = S[2,1,4,1].coproduct()
    sage: def linear_morphism(x, y):
    ....:     return x.duality_pairing(M[1,3]) * y.duality_pairing(M[2,1,1])
    sage: X.apply_multilinear_morphism(linear_morphism, codomain=ZZ)
    1

Similarly, if `H` is an element of `QSym` and `g` and `h` are elements of
`NCSF`, then

.. MATH::

    [ H, g \cdot h ] = [ \Delta(H), g \otimes h ].

For example, the coefficient of ``R[2,3,1]`` in ``R[2,1]*R[2,1]`` is computed
with the duality pairing by the following command::

    sage: (R[2,1]*R[2,1]).duality_pairing(F[2,3,1])
    1
    sage: R[2,1]*R[2,1]
    R[2, 1, 2, 1] + R[2, 3, 1]

This coefficient should then be equal to the coefficient of ``F[2,1] # F[2,1]``
in ``F[2,3,1].coproduct()``::

    sage: F[2,3,1].coproduct()
    F[] # F[2, 3, 1] + ... + F[2, 1] # F[2, 1]  + ... + F[2, 3, 1] # F[]

This can also be computed by the duality pairing on the tensor space,
as above::

    sage: X = F[2,3,1].coproduct()
    sage: def linear_morphism(x, y):
    ....:     return x.duality_pairing(R[2,1]) * y.duality_pairing(R[2,1])
    sage: X.apply_multilinear_morphism(linear_morphism, codomain=ZZ)
    1

.. rubric:: The Operation Adjoint to Multiplication by a Non-Commutative Symmetric Function

Let `g \in NCSF` and consider the linear endomorphism of `NCSF` defined by
left (respectively, right) multiplication by `g`. Since there is a duality
between `QSym` and `NCSF`, this linear transformation induces an operator
`g^\perp` on `QSym` satisfying

.. MATH::

    [ g^\perp(H), h ] = [ H, g \cdot h ].

for any non-commutative symmetric function `h`.

This is implemented by the method
:meth:`~sage.combinat.ncsf_qsym.generic_basis_code.BasesOfQSymOrNCSF.ElementMethods.skew_by()`.
Explicitly, if ``H`` is a quasisymmetric function and ``g``
a non-commutative symmetric function, then ``H.skew_by(g)`` and
``H.skew_by(g, side='right')`` are expressions that satisfy,
for any non-commutative symmetric function ``h``, the following
identities::

    H.skew_by(g).duality_pairing(h) == H.duality_pairing(g*h)
    H.skew_by(g, side='right').duality_pairing(h) == H.duality_pairing(h*g)

For example, ``M[J].skew_by(S[I])`` is `0` unless the composition `J`
begins with `I` and ``M(J).skew_by(S(I), side='right')`` is `0` unless
the composition `J` ends with `I`::

    sage: M[3,2,2].skew_by(S[3])
    M[2, 2]
    sage: M[3,2,2].skew_by(S[2])
    0
    sage: M[3,2,2].coproduct().apply_multilinear_morphism( lambda x,y: x.duality_pairing(S[3])*y )
    M[2, 2]
    sage: M[3,2,2].skew_by(S[3], side='right')
    0
    sage: M[3,2,2].skew_by(S[2], side='right')
    M[3, 2]

.. rubric:: The antipode

The antipode sends the Fundamental basis element indexed by the
composition `I` to `-1` to the size of `I` times the Fundamental
basis element indexed by the conjugate composition to `I`::

    sage: F[3,2,2].antipode()
    -F[1, 2, 2, 1, 1]
    sage: Composition([3,2,2]).conjugate()
    [1, 2, 2, 1, 1]
    sage: M[3,2,2].antipode()
    -M[2, 2, 3] - M[2, 5] - M[4, 3] - M[7]

We demonstrate here the defining relation of the antipode::

    sage: X = F[3,2,2].coproduct()
    sage: X.apply_multilinear_morphism(lambda x,y: x*y.antipode())
    0
    sage: X.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
    0

REFERENCES:

.. [HHL05] *A combinatorial formula for Macdonald polynomials*.
   Haiman, Haglund, and Loehr.
   J. Amer. Math. Soc. 18 (2005), no. 3, 735-761.

.. [LW12] *Quasisymmetric expansions of Schur-function plethysms*.
   Loehr and Warrington.
   Proc. Amer. Math. Soc. 140 (2012), no. 4, 1159-1171.

.. [KT97] *Noncommutative symmetric functions IV: Quantum linear groups and
   Hecke algebras at* `q = 0`.
   Krob and Thibon.
   Journal of Algebraic Combinatorics 6 (1997), 339-376.
"""
