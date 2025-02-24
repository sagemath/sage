# sage.doctest: needs sage.combinat sage.modules
r"""
Quasisymmetric functions

REFERENCES:

.. [Ges] \I. Gessel, *Multipartite P-partitions and inner products of skew Schur
   functions*, Contemp. Math. **34** (1984), 289-301.
   http://people.brandeis.edu/~gessel/homepage/papers/multipartite.pdf

.. [MR] \C. Malvenuto and C. Reutenauer, *Duality between quasi-symmetric
   functions and the Solomon descent algebra*, J. Algebra **177** (1995),
   no. 3, 967-982. http://www.mat.uniroma1.it/people/malvenuto/Duality.pdf

.. [GriRei18]_

.. [Mal1993] Claudia Malvenuto, *Produits et coproduits des fonctions
   quasi-symétriques et de l'algèbre des descentes*,
   thesis, November 1993.
   http://www1.mat.uniroma1.it/people/malvenuto/Thesis.pdf

.. [Haz2004] Michiel Hazewinkel, *Explicit polynomial generators for the
   ring of quasisymmetric functions over the integers*.
   :arxiv:`math/0410366v1`

.. [Rad1979] David E. Radford, *A natural ring basis for the shuffle algebra
   and an application to group schemes*, J. Algebra **58** (1979), 432-454.

.. [NCSF1] Israel Gelfand, D. Krob, Alain Lascoux, B. Leclerc,
   V. S. Retakh, J.-Y. Thibon,
   *Noncommutative symmetric functions*.
   :arxiv:`hep-th/9407124v1`

.. [NCSF2] \D. Krob, B. Leclerc, J.-Y. Thibon,
   *Noncommutative symmetric functions II: Transformations of alphabets*.
   http://www-igm.univ-mlv.fr/~jyt/ARTICLES/NCSF2.ps

.. [HLNT09] \F. Hivert, J.-G. Luque, J.-C. Novelli, J.-Y. Thibon,
   *The (1-E)-transform in combinatorial Hopf algebras*.
   :arxiv:`math/0912.0184v2`

.. [LMvW13] Kurt Luoto, Stefan Mykytiuk and Stephanie van Willigenburg,
   *An introduction to quasisymmetric Schur functions -- Hopf algebras,
   quasisymmetric functions, and Young composition tableaux*,
   May 23, 2013, Springer.
   http://www.math.ubc.ca/%7Esteph/papers/QuasiSchurBook.pdf

.. [BBSSZ2012] Chris Berg, Nantel Bergeron, Franco Saliola,
   Luis Serrano, Mike Zabrocki,
   *A lift of the Schur and Hall-Littlewood bases to
   non-commutative symmetric functions*,
   :arxiv:`1208.5191v3`.

.. [Hoff2015] Michael Hoffman.
   *Quasi-symmetric functions and mod* `p` *multiple harmonic sums*.
   Kyushu J. Math. **69** (2015), pp. 345-366.
   :doi:`10.2206/kyushujm.69.345`, :arxiv:`math/0401319v3`.

.. [BDHMN2017] Cristina Ballantine, Zajj Daugherty, Angela Hicks, Sarah Mason,
   Elizabeth Niese. *Quasisymmetric power sums*. :arxiv:`1710.11613`.

.. [AHM2018] Edward Allen, Joshua Hallam, Sarah Mason, *Dual Immaculate
   Quasisymmetric Functions Expand Positively into Young Quasisymmetric
   Schur Functions*. :arxiv:`1606.03519`

.. [SW2010] John Shareshian and Michelle Wachs.
   *Eulerian quasisymmetric functions*. (2010).
   :arxiv:`0812.0764v2`

AUTHOR:

- Jason Bandlow
- Franco Saliola
- Chris Berg
- Darij Grinberg
"""
# ****************************************************************************
#       Copyright (C) 2010 Jason Bandlow <jbandlow@gmail.com>,
#                     2012 Franco Saliola <saliola@gmail.com>,
#                     2012 Chris Berg <chrisjamesberg@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.rings import Rings
from sage.categories.fields import Fields
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.permutation import Permutations
from sage.combinat.composition import Composition, Compositions
from sage.combinat.partition import Partitions, _Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.ncsf_qsym.generic_basis_code import BasesOfQSymOrNCSF
from sage.combinat.ncsf_qsym.combinatorics import (number_of_fCT, number_of_SSRCT,
                   compositions_order, coeff_pi, coeff_lp, coeff_sp, coeff_ell)
from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
from sage.combinat.words.word import Word
from sage.combinat.tableau import StandardTableaux
from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.superseded import deprecated_function_alias
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.matrix.matrix_space', 'MatrixSpace')


class QuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    .. rubric:: The Hopf algebra of quasisymmetric functions.

    Let `R` be a commutative ring with unity.
    The `R`-algebra of quasi-symmetric functions may be realized as an
    `R`-subalgebra of the ring of power series in countably many
    variables `R[[x_1, x_2, x_3, \ldots]]`. It consists of those
    formal power series `p` which are degree-bounded (i. e., the degrees
    of all monomials occurring with nonzero coefficient in `p` are bounded
    from above, although the bound can depend on `p`) and satisfy the
    following condition: For every tuple `(a_1, a_2, \ldots, a_m)` of
    positive integers, the coefficient of the monomial
    `x_{i_1}^{a_1} x_{i_2}^{a_2} \cdots x_{i_m}^{a_m}` in `p` is the same
    for all strictly increasing sequences `(i_1 < i_2 < \cdots < i_m)` of
    positive integers. (In other words, the coefficient of a monomial in `p`
    depends only on the sequence of nonzero exponents in the monomial. If
    "sequence" were to be replaced by "multiset" here, we would obtain
    the definition of a symmetric function.)

    The `R`-algebra of quasi-symmetric functions is commonly called
    `\mathrm{QSym}_R` or occasionally just `\mathrm{QSym}` (when
    `R` is clear from the context or `\ZZ` or `\QQ`). It is graded by
    the total degree of the power series. Its homogeneous elements of degree
    `k` form a free `R`-submodule of rank equal to the number of
    compositions of `k` (that is, `2^{k-1}` if `k \geq 1`, else `1`).

    The two classical bases of `\mathrm{QSym}`, the monomial basis
    `(M_I)_I` and the fundamental basis `(F_I)_I`, are indexed by
    compositions `I = (I_1, I_2, \cdots, I_\ell )` and defined by the
    formulas:

    .. MATH::

        M_I = \sum_{1 \leq i_1 < i_2 < \cdots < i_\ell} x_{i_1}^{I_1}
        x_{i_2}^{I_2} \cdots x_{i_\ell}^{I_\ell}

    and

    .. MATH::

        F_I = \sum_{(j_1, j_2, \ldots, j_n)} x_{j_1} x_{j_2} \cdots
        x_{j_n}

    where in the second equation the sum runs over all weakly increasing
    `n`-tuples `(j_1, j_2, \ldots, j_n)` of positive integers
    (where `n` is the size of `I`) which increase strictly from `j_r`
    to `j_{r+1}` if `r` is a descent of the composition `I`.

    These bases are related by the formula

    `F_I = \sum_{J \leq I} M_J`

    where the inequality `J \leq I` indicates that `J` is finer than `I`.

    The `R`-algebra of quasi-symmetric functions is a Hopf algebra,
    with the coproduct satisfying

    .. MATH::

        \Delta M_I = \sum_{k=0}^{\ell} M_{(I_1, I_2, \cdots, I_k)} \otimes
        M_{(I_{k+1}, I_{k+2}, \cdots , I_{\ell})}

    for every composition `I = (I_1, I_2, \cdots , I_\ell )`.

    It is possible to define an `R`-algebra of quasi-symmetric
    functions in a finite number of variables as well (but it is not
    a bialgebra). These quasi-symmetric functions are actual polynomials
    then, not just power series.

    Chapter 5 of [GriRei18]_ and Section 11 of [HazWitt1]_ are devoted
    to quasi-symmetric functions, as are Malvenuto's thesis [Mal1993]_
    and part of Chapter 7 of [Sta-EC2]_.

    .. rubric:: The implementation of the quasi-symmetric function Hopf algebra

    We realize the `R`-algebra of quasi-symmetric functions in Sage as
    a graded Hopf algebra with basis elements indexed by compositions::

        sage: QSym = QuasiSymmetricFunctions(QQ)
        sage: QSym.category()
        Join of Category of Hopf algebras over Rational Field
            and Category of graded algebras over Rational Field
            and Category of commutative algebras over Rational Field
            and Category of monoids with realizations
            and Category of graded coalgebras over Rational Field
            and Category of coalgebras over Rational Field with realizations

    The most standard two bases for this `R`-algebra are the monomial and
    fundamental bases, and are accessible by the ``M()`` and ``F()`` methods::

        sage: M = QSym.M()
        sage: F = QSym.F()
        sage: M(F[2,1,2])
        M[1, 1, 1, 1, 1] + M[1, 1, 1, 2] + M[2, 1, 1, 1] + M[2, 1, 2]
        sage: F(M[2,1,2])
        F[1, 1, 1, 1, 1] - F[1, 1, 1, 2] - F[2, 1, 1, 1] + F[2, 1, 2]

    The product on this space is commutative and is inherited from the product
    on the realization within the ring of power series::

        sage: M[3]*M[1,1] == M[1,1]*M[3]
        True
        sage: M[3]*M[1,1]
        M[1, 1, 3] + M[1, 3, 1] + M[1, 4] + M[3, 1, 1] + M[4, 1]
        sage: F[3]*F[1,1]
        F[1, 1, 3] + F[1, 2, 2] + F[1, 3, 1] + F[1, 4] + F[2, 1, 2]
         + F[2, 2, 1] + F[2, 3] + F[3, 1, 1] + F[3, 2] + F[4, 1]
        sage: M[3]*F[2]
        M[1, 1, 3] + M[1, 3, 1] + M[1, 4] + M[2, 3] + M[3, 1, 1] + M[3, 2]
         + M[4, 1] + M[5]
        sage: F[2]*M[3]
        F[1, 1, 1, 2] - F[1, 2, 2] + F[2, 1, 1, 1] - F[2, 1, 2] - F[2, 2, 1]
         + F[5]

    There is a coproduct on `\mathrm{QSym}` as well, which in the Monomial
    basis acts by cutting the composition into a left half and a right
    half. The coproduct is not co-commutative::

        sage: M[1,3,1].coproduct()
        M[] # M[1, 3, 1] + M[1] # M[3, 1] + M[1, 3] # M[1] + M[1, 3, 1] # M[]
        sage: F[1,3,1].coproduct()
        F[] # F[1, 3, 1] + F[1] # F[3, 1] + F[1, 1] # F[2, 1]
         + F[1, 2] # F[1, 1] + F[1, 3] # F[1] + F[1, 3, 1] # F[]

    .. rubric:: The duality pairing with non-commutative symmetric functions

    These two operations endow the quasi-symmetric functions
    `\mathrm{QSym}` with the structure of a Hopf algebra. It is the graded
    dual Hopf algebra of the non-commutative symmetric functions `NCSF`.
    Under this duality, the Monomial basis of `\mathrm{QSym}` is dual to
    the Complete basis of `NCSF`, and the Fundamental basis of
    `\mathrm{QSym}` is dual to the Ribbon basis of `NCSF` (see [MR]_).

    ::

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

    Let `H` and `G` be elements of `\mathrm{QSym}`, and `h` an element of
    `NCSF`. Then, if we represent the duality pairing with the
    mathematical notation `[ \cdot, \cdot ]`,

    `[H G, h] = [H \otimes G, \Delta(h)]~.`

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
    coefficient, but currently the method ``duality_pairing`` is not defined on
    the tensor squared space. However, we can extend this functionality by
    applying a linear morphism to the terms in the coproduct, as follows::

        sage: X = S[2,1,4,1].coproduct()
        sage: def linear_morphism(x, y):
        ....:   return x.duality_pairing(M[1,3]) * y.duality_pairing(M[2,1,1])
        sage: X.apply_multilinear_morphism(linear_morphism, codomain=ZZ)
        1

    Similarly, if `H` is an element of `\mathrm{QSym}` and `g` and `h` are
    elements of `NCSF`, then

    .. MATH::

        [ H, g h ] = [ \Delta(H), g \otimes h ].

    For example, the coefficient of ``R[2,3,1]`` in ``R[2,1]*R[2,1]`` is
    computed with the duality pairing by the following command::

        sage: (R[2,1]*R[2,1]).duality_pairing(F[2,3,1])
        1
        sage: R[2,1]*R[2,1]
        R[2, 1, 2, 1] + R[2, 3, 1]

    This coefficient should then be equal to the coefficient of
    ``F[2,1] # F[2,1]`` in ``F[2,3,1].coproduct()``::

        sage: F[2,3,1].coproduct()
        F[] # F[2, 3, 1] + ... + F[2, 1] # F[2, 1]  + ... + F[2, 3, 1] # F[]

    This can also be computed by the duality pairing on the tensor space,
    as above::

        sage: X = F[2,3,1].coproduct()
        sage: def linear_morphism(x, y):
        ....:   return x.duality_pairing(R[2,1]) * y.duality_pairing(R[2,1])
        sage: X.apply_multilinear_morphism(linear_morphism, codomain=ZZ)
        1

    .. rubric:: The operation dual to multiplication by a non-commutative symmetric function

    Let `g \in NCSF` and consider the linear endomorphism of `NCSF` defined by
    left (respectively, right) multiplication by `g`. Since there is a duality
    between `\mathrm{QSym}` and `NCSF`, this linear transformation induces an
    operator `g^{\perp}` on `\mathrm{QSym}` satisfying

    .. MATH::

        [ g^{\perp}(H), h ] = [ H, gh ].

    for any non-commutative symmetric function `h`.

    This is implemented by the method
    :meth:`~sage.combinat.ncsf_qsym.generic_basis_code.BasesOfQSymOrNCSF.ElementMethods.skew_by`.
    Explicitly, if ``H`` is a quasi-symmetric function and ``g``
    a non-commutative symmetric function, then ``H.skew_by(g)`` and
    ``H.skew_by(g, side='right')`` are expressions that satisfy,
    for any non-commutative symmetric function ``h``, the following
    equalities::

        H.skew_by(g).duality_pairing(h) == H.duality_pairing(g*h)
        H.skew_by(g, side='right').duality_pairing(h) == H.duality_pairing(h*g)

    For example, ``M[J].skew_by(S[I])`` is `0` unless the composition ``J``
    begins with ``I`` and ``M(J).skew_by(S(I), side='right')`` is `0` unless
    the composition ``J`` ends with ``I``. For example::

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

    .. rubric:: The counit

    The counit is defined by sending all elements of positive degree to zero::

        sage: M[3].degree(), M[3,1,2].degree(), M.one().degree()
        (3, 6, 0)
        sage: M[3].counit()
        0
        sage: M[3,1,2].counit()
        0
        sage: M.one().counit()
        1
        sage: (M[3] - 2*M[3,1,2] + 7).counit()
        7
        sage: (F[3] - 2*F[3,1,2] + 7).counit()
        7

    .. rubric:: The antipode

    The antipode sends the Fundamental basis element indexed by the
    composition `I` to `(-1)^{|I|}` times the Fundamental
    basis element indexed by the conjugate composition to `I`
    (where `|I|` stands for the size of `I`, that is, the sum of all
    entries of `I`).
    ::

        sage: F[3,2,2].antipode()
        -F[1, 2, 2, 1, 1]
        sage: Composition([3,2,2]).conjugate()
        [1, 2, 2, 1, 1]

    The antipodes of the Monomial quasisymmetric functions can also be
    computed easily: Every composition `I` satisfies

    .. MATH::

        \omega(M_I) = (-1)^{\ell(I)} \sum M_J,

    where the sum ranges over all compositions `J` of `|I|`
    which are coarser than the reversed composition `I^r` of
    `I`. Here, `\ell(I)` denotes the length of the composition `I`
    (that is, the number of its parts). ::

        sage: M[3,2,1].antipode()
        -M[1, 2, 3] - M[1, 5] - M[3, 3] - M[6]
        sage: M[3,2,2].antipode()
        -M[2, 2, 3] - M[2, 5] - M[4, 3] - M[7]

    We demonstrate here the defining relation of the antipode::

        sage: X = F[3,2,2].coproduct()
        sage: X.apply_multilinear_morphism(lambda x,y: x*y.antipode())
        0
        sage: X.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
        0

    .. rubric:: The relation with symmetric functions

    The quasi-symmetric functions are a ring which contain the
    symmetric functions as a subring.  The Monomial quasi-symmetric
    functions are related to the monomial symmetric functions by

    .. MATH::

        m_\lambda = \sum_{\operatorname{sort}(I) = \lambda} M_I

    (where `\operatorname{sort}(I)` denotes the result of sorting
    the entries of `I` in decreasing order).

    There are methods to test if an expression in the quasi-symmetric
    functions is a symmetric function and, if it is, send it to an
    expression in the symmetric functions::

        sage: f = M[1,1,2] + M[1,2,1]
        sage: f.is_symmetric()
        False
        sage: g = M[3,1] + M[1,3]
        sage: g.is_symmetric()
        True
        sage: g.to_symmetric_function()
        m[3, 1]

    The expansion of the Schur function in terms of the Fundamental quasi-symmetric
    functions is due to [Ges]_.  There is one term in the expansion for each standard
    tableau of shape equal to the partition indexing the Schur function.
    ::

        sage: f = F[3,2] + F[2,2,1] + F[2,3] + F[1,3,1] + F[1,2,2]
        sage: f.is_symmetric()
        True
        sage: f.to_symmetric_function()
        5*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + 2*m[2, 2, 1] + m[3, 1, 1] + m[3, 2]
        sage: s = SymmetricFunctions(QQ).s()
        sage: s(f.to_symmetric_function())
        s[3, 2]

    It is also possible to convert a symmetric function to a
    quasi-symmetric function::

        sage: m = SymmetricFunctions(QQ).m()
        sage: M( m[3,1,1] )
        M[1, 1, 3] + M[1, 3, 1] + M[3, 1, 1]
        sage: F( s[2,2,1] )
        F[1, 1, 2, 1] + F[1, 2, 1, 1] + F[1, 2, 2] + F[2, 1, 2] + F[2, 2, 1]

    It is possible to experiment with the quasi-symmetric function expansion of other
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

    .. rubric:: The map to the ring of polynomials

    The quasi-symmetric functions can be seen as an inverse limit
    of a subring of a polynomial ring as the number of variables
    increases. Indeed, there exists a projection from the
    quasi-symmetric functions onto the polynomial ring
    `R[x_1, x_2, \ldots, x_n]`. This projection is defined by
    sending the variables `x_{n+1}, x_{n+2}, \cdots` to `0`, while
    the remaining `n` variables remain fixed. Note that this
    projection sends `M_I` to `0` if the length of the composition
    `I` is higher than `n`. ::

        sage: M[1,3,1].expand(4)
        x0*x1^3*x2 + x0*x1^3*x3 + x0*x2^3*x3 + x1*x2^3*x3
        sage: F[1,3,1].expand(4)
        x0*x1^3*x2 + x0*x1^3*x3 + x0*x1^2*x2*x3 + x0*x1*x2^2*x3 + x0*x2^3*x3 + x1*x2^3*x3
        sage: M[1,3,1].expand(2)
        0

    TESTS::

        sage: QSym = QuasiSymmetricFunctions(QQ); QSym
        Quasisymmetric functions over the Rational Field
        sage: QSym.base_ring()
        Rational Field
        sage: algebras.QSym(QQ) is QSym
        True
    """

    def __init__(self, R):
        """
        The Hopf algebra of quasi-symmetric functions.
        See ``QuasiSymmetricFunctions`` for full documentation.

        EXAMPLES::

            sage: QuasiSymmetricFunctions(QQ)
            Quasisymmetric functions over the Rational Field
            sage: QSym1 = QuasiSymmetricFunctions(FiniteField(23))
            sage: QSym2 = QuasiSymmetricFunctions(Integers(23))
            sage: TestSuite(QuasiSymmetricFunctions(QQ)).run()
        """
        # change the line below to assert R in Rings() once MRO issues from #15536, #15475 are resolved
        assert R in Fields() or R in Rings()  # side effect of this statement assures MRO exists for R
        self._base = R # Won't be needed once CategoryObject won't override base_ring
        category = GradedHopfAlgebras(R).Commutative()
        self._category = category
        Parent.__init__(self, category=category.WithRealizations())

        # Bases
        Monomial = self.Monomial()
        Fundamental = self.Fundamental()
        dualImmaculate = self.dualImmaculate()
        QS = self.Quasisymmetric_Schur()

        # Change of bases
        Fundamental.module_morphism(Monomial.sum_of_finer_compositions,
                                    codomain=Monomial, category=category
                                    ).register_as_coercion()
        Monomial   .module_morphism(Fundamental.alternating_sum_of_finer_compositions,
                                    codomain=Fundamental, category=category
                                    ).register_as_coercion()
        #This changes dualImmaculate into Monomial
        dualImmaculate.module_morphism(dualImmaculate._to_Monomial_on_basis,
                                          codomain=Monomial, category=category
                                          ).register_as_coercion()
        #This changes Monomial into dualImmaculate
        Monomial.module_morphism(dualImmaculate._from_Monomial_on_basis,
                                          codomain=dualImmaculate, category=category
                                          ).register_as_coercion()
        #This changes Quasisymmetric Schur into Monomial
        QS         .module_morphism(QS._to_monomial_on_basis,
                                    codomain=Monomial, category=category
                                    ).register_as_coercion()
        #This changes Monomial into Quasisymmetric Schur
        Monomial.module_morphism(QS._from_monomial_on_basis,
                                    codomain=QS, category=category
                                    ).register_as_coercion()

        # Embedding of Sym into QSym in the monomial bases
        Sym = SymmetricFunctions(self.base_ring())
        Sym_m_to_M = Sym.m().module_morphism(Monomial.sum_of_partition_rearrangements,
                                           triangular='upper', inverse_on_support=Monomial._comp_to_par,
                                           codomain=Monomial, category=category)
        Sym_m_to_M.register_as_coercion()
        self.to_symmetric_function = Sym_m_to_M.section()

        Sym_s_to_F = Sym.s().module_morphism(Fundamental._from_schur_on_basis,
                                             unitriangular='upper',
                                             codomain=Fundamental, category=category)
        Sym_s_to_F.register_as_coercion()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M = QuasiSymmetricFunctions(ZZ).M()
            sage: M._repr_()
            'Quasisymmetric functions over the Integer Ring in the Monomial basis'
        """
        return "Quasisymmetric functions over the %s" % self.base_ring()

    def a_realization(self):
        r"""
        Return the realization of the Monomial basis of the ring of quasi-symmetric functions.

        OUTPUT: the Monomial basis of quasi-symmetric functions

        EXAMPLES::

            sage: QuasiSymmetricFunctions(QQ).a_realization()
            Quasisymmetric functions over the Rational Field in the Monomial basis
        """
        return self.Monomial()

    _shorthands = tuple(['M', 'F', 'E', 'dI', 'QS', 'YQS', 'phi', 'psi'])

    def dual(self):
        r"""
        Return the dual Hopf algebra of the quasi-symmetric functions, which is the
        non-commutative symmetric functions.

        OUTPUT: the non-commutative symmetric functions

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QSym.dual()
            Non-Commutative Symmetric Functions over the Rational Field
        """
        return NonCommutativeSymmetricFunctions(self.base_ring())

    def from_polynomial(self, f, check=True):
        """
        Return the quasi-symmetric function in the Monomial basis
        corresponding to the quasi-symmetric polynomial ``f``.

        INPUT:

        - ``f`` -- a polynomial in finitely many variables over the same base
          ring as ``self``. It is assumed that this polynomial is
          quasi-symmetric.
        - ``check`` -- boolean (default: ``True``); checks whether the
          polynomial is indeed quasi-symmetric

        OUTPUT: quasi-symmetric function in the Monomial basis

        EXAMPLES::

            sage: P = PolynomialRing(QQ, 'x', 3)
            sage: x = P.gens()
            sage: f = x[0] + x[1] + x[2]
            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QSym.from_polynomial(f)
            M[1]

        Beware of setting ``check=False``::

            sage: f = x[0] + 2*x[1] + x[2]
            sage: QSym.from_polynomial(f, check=True)
            Traceback (most recent call last):
            ...
            ValueError: x0 + 2*x1 + x2 is not a quasi-symmetric polynomial
            sage: QSym.from_polynomial(f, check=False)
            M[1]

        To expand the quasi-symmetric function in a basis other than the
        Monomial basis, the following shorthands are provided::

            sage: M = QSym.Monomial()
            sage: f = x[0]**2+x[1]**2+x[2]**2
            sage: g = M.from_polynomial(f); g
            M[2]
            sage: F = QSym.Fundamental()
            sage: F(g)
            -F[1, 1] + F[2]
            sage: F.from_polynomial(f)
            -F[1, 1] + F[2]
        """
        assert self.base_ring() == f.base_ring()
        exponent_coefficient = f.monomial_coefficients()
        z = {}
        for e, c in exponent_coefficient.items():
            I = Compositions()([ei for ei in e if ei])
            if I not in z:
                z[I] = c
        out = self.Monomial()._from_dict(z)
        if check and out.expand(f.parent().ngens(), f.parent().variable_names()) != f:
            raise ValueError("%s is not a quasi-symmetric polynomial" % f)
        return out

    class Bases(Category_realization_of_parent):
        r"""
        Category of bases of quasi-symmetric functions.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QSym.Bases()
            Category of bases of Quasisymmetric functions over the Rational Field
        """

        def super_categories(self):
            r"""
            Return the super categories of bases of the Quasi-symmetric functions.

            OUTPUT: list of categories

            TESTS::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QSym.Bases().super_categories()
                [Category of commutative bases of Non-Commutative Symmetric Functions or Quasisymmetric functions over the Rational Field]
            """
            return [BasesOfQSymOrNCSF(self.base()).Commutative()]

        class ParentMethods:
            r"""
            Methods common to all bases of ``QuasiSymmetricFunctions``.
            """

            def from_polynomial(self, f, check=True):
                """
                The quasi-symmetric function expanded in this basis
                corresponding to the quasi-symmetric polynomial ``f``.

                This is a default implementation that computes
                the expansion in the Monomial basis and converts
                to this basis.

                INPUT:

                - ``f`` -- a polynomial in finitely many variables over the same base
                  ring as ``self``. It is assumed that this polynomial is
                  quasi-symmetric.
                - ``check`` -- boolean (default: ``True``); checks whether the
                  polynomial is indeed quasi-symmetric

                OUTPUT: quasi-symmetric function

                EXAMPLES::

                    sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                    sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                    sage: P = PolynomialRing(QQ, 'x', 3)
                    sage: x = P.gens()
                    sage: f = x[0] + x[1] + x[2]
                    sage: M.from_polynomial(f)
                    M[1]
                    sage: F.from_polynomial(f)
                    F[1]
                    sage: f = x[0]**2+x[1]**2+x[2]**2
                    sage: M.from_polynomial(f)
                    M[2]
                    sage: F.from_polynomial(f)
                    -F[1, 1] + F[2]

                If the polynomial is not quasi-symmetric, an error
                is raised::

                    sage: f = x[0]^2+x[1]
                    sage: M.from_polynomial(f)
                    Traceback (most recent call last):
                    ...
                    ValueError: x0^2 + x1 is not a quasi-symmetric polynomial
                    sage: F.from_polynomial(f)
                    Traceback (most recent call last):
                    ...
                    ValueError: x0^2 + x1 is not a quasi-symmetric polynomial

                TESTS:

                We convert some quasi-symmetric functions to quasi-symmetric
                polynomials and back::

                    sage: f = (M[1,2] + M[1,1]).expand(3); f
                    x0*x1^2 + x0*x2^2 + x1*x2^2 + x0*x1 + x0*x2 + x1*x2
                    sage: M.from_polynomial(f)
                    M[1, 1] + M[1, 2]
                    sage: f = (2*M[2,1]+M[1,1]+3*M[3]).expand(3)
                    sage: M.from_polynomial(f)
                    M[1, 1] + 2*M[2, 1] + 3*M[3]
                    sage: f = (F[1,2] + F[1,1]).expand(3); f
                    x0*x1^2 + x0*x1*x2 + x0*x2^2 + x1*x2^2 + x0*x1 + x0*x2 + x1*x2
                    sage: F.from_polynomial(f)
                    F[1, 1] + F[1, 2]
                    sage: f = (2*F[2,1]+F[1,1]+3*F[3]).expand(3)
                    sage: F.from_polynomial(f)
                    F[1, 1] + 2*F[2, 1] + 3*F[3]
                """
                g = self.realization_of().from_polynomial(f, check=check)
                return self(g)

            def Eulerian(self, n, j, k=None):
                """
                Return the Eulerian (quasi)symmetric function `Q_{n,j}` in
                terms of ``self``.

                INPUT:

                - ``n`` -- the value `n` or a partition
                - ``j`` -- the number of excedances
                - ``k`` -- (optional) if specified, determines the number of
                  fixed points of the permutation

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.M()
                    sage: M.Eulerian(3, 1)
                    4*M[1, 1, 1] + 3*M[1, 2] + 3*M[2, 1] + 2*M[3]
                    sage: M.Eulerian(4, 1, 2)
                    6*M[1, 1, 1, 1] + 4*M[1, 1, 2] + 4*M[1, 2, 1]
                     + 2*M[1, 3] + 4*M[2, 1, 1] + 3*M[2, 2] + 2*M[3, 1] + M[4]
                    sage: QS = QSym.QS()
                    sage: QS.Eulerian(4, 2)
                    2*QS[1, 3] + QS[2, 2] + 2*QS[3, 1] + 3*QS[4]
                    sage: QS.Eulerian([2, 2, 1], 2)
                    QS[1, 2, 2] + QS[1, 4] + QS[2, 1, 2] + QS[2, 2, 1]
                     + QS[2, 3] + QS[3, 2] + QS[4, 1] + QS[5]
                    sage: dI = QSym.dI()
                    sage: dI.Eulerian(5, 2)
                    -dI[1, 3, 1] - 5*dI[1, 4] + dI[2, 2, 1] + dI[3, 1, 1]
                     + 5*dI[3, 2] + 6*dI[4, 1] + 6*dI[5]
                """
                F = self.realization_of().F()
                if n in _Partitions:
                    n = _Partitions(n)
                return self(F.Eulerian(n, j, k))

        class ElementMethods:
            r"""
            Methods common to all elements of ``QuasiSymmetricFunctions``.
            """

            def internal_coproduct(self):
                r"""
                Return the inner coproduct of ``self`` in the basis of ``self``.

                The inner coproduct (also known as the Kronecker coproduct,
                or as the second comultiplication on the `R`-algebra of
                quasi-symmetric functions) is an `R`-algebra homomorphism
                `\Delta^{\times}` from the `R`-algebra of quasi-symmetric
                functions to the tensor square (over `R`) of quasi-symmetric
                functions. It can be defined in the following two ways:

                #. If `I` is a composition, then a `(0, I)`-matrix will mean a
                   matrix whose entries are nonnegative integers such that no
                   row and no column of this matrix is zero, and such that if
                   all the nonzero entries of the matrix are read (row by row,
                   starting at the topmost row, reading every row from left to
                   right), then the reading word obtained is `I`. If `A` is
                   a `(0, I)`-matrix, then `\mathrm{row}(A)` will denote the
                   vector of row sums of `A` (regarded as a composition), and
                   `\mathrm{column}(A)` will denote the vector of column sums
                   of `A` (regarded as a composition).

                   For every composition `I`, the internal coproduct
                   `\Delta^{\times}(M_I)` of the `I`-th monomial quasisymmetric
                   function `M_I` is the sum

                   .. MATH::

                       \sum_{A \hbox{ is a } (0, I) \text{-matrix}}
                       M_{\mathrm{row}(A)} \otimes M_{\mathrm{column}(A)}.

                   See Section 11.39 of [HazWitt1]_.

                #. For every permutation `w`, let `C(w)` denote the descent
                   composition of `w`. Then, for any composition `I` of size
                   `n`, the internal coproduct `\Delta^{\times}(F_I)` of the
                   `I`-th fundamental quasisymmetric function `F_I` is the sum

                   .. MATH::

                       \sum_{\substack{\sigma \in S_n,\\ \tau \in S_n,\\
                       \tau \sigma = \pi}} F_{C(\sigma)} \otimes F_{C(\tau)},

                   where `\pi` is any permutation in `S_n` having descent
                   composition `I` and where permutations act from the left and
                   multiply accordingly, so `\tau \sigma` means first applying
                   `\sigma` and then `\tau`. See Theorem 4.23 in [Mal1993]_,
                   but beware of the notations which are apparently different
                   from those in [HazWitt1]_.

                The restriction of the internal coproduct to the
                `R`-algebra of symmetric functions is the well-known
                internal coproduct on the symmetric functions.

                The method :meth:`kronecker_coproduct` is a synonym of this one.

                EXAMPLES:

                Let us compute the internal coproduct of `M_{21}` (which is
                short for `M_{[2, 1]}`). The `(0, [2,1])`-matrices are

                .. MATH::

                    \begin{bmatrix} 2 & 1 \end{bmatrix},
                    \begin{bmatrix} 2 \\ 1 \end{bmatrix},
                    \begin{bmatrix} 2 & 0 \\ 0 & 1 \end{bmatrix}, \hbox{ and }
                    \begin{bmatrix} 0 & 2 \\ 1 & 0 \end{bmatrix}

                so

                .. MATH::

                    \Delta^\times(M_{21}) = M_{3} \otimes M_{21} +
                    M_{21} \otimes M_3 + M_{21} \otimes M_{21} +
                    M_{21} \otimes M_{12}.

                This is confirmed by the following Sage computation
                (incidentally demonstrating the non-cocommutativity of
                the internal coproduct)::

                    sage: M = QuasiSymmetricFunctions(ZZ).M()
                    sage: a = M([2,1])
                    sage: a.internal_coproduct()
                    M[2, 1] # M[1, 2] + M[2, 1] # M[2, 1] + M[2, 1] # M[3] + M[3] # M[2, 1]

                Further examples::

                    sage: all( M([i]).internal_coproduct() == tensor([M([i]), M([i])])
                    ....:      for i in range(1, 4) )
                    True

                    sage: M([1, 2]).internal_coproduct()
                    M[1, 2] # M[1, 2] + M[1, 2] # M[2, 1] + M[1, 2] # M[3] + M[3] # M[1, 2]

                The definition of `\Delta^{\times}(M_I)` in terms of
                `(0, I)`-matrices is not suitable for computation in
                cases where the length of `I` is large, but we can use
                it as a doctest. Here is a naive implementation::

                    sage: def naive_internal_coproduct_on_M(I):
                    ....:     # INPUT: composition I
                    ....:     #        (not quasi-symmetric function)
                    ....:     # OUTPUT: interior coproduct of M_I
                    ....:     M = QuasiSymmetricFunctions(ZZ).M()
                    ....:     M2 = M.tensor(M)
                    ....:     res = M2.zero()
                    ....:     l = len(I)
                    ....:     n = I.size()
                    ....:     for S in Subsets(range(l**2), l):
                    ....:         M_list = sorted(S)
                    ....:         row_M = [sum([I[M_list.index(l * i + j)]
                    ....:                       for j in range(l) if
                    ....:                       l * i + j in S])
                    ....:                  for i in range(l)]
                    ....:         col_M = [sum([I[M_list.index(l * i + j)]
                    ....:                       for i in range(l) if
                    ....:                       l * i + j in S])
                    ....:                  for j in range(l)]
                    ....:         if 0 in row_M:
                    ....:             first_zero = row_M.index(0)
                    ....:             row_M = row_M[:first_zero]
                    ....:             if sum(row_M) != n:
                    ....:                 continue
                    ....:         if 0 in col_M:
                    ....:             first_zero = col_M.index(0)
                    ....:             col_M = col_M[:first_zero]
                    ....:             if sum(col_M) != n:
                    ....:                 continue
                    ....:         res += tensor([M(Compositions(n)(row_M)),
                    ....:                        M(Compositions(n)(col_M))])
                    ....:     return res
                    sage: all( naive_internal_coproduct_on_M(I)
                    ....:      == M(I).internal_coproduct()
                    ....:      for I in Compositions(3) )
                    True

                TESTS:

                Border cases::

                    sage: M = QuasiSymmetricFunctions(ZZ).M()
                    sage: F = QuasiSymmetricFunctions(ZZ).F()
                    sage: M([]).internal_coproduct()
                    M[] # M[]
                    sage: F([]).internal_coproduct()
                    F[] # F[]

                The implementations on the ``F`` and ``M`` bases agree
                with each other::

                    sage: M = QuasiSymmetricFunctions(ZZ).M()
                    sage: F = QuasiSymmetricFunctions(ZZ).F()
                    sage: def int_copr_on_F_via_M(I):
                    ....:     result = tensor([F.zero(), F.zero()])
                    ....:     w = M(F(I)).internal_coproduct()
                    ....:     for lam, a in w:
                    ....:         (U, V) = lam
                    ....:         result += a * tensor([F(M(U)), F(M(V))])
                    ....:     return result
                    sage: all( int_copr_on_F_via_M(I) == F(I).internal_coproduct()
                    ....:      for I in Compositions(3) )
                    True
                    sage: all( int_copr_on_F_via_M(I) == F(I).internal_coproduct()
                    ....:      for I in Compositions(4) )
                    True

                Restricting to the subring of symmetric functions gives the
                standard internal coproduct on the latter::

                    sage: M = QuasiSymmetricFunctions(ZZ).M()
                    sage: e = SymmetricFunctions(ZZ).e()
                    sage: def int_copr_of_e_in_M(mu):
                    ....:     result = tensor([M.zero(), M.zero()])
                    ....:     w = e(mu).internal_coproduct()
                    ....:     for lam, a in w:
                    ....:         (nu, kappa) = lam
                    ....:         result += a * tensor([M(e(nu)), M(e(kappa))])
                    ....:     return result
                    sage: all( int_copr_of_e_in_M(mu) == M(e(mu)).internal_coproduct()
                    ....:      for mu in Partitions(3) )
                    True
                    sage: all( int_copr_of_e_in_M(mu) == M(e(mu)).internal_coproduct()
                    ....:      for mu in Partitions(4) )
                    True

                .. TODO::

                    Implement this directly on the monomial basis maybe?
                    The `(0, I)`-matrices are a pain to generate from their
                    definition, but maybe there is a good algorithm.
                    If so, the above "further examples" should be moved
                    to the M-method.
                """
                # Coerce to F basis and back, doing the actual computation in that basis.
                parent = self.parent()
                F = parent.realization_of().F()
                from sage.categories.tensor import tensor
                result = tensor([parent.zero(), parent.zero()])
                for lam, a in F(self).internal_coproduct():
                    (I, J) = lam
                    result += a * tensor([parent(F(I)), parent(F(J))])
                return result

            kronecker_coproduct = internal_coproduct

            def adams_operator(self, n):
                r"""
                Return the image of the quasi-symmetric function ``self``
                under the `n`-th Adams operator.

                The `n`-th Adams operator `\mathbf{f}_n` is defined to be
                the map from the `R`-algebra of quasi-symmetric functions
                to itself that sends every symmetric function
                `P(x_1, x_2, x_3, \ldots)` to
                `P(x_1^n, x_2^n, x_3^n, \ldots)`. This operator `\mathbf{f}_n`
                is a Hopf algebra endomorphism, and satisfies

                .. MATH::

                    f_n M_{(i_1, i_2, i_3, \ldots)} =
                    M_{(ni_1, ni_2, ni_3, \ldots)}

                for every composition `(i_1, i_2, i_3, \ldots)`
                (where `M` means the monomial basis).

                The `n`-th Adams operator is also called the `n`-th
                Frobenius endomorphism. It is not related to the Frobenius map
                which connects the ring of symmetric functions with the
                representation theory of the symmetric group.

                The `n`-th Adams operator is the `n`-th Adams operator
                of the `\Lambda`-ring of quasi-symmetric functions over the
                integers.

                The restriction of the `n`-th Adams operator to the
                subring formed by all symmetric functions is, not
                unexpectedly, the `n`-th Adams operator of the ring of
                symmetric functions.

                .. SEEALSO::

                    :meth:`Symmetric functions plethysm
                    <sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.plethysm>`

                INPUT:

                - ``n`` -- positive integer

                OUTPUT:

                The result of applying the `n`-th Adams operator (on the
                ring of quasi-symmetric functions) to ``self``.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: M = QSym.M()
                    sage: F = QSym.F()
                    sage: M[3,2].adams_operator(2)
                    M[6, 4]
                    sage: (M[2,1] - 2*M[3]).adams_operator(4)
                    M[8, 4] - 2*M[12]
                    sage: M([]).adams_operator(3)
                    M[]
                    sage: F[1,1].adams_operator(2)
                    F[1, 1, 1, 1] - F[1, 1, 2] - F[2, 1, 1] + F[2, 2]

                The Adams endomorphisms are multiplicative::

                    sage: all( all( M(I).adams_operator(3) * M(J).adams_operator(3)
                    ....:           == (M(I) * M(J)).adams_operator(3)
                    ....:           for I in Compositions(3) )
                    ....:      for J in Compositions(2) )
                    True

                Being Hopf algebra endomorphisms, the Adams operators
                commute with the antipode::

                    sage: all( M(I).adams_operator(4).antipode()
                    ....:      == M(I).antipode().adams_operator(4)
                    ....:      for I in Compositions(3) )
                    True

                The restriction of the Adams operators to the subring
                of symmetric functions are the Adams operators of
                the latter::

                    sage: e = SymmetricFunctions(ZZ).e()
                    sage: all( M(e(lam)).adams_operator(3)
                    ....:      == M(e(lam).adams_operator(3))
                    ....:      for lam in Partitions(3) )
                    True
                """
                # Convert to the monomial basis, there apply componentwise,
                # then convert back.
                parent = self.parent()
                M = parent.realization_of().M()
                C = parent._indices
                dct = {C([n * i for i in I]): coeff
                       for (I, coeff) in M(self)}
                result_in_M_basis = M._from_dict(dct)
                return parent(result_in_M_basis)

            frobenius = deprecated_function_alias(36396, adams_operator)

            def star_involution(self):
                r"""
                Return the image of the quasisymmetric function ``self`` under
                the star involution.

                The star involution is defined as the linear map
                `QSym \to QSym` which, for every composition `I`, sends the
                monomial quasisymmetric function `M_I` to `M_{I^r}`. Here, if
                `I` is a composition, we denote by `I^r` the reversed
                composition of `I`. Denoting by `f^{\ast}` the image of an
                element `f \in QSym` under the star involution, it can be shown
                that every composition `I` satisfies

                .. MATH::

                    (M_I)^{\ast} = M_{I^r}, \quad (F_I)^{\ast} = F_{I^r},

                where `F_I` denotes the fundamental quasisymmetric function
                corresponding to the composition `I`. The star involution is an
                involution, an algebra automorphism and a coalgebra
                anti-automorphism of `QSym`. It also is an automorphism of the
                graded vector space `QSym`, and is the identity on the subspace
                `Sym` of `QSym`. It is adjoint to the star involution on `NCSF`
                by the standard adjunction between `NCSF` and `QSym`.

                The star involution has been denoted by `\rho` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`star involution on NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: M = QSym.M()
                    sage: M[3,2].star_involution()
                    M[2, 3]
                    sage: M[6,3].star_involution()
                    M[3, 6]
                    sage: (M[9,1] - M[6,2] + 2*M[6,4] - 3*M[3] + 4*M[[]]).star_involution()
                    4*M[] + M[1, 9] - M[2, 6] - 3*M[3] + 2*M[4, 6]
                    sage: (M[3,3] - 2*M[2]).star_involution()
                    -2*M[2] + M[3, 3]
                    sage: M([4,2]).star_involution()
                    M[2, 4]
                    sage: dI = QSym.dI()
                    sage: dI([1,2]).star_involution()
                    -dI[1, 2] + dI[2, 1]
                    sage: dI.zero().star_involution()
                    0

                The star involution commutes with the antipode::

                    sage: all( M(I).star_involution().antipode()
                    ....:      == M(I).antipode().star_involution()
                    ....:      for I in Compositions(4) )
                    True

                The star involution is the identity on `Sym`::

                    sage: Sym = SymmetricFunctions(ZZ)
                    sage: e = Sym.e()
                    sage: all( M(e(lam)).star_involution() == M(e(lam))
                    ....:      for lam in Partitions(4) )
                    True
                """
                # Convert to the homogeneous basis, there apply the star
                # involution componentwise, then convert back.
                parent = self.parent()
                M = parent.realization_of().M()
                dct = {I.reversed(): coeff for (I, coeff) in M(self)}
                return parent(M._from_dict(dct))

            def omega_involution(self):
                r"""
                Return the image of the quasisymmetric function
                ``self`` under the omega involution.

                The omega involution is defined as the linear map
                `QSym \to QSym` which, for every composition `I`, sends
                the fundamental quasisymmetric function `F_I` to
                `F_{I^t}`, where `I^t` denotes the conjugate
                (:meth:`~sage.combinat.composition.Composition.conjugate`)
                of the composition `I`. This map is commonly denoted by
                `\omega`. It is an algebra homomorphism and a coalgebra
                antihomomorphism; it also is an involution and an
                automorphism of the graded vector space `QSym`. Also,
                every composition `I` satisfies

                .. MATH::

                    \omega(M_I) = (-1)^{|I|-\ell(I)} \sum M_J,

                where the sum ranges over all compositions `J` of `|I|`
                which are coarser than the reversed composition `I^r` of
                `I`. Here, `\ell(I)` denotes the length of the composition
                `I` (that is, the number of parts of `I`).

                If `f` is a homogeneous element of `NCSF` of degree `n`,
                then

                .. MATH::

                    \omega(f) = (-1)^n S(f),

                where `S` denotes the antipode of `QSym`.

                The restriction of `\omega` to the ring of symmetric
                functions (which is a subring of `QSym`) is precisely the
                omega involution
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.omega`)
                of said ring.

                The omega involution on `QSym` is adjoint to the omega
                involution on `NCSF` by the standard adjunction between `NCSF`
                and `QSym`.

                The omega involution has been denoted by `\omega` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`omega involution on NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.omega_involution>`,
                    :meth:`psi involution on QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution on QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: F = QSym.F()
                    sage: F[3,2].omega_involution()
                    F[1, 2, 1, 1]
                    sage: F[6,3].omega_involution()
                    F[1, 1, 2, 1, 1, 1, 1, 1]
                    sage: (F[9,1] - F[8,2] + 2*F[2,4] - 3*F[3] + 4*F[[]]).omega_involution()
                    4*F[] - 3*F[1, 1, 1] + 2*F[1, 1, 1, 2, 1] - F[1, 2, 1, 1, 1, 1, 1, 1, 1] + F[2, 1, 1, 1, 1, 1, 1, 1, 1]
                    sage: (F[3,3] - 2*F[2]).omega_involution()
                    -2*F[1, 1] + F[1, 1, 2, 1, 1]
                    sage: F([2,1,1]).omega_involution()
                    F[3, 1]
                    sage: M = QSym.M()
                    sage: M([2,1]).omega_involution()
                    -M[1, 2] - M[3]
                    sage: M.zero().omega_involution()
                    0

                Testing the fact that the restriction of `\omega` to `Sym`
                is the omega automorphism of `Sym`::

                    sage: Sym = SymmetricFunctions(ZZ)
                    sage: e = Sym.e()
                    sage: all( F(e[lam]).omega_involution() == F(e[lam].omega())
                    ....:      for lam in Partitions(4) )
                    True
                """
                # Use the `\omega(f) = (-1)^n S(f)` formula.
                return self.antipode().degree_negation()

            def psi_involution(self):
                r"""
                Return the image of the quasisymmetric function ``self``
                under the involution `\psi`.

                The involution `\psi` is defined as the linear map
                `QSym \to QSym` which, for every composition `I`, sends the
                fundamental quasisymmetric function `F_I` to `F_{I^c}`, where
                `I^c` denotes the complement of the composition `I`.
                The map `\psi` is an involution and a graded Hopf algebra
                automorphism of `QSym`. Its restriction to the ring of
                symmetric functions coincides with the omega automorphism of
                the latter ring.

                The involution `\psi` of `QSym` is adjoint to the involution
                `\psi` of `NCSF` by the standard adjunction between `NCSF` and
                `QSym`.

                The involution `\psi` has been denoted by `\psi` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`psi involution on NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution on QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: F = QSym.F()
                    sage: F[3,2].psi_involution()
                    F[1, 1, 2, 1]
                    sage: F[6,3].psi_involution()
                    F[1, 1, 1, 1, 1, 2, 1, 1]
                    sage: (F[9,1] - F[8,2] + 2*F[2,4] - 3*F[3] + 4*F[[]]).psi_involution()
                    4*F[] - 3*F[1, 1, 1] + F[1, 1, 1, 1, 1, 1, 1, 1, 2] - F[1, 1, 1, 1, 1, 1, 1, 2, 1] + 2*F[1, 2, 1, 1, 1]
                    sage: (F[3,3] - 2*F[2]).psi_involution()
                    -2*F[1, 1] + F[1, 1, 2, 1, 1]
                    sage: F([2,1,1]).psi_involution()
                    F[1, 3]
                    sage: M = QSym.M()
                    sage: M([2,1]).psi_involution()
                    -M[2, 1] - M[3]
                    sage: M.zero().psi_involution()
                    0

                The involution `\psi` commutes with the antipode::

                    sage: all( F(I).psi_involution().antipode()
                    ....:      == F(I).antipode().psi_involution()
                    ....:      for I in Compositions(4) )
                    True

                Testing the fact that the restriction of `\psi` to `Sym`
                is the omega automorphism of `Sym`::

                    sage: Sym = SymmetricFunctions(ZZ)
                    sage: e = Sym.e()
                    sage: all( F(e[lam]).psi_involution() == F(e[lam].omega())
                    ....:      for lam in Partitions(4) )
                    True
                """
                # Convert to the fundamental basis, there apply the psi
                # involution componentwise, then convert back.
                parent = self.parent()
                F = parent.realization_of().F()
                dct = {I.complement(): coeff for (I, coeff) in F(self)}
                return parent(F._from_dict(dct))

            def dendriform_less(self, other):
                r"""
                Return the result of applying the dendriform smaller
                operation to the two quasi-symmetric functions ``self``
                and ``other``.

                The dendriform smaller operation is a binary operation,
                denoted by `\prec` and written infix, on the ring of
                quasi-symmetric functions. It can be defined as a
                restriction of a binary operation (denoted by `\prec`
                and written infix as well) on the ring of formal power
                series `R[[x_1, x_2, x_3, \ldots]]`, which is defined
                as follows: If `m` and `n` are two monomials in
                `x_1, x_2, x_3, \ldots`, then we let `m \prec n` be
                the product `mn` if the smallest positive integer `i`
                for which `x_i` occurs in `m` is smaller than the
                smallest positive integer `j` for which `x_j` occurs
                in `n` (this is understood to be false when `m = 1`,
                and true when `m \neq 1` and `n = 1`), and we let
                `m \prec n` be `0` otherwise. Having thus defined
                `\prec` on monomials, we extend `\prec` to a binary
                operation on `R[[x_1, x_2, x_3, \ldots]]` by requiring
                it to be continuous (in both inputs) and `R`-bilinear.
                It is easily seen that `QSym \prec QSym \subseteq
                QSym`, so that `\prec` restricts to a binary operation
                on `QSym`.

                .. SEEALSO::

                    :meth:`dendriform_leq`

                INPUT:

                - ``other`` -- a quasi-symmetric function over the
                  same base ring as ``self``

                OUTPUT:

                The quasi-symmetric function ``self`` `\prec`
                ``other``, written in the basis of ``self``.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.M()
                    sage: M[2, 1].dendriform_less(M[1, 2])
                    2*M[2, 1, 1, 2] + M[2, 1, 2, 1] + M[2, 1, 3] + M[2, 2, 2]
                    sage: F = QSym.F()
                    sage: F[2, 1].dendriform_less(F[1, 2])
                    F[1, 1, 2, 1, 1] + F[1, 1, 2, 2] + F[1, 1, 3, 1]
                     + F[1, 2, 1, 2] + F[1, 2, 2, 1] + F[1, 2, 3]
                     + F[2, 1, 1, 2] + F[2, 1, 2, 1] + F[2, 1, 3] + F[2, 2, 2]

                The operation `\prec` can be used to recursively
                construct the dual immaculate basis: For every positive
                integer `m` and every composition `I`, the dual
                immaculate function `\operatorname{dI}_{[m, I]}` of the
                composition `[m, I]` (this composition is `I` with `m`
                prepended to it) is `F_{[m]} \prec \operatorname{dI}_I`. ::

                    sage: dI = QSym.dI()
                    sage: dI(F[2]).dendriform_less(dI[1, 2])
                    dI[2, 1, 2]

                We check with the identity element::

                    sage: M.one().dendriform_less(M[1,2])
                    0
                    sage: M[1,2].dendriform_less(M.one())
                    M[1, 2]

                The operation `\prec` is not symmetric, nor if
                `a \prec b \neq 0`, then `b \prec a = 0` (as it would be
                for a single monomial)::

                    sage: M[1,2,1].dendriform_less(M[1,2])
                    M[1, 1, 2, 1, 2] + 2*M[1, 1, 2, 2, 1] + M[1, 1, 2, 3]
                     + M[1, 1, 4, 1] + 2*M[1, 2, 1, 1, 2] + M[1, 2, 1, 2, 1]
                     + M[1, 2, 1, 3] + M[1, 2, 2, 2] + M[1, 3, 1, 2]
                     + M[1, 3, 2, 1] + M[1, 3, 3]
                    sage: M[1,2].dendriform_less(M[1,2,1])
                    M[1, 1, 2, 1, 2] + 2*M[1, 1, 2, 2, 1] + M[1, 1, 2, 3]
                     + M[1, 1, 4, 1] + M[1, 2, 1, 2, 1] + M[1, 3, 2, 1]
                """
                # Convert to the monomial basis, there do restricted
                # shuffle product, then convert back to self.parent().
                P = self.parent()
                M = P.realization_of().M()
                a = M(self)
                b = M(other)
                res = M.zero()
                for I, I_coeff in a:
                    if not I:
                        continue
                    i_head = I[0]
                    I_tail = Composition(I[1:])
                    for J, J_coeff in b:
                        shufpro = I_tail.shuffle_product(J, overlap=True)
                        res += J_coeff * M.sum_of_monomials(Composition([i_head] + list(K))
                                                            for K in shufpro)
                return P(res)

            def dendriform_leq(self, other):
                r"""
                Return the result of applying the dendriform
                smaller-or-equal operation to the two quasi-symmetric
                functions ``self`` and ``other``.

                The dendriform smaller-or-equal operation is a binary
                operation, denoted by `\preceq` and written infix, on
                the ring of quasi-symmetric functions. It can be
                defined as a restriction of a binary operation
                (denoted by `\preceq` and written infix as well) on
                the ring of formal power series
                `R[[x_1, x_2, x_3, \ldots]]`, which is defined as
                follows: If `m` and `n` are two monomials in
                `x_1, x_2, x_3, \ldots`, then we let `m \preceq n` be
                the product `mn` if the smallest positive integer `i`
                for which `x_i` occurs in `m` is smaller or equal to
                the smallest positive integer `j` for which `x_j`
                occurs in `n` (this is understood to be false when
                `m = 1` and `n \neq 1`, and true when `n = 1`), and we
                let `m \preceq n` be `0` otherwise. Having thus
                defined `\preceq` on monomials, we extend `\preceq` to
                a binary operation on `R[[x_1, x_2, x_3, \ldots]]` by
                requiring it to be continuous (in both inputs) and
                `R`-bilinear. It is easily seen that
                `QSym \preceq QSym \subseteq QSym`, so that `\preceq`
                restricts to a binary operation on `QSym`.

                This operation `\preceq` is related to the dendriform
                smaller relation `\prec` (:meth:`dendriform_less`).
                Namely, if we define a binary operation `\succ` on
                `QSym` by `a \succ b = b \prec a`, then
                `(QSym, \preceq, \succ)` is a dendriform `R`-algebra.
                Thus, any `a, b \in QSym` satisfy
                `a \preceq b = ab - b \prec a`.

                .. SEEALSO::

                    :meth:`dendriform_less`

                INPUT:

                - ``other`` -- a quasi-symmetric function over the
                  same base ring as ``self``

                OUTPUT:

                The quasi-symmetric function ``self`` `\preceq`
                ``other``, written in the basis of ``self``.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.M()
                    sage: M[2, 1].dendriform_leq(M[1, 2])
                    2*M[2, 1, 1, 2] + M[2, 1, 2, 1] + M[2, 1, 3] + M[2, 2, 2]
                     + M[3, 1, 2] + M[3, 2, 1] + M[3, 3]
                    sage: F = QSym.F()
                    sage: F[2, 1].dendriform_leq(F[1, 2])
                    F[2, 1, 1, 2] + F[2, 1, 2, 1] + F[2, 1, 3] + F[2, 2, 1, 1]
                     + 2*F[2, 2, 2] + F[2, 3, 1] + F[3, 1, 2] + F[3, 2, 1] + F[3, 3]
                """
                # This might be somewhat slow...
                P = self.parent()
                return self * P(other) - P(other.dendriform_less(self))

            def expand(self, n, alphabet='x'):
                r"""
                Expand the quasi-symmetric function into ``n`` variables in
                an alphabet, which by default is ``'x'``.

                INPUT:

                - ``n`` -- nonnegative integer; the number of variables
                  in the expansion
                - ``alphabet`` -- (default: ``'x'``) the alphabet in
                  which ``self`` is to be expanded

                OUTPUT: an expansion of ``self`` into the ``n`` variables
                specified by ``alphabet``

                EXAMPLES::

                    sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                    sage: F[3].expand(3)
                    x0^3 + x0^2*x1 + x0*x1^2 + x1^3 + x0^2*x2 + x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x2^3
                    sage: F[2,1].expand(3)
                    x0^2*x1 + x0^2*x2 + x0*x1*x2 + x1^2*x2

                One can use a different set of variable by adding an optional
                argument ``alphabet=...`` ::

                    sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                    sage: F[3].expand(2,alphabet='y')
                    y0^3 + y0^2*y1 + y0*y1^2 + y1^3

                TESTS::

                    sage: (3*F([])).expand(2)
                    3
                    sage: F[4,2].expand(0)
                    0
                    sage: F([]).expand(0)
                    1
                """
                M = self.parent().realization_of().Monomial()
                return M(self).expand(n,alphabet)

            def is_symmetric(self):
                r"""
                Return ``True`` if ``self`` is an element of the symmetric
                functions.

                This is being tested by looking at the expansion in
                the Monomial basis and checking if the coefficients are
                the same if the indexing compositions are permutations
                of each other.

                OUTPUT:

                - ``True`` if ``self`` is symmetric.
                  ``False`` if ``self`` is not symmetric.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: F = QSym.Fundamental()
                    sage: (F[3,2] + F[2,3]).is_symmetric()
                    False
                    sage: (F[1, 1, 1, 2] + F[1, 1, 3] + F[1, 3, 1] + F[2, 1, 1, 1] + F[3, 1, 1]).is_symmetric()
                    True
                    sage: F([]).is_symmetric()
                    True
                """
                M = self.parent().realization_of().Monomial()
                return M(self).is_symmetric()

            def to_symmetric_function(self):
                r"""
                Convert a quasi-symmetric function to a symmetric function.

                OUTPUT:

                - If ``self`` is a symmetric function, then return the expansion
                  in the monomial basis.  Otherwise raise an error.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: F = QSym.Fundamental()
                    sage: (F[3,2] + F[2,3]).to_symmetric_function()
                    Traceback (most recent call last):
                    ...
                    ValueError: F[2, 3] + F[3, 2] is not a symmetric function
                    sage: m = SymmetricFunctions(QQ).m()
                    sage: s = SymmetricFunctions(QQ).s()
                    sage: F(s[3,1,1]).to_symmetric_function()
                    6*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + m[2, 2, 1] + m[3, 1, 1]
                    sage: m(s[3,1,1])
                    6*m[1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1] + m[2, 2, 1] + m[3, 1, 1]
                """
                if self.is_symmetric():
                    M = self.parent().realization_of().Monomial()
                    return M( self ).to_symmetric_function()
                else:
                    raise ValueError("%s is not a symmetric function" % self)

    class Monomial(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric function in the Monomial basis.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: M = QSym.M()
            sage: F = QSym.F()
            sage: M(F[2,2])
            M[1, 1, 1, 1] + M[1, 1, 2] + M[2, 1, 1] + M[2, 2]
            sage: m = SymmetricFunctions(QQ).m()
            sage: M(m[3,1,1])
            M[1, 1, 3] + M[1, 3, 1] + M[3, 1, 1]
            sage: (1+M[1])^3
            M[] + 3*M[1] + 6*M[1, 1] + 6*M[1, 1, 1] + 3*M[1, 2] + 3*M[2] + 3*M[2, 1] + M[3]
            sage: M[1,2,1].coproduct()
            M[] # M[1, 2, 1] + M[1] # M[2, 1] + M[1, 2] # M[1] + M[1, 2, 1] # M[]

        The following is an alias for this basis::

            sage: QSym.Monomial()
            Quasisymmetric functions over the Rational Field in the Monomial basis

        TESTS::

            sage: M(F([]))
            M[]
            sage: M(F(0))
            0
            sage: M(m([]))
            M[]
        """

        def __init__(self, QSym):
            """
            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).Monomial(); M
                Quasisymmetric functions over the Rational Field in the Monomial basis
                sage: TestSuite(M).run()
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='M', bracket=False,
                                             category=QSym.Bases())

        def dual(self):
            r"""
            Return the dual basis to the Monomial basis. This is the complete basis of the
            non-commutative symmetric functions.

            OUTPUT: the complete basis of the non-commutative symmetric functions

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M.dual()
                Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
            """
            return self.realization_of().dual().Complete()

        def product_on_basis(self, I, J):
            """
            The product on Monomial basis elements.

            The product of the basis elements indexed by two compositions
            `I` and `J` is the sum of the basis elements indexed by
            compositions in the stuffle product (also called the
            overlapping shuffle product) of `I` and `J`.

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The product of the Monomial quasi-symmetric functions indexed by ``I`` and
              ``J``, expressed in the Monomial basis.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: c1 = Composition([2])
                sage: c2 = Composition([1,3])
                sage: M.product_on_basis(c1, c2)
                M[1, 2, 3] + M[1, 3, 2] + M[1, 5] + M[2, 1, 3] + M[3, 3]
                sage: M.product_on_basis(c1, Composition([]))
                M[2]
            """
            return self.sum_of_monomials(I.shuffle_product(J, overlap=True))

        def antipode_on_basis(self, compo):
            r"""
            Return the result of the antipode applied to a quasi-symmetric Monomial basis
            element.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The result of the antipode applied to the composition ``compo``, expressed
              in the Monomial basis.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M.antipode_on_basis(Composition([2,1]))
                M[1, 2] + M[3]
                sage: M.antipode_on_basis(Composition([]))
                M[]
            """
            return (-1)**(len(compo)) * self.sum_of_fatter_compositions(compo.reversed())

        def coproduct_on_basis(self, compo):
            r"""
            Return the coproduct of a Monomial basis element.

            Combinatorial rule: deconcatenation.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The coproduct applied to the Monomial quasi-symmetric function indexed by
              ``compo``, expressed in the Monomial basis.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: M[4,2,3].coproduct()
                M[] # M[4, 2, 3] + M[4] # M[2, 3] + M[4, 2] # M[3] + M[4, 2, 3] # M[]
                sage: M.coproduct_on_basis(Composition([]))
                M[] # M[]
            """
            return self.tensor_square().sum_of_monomials((self._indices(compo[:i]),
                                                          self._indices(compo[i:]))
                                                         for i in range(len(compo)+1))

        def lambda_of_monomial(self, I, n):
            r"""
            Return the image of the monomial quasi-symmetric function
            `M_I` under the lambda-map `\lambda^n`, expanded in the
            monomial basis.

            The ring of quasi-symmetric functions over the integers,
            `\mathrm{QSym}_{\ZZ}` (and more generally, the ring of
            quasi-symmetric functions over any binomial ring) becomes
            a `\lambda`-ring (with the `\lambda`-structure inherited
            from the ring of formal power series, so that
            `\lambda^i(x_j)` is `x_j` if `i = 1` and `0` if `i > 1`).

            The Adams operations of this `\lambda`-ring are the
            Adams endomorphisms `\mathbf{f}_n` (see
            :meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.adams_operator`
            for their definition). Using these endomorphisms, the
            `\lambda`-operations can be explicitly computed via the formula

            .. MATH::

                \exp \left(- \sum_{n=1}^{\infty}
                \frac{1}{n} \mathbf{f}_n(x) t^n \right)
                = \sum_{j=0}^{\infty} (-1)^j \lambda^j(x) t^j

            in the ring of formal power series in a variable `t` over
            the ring of quasi-symmetric functions. In particular,
            every composition `I = (I_1, I_2, \cdots, I_\ell )` satisfies

            .. MATH::

                \exp \left(- \sum_{n=1}^{\infty}
                \frac{1}{n} M_{(nI_1, nI_2, \cdots, nI_\ell )} t^n \right)
                = \sum_{j=0}^{\infty} (-1)^j \lambda^j(M_I) t^j

            (corrected version of Remark 2.4 in [Haz2004]_).

            The quasi-symmetric functions `\lambda^i(M_I)` with `n`
            ranging over the positive integers and `I` ranging over
            the reduced Lyndon compositions (i. e., compositions
            which are Lyndon words and have the gcd of their entries
            equal to `1`) form a set of free polynomial generators
            for `\mathrm{QSym}`. See [GriRei18]_, Chapter 6, for
            the proof, and [Haz2004]_ for a major part of it.

            INPUT:

            - ``I`` -- composition
            - ``n`` -- nonnegative integer

            OUTPUT:

            The quasi-symmetric function `\lambda^n(M_I)`, expanded in
            the monomial basis over the ground ring of ``self``.

            EXAMPLES::

                sage: M = QuasiSymmetricFunctions(CyclotomicField()).Monomial()
                sage: M.lambda_of_monomial([1, 2], 2)
                2*M[1, 1, 2, 2] + M[1, 1, 4] + M[1, 2, 1, 2] + M[1, 3, 2] + M[2, 2, 2]
                sage: M.lambda_of_monomial([1, 1], 2)
                3*M[1, 1, 1, 1] + M[1, 1, 2] + M[1, 2, 1] + M[2, 1, 1]
                sage: M = QuasiSymmetricFunctions(Integers(19)).Monomial()
                sage: M.lambda_of_monomial([1, 2], 3)
                6*M[1, 1, 1, 2, 2, 2] + 3*M[1, 1, 1, 2, 4] + 3*M[1, 1, 1, 4, 2]
                + M[1, 1, 1, 6] + 4*M[1, 1, 2, 1, 2, 2] + 2*M[1, 1, 2, 1, 4]
                + 2*M[1, 1, 2, 2, 1, 2] + 2*M[1, 1, 2, 3, 2] + 4*M[1, 1, 3, 2, 2]
                + 2*M[1, 1, 3, 4] + M[1, 1, 4, 1, 2] + M[1, 1, 5, 2]
                + 2*M[1, 2, 1, 1, 2, 2] + M[1, 2, 1, 1, 4] + M[1, 2, 1, 2, 1, 2]
                + M[1, 2, 1, 3, 2] + 4*M[1, 2, 2, 2, 2] + M[1, 2, 2, 4] + M[1, 2, 4, 2]
                + 2*M[1, 3, 1, 2, 2] + M[1, 3, 1, 4] + M[1, 3, 2, 1, 2] + M[1, 3, 3, 2]
                + M[1, 4, 2, 2] + 3*M[2, 1, 2, 2, 2] + M[2, 1, 2, 4] + M[2, 1, 4, 2]
                + 2*M[2, 2, 1, 2, 2] + M[2, 2, 1, 4] + M[2, 2, 2, 1, 2] + M[2, 2, 3, 2]
                + 2*M[2, 3, 2, 2] + M[2, 3, 4] + M[3, 2, 2, 2]

            The map `\lambda^0` sends everything to `1`::

                sage: M = QuasiSymmetricFunctions(ZZ).Monomial()
                sage: all( M.lambda_of_monomial(I, 0) == M.one()
                ....:      for I in Compositions(3) )
                True

            The map `\lambda^1` is the identity map::

                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: all( M.lambda_of_monomial(I, 1) == M(I)
                ....:      for I in Compositions(3) )
                True
                sage: M = QuasiSymmetricFunctions(Integers(5)).Monomial()
                sage: all( M.lambda_of_monomial(I, 1) == M(I)
                ....:      for I in Compositions(3) )
                True
                sage: M = QuasiSymmetricFunctions(ZZ).Monomial()
                sage: all( M.lambda_of_monomial(I, 1) == M(I)
                ....:      for I in Compositions(3) )
                True
            """
            # The following algorithm is a rewriting of the formula in the docstring.
            # We are working over QQ because there are denominators which don't
            # immediately cancel.
            from sage.rings.integer_ring import ZZ
            from sage.rings.rational_field import QQ
            QQM = QuasiSymmetricFunctions(QQ).M()
            QQ_result = QQM.zero()
            for lam in Partitions(n):
                coeff = QQ((-1) ** len(lam)) / lam.centralizer_size()
                QQ_result += coeff * QQM.prod([QQM(self._indices([k * i for i in I]))
                                               for k in lam])
            QQ_result *= (-1) ** n
            # QQ_result is now \lambda^n(M_I) over QQ.
            result = self.sum_of_terms([(J, ZZ(coeff)) for (J, coeff) in QQ_result],
                                        distinct=True)
            return result

        class Element(CombinatorialFreeModule.Element):
            r"""
            Element methods of the ``Monomial`` basis of ``QuasiSymmetricFunctions``.
            """

            def psi_involution(self):
                r"""
                Return the image of the quasisymmetric function ``self``
                under the involution `\psi`.

                The involution `\psi` is defined as the linear map
                `QSym \to QSym` which, for every composition `I`, sends the
                fundamental quasisymmetric function `F_I` to `F_{I^c}`, where
                `I^c` denotes the complement of the composition `I`.
                The map `\psi` is an involution and a graded Hopf algebra
                automorphism of `QSym`. Its restriction to the ring of
                symmetric functions coincides with the omega automorphism of
                the latter ring.

                The involution `\psi` of `QSym` is adjoint to the involution
                `\psi` of `NCSF` by the standard adjunction between `NCSF` and
                `QSym`.

                The involution `\psi` has been denoted by `\psi` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`psi involution on QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`psi involution on NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution on QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: M = QSym.M()
                    sage: M[3,2].psi_involution()
                    -M[3, 2] - M[5]
                    sage: M[3,1].psi_involution()
                    M[3, 1] + M[4]
                    sage: M[3,1,1].psi_involution()
                    M[3, 1, 1] + M[3, 2] + M[4, 1] + M[5]
                    sage: M[1,1,1].psi_involution()
                    M[1, 1, 1] + M[1, 2] + M[2, 1] + M[3]
                    sage: M[[]].psi_involution()
                    M[]
                    sage: M(0).psi_involution()
                    0
                    sage: (2*M[[]] - M[3,1] + 4*M[2]).psi_involution()
                    2*M[] - 4*M[2] - M[3, 1] - M[4]

                This particular implementation is tailored to the monomial
                basis. It is semantically equivalent to the generic
                implementation it overshadows::

                    sage: F = QSym.F()
                    sage: all( F(M[I].psi_involution()) == F(M[I]).psi_involution()
                    ....:      for I in Compositions(3) )
                    True

                    sage: F = QSym.F()
                    sage: all( F(M[I].psi_involution()) == F(M[I]).psi_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                return parent.sum( (-1) ** (I.size() - len(I)) * coeff
                                   * parent.sum_of_fatter_compositions(I)
                                   for I, coeff in
                                   self._monomial_coefficients.items() )

            def expand(self, n, alphabet='x'):
                r"""
                Expand the quasi-symmetric function written in the monomial basis in
                `n` variables.

                INPUT:

                - ``n`` -- integer
                - ``alphabet`` -- string (default: ``'x'``)

                OUTPUT: the quasi-symmetric function ``self`` expressed in the
                `n` variables described by ``alphabet``

                .. TODO:: accept an *alphabet* as input

                EXAMPLES::

                    sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                    sage: M[4,2].expand(3)
                    x0^4*x1^2 + x0^4*x2^2 + x1^4*x2^2

                One can use a different set of variables by using the
                optional argument ``alphabet``::

                    sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                    sage: M[2,1,1].expand(4,alphabet='y')
                    y0^2*y1*y2 + y0^2*y1*y3 + y0^2*y2*y3 + y1^2*y2*y3

                TESTS::

                    sage: (3*M([])).expand(2)
                    3
                    sage: M[4,2].expand(0)
                    0
                    sage: M([]).expand(0)
                    1
                """
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                M = self.parent()
                P = PolynomialRing(M.base_ring(), n, alphabet)
                x = P.gens()

                def on_basis(comp, i):
                    if not comp:
                        return P.one()
                    elif len(comp) > i:
                        return P.zero()
                    else:
                        return x[i-1]**comp[-1] * on_basis(comp[:-1], i-1) + \
                                                  on_basis(comp, i-1)
                return M._apply_module_morphism(self, lambda comp: on_basis(comp,n),
                                                codomain=P)

            def is_symmetric( self ):
                r"""
                Determine if a quasi-symmetric function, written in the Monomial basis,
                is symmetric.

                This is being tested by looking at the expansion in the Monomial
                basis and checking if the coefficients are the same if the indexing
                compositions are permutations of each other.

                OUTPUT:

                - ``True`` if ``self`` is an element of the symmetric
                  functions and ``False`` otherwise.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.Monomial()
                    sage: (M[3,2] + M[2,3] + M[4,1]).is_symmetric()
                    False
                    sage: (M[3,2] + M[2,3]).is_symmetric()
                    True
                    sage: (M[1,2,1] + M[1,1,2]).is_symmetric()
                    False
                    sage: (M[1,2,1] + M[1,1,2] + M[2,1,1]).is_symmetric()
                    True
                """
                # We must check that every rearrangement of a composition
                # that appears in self appears with the same coefficient.
                # We use a dictionary to keep track of the coefficient
                # and how many rearrangements of the composition we've seen.
                from sage.combinat.permutation import Permutations_mset
                d = {}
                for (I, coeff) in self:
                    partition = I.to_partition()
                    if partition not in d:
                        d[partition] = [coeff, 1]
                    else:
                        if d[partition][0] != coeff:
                            return False
                        else:
                            d[partition][1] += 1
                # make sure we've seen each rearrangement of the composition
                return all(d[partition][1] == Permutations_mset(partition).cardinality()
                            for partition in d)

            def to_symmetric_function(self):
                r"""
                Take a quasi-symmetric function, expressed in the monomial basis, and
                return its symmetric realization, when possible, expressed in the
                monomial basis of symmetric functions.

                OUTPUT:

                - If ``self`` is a symmetric function, then the expansion
                  in the monomial basis of the symmetric functions is returned.
                  Otherwise an error is raised.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.Monomial()
                    sage: (M[3,2] + M[2,3] + M[4,1]).to_symmetric_function()
                    Traceback (most recent call last):
                    ...
                    ValueError: M[2, 3] + M[3, 2] + M[4, 1] is not a symmetric function
                    sage: (M[3,2] + M[2,3] + 2*M[4,1] + 2*M[1,4]).to_symmetric_function()
                    m[3, 2] + 2*m[4, 1]
                    sage: m = SymmetricFunctions(QQ).m()
                    sage: M(m[3,1,1]).to_symmetric_function()
                    m[3, 1, 1]
                    sage: (M(m[2,1])*M(m[2,1])).to_symmetric_function()-m[2,1]*m[2,1]
                    0

                TESTS::

                    sage: (M(0)).to_symmetric_function()
                    0
                    sage: (M([])).to_symmetric_function()
                    m[]
                    sage: (2*M([])).to_symmetric_function()
                    2*m[]

                We check that the result is indexed by partitions::

                    sage: M([]).to_symmetric_function().leading_support().parent()
                    Partitions
                """
                m = SymmetricFunctions(self.parent().base_ring()).monomial()
                if self.is_symmetric():
                    return m._from_dict({_Partitions(list(I)): coeff
                                         for I, coeff in self
                        if list(I) in _Partitions}, remove_zeros=False)
                else:
                    raise ValueError("%s is not a symmetric function" % self)

    M = Monomial

    class Fundamental(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric functions in the Fundamental basis.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: F = QSym.F()
            sage: M = QSym.M()
            sage: F(M[2,2])
            F[1, 1, 1, 1] - F[1, 1, 2] - F[2, 1, 1] + F[2, 2]
            sage: s = SymmetricFunctions(QQ).s()
            sage: F(s[3,2])
            F[1, 2, 2] + F[1, 3, 1] + F[2, 2, 1] + F[2, 3] + F[3, 2]
            sage: (1+F[1])^3
            F[] + 3*F[1] + 3*F[1, 1] + F[1, 1, 1] + 2*F[1, 2] + 3*F[2] + 2*F[2, 1] + F[3]
            sage: F[1,2,1].coproduct()
            F[] # F[1, 2, 1] + F[1] # F[2, 1] + F[1, 1] # F[1, 1] + F[1, 2] # F[1] + F[1, 2, 1] # F[]

        The following is an alias for this basis::

            sage: QSym.Fundamental()
            Quasisymmetric functions over the Rational Field in the Fundamental basis

        TESTS::

            sage: F(M([]))
            F[]
            sage: F(M(0))
            0
            sage: F(s([]))
            F[]
            sage: F(s(0))
            0
        """

        def __init__(self, QSym):
            """
            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).Fundamental(); F
                Quasisymmetric functions over the Rational Field in the Fundamental basis
                sage: TestSuite(F).run()
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='F', bracket=False,
                                             category=QSym.Bases())

        def _from_schur_on_basis(self, la):
            r"""
            Map the Schur symmetric function indexed by ``la`` to the
            Fundamental basis.

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: F(s[2,2]-s[3,1])                                          # indirect doctest
                F[1, 2, 1] - F[1, 3] - F[3, 1]

                sage: F._from_schur_on_basis(Partition([]))
                F[]

                sage: F._from_schur_on_basis(Partition([2,2,2]))
                F[1, 1, 2, 1, 1] + F[1, 2, 1, 2] + F[1, 2, 2, 1] + F[2, 1, 2, 1] + F[2, 2, 2]
            """
            C = self._indices
            res = 0
            n = la.size()
            for T in StandardTableaux(la):
                des = T.standard_descents()
                comp = C.from_descents([d-1 for d in des], n)
                res += self.monomial(comp)
            return res

        def dual(self):
            r"""
            Return the dual basis to the Fundamental basis. This is the ribbon
            basis of the non-commutative symmetric functions.

            OUTPUT: the ribbon basis of the non-commutative symmetric functions

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.dual()
                Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
            """
            return self.realization_of().dual().Ribbon()

        def antipode_on_basis(self, compo):
            r"""
            Return the antipode to a Fundamental quasi-symmetric basis element.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The result of the antipode applied to the quasi-symmetric
              Fundamental basis element indexed by ``compo``.

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.antipode_on_basis(Composition([2,1]))
                -F[2, 1]
            """
            return (-1)**(compo.size()) * self.monomial(compo.conjugate())

        def coproduct_on_basis(self, compo):
            r"""
            Return the coproduct to a Fundamental quasi-symmetric basis element.

            Combinatorial rule: quasi-deconcatenation.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The application of the coproduct to the Fundamental quasi-symmetric
              function indexed by the composition ``compo``.

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: F[4].coproduct()
                F[] # F[4] + F[1] # F[3] + F[2] # F[2] + F[3] # F[1] + F[4] # F[]
                sage: F[2,1,3].coproduct()
                F[] # F[2, 1, 3] + F[1] # F[1, 1, 3] + F[2] # F[1, 3] + F[2, 1] # F[3] + F[2, 1, 1] # F[2] + F[2, 1, 2] # F[1] + F[2, 1, 3] # F[]

            TESTS::

                sage: F.coproduct_on_basis(Composition([2,1,3]))
                F[] # F[2, 1, 3] + F[1] # F[1, 1, 3] + F[2] # F[1, 3] + F[2, 1] # F[3] + F[2, 1, 1] # F[2] + F[2, 1, 2] # F[1] + F[2, 1, 3] # F[]
                sage: F.one().coproduct()          # generic for graded / graded connected
                F[] # F[]
            """
            T = self.tensor_square()
            C = Composition
            return T.sum_of_monomials( ( C(compo[:i]), C(compo[i:]) ) for i in range(len(compo)+1) ) + \
                   T.sum_of_monomials( ( C(compo[:i]+[j]), C([compo[i]-j]+compo[i+1:]) )
                                       for i in range(len(compo))
                                       for j in range(1, compo[i]) )

        @cached_method
        def Eulerian(self, n, j, k=None):
            r"""
            Return the Eulerian (quasi)symmetric function `Q_{n,j}` (with
            `n` either an integer or a partition) defined in [SW2010]_ in
            terms of the fundamental quasisymmetric functions.
            Or, if the optional argument ``k`` is specified, return the
            function `Q_{n,j,k}` defined ibidem.

            If `n` and `j` are nonnegative integers, then the Eulerian
            quasisymmetric function `Q_{n,j}` is defined as

            .. MATH::

                Q_{n,j} := \sum_{\sigma} F_{\mathrm{Dex}(\sigma)},

            where we sum over all permutations `\sigma \in S_n` such that
            the number of excedances of `\sigma` is `j`, and where
            `\mathrm{Dex}(\sigma)` is a composition of `n` defined as follows:
            Let `S` be the set of all `i \in \{ 1, 2, \ldots, n-1 \}` such
            that either `\sigma_i > \sigma_{i+1} > i+1` or
            `i \geq \sigma_i > \sigma_{i+1}` or
            `\sigma_{i+1} > i + 1 > \sigma_i`. Then,
            `\mathrm{Dex}(\sigma)` is set to be the composition of `n` whose
            descent set is `S`.

            Here, an excedance of a permutation `\sigma \in S_n` means an
            element `i \in \{ 1, 2, \ldots, n-1 \}` satisfying `\sigma_i > i`.

            Similarly we can define a quasisymmetric function `Q_{\lambda, j}`
            for every partition `\lambda` and every nonnegative integer `j`.
            This differs from `Q_{n,j}` only in that the sum is restricted to
            all permutations `\sigma \in S_n` whose cycle type is `\lambda`
            (where `n = |\lambda|`, and where we still require the number of
            excedances to be `j`). The method at hand allows computing these
            functions by passing `\lambda` as the ``n`` parameter.

            Analogously we can define a quasisymmetric function `Q_{n,j,k}` for
            any nonnegative integers `n`, `j` and `k` by restricting the sum to
            all permutations `\sigma \in S_n` that have exactly `k` fixed
            points (and `j` excedances). This can be obtained by specifying the
            optional ``k`` argument in this method.

            All three versions of Eulerian quasisymmetric functions
            (`Q_{n,j}`, `Q_{\lambda,j}` and `Q_{n,j,k}`) are actually
            symmetric functions. See
            :meth:`~sage.combinat.sf.SymmetricFunctionsBases.ParentMethods.Eulerian`.

            INPUT:

            - ``n`` -- the nonnegative integer `n` or a partition
            - ``j`` -- the number of excedances
            - ``k`` -- (optional) if specified, determines the number of fixed
              points of the permutations which are being summed over

            EXAMPLES::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.Eulerian(3, 1)
                F[1, 2] + F[2, 1] + 2*F[3]
                sage: F.Eulerian(4, 2)
                F[1, 2, 1] + 2*F[1, 3] + 3*F[2, 2] + 2*F[3, 1] + 3*F[4]
                sage: F.Eulerian(5, 2)
                F[1, 1, 2, 1] + F[1, 1, 3] + F[1, 2, 1, 1] + 7*F[1, 2, 2] + 6*F[1, 3, 1] + 6*F[1, 4] + 2*F[2, 1, 2] + 7*F[2, 2, 1] + 11*F[2, 3] + F[3, 1, 1] + 11*F[3, 2] + 6*F[4, 1] + 6*F[5]
                sage: F.Eulerian(4, 0)
                F[4]
                sage: F.Eulerian(4, 3)
                F[4]
                sage: F.Eulerian(4, 1, 2)
                F[1, 2, 1] + F[1, 3] + 2*F[2, 2] + F[3, 1] + F[4]
                sage: F.Eulerian(Partition([2, 2, 1]), 2)
                F[1, 1, 2, 1] + F[1, 2, 1, 1] + 2*F[1, 2, 2] + F[1, 3, 1]
                 + F[1, 4] + F[2, 1, 2] + 2*F[2, 2, 1] + 2*F[2, 3]
                 + 2*F[3, 2] + F[4, 1] + F[5]
                sage: F.Eulerian(0, 0)
                F[]
                sage: F.Eulerian(0, 1)
                0
                sage: F.Eulerian(1, 0)
                F[1]
                sage: F.Eulerian(1, 1)
                0

            TESTS::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.Eulerian(Partition([3, 1]), 1, 1)
                Traceback (most recent call last):
                ...
                ValueError: invalid input, k cannot be specified
            """
            from sage.combinat.partition import _Partitions
            if n in _Partitions:
                if k is not None:
                    raise ValueError("invalid input, k cannot be specified")
                la = _Partitions(n)
                n = sum(la)
            else:
                la = None

            if n == 0:
                if k is not None and k != 0:
                    return self.zero()
                if j != 0:
                    return self.zero()
                return self.one()

            monomials = []
            for p in Permutations(n):
                dex = []
                exc = 0
                for i in range(n-1):
                    if p[i] > i + 1:
                        exc += 1
                    if (p[i] > p[i+1] or (p[i] <= i+1 and p[i+1] > i+2)) \
                            and not (p[i] > i+1 and p[i+1] <= i+2):
                        dex.append(i)

                if exc != j:
                    continue
                if k is not None and p.number_of_fixed_points() != k:
                    continue
                if la is not None and p.cycle_type() != la:
                    continue

                # Converting to a composition
                d = [-1] + dex + [n-1]
                monomials.append(Compositions()( [d[i+1]-d[i] for i in range(len(d)-1)] ))

            return self.sum_of_monomials(monomials)

        class Element(CombinatorialFreeModule.Element):

            def internal_coproduct(self):
                r"""
                Return the inner coproduct of ``self`` in the Fundamental basis.

                The inner coproduct (also known as the Kronecker coproduct,
                or as the second comultiplication on the `R`-algebra of
                quasi-symmetric functions) is an `R`-algebra homomorphism
                `\Delta^{\times}` from the `R`-algebra of quasi-symmetric
                functions to the tensor square (over `R`) of quasi-symmetric
                functions. It can be defined in the following two ways:

                #. If `I` is a composition, then a `(0, I)`-matrix will mean a
                   matrix whose entries are nonnegative integers such that no
                   row and no column of this matrix is zero, and such that if
                   all the nonzero entries of the matrix are read (row by row,
                   starting at the topmost row, reading every row from left to
                   right), then the reading word obtained is `I`. If `A` is
                   a `(0, I)`-matrix, then `\mathrm{row}(A)` will denote the
                   vector of row sums of `A` (regarded as a composition), and
                   `\mathrm{column}(A)` will denote the vector of column sums
                   of `A` (regarded as a composition).

                   For every composition `I`, the internal coproduct
                   `\Delta^{\times}(M_I)` of the `I`-th monomial quasisymmetric
                   function `M_I` is the sum

                   .. MATH::

                       \sum_{A \hbox{ is a } (0, I) \text{-matrix}}
                       M_{\mathrm{row}(A)} \otimes M_{\mathrm{column}(A)}.

                   See Section 11.39 of [HazWitt1]_.

                #. For every permutation `w`, let `C(w)` denote the descent
                   composition of `w`. Then, for any composition `I` of size
                   `n`, the internal coproduct `\Delta^{\times}(F_I)` of the
                   `I`-th fundamental quasisymmetric function `F_I` is the sum

                   .. MATH::

                       \sum_{\substack{\sigma \in S_n,\\ \tau \in S_n,\\
                       \tau \sigma = \pi}} F_{C(\sigma)} \otimes F_{C(\tau)},

                   where `\pi` is any permutation in `S_n` having descent
                   composition `I` and where permutations act from the left and
                   multiply accordingly, so `\tau \sigma` means first applying
                   `\sigma` and then `\tau`. See Theorem 4.23 in [Mal1993]_,
                   but beware of the notations which are apparently different
                   from those in [HazWitt1]_.

                The restriction of the internal coproduct to the
                `R`-algebra of symmetric functions is the well-known
                internal coproduct on the symmetric functions.

                The method :meth:`kronecker_coproduct` is a synonym of this one.

                EXAMPLES:

                Let us compute the internal coproduct of `M_{21}` (which is
                short for `M_{[2, 1]}`). The `(0, [2,1])`-matrices are

                .. MATH::

                    \begin{bmatrix}2& 1\end{bmatrix},
                    \begin{bmatrix}2\\1\end{bmatrix},
                    \begin{bmatrix}2& 0\\0& 1\end{bmatrix}, \hbox{ and }
                    \begin{bmatrix}0&2\\1&0\end{bmatrix}

                so

                .. MATH::

                    \Delta^\times(M_{21}) = M_{3} \otimes M_{21} +
                    M_{21} \otimes M_3 + M_{21} \otimes M_{21} +
                    M_{21} \otimes M_{12}.

                This is confirmed by the following Sage computation (incidentally
                demonstrating the non-cocommutativity of the internal
                coproduct)::

                    sage: M = QuasiSymmetricFunctions(ZZ).M()
                    sage: a = M([2,1])
                    sage: a.internal_coproduct()
                    M[2, 1] # M[1, 2] + M[2, 1] # M[2, 1] + M[2, 1] # M[3] + M[3] # M[2, 1]

                Some examples on the Fundamental basis::

                    sage: F = QuasiSymmetricFunctions(ZZ).F()
                    sage: F([1,1]).internal_coproduct()
                    F[1, 1] # F[2] + F[2] # F[1, 1]
                    sage: F([2]).internal_coproduct()
                    F[1, 1] # F[1, 1] + F[2] # F[2]
                    sage: F([3]).internal_coproduct()
                    F[1, 1, 1] # F[1, 1, 1] + F[1, 2] # F[1, 2] + F[1, 2] # F[2, 1]
                     + F[2, 1] # F[1, 2] + F[2, 1] # F[2, 1] + F[3] # F[3]
                    sage: F([1,2]).internal_coproduct()
                    F[1, 1, 1] # F[1, 2] + F[1, 2] # F[2, 1] + F[1, 2] # F[3]
                     + F[2, 1] # F[1, 1, 1] + F[2, 1] # F[2, 1] + F[3] # F[1, 2]
                """
                F = self.parent()
                F2 = F.tensor(F)
                result = F2.zero()
                from sage.categories.tensor import tensor
                from sage.combinat.permutation import Permutation
                for I, a in self:
                    # We must add a * \Delta^\times(F_I) to result.
                    from sage.combinat.permutation import descents_composition_last
                    pi = descents_composition_last(I)
                    n = I.size()
                    for sigma in Permutations(n):
                        sigma_inverse = sigma.inverse()
                        # If the __mul__ of permutations wasn't such a mess,
                        # the next line could be as simple as
                        # tau = pi * sigma_inverse.
                        tau = Permutation([pi(i) for i in sigma_inverse])
                        result += a * tensor([F(sigma.descents_composition()),
                                              F(tau.descents_composition())])
                return result

            kronecker_coproduct = internal_coproduct

            def star_involution(self):
                r"""
                Return the image of the quasisymmetric function ``self`` under
                the star involution.

                The star involution is defined as the linear map
                `QSym \to QSym` which, for every composition `I`, sends the
                monomial quasisymmetric function `M_I` to `M_{I^r}`. Here, if
                `I` is a composition, we denote by `I^r` the reversed
                composition of `I`. Denoting by `f^{\ast}` the image of an
                element `f \in QSym` under the star involution, it can be shown
                that every composition `I` satisfies

                .. MATH::

                    (M_I)^{\ast} = M_{I^r}, \quad (F_I)^{\ast} = F_{I^r},

                where `F_I` denotes the fundamental quasisymmetric function
                corresponding to the composition `I`. The star involution is an
                involution, an algebra automorphism and a coalgebra
                anti-automorphism of `QSym`. It also is an automorphism of the
                graded vector space `QSym`, and is the identity on the subspace
                `Sym` of `QSym`. It is adjoint to the star involution on `NCSF`
                by the standard adjunction between `NCSF` and `QSym`.

                The star involution has been denoted by `\rho` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`star involution on QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`,
                    :meth:`star involution on NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: F = QSym.F()
                    sage: F[3,1].star_involution()
                    F[1, 3]
                    sage: F[5,3].star_involution()
                    F[3, 5]
                    sage: (F[9,1] - F[6,2] + 2*F[6,4] - 3*F[3] + 4*F[[]]).star_involution()
                    4*F[] + F[1, 9] - F[2, 6] - 3*F[3] + 2*F[4, 6]
                    sage: (F[3,3] - 2*F[2]).star_involution()
                    -2*F[2] + F[3, 3]
                    sage: F([4,2]).star_involution()
                    F[2, 4]
                    sage: dI = QSym.dI()
                    sage: dI([1,2]).star_involution()
                    -dI[1, 2] + dI[2, 1]
                    sage: dI.zero().star_involution()
                    0
                """
                parent = self.parent()
                dct = {I.reversed(): coeff for (I, coeff) in self}
                return parent._from_dict(dct)

    F = Fundamental

    class Essential(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric functions in the Essential basis.

        The Essential quasi-symmetric functions are defined by

        .. MATH::

            E_I = \sum_{J \geq I} M_J = \sum_{i_1 \leq \cdots \leq i_k}
            x_{i_1}^{I_1} \cdots x_{i_k}^{I_k},

        where `I = (I_1, \ldots, I_k)`.

        .. NOTE::

            Our convention of `\leq` and `\geq` of compositions is
            opposite that of [Hoff2015]_.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: E = QSym.E()
            sage: M = QSym.M()
            sage: E(M[2,2])
            E[2, 2] - E[4]
            sage: s = SymmetricFunctions(QQ).s()
            sage: E(s[3,2])
            5*E[1, 1, 1, 1, 1] - 2*E[1, 1, 1, 2] - 2*E[1, 1, 2, 1]
             - 2*E[1, 2, 1, 1] + E[1, 2, 2] - 2*E[2, 1, 1, 1]
             + E[2, 1, 2] + E[2, 2, 1]
            sage: (1 + E[1])^3
            E[] + 3*E[1] + 6*E[1, 1] + 6*E[1, 1, 1] - 3*E[1, 2]
             - 3*E[2] - 3*E[2, 1] + E[3]
            sage: E[1,2,1].coproduct()
            E[] # E[1, 2, 1] + E[1] # E[2, 1] + E[1, 2] # E[1] + E[1, 2, 1] # E[]

        The following is an alias for this basis::

            sage: QSym.Essential()
            Quasisymmetric functions over the Rational Field in the Essential basis

        TESTS::

            sage: E(M([]))
            E[]
            sage: E(M(0))
            0
            sage: E(s([]))
            E[]
            sage: E(s(0))
            0
        """

        def __init__(self, QSym):
            """
            EXAMPLES::

                sage: E = QuasiSymmetricFunctions(QQ).Essential()
                sage: TestSuite(E).run()

            TESTS::

                sage: E = QuasiSymmetricFunctions(QQ).E()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: all(E(M(E[c])) == E[c] for n in range(5)
                ....:     for c in Compositions(n))
                True
                sage: all(M(E(M[c])) == M[c] for n in range(5)
                ....:     for c in Compositions(n))
                True
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='E', bracket=False,
                                             category=QSym.Bases())

            M = QSym.M()
            category = self.realization_of()._category
            # This changes Monomial into Essential
            M.module_morphism(self.alternating_sum_of_fatter_compositions,
                              codomain=self, category=category
                              ).register_as_coercion()
            # This changes Essential into Monomial
            self.module_morphism(M.sum_of_fatter_compositions,
                                 codomain=M, category=category
                                 ).register_as_coercion()

        def antipode_on_basis(self, compo):
            r"""
            Return the result of the antipode applied to a quasi-symmetric
            Essential basis element.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The result of the antipode applied to the composition ``compo``,
              expressed in the Essential basis.

            EXAMPLES::

                sage: E = QuasiSymmetricFunctions(QQ).E()
                sage: E.antipode_on_basis(Composition([2,1]))
                E[1, 2] - E[3]
                sage: E.antipode_on_basis(Composition([]))
                E[]

            TESTS::

                sage: E = QuasiSymmetricFunctions(QQ).E()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: all(E(M(E[c]).antipode()) == E[c].antipode()
                ....:     for n in range(5) for c in Compositions(n))
                True

                sage: all((-1)**len(I) * E[I] == M[I].star_involution().antipode()
                ....:     for k in [3,4] for I in Compositions(k))
                True
            """
            return (-1)**len(compo) * self.alternating_sum_of_fatter_compositions(compo.reversed())

        def coproduct_on_basis(self, compo):
            r"""
            Return the coproduct of a Essential basis element.

            Combinatorial rule: deconcatenation.

            INPUT:

            - ``compo`` -- composition

            OUTPUT:

            - The coproduct applied to the Essential quasi-symmetric function
              indexed by ``compo``, expressed in the Essential basis.

            EXAMPLES::

                sage: E = QuasiSymmetricFunctions(QQ).Essential()
                sage: E[4,2,3].coproduct()
                E[] # E[4, 2, 3] + E[4] # E[2, 3] + E[4, 2] # E[3] + E[4, 2, 3] # E[]
                sage: E.coproduct_on_basis(Composition([]))
                E[] # E[]
            """
            return self.tensor_square().sum_of_monomials((self._indices(compo[:i]),
                                                          self._indices(compo[i:]))
                                                         for i in range(len(compo)+1))

        def product_on_basis(self, I, J):
            r"""
            The product on Essential basis elements.

            The product of the basis elements indexed by two compositions
            `I` and `J` is the sum of the basis elements indexed by
            compositions `K` in the stuffle product (also called the
            overlapping shuffle product) of `I` and `J` with a
            coefficient of `(-1)^{\ell(I) + \ell(J) - \ell(K)}`,
            where `\ell(C)` is the length of the composition `C`.

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The product of the Essential quasi-symmetric functions indexed
              by ``I`` and ``J``, expressed in the Essential basis.

            EXAMPLES::

                sage: E = QuasiSymmetricFunctions(QQ).E()
                sage: c1 = Composition([2])
                sage: c2 = Composition([1,3])
                sage: E.product_on_basis(c1, c2)
                E[1, 2, 3] + E[1, 3, 2] - E[1, 5] + E[2, 1, 3] - E[3, 3]
                sage: E.product_on_basis(c1, Composition([]))
                E[2]
                sage: E.product_on_basis(c1, Composition([3]))
                E[2, 3] + E[3, 2] - E[5]

            TESTS::

                sage: E = QuasiSymmetricFunctions(QQ).E()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: all(E(M(E[cp])*M(E[c])) == E[cp]*E[c]  # long time
                ....:     for c in Compositions(3) for cp in Compositions(5))
                True
            """
            n = len(I) + len(J)
            return self.sum_of_terms((K, (-1)**(n - len(K)))
                                     for K in I.shuffle_product(J, overlap=True))

    E = Essential

    class Quasisymmetric_Schur(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric function in the Quasisymmetric
        Schur basis.

        The basis of Quasisymmetric Schur functions is defined in [QSCHUR]_
        and in Definition 5.1.1 of [LMvW13]_.
        Don't mistake them for the completely unrelated quasi-Schur
        functions of [NCSF1]_!

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: QS = QSym.QS()
            sage: F = QSym.F()
            sage: M = QSym.M()
            sage: F(QS[1,2])
            F[1, 2]
            sage: M(QS[1,2])
            M[1, 1, 1] + M[1, 2]
            sage: s = SymmetricFunctions(QQ).s()
            sage: QS(s[2,1,1])
            QS[1, 1, 2] + QS[1, 2, 1] + QS[2, 1, 1]
        """

        def __init__(self, QSym):
            r"""
            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: F = QSym.Fundamental()
                sage: QS = QSym.Quasisymmetric_Schur()
                sage: QS(F(QS.an_element())) == QS.an_element()
                True
                sage: F(QS(F.an_element())) == F.an_element()
                True
                sage: M = QSym.Monomial()
                sage: QS(M(QS.an_element())) == QS.an_element()
                True
                sage: M(QS(M.an_element())) == M.an_element()
                True
                sage: TestSuite(QS).run() # long time
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='QS', bracket=False,
                                             category=QSym.Bases())

        def _realization_name(self):
            r"""
            Return a nicer name for ``self`` than what is inherited
            from :mod:`sage.categories.sets_cat`.

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._realization_name()
                'Quasisymmetric Schur'
            """
            return "Quasisymmetric Schur"

        @cached_method
        def _from_monomial_transition_matrix(self, n):
            r"""
            A matrix representing the transition coefficients to
            the complete basis and the ordering of compositions.

            INPUT:

            - ``n`` -- integer

            OUTPUT: a pair of a square matrix and the ordered list of compositions

            EXAMPLES::

                sage: QS = QuasiSymmetricFunctions(QQ).QS()
                sage: QS._from_monomial_transition_matrix(4)[0]
                [ 1 -1 -1  0  1  1  1 -1]
                [ 0  1  0  0 -1 -1  0  1]
                [ 0  0  1 -1  0  0 -1  1]
                [ 0  0  0  1 -1 -1 -1  1]
                [ 0  0  0  0  1  0  0 -1]
                [ 0  0  0  0  0  1  0 -1]
                [ 0  0  0  0  0  0  1 -1]
                [ 0  0  0  0  0  0  0  1]
            """
            if n == 0:
                return (matrix([[]]), [])
            CO = compositions_order(n)
            # ZZ is faster than over QQ for inverting a matrix
            from sage.rings.integer_ring import ZZ
            MS = MatrixSpace(ZZ, len(CO))
            M = MS([[number_of_SSRCT(al, be) for al in CO] for be in CO])
            return (M.inverse_of_unit(), CO)

        @cached_method
        def _from_monomial_on_basis(self, comp):
            r"""
            Map the Monomial quasi-symmetric function indexed by
            ``comp`` to the Quasisymmetric Schur basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT: a quasi-symmetric function in the Quasisymmetric Schur basis

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: M = QSym.M()
                sage: QS._from_monomial_on_basis(Composition([1,3,1]))
                QS[1, 1, 1, 1, 1] - QS[1, 1, 2, 1] + QS[1, 3, 1] - QS[2, 2, 1]
                sage: QS._from_monomial_on_basis(Composition([2]))
                -QS[1, 1] + QS[2]
            """
            comp = Composition(comp)
            if not comp._list:
                return self.one()
            T, comps = self._from_monomial_transition_matrix(comp.size())
            i = comps.index(comp)
            return self._from_dict({c: T[i,j] for j,c in enumerate(comps)
                                    if T[i,j] != 0},
                                   remove_zeros=False)

        @cached_method
        def _to_monomial_on_basis(self, comp_shape):
            r"""
            Expand the quasi-symmetric Schur function in the Monomial basis.

            The expansion of the quasi-symmetric Schur function indexed
            by ``comp_shape`` has coefficients which are given by the method
            :meth:`~sage.combinat.ncsf_qsym.combinatorics.number_of_SSRCT`.

            INPUT:

            - ``comp_shape`` -- a composition

            OUTPUT: a quasi-symmetric function in the Monomial basis

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._to_monomial_on_basis(Composition([1,3,1]))
                2*M[1, 1, 1, 1, 1] + 2*M[1, 1, 2, 1] + M[1, 2, 1, 1] + M[1, 3, 1] + M[2, 1, 1, 1] + M[2, 2, 1]
            """
            M = self.realization_of().Monomial()
            if not comp_shape:
                return M([])
            return M.sum_of_terms(((comp_content,
                    number_of_SSRCT(comp_content, comp_shape))
                    for comp_content in Compositions(sum(comp_shape))),
                    distinct=True )

        def dual(self):
            r"""
            The dual basis to the Quasisymmetric Schur basis.

            The dual basis to the Quasisymmetric Schur basis is
            implemented as dual.

            OUTPUT:

            - the dual Quasisymmetric Schur basis of the
              non-commutative symmetric functions

            EXAMPLES::

                sage: QS = QuasiSymmetricFunctions(QQ).Quasisymmetric_Schur()
                sage: QS.dual()
                Non-Commutative Symmetric Functions over the Rational Field
                 in the dual Quasisymmetric-Schur basis
            """
            return self.realization_of().dual().dualQuasisymmetric_Schur()

    QS = Quasisymmetric_Schur

    class Young_Quasisymmetric_Schur(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric functions in the Young
        Quasisymmetric Schur basis.

        The basis of Young Quasisymmetric Schur functions is from
        Definition 5.2.1 of [LMvW13]_.

        This basis is related to the Quasisymmetric Schur basis ``QS`` by
        ``QS(alpha.reversed()) == YQS(alpha).star_involution()`` .

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: YQS = QSym.YQS()
            sage: F = QSym.F()
            sage: QS = QSym.QS()
            sage: F(YQS[1,2])
            F[1, 2]
            sage: all(QS(al.reversed())==YQS(al).star_involution() for al in Compositions(5))
            True
            sage: s = SymmetricFunctions(QQ).s()
            sage: YQS(s[2,1,1])
            YQS[1, 1, 2] + YQS[1, 2, 1] + YQS[2, 1, 1]
        """

        def __init__(self, QSym):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: F = QSym.Fundamental()
                sage: YQS = QSym.Young_Quasisymmetric_Schur()
                sage: YQS(F(YQS.an_element())) == YQS.an_element()
                True
                sage: F(YQS(F.an_element())) == F.an_element()
                True
                sage: M = QSym.Monomial()
                sage: YQS(M(YQS.an_element())) == YQS.an_element()
                True
                sage: M(YQS(M.an_element())) == M.an_element()
                True
                sage: TestSuite(YQS).run() # long time
            """
            self._QS = QSym.QS()
            self._M = QSym.M()
            CombinatorialFreeModule.__init__(self, QSym.base_ring(),
                                             Compositions(), prefix='YQS',
                                             bracket=False, category=QSym.Bases())
            self.module_morphism(self._to_monomial_on_basis,
                                 codomain=self._M, category=QSym.Bases()).register_as_coercion()
            self._M.module_morphism(self._from_monomial_on_basis,
                                    codomain=self, category=QSym.Bases()).register_as_coercion()

        def _realization_name(self):
            r"""
            Return a nicer name for ``self`` than what is inherited
            from :mod:`sage.categories.sets_cat`.

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: QS = QSym.QS()
                sage: QS._realization_name()
                'Quasisymmetric Schur'
            """
            return "Young Quasisymmetric Schur"

        @cached_method
        def _to_monomial_on_basis(self, comp):
            r"""
            Expand the Young quasi-symmetric Schur function in the
            Monomial basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT: a quasi-symmetric function in the Monomial basis

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: YQS = QSym.QS()
                sage: YQS._to_monomial_on_basis(Composition([2]))
                M[1, 1] + M[2]
                sage: YQS._to_monomial_on_basis(Composition([1,3,1]))
                2*M[1, 1, 1, 1, 1] + 2*M[1, 1, 2, 1] + M[1, 2, 1, 1] + M[1, 3, 1] + M[2, 1, 1, 1] + M[2, 2, 1]
            """
            return self._M(self._QS.monomial(comp.reversed())).star_involution()

        @cached_method
        def _from_monomial_on_basis(self, comp):
            r"""
            The expansion of the quasi-symmetric Schur function indexed
            by ``comp_shape`` has coefficients which are given by the method
            :meth:`~sage.combinat.ncsf_qsym.combinatorics.number_of_SSRCT`.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT: a quasi-symmetric function in the Young Quasisymmetric Schur basis

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: YQS = QSym.YQS()
                sage: YQS._from_monomial_on_basis(Composition([2]))
                -YQS[1, 1] + YQS[2]
                sage: YQS._from_monomial_on_basis(Composition([1,3,1]))
                YQS[1, 1, 1, 1, 1] - YQS[1, 2, 1, 1] - YQS[1, 2, 2] + YQS[1, 3, 1]
            """
            elt = self._QS(self._M.monomial(comp.reversed()))
            return self._from_dict({al.reversed(): c for al,c in elt},
                                   coerce=False, remove_zeros=False)

    YQS = Young_Quasisymmetric_Schur

    class dualImmaculate(CombinatorialFreeModule, BindableClass):
        def __init__(self, QSym):
            r"""
            The dual immaculate basis of the quasi-symmetric functions.

            This basis first appears in [BBSSZ2012]_.

            EXAMPLES::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: dI = QSym.dI()
                sage: dI([1,3,2])*dI([1])  # long time (6s on sage.math, 2013)
                dI[1, 1, 3, 2] + dI[2, 3, 2]
                sage: dI([1,3])*dI([1,1])
                dI[1, 1, 1, 3] + dI[1, 1, 4] + dI[1, 2, 3] - dI[1, 3, 2] - dI[1, 4, 1] - dI[1, 5] + dI[2, 3, 1] + dI[2, 4]
                sage: dI([3,1])*dI([2,1])  # long time (7s on sage.math, 2013)
                dI[1, 1, 5] - dI[1, 4, 1, 1] - dI[1, 4, 2] - 2*dI[1, 5, 1] - dI[1, 6] - dI[2, 4, 1] - dI[2, 5] - dI[3, 1, 3] + dI[3, 2, 1, 1] + dI[3, 2, 2] + dI[3, 3, 1] + dI[4, 1, 1, 1] + 2*dI[4, 2, 1] + dI[4, 3] + dI[5, 1, 1] + dI[5, 2]
                sage: F = QSym.F()
                sage: dI(F[1,3,1])
                -dI[1, 1, 1, 2] + dI[1, 1, 2, 1] - dI[1, 2, 2] + dI[1, 3, 1]
                sage: F(dI(F([2,1,3])))
                F[2, 1, 3]
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='dI', bracket=False,
                                             category=QSym.Bases())

        def _to_Monomial_on_basis(self, J):
            r"""
            Expand the dual immaculate function labelled by a composition
            ``J`` in the quasi-symmetric monomial basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT: a quasi-symmetric function in the monomial basis

            EXAMPLES::

                sage: dI = QuasiSymmetricFunctions(QQ).dI()
                sage: dI._to_Monomial_on_basis(Composition([1,3]))
                M[1, 1, 1, 1] + M[1, 1, 2] + M[1, 2, 1] + M[1, 3]
                sage: dI._to_Monomial_on_basis(Composition([]))
                M[]
                sage: dI._to_Monomial_on_basis(Composition([2,1,2]))
                4*M[1, 1, 1, 1, 1] + 3*M[1, 1, 1, 2] + 2*M[1, 1, 2, 1] + M[1, 1, 3] + M[1, 2, 1, 1] + M[1, 2, 2] + M[2, 1, 1, 1] + M[2, 1, 2]
            """
            M = self.realization_of().Monomial()
            if not J._list:
                return M([])
            C = Compositions()
            C_size = Compositions(J.size())
            return M.sum_of_terms( ( (C(I), number_of_fCT(C(I), J)) for I in C_size ),
                                   distinct=True )

        @cached_method
        def _matrix_monomial_to_dual_immaculate(self, n):
            r"""
            This function caches the change of basis matrix from the
            quasisymmetric monomial basis to the dual immaculate basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT:

            - A list. Each entry in the list is a row in the
              change-of-basis matrix.

            EXAMPLES::

                sage: dI = QuasiSymmetricFunctions(QQ).dI()
                sage: dI._matrix_monomial_to_dual_immaculate(3)
                [[1, -1, -1, 1], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1]]
                sage: dI._matrix_monomial_to_dual_immaculate(0)
                [[1]]
            """
            N = NonCommutativeSymmetricFunctions(self.base_ring())
            I = N.I()
            S = N.S()
            mat = []
            C = Compositions()
            C_n = Compositions(n)
            for alp in C_n:
                row = []
                expansion = S(I(C(alp)))
                for bet in C_n:
                    row.append(expansion.coefficient(C(bet)))
                mat.append(row)
            return mat

        def _from_Monomial_on_basis(self, J):
            r"""
            Expand the monomial quasi-symmetric function labelled by the
            composition ``J`` in the dual immaculate basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT: a quasi-symmetric function in the dual immaculate basis

            EXAMPLES::

                sage: dI = QuasiSymmetricFunctions(QQ).dI()
                sage: dI._from_Monomial_on_basis(Composition([]))
                dI[]
                sage: dI._from_Monomial_on_basis(Composition([2,1]))
                -dI[1, 1, 1] - dI[1, 2] + dI[2, 1]
                sage: dI._from_Monomial_on_basis(Composition([3,1,2]))
                -dI[1, 1, 1, 1, 1, 1] + dI[1, 1, 1, 1, 2] + dI[1, 1, 1, 3] - dI[1, 1, 4] - dI[1, 2, 1, 1, 1] + dI[1, 2, 3] + dI[2, 1, 1, 1, 1] - dI[2, 1, 1, 2] + dI[2, 2, 1, 1] - dI[2, 2, 2] - dI[3, 1, 1, 1] + dI[3, 1, 2]
            """
            n = J.size()
            C = Compositions()
            C_n = Compositions(n)
            mat = self._matrix_monomial_to_dual_immaculate(n)
            column = C_n.list().index(J)
            return self.sum_of_terms( ( (C(I), mat[C_n.list().index(I)][column])
                                        for I in C_n ),
                                      distinct=True )

    dI = dualImmaculate

    class HazewinkelLambda(CombinatorialFreeModule, BindableClass):
        r"""
        The Hazewinkel lambda basis of the quasi-symmetric functions.

        This basis goes back to [Haz2004]_, albeit it is indexed in a
        different way here. It is a multiplicative basis in a weak
        sense of this word (the product of any two basis elements is a
        basis element, but of course not the one obtained by
        concatenating the indexing compositions).

        In [Haz2004]_, Hazewinkel showed that the `\mathbf{k}`-algebra
        `\mathrm{QSym}` is a polynomial algebra. (The proof is correct
        but rests upon an unproven claim that the lexicographically
        largest term of the `n`-th shuffle power of a Lyndon word is
        the `n`-fold concatenation of this Lyndon word with
        itself, occurring `n!` times in that shuffle power. But this
        can be deduced from Section 2 of [Rad1979]_. See also
        Chapter 6 of [GriRei18]_, specifically Theorem 6.5.13, for a
        complete proof.) More precisely, he showed that
        `\mathrm{QSym}` is generated, as a free commutative
        `\mathbf{k}`-algebra, by the elements `\lambda^n(M_I)`, where
        `n` ranges over the positive integers, and `I` ranges over
        all compositions which are Lyndon words and whose entries
        have gcd `1`. Here, `\lambda^n` denotes the `n`-th lambda
        operation as explained in
        :meth:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Monomial.lambda_of_monomial`.

        Thus, products of these generators form a `\mathbf{k}`-module
        basis of `\mathrm{QSym}`. We index this basis by compositions
        here. More precisely, we define the Hazewinkel lambda basis
        `(\mathrm{HWL}_I)_I` (with `I` ranging over all compositions)
        as follows:

        Let `I` be a composition. Let `I = I_1 I_2 \ldots I_k` be the
        Chen-Fox-Lyndon factorization of `I` (see
        :meth:`~sage.combinat.words.finite_word.FiniteWord_class.lyndon_factorization`).
        For every `j \in \{1, 2, \ldots , k\}`, let `g_j` be the
        gcd of the entries of the Lyndon word `I_j`, and let `J_j` be
        the result of dividing the entries of `I_j` by this gcd. Then,
        `\mathrm{HWL}_I` is defined to be
        `\prod_{j=1}^{k} \lambda^{g_j} (M_{J_j})`.

        .. TODO::

            The conversion from the M basis to the HWL basis is
            currently implemented in the naive way (inverting the
            base-change matrix in the other direction). This matrix
            is not triangular (not even after any permutations of
            the bases), and there could very well be a faster method
            (the one given by Hazewinkel?).

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(ZZ)
            sage: HWL = QSym.HazewinkelLambda()
            sage: M = QSym.M()
            sage: M(HWL([2]))
            M[1, 1]
            sage: M(HWL([1,1]))
            2*M[1, 1] + M[2]
            sage: M(HWL([1,2]))
            M[1, 2]
            sage: M(HWL([2,1]))
            3*M[1, 1, 1] + M[1, 2] + M[2, 1]
            sage: M(HWL(Composition([])))
            M[]
            sage: HWL(M([1,1]))
            HWL[2]
            sage: HWL(M(Composition([2])))
            HWL[1, 1] - 2*HWL[2]
            sage: HWL(M([1]))
            HWL[1]

        TESTS:

        Transforming from the M-basis into the HWL-basis and back
        returns us to where we started::

            sage: all( M(HWL(M[I])) == M[I] for I in Compositions(3) )
            True
            sage: all( HWL(M(HWL[I])) == HWL[I] for I in Compositions(4) )
            True

        Checking the HWL basis elements corresponding to Lyndon
        words::

            sage: all( M(HWL[Composition(I)])
            ....:      == M.lambda_of_monomial([i // gcd(I) for i in I], gcd(I))
            ....:      for I in LyndonWords(e=3, k=2) )
            True
        """

        def __init__(self, QSym):
            r"""
            TESTS::

                sage: HWL = QuasiSymmetricFunctions(QQ).HazewinkelLambda()
                sage: TestSuite(HWL).run()
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='HWL', bracket=False,
                                             category=QSym.Bases())

        def __init_extra__(self):
            """
            Set up caches for the transition maps to and from the monomial
            basis, and register them as coercions.

            TESTS::

                sage: HWL = QuasiSymmetricFunctions(QQ).HazewinkelLambda()
                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: M2HWL = HWL.coerce_map_from(M); M2HWL
                Generic morphism:
                  From: Quasisymmetric functions over the Rational Field in the Monomial basis
                  To:   Quasisymmetric functions over the Rational Field in the HazewinkelLambda basis
                sage: HWL2M = M.coerce_map_from(HWL); HWL2M
                Generic morphism:
                  From: Quasisymmetric functions over the Rational Field in the HazewinkelLambda basis
                  To:   Quasisymmetric functions over the Rational Field in the Monomial basis
                sage: HWL2M(HWL[2])
                M[1, 1]
                sage: M2HWL(M[2])
                HWL[1, 1] - 2*HWL[2]
            """
            M = self.realization_of().M()
            category = self.realization_of()._category
            # This changes Monomial into Hazewinkel Lambda
            M.module_morphism(self._from_Monomial_on_basis,
                                              codomain=self, category=category
                                              ).register_as_coercion()
            # This changes Hazewinkel Lambda into Monomial
            self.module_morphism(self._to_Monomial_on_basis,
                                              codomain=M, category=category
                                              ).register_as_coercion()

            # cache for the coordinates of the elements
            # of the monomial basis with respect to the HWL basis
            self._M_to_self_cache = {}
            # cache for the coordinates of the elements
            # of the HWL basis with respect to the monomial basis
            self._M_from_self_cache = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the monomial basis with respect to the HWL basis
            self._M_transition_matrices = {}
            # cache for transition matrices which contain the coordinates of
            # the elements of the HWL basis with respect to the monomial basis
            self._M_inverse_transition_matrices = {}

        def _precompute_cache(self, n, to_self_cache, from_self_cache, transition_matrices, inverse_transition_matrices, from_self_gen_function):
            r"""
            Compute the transition matrices between ``self`` and the
            monomial basis in the homogeneous components of degree `n`.
            The results are not returned, but rather stored in the caches.

            This assumes that the transition matrices in all degrees smaller
            than `n` have already been computed and cached!

            INPUT:

            - ``n`` -- nonnegative integer
            - ``to_self_cache`` -- a cache which stores the coordinates of
              the elements of the monomial basis with respect to the
              basis ``self``
            - ``from_self_cache`` -- a cache which stores the coordinates
              of the elements of ``self`` with respect to the monomial
              basis
            - ``transition_matrices`` -- a cache for transition matrices
              which contain the coordinates of the elements of the monomial
              basis with respect to ``self``
            - ``inverse_transition_matrices`` -- a cache for transition
              matrices which contain the coordinates of the elements of
              ``self`` with respect to the monomial basis
            - ``from_self_gen_function`` -- a function which takes a
              Lyndon word `I` and returns the Hazewinkel Lambda basis
              element indexed by `I` expanded with respect to the
              monomial basis (as an element of the monomial basis, not as
              a dictionary)

            EXAMPLES:

            The examples below demonstrate how the caches are built
            step by step using the ``_precompute_cache`` method. In order
            not to influence the outcome of other doctests, we make sure
            not to use the caches internally used by this class, but
            rather to create new caches. This allows us to compute the
            transition matrices for a slight variation of the Hazewinkel
            Lambda basis, namely the basis whose `I`-th element is given
            by the simple formula `\prod_{j=1}^{k} M_{I_j}` instead of
            `\prod_{j=1}^{k} \lambda^{g_j} (M_{J_j})`. We will see that
            this is only a `\QQ`-basis rather than a `\ZZ`-basis (a
            reason why the Ditters conjecture took so long to prove)::

                sage: QSym = QuasiSymmetricFunctions(QQ)
                sage: HWL = QSym.HazewinkelLambda()
                sage: M = QSym.M()
                sage: toy_to_self_cache = {}
                sage: toy_from_self_cache = {}
                sage: toy_transition_matrices = {}
                sage: toy_inverse_transition_matrices = {}
                sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
                sage: l(toy_to_self_cache)
                []
                sage: def toy_gen_function(I):
                ....:     return M[I]
                sage: HWL._precompute_cache(0, toy_to_self_cache,
                ....:                          toy_from_self_cache,
                ....:                          toy_transition_matrices,
                ....:                          toy_inverse_transition_matrices,
                ....:                          toy_gen_function)
                sage: l(toy_to_self_cache)
                [([], [([], 1)])]
                sage: HWL._precompute_cache(1, toy_to_self_cache,
                ....:                          toy_from_self_cache,
                ....:                          toy_transition_matrices,
                ....:                          toy_inverse_transition_matrices,
                ....:                          toy_gen_function)
                sage: l(toy_to_self_cache)
                [([], [([], 1)]), ([1], [([1], 1)])]
                sage: HWL._precompute_cache(2, toy_to_self_cache,
                ....:                          toy_from_self_cache,
                ....:                          toy_transition_matrices,
                ....:                          toy_inverse_transition_matrices,
                ....:                          toy_gen_function)
                sage: l(toy_to_self_cache)
                [([], [([], 1)]), ([1], [([1], 1)]), ([1, 1], [([1, 1], 1/2), ([2], -1/2)]), ([2], [([2], 1)])]
                sage: toy_transition_matrices[2]
                [ 1/2 -1/2]
                [   0    1]
                sage: toy_inverse_transition_matrices[2]
                [2 1]
                [0 1]
                sage: sorted(toy_transition_matrices)
                [0, 1, 2]

            As we see from the fractions in the transition matrices, this
            is only a basis over `\QQ`, not over `\ZZ`.

            Let us try another variation on the definition of the basis:
            `\prod_{j=1}^{k} \lambda^{g_j} (M_{J_j})` will now be replaced
            by `\prod_{j=1}^{k} M_{J_j^{g_j}}`, where `J_g` means the
            `g`-fold concatenation of the composition `J` with itself::

                sage: toy_to_self_cache = {}
                sage: toy_from_self_cache = {}
                sage: toy_transition_matrices = {}
                sage: toy_inverse_transition_matrices = {}
                sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
                sage: l(toy_to_self_cache)
                []
                sage: def toy_gen_function(I):
                ....:     xs = [i // gcd(I) for i in I] * gcd(I)
                ....:     return M[xs]
                sage: HWL._precompute_cache(0, toy_to_self_cache,
                ....:                          toy_from_self_cache,
                ....:                          toy_transition_matrices,
                ....:                          toy_inverse_transition_matrices,
                ....:                          toy_gen_function)
                sage: l(toy_to_self_cache)
                [([], [([], 1)])]
                sage: HWL._precompute_cache(1, toy_to_self_cache,
                ....:                          toy_from_self_cache,
                ....:                          toy_transition_matrices,
                ....:                          toy_inverse_transition_matrices,
                ....:                          toy_gen_function)
                sage: l(toy_to_self_cache)
                [([], [([], 1)]), ([1], [([1], 1)])]
                sage: HWL._precompute_cache(2, toy_to_self_cache,
                ....:                          toy_from_self_cache,
                ....:                          toy_transition_matrices,
                ....:                          toy_inverse_transition_matrices,
                ....:                          toy_gen_function)
                sage: l(toy_to_self_cache)
                [([], [([], 1)]), ([1], [([1], 1)]), ([1, 1], [([2], 1)]), ([2], [([1, 1], 1), ([2], -2)])]

            This appears to form another `\ZZ`-basis of
            `\mathrm{QSym}`, but the appearance is deceiving: it
            fails to span the degree-`9` part of `\mathrm{QSym}`.
            (The corresponding computation is not tested as it takes
            a few minutes.) We have not checked if it spans
            `\mathrm{QSym}` over `\QQ`.
            """
            # Much of this code is adapted from sage/combinat/sf/dual.py
            base_ring = self.base_ring()
            zero = base_ring.zero()

            # Handle the n == 0 case separately
            if n == 0:
                part = self._indices([])
                one = base_ring.one()
                to_self_cache[ part ] = { part: one }
                from_self_cache[ part ] = { part: one }
                transition_matrices[n] = matrix(base_ring, [[one]])
                inverse_transition_matrices[n] = matrix(base_ring, [[one]])
                return

            compositions_n = Compositions(n).list()
            len_compositions_n = 2 ** (n-1)     # since n > 0 by now.
            M = self.realization_of().M()

            # The monomial basis will be called M from now on.

            # This contains the data for the transition matrix from the
            # monomial basis M to the Hazewinkel lambda basis self.
            transition_matrix_n = matrix(base_ring, len_compositions_n, len_compositions_n)

            # This first section calculates how the basis elements of the
            # Hazewinkel lambda basis self decompose in the monomial
            # basis M.

            # For every composition I of size n, expand self[I] in terms
            # of the monomial basis M.
            i = 0
            for I in compositions_n:
                # M_coeffs will be M(self[I])._monomial_coefficients
                M_coeffs = {}

                self_I_in_M_basis = M.prod([from_self_gen_function(self._indices(list(J)))
                                            for J in Word(I).lyndon_factorization()])

                j = 0

                for J in compositions_n:
                    if J in self_I_in_M_basis._monomial_coefficients:
                        sp = self_I_in_M_basis._monomial_coefficients[J]
                        M_coeffs[J] = sp
                        transition_matrix_n[i,j] = sp

                    j += 1

                from_self_cache[I] = M_coeffs
                i += 1

            # Save the transition matrix
            inverse_transition_matrices[n] = transition_matrix_n

            # This second section calculates how the basis elements of the
            # monomial basis M expand in terms of the basis B. We do this
            # by inverting the above transition matrix.
            #
            # TODO: The way given in Hazewinkel's [Haz2004]_ paper might
            # be faster.
            inverse_transition = (~transition_matrix_n).change_ring(base_ring)
            # Note that we don't simply write
            # "inverse_transition = ~transition_matrix_n" because that
            # tends to cast the entries of the matrix into a quotient
            # field even if this is unnecessary.
            # MAYBE use .inverse_of_unit() ?

            # TODO: This still looks fragile when the base ring is weird!
            # Possibly work over ZZ in this method?

            for i in range(len_compositions_n):
                self_coeffs = {}
                for j in range(len_compositions_n):
                    if inverse_transition[i,j] != zero:
                        self_coeffs[ compositions_n[j] ] = inverse_transition[i,j]

                to_self_cache[ compositions_n[i] ] = self_coeffs

            transition_matrices[n] = inverse_transition

        def _precompute_M(self, n):
            """
            Compute the transition matrices between ``self`` and the
            monomial basis in the homogeneous components of degree `n`
            (and in those of smaller degree, if not already computed).
            The result is not returned, but rather stored in the cache.

            INPUT:

            - ``n`` -- nonnegative integer

            EXAMPLES:

            The examples below demonstrate how the caches of ``self`` are
            built step by step using the ``_precompute_M`` method. This
            demonstration relies on an untouched Hazewinkel Lambda basis
            (because any computations might have filled the cache already).
            We obtain such a basis by choosing a ground ring unlikely to
            appear elsewhere::

                sage: QSym = QuasiSymmetricFunctions(ZZ['hell', 'yeah'])
                sage: HWL = QSym.HazewinkelLambda()
                sage: M = QSym.M()
                sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
                sage: l(HWL._M_to_self_cache)
                []
                sage: HWL._precompute_M(0)
                sage: l(HWL._M_to_self_cache)
                [([], [([], 1)])]
                sage: HWL._precompute_M(1)
                sage: l(HWL._M_to_self_cache)
                [([], [([], 1)]), ([1], [([1], 1)])]
                sage: HWL._precompute_M(2)
                sage: l(HWL._M_to_self_cache)
                [([], [([], 1)]),
                 ([1], [([1], 1)]),
                 ([1, 1], [([2], 1)]),
                 ([2], [([1, 1], 1), ([2], -2)])]
                sage: HWL._M_transition_matrices[2]
                [ 0  1]
                [ 1 -2]
                sage: HWL._M_inverse_transition_matrices[2]
                [2 1]
                [1 0]
                sage: sorted(HWL._M_transition_matrices)
                [0, 1, 2]

            We did not have to call ``HWL._precompute_M(0)``,
            ``HWL._precompute_M(1)`` and ``HWL._precompute_M(2)``
            in this order; it would be enough to just call
            ``HWL._precompute_M(2)``::

                sage: QSym = QuasiSymmetricFunctions(ZZ['lol', 'wut'])
                sage: HWL = QSym.HazewinkelLambda()
                sage: M = QSym.M()
                sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
                sage: l(HWL._M_to_self_cache)
                []
                sage: HWL._precompute_M(2)
                sage: l(HWL._M_to_self_cache)
                [([], [([], 1)]),
                 ([1], [([1], 1)]),
                 ([1, 1], [([2], 1)]),
                 ([2], [([1, 1], 1), ([2], -2)])]
                sage: HWL._precompute_M(1)
                sage: l(HWL._M_to_self_cache)
                [([], [([], 1)]),
                 ([1], [([1], 1)]),
                 ([1, 1], [([2], 1)]),
                 ([2], [([1, 1], 1), ([2], -2)])]
            """
            l = len(self._M_transition_matrices)
            M = self.realization_of().M()
            if l <= n:
                from sage.misc.cachefunc import cached_function
                from sage.arith.misc import gcd

                @cached_function
                def monolambda(I):
                    # expansion of self[I] in monomial basis,
                    # where I is a composition which is a Lyndon word
                    g = gcd(I)
                    I_reduced = [i // g for i in I]
                    return M.lambda_of_monomial(I_reduced, g)
                for i in range(l, n + 1):
                    self._precompute_cache(i, self._M_to_self_cache,
                                           self._M_from_self_cache,
                                           self._M_transition_matrices,
                                           self._M_inverse_transition_matrices,
                                           monolambda)

        def _to_Monomial_on_basis(self, J):
            r"""
            Expand the Hazewinkel Lambda basis element labelled by a
            composition ``J`` in the quasi-symmetric Monomial basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT: a quasi-symmetric function in the Monomial basis

            EXAMPLES::

                sage: HWL = QuasiSymmetricFunctions(QQ).HazewinkelLambda()
                sage: J = Composition([1, 2, 1])
                sage: HWL._to_Monomial_on_basis(J)
                2*M[1, 1, 2] + M[1, 2, 1] + M[1, 3] + M[2, 2]
            """
            n = sum(J)
            self._precompute_M(n)
            return self.realization_of().M()._from_dict(self._M_from_self_cache[J])

        def _from_Monomial_on_basis(self, J):
            r"""
            Expand the Monomial quasi-symmetric function labelled by the
            composition ``J`` in the Hazewinkel Lambda basis.

            INPUT:

            - ``J`` -- a composition

            OUTPUT: a quasi-symmetric function in the Hazewinkel lambda basis

            EXAMPLES::

                sage: HWL = QuasiSymmetricFunctions(QQ).HazewinkelLambda()
                sage: J = Composition([1, 2, 1])
                sage: HWL._from_Monomial_on_basis(J)
                -2*HWL[1, 1, 2] + HWL[1, 2, 1] - HWL[1, 3] - HWL[2, 2] + 2*HWL[3, 1] - 2*HWL[4]
            """
            n = sum(J)
            self._precompute_M(n)
            return self._from_dict(self._M_to_self_cache[J])

        def product_on_basis(self, I, J):
            """
            The product on Hazewinkel Lambda basis elements.

            The product of the basis elements indexed by two compositions
            `I` and `J` is the basis element obtained by concatenating the
            Lyndon factorizations of the words `I` and `J`, then reordering
            the Lyndon factors in nonincreasing order, and finally
            concatenating them in this order (giving a new composition).

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The product of the Hazewinkel Lambda quasi-symmetric
              functions indexed by ``I`` and ``J``, expressed in the
              Hazewinkel Lambda basis.

            EXAMPLES::

                sage: HWL = QuasiSymmetricFunctions(QQ).HazewinkelLambda()
                sage: c1 = Composition([1, 2, 1])
                sage: c2 = Composition([2, 1, 3, 2])
                sage: HWL.product_on_basis(c1, c2)
                HWL[2, 1, 3, 2, 1, 2, 1]
                sage: HWL.product_on_basis(c1, Composition([]))
                HWL[1, 2, 1]
                sage: HWL.product_on_basis(Composition([]), Composition([]))
                HWL[]

            TESTS::

                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: all( all( M(HWL[I] * HWL[J]) == M(HWL[I]) * M(HWL[J])  # long time
                ....:           for I in Compositions(3) )
                ....:      for J in Compositions(3) )
                True
            """
            from sage.misc.flatten import flatten
            I_factors = [list(i) for i in Word(I).lyndon_factorization()]
            J_factors = [list(j) for j in Word(J).lyndon_factorization()]
            # This uses the convenient fact that comparison of lists in
            # Python works exactly as one wants in Lyndon word theory:
            # [a_1, a_2, ..., a_n] > [b_1, b_2, ..., b_m] if and only if
            # either some i satisfies a_i > b_i and (a_j == b_j for all
            # j < i), or we have n > m and all i <= m satisfy a_i == b_i.
            new_factors = sorted(I_factors + J_factors, reverse=True)
            return self.monomial(self._indices(flatten(new_factors)))

    class psi(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric functions in the `\psi` basis.

        The `\psi` basis is defined as a rescaled Hopf dual of the `\Psi`
        basis of the non-commutative symmetric functions (see Section 3.1
        of [BDHMN2017]_), where the pairing is

        .. MATH::

            (\psi_I, \Psi_J) = z_I \delta_{I,J},

        where `z_I = 1^{m_1} m_1! 2^{m_2} m_2! \cdots` with `m_i` being the
        multiplicity of `i` in the composition `I`. Therefore, we call these
        the *quasi-symmetric power sums of the first kind*.

        Using the duality, we can directly define the `\psi` basis by

        .. MATH::

            \psi_I = \sum_{J \succ I} z_I / \pi_{I,J} M_J,

        where `\pi_{I,J}` is as defined in [NCSF]_.

        The `\psi`-basis is well-defined only when the base ring is a
        `\QQ`-algebra.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: psi = QSym.psi(); psi
            Quasisymmetric functions over the Rational Field in the psi basis
            sage: psi.an_element()
            2*psi[] + 2*psi[1] + 3*psi[1, 1]
            sage: p = SymmetricFunctions(QQ).p()
            sage: psi(p[2,2,1])
            psi[1, 2, 2] + psi[2, 1, 2] + psi[2, 2, 1]
            sage: all(sum(psi(list(al)) for al in Permutations(la))==psi(p(la)) for la in Partitions(6))
            True
            sage: p = SymmetricFunctions(QQ).p()
            sage: psi(p[3,2,2])
            psi[2, 2, 3] + psi[2, 3, 2] + psi[3, 2, 2]

        Checking the equivalent definition of `\psi_n`::

            sage: def test_psi(n):
            ....:     psi = QuasiSymmetricFunctions(QQ).psi()
            ....:     Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
            ....:     M = matrix([[psi[I].duality_pairing(Psi[J])
            ....:                  for I in Compositions(n)]
            ....:                 for J in Compositions(n)])
            ....:     def z(J): return J.to_partition().centralizer_size()
            ....:     return M == matrix.diagonal([z(I) for I in Compositions(n)])
            sage: all(test_psi(k) for k in range(1,5))
            True
        """

        def __init__(self, QSym):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: psi = QuasiSymmetricFunctions(QQ).psi()
                sage: TestSuite(psi).run()

            TESTS::

                sage: psi = QuasiSymmetricFunctions(QQ).psi()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: all(psi(M(psi[c])) == psi[c] for n in range(5)
                ....:     for c in Compositions(n))
                True
                sage: all(M(psi(M[c])) == M[c] for n in range(5)
                ....:     for c in Compositions(n))
                True
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='psi', bracket=False,
                                             category=QSym.Bases())

            category = self.realization_of()._category
            Monomial = self.realization_of().Monomial()
            self.module_morphism(self._to_Monomial_on_basis,
                                 codomain=Monomial, category=category
                                 ).register_as_coercion()
            Monomial.module_morphism(self._from_Monomial_on_basis,
                                     codomain=self, category=category
                                     ).register_as_coercion()

        def _from_Monomial_on_basis(self, I):
            r"""
            Expand a Monomial basis element indexed by ``I`` in the
            `\psi` basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT: a quasi-symmetric function in the `\psi` basis

            TESTS::

                sage: psi = QuasiSymmetricFunctions(QQ).psi()
                sage: I = Composition([2, 3, 2])
                sage: psi._from_Monomial_on_basis(I)
                1/2*psi[2, 3, 2] - 2/5*psi[2, 5] - 3/5*psi[5, 2] + 2/7*psi[7]
            """
            R = self.base_ring()
            minus_one = -R.one()

            def z(J):
                return R(J.to_partition().centralizer_size())
            return self._from_dict({J: minus_one**(len(I)-len(J)) / z(J) * coeff_lp(I, J)
                                    for J in I.fatter()})

        def _to_Monomial_on_basis(self, I):
            r"""
            Expand a `\psi` basis element indexed by ``I`` in the
            Monomial basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT: a quasi-symmetric function in the Monomial basis

            TESTS::

                sage: psi = QuasiSymmetricFunctions(QQ).psi()
                sage: I = Composition([2, 3, 2])
                sage: psi._to_Monomial_on_basis(I)
                2*M[2, 3, 2] + 4/5*M[2, 5] + 6/5*M[5, 2] + 12/35*M[7]
            """
            R = self.base_ring()
            z = R(I.to_partition().centralizer_size())
            Monomial = self.realization_of().Monomial()
            return Monomial._from_dict({J: z / coeff_pi(I,J) for J in I.fatter()})

    class phi(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of quasi-symmetric functions in the `\phi` basis.

        The `\phi` basis is defined as a rescaled Hopf dual of the `\Phi`
        basis of the non-commutative symmetric functions (see Section 3.1
        of [BDHMN2017]_), where the pairing is

        .. MATH::

            (\phi_I, \Phi_J) = z_I \delta_{I,J},

        where `z_I = 1^{m_1} m_1! 2^{m_2} m_2! \cdots` with `m_i` being the
        multiplicity of `i` in the composition `I`. Therefore, we call these
        the *quasi-symmetric power sums of the second kind*.

        Using the duality, we can directly define the `\phi` basis by

        .. MATH::

            \phi_I = \sum_{J \succ I} z_I / sp_{I,J} M_J,

        where `sp_{I,J}` is as defined in [NCSF]_.

        The `\phi`-basis is well-defined only when the base ring is a
        `\QQ`-algebra.

        EXAMPLES::

            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: phi = QSym.phi(); phi
            Quasisymmetric functions over the Rational Field in the phi basis
            sage: phi.an_element()
            2*phi[] + 2*phi[1] + 3*phi[1, 1]
            sage: p = SymmetricFunctions(QQ).p()
            sage: phi(p[2,2,1])
            phi[1, 2, 2] + phi[2, 1, 2] + phi[2, 2, 1]
            sage: all(sum(phi(list(al)) for al in Permutations(la))==phi(p(la)) for la in Partitions(6))
            True
            sage: p = SymmetricFunctions(QQ).p()
            sage: phi(p[3,2,2])
            phi[2, 2, 3] + phi[2, 3, 2] + phi[3, 2, 2]

        Checking the equivalent definition of `\phi_n`::

            sage: def test_phi(n):
            ....:     phi = QuasiSymmetricFunctions(QQ).phi()
            ....:     Phi = NonCommutativeSymmetricFunctions(QQ).Phi()
            ....:     M = matrix([[phi[I].duality_pairing(Phi[J])
            ....:                  for I in Compositions(n)]
            ....:                 for J in Compositions(n)])
            ....:     def z(J): return J.to_partition().centralizer_size()
            ....:     return M == matrix.diagonal([z(I) for I in Compositions(n)])
            sage: all(test_phi(k) for k in range(1,5))
            True
        """

        def __init__(self, QSym):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: phi = QuasiSymmetricFunctions(QQ).phi()
                sage: TestSuite(phi).run()

            TESTS::

                sage: phi = QuasiSymmetricFunctions(QQ).phi()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: all(phi(M(phi[c])) == phi[c] for n in range(5)
                ....:     for c in Compositions(n))
                True
                sage: all(M(phi(M[c])) == M[c] for n in range(5)
                ....:     for c in Compositions(n))
                True
            """
            CombinatorialFreeModule.__init__(self, QSym.base_ring(), Compositions(),
                                             prefix='phi', bracket=False,
                                             category=QSym.Bases())

            category = self.realization_of()._category
            Monomial = self.realization_of().Monomial()
            self.module_morphism(self._to_Monomial_on_basis,
                                 codomain=Monomial, category=category
                                 ).register_as_coercion()
            Monomial.module_morphism(self._from_Monomial_on_basis,
                                     codomain=self, category=category
                                     ).register_as_coercion()

        def _from_Monomial_on_basis(self, I):
            r"""
            Expand a Monomial basis element indexed by ``I`` in the
            `\phi` basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT: a quasi-symmetric function in the `\psi` basis

            TESTS::

                sage: phi = QuasiSymmetricFunctions(QQ).phi()
                sage: I = Composition([3, 2, 2])
                sage: phi._from_Monomial_on_basis(I)
                1/2*phi[3, 2, 2] - 1/2*phi[3, 4] - 1/2*phi[5, 2] + 1/3*phi[7]
            """
            R = self.base_ring()
            minus_one = -R.one()

            def z(J):
                return R(J.to_partition().centralizer_size())
            return self._from_dict({J: minus_one**(len(I)-len(J)) * R.prod(J) / (coeff_ell(I, J) * z(J))
                                    for J in I.fatter()})

        def _to_Monomial_on_basis(self, I):
            r"""
            Expand a `\phi` basis element indexed by ``I`` in the
            Monomial basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT: a quasi-symmetric function in the Monomial basis

            TESTS::

                sage: phi = QuasiSymmetricFunctions(QQ).phi()
                sage: I = Composition([3, 2, 2])
                sage: phi._to_Monomial_on_basis(I)
                2*M[3, 2, 2] + M[3, 4] + M[5, 2] + 1/3*M[7]
            """
            R = self.base_ring()
            z = R(I.to_partition().centralizer_size())
            Monomial = self.realization_of().Monomial()
            return Monomial._from_dict({J: z / coeff_sp(I,J) for J in I.fatter()})
