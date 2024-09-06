r"""
Free modules of finite rank

The class :class:`FiniteRankFreeModule` implements free modules of finite rank
over a commutative ring.

A *free module of finite rank* over a commutative ring `R` is a module `M` over
`R` that admits a *finite basis*, i.e. a finite family of linearly independent
generators. Since `R` is commutative, it has the invariant basis number
property, so that the rank of the free module `M` is defined uniquely, as the
cardinality of any basis of `M`.

No distinguished basis of `M` is assumed. On the contrary, many bases can be
introduced on the free module along with change-of-basis rules (as module
automorphisms). Each
module element has then various representations over the various bases.

.. NOTE::

    The class :class:`FiniteRankFreeModule` does not inherit from
    class :class:`~sage.modules.free_module.FreeModule_generic`
    nor from class
    :class:`~sage.combinat.free_module.CombinatorialFreeModule`, since
    both classes deal with modules with a *distinguished basis* (see
    details :ref:`below <diff-FreeModule>`). Accordingly, the class
    :class:`FiniteRankFreeModule` inherits directly from the generic class
    :class:`~sage.structure.parent.Parent` with the category set to
    :class:`~sage.categories.modules.Modules` (and not to
    :class:`~sage.categories.modules_with_basis.ModulesWithBasis`).

.. TODO::

    - implement submodules
    - create a FreeModules category (cf. the *TODO* statement in the
      documentation of :class:`~sage.categories.modules.Modules`: *Implement
      a ``FreeModules(R)`` category, when so prompted by a concrete use case*)

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): category set to ``Modules(ring).FiniteDimensional()``
  (:issue:`20770`)
- Michael Jung (2019): improve treatment of the zero element
- Eric Gourgoulhon (2021): unicode symbols for tensor and exterior products
- Matthias Koeppe (2022): ``FiniteRankFreeModule_abstract``, symmetric powers

REFERENCES:

- Chap. 10 of R. Godement : *Algebra* [God1968]_
- Chap. 3 of S. Lang : *Algebra* [Lan2002]_

EXAMPLES:

Let us define a free module of rank 2 over `\ZZ`::

    sage: M = FiniteRankFreeModule(ZZ, 2, name='M') ; M
    Rank-2 free module M over the Integer Ring
    sage: M.category()
    Category of finite dimensional modules over Integer Ring

We introduce a first basis on ``M``::

    sage: e = M.basis('e') ; e
    Basis (e_0,e_1) on the Rank-2 free module M over the Integer Ring

The elements of the basis are of course module elements::

    sage: e[0]
    Element e_0 of the Rank-2 free module M over the Integer Ring
    sage: e[1]
    Element e_1 of the Rank-2 free module M over the Integer Ring
    sage: e[0].parent()
    Rank-2 free module M over the Integer Ring

We define a module element by its components w.r.t. basis ``e``::

    sage: u = M([2,-3], basis=e, name='u')
    sage: u.display(e)
    u = 2 e_0 - 3 e_1

Module elements can be also be created by arithmetic expressions::

    sage: v = -2*u + 4*e[0] ; v
    Element of the Rank-2 free module M over the Integer Ring
    sage: v.display(e)
    6 e_1
    sage: u == 2*e[0] - 3*e[1]
    True

We define a second basis on ``M`` from a family of linearly independent
elements::

    sage: f = M.basis('f', from_family=(e[0]-e[1], -2*e[0]+3*e[1])) ; f
    Basis (f_0,f_1) on the Rank-2 free module M over the Integer Ring
    sage: f[0].display(e)
    f_0 = e_0 - e_1
    sage: f[1].display(e)
    f_1 = -2 e_0 + 3 e_1

We may of course express the elements of basis ``e`` in terms of basis ``f``::

    sage: e[0].display(f)
    e_0 = 3 f_0 + f_1
    sage: e[1].display(f)
    e_1 = 2 f_0 + f_1

as well as any module element::

    sage: u.display(f)
    u = -f_1
    sage: v.display(f)
    12 f_0 + 6 f_1

The two bases are related by a module automorphism::

    sage: a = M.change_of_basis(e,f) ; a
    Automorphism of the Rank-2 free module M over the Integer Ring
    sage: a.parent()
    General linear group of the Rank-2 free module M over the Integer Ring
    sage: a.matrix(e)
    [ 1 -2]
    [-1  3]

Let us check that basis ``f`` is indeed the image of basis ``e`` by ``a``::

    sage: f[0] == a(e[0])
    True
    sage: f[1] == a(e[1])
    True

The reverse change of basis is of course the inverse automorphism::

    sage: M.change_of_basis(f,e) == a^(-1)
    True

We introduce a new module element via its components w.r.t. basis ``f``::

    sage: v = M([2,4], basis=f, name='v')
    sage: v.display(f)
    v = 2 f_0 + 4 f_1

The sum of the two module elements ``u`` and ``v`` can be performed even if
they have been defined on different bases, thanks to the known relation
between the two bases::

    sage: s = u + v ; s
    Element u+v of the Rank-2 free module M over the Integer Ring

We can display the result in either basis::

    sage: s.display(e)
    u+v = -4 e_0 + 7 e_1
    sage: s.display(f)
    u+v = 2 f_0 + 3 f_1

Tensor products of elements are implemented::

    sage: t = u*v ; t
    Type-(2,0) tensor u⊗v on the Rank-2 free module M over the Integer Ring
    sage: t.parent()
    Free module of type-(2,0) tensors on the
     Rank-2 free module M over the Integer Ring
    sage: t.display(e)
    u⊗v = -12 e_0⊗e_0 + 20 e_0⊗e_1 + 18 e_1⊗e_0 - 30 e_1⊗e_1
    sage: t.display(f)
    u⊗v = -2 f_1⊗f_0 - 4 f_1⊗f_1

We can access to tensor components w.r.t. to a given basis via the square
bracket operator::

    sage: t[e,0,1]
    20
    sage: t[f,1,0]
    -2
    sage: u[e,0]
    2
    sage: u[e,:]
    [2, -3]
    sage: u[f,:]
    [0, -1]

The parent of the automorphism ``a`` is the group `\mathrm{GL}(M)`, but
``a`` can also be considered as a tensor of type `(1,1)` on ``M``::

    sage: a.parent()
    General linear group of the Rank-2 free module M over the Integer Ring
    sage: a.tensor_type()
    (1, 1)
    sage: a.display(e)
    e_0⊗e^0 - 2 e_0⊗e^1 - e_1⊗e^0 + 3 e_1⊗e^1
    sage: a.display(f)
    f_0⊗f^0 - 2 f_0⊗f^1 - f_1⊗f^0 + 3 f_1⊗f^1

As such, we can form its tensor product with ``t``, yielding a tensor of
type `(3,1)`::

    sage: t*a
    Type-(3,1) tensor on the Rank-2 free module M over the Integer Ring
    sage: (t*a).display(e)
    -12 e_0⊗e_0⊗e_0⊗e^0 + 24 e_0⊗e_0⊗e_0⊗e^1 + 12 e_0⊗e_0⊗e_1⊗e^0
     - 36 e_0⊗e_0⊗e_1⊗e^1 + 20 e_0⊗e_1⊗e_0⊗e^0 - 40 e_0⊗e_1⊗e_0⊗e^1
     - 20 e_0⊗e_1⊗e_1⊗e^0 + 60 e_0⊗e_1⊗e_1⊗e^1 + 18 e_1⊗e_0⊗e_0⊗e^0
     - 36 e_1⊗e_0⊗e_0⊗e^1 - 18 e_1⊗e_0⊗e_1⊗e^0 + 54 e_1⊗e_0⊗e_1⊗e^1
     - 30 e_1⊗e_1⊗e_0⊗e^0 + 60 e_1⊗e_1⊗e_0⊗e^1 + 30 e_1⊗e_1⊗e_1⊗e^0
     - 90 e_1⊗e_1⊗e_1⊗e^1

The parent of `t\otimes a` is itself a free module of finite rank over `\ZZ`::

        sage: T = (t*a).parent() ; T
        Free module of type-(3,1) tensors on the Rank-2 free module M over the
         Integer Ring
        sage: T.base_ring()
        Integer Ring
        sage: T.rank()
        16

.. _diff-FreeModule:

.. RUBRIC:: Differences between ``FiniteRankFreeModule`` and ``FreeModule``
  (or ``VectorSpace``)

To illustrate the differences, let us create two free modules of rank 3 over
`\ZZ`, one with ``FiniteRankFreeModule`` and the other one with
``FreeModule``::

    sage: M = FiniteRankFreeModule(ZZ, 3, name='M') ; M
    Rank-3 free module M over the Integer Ring
    sage: N = FreeModule(ZZ, 3) ; N
    Ambient free module of rank 3 over the principal ideal domain Integer Ring

The main difference is that ``FreeModule`` returns a free module with a
distinguished basis, while ``FiniteRankFreeModule`` does not::

    sage: N.basis()
    [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]
    sage: M.bases()
    []
    sage: M.print_bases()
    No basis has been defined on the Rank-3 free module M over the Integer Ring

This is also revealed by the category of each module::

    sage: M.category()
    Category of finite dimensional modules over Integer Ring
    sage: N.category()
    Category of finite dimensional modules with basis over
     (Dedekind domains and euclidean domains and noetherian rings
     and infinite enumerated sets and metric spaces)

In other words, the module created by ``FreeModule`` is actually `\ZZ^3`,
while, in the absence of any distinguished basis, no *canonical* isomorphism
relates the module created by ``FiniteRankFreeModule`` to `\ZZ^3`::

    sage: N is ZZ^3
    True
    sage: M is ZZ^3
    False
    sage: M == ZZ^3
    False

Because it is `\ZZ^3`, ``N`` is unique, while there may be various modules
of the same rank over the same ring created by ``FiniteRankFreeModule``;
they are then distinguished by their names (actually by the complete
sequence of arguments of ``FiniteRankFreeModule``)::

    sage: N1 = FreeModule(ZZ, 3) ; N1
    Ambient free module of rank 3 over the principal ideal domain Integer Ring
    sage: N1 is N  # FreeModule(ZZ, 3) is unique
    True
    sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M_1') ; M1
    Rank-3 free module M_1 over the Integer Ring
    sage: M1 is M  # M1 and M are different rank-3 modules over ZZ
    False
    sage: M1b = FiniteRankFreeModule(ZZ, 3, name='M_1') ; M1b
    Rank-3 free module M_1 over the Integer Ring
    sage: M1b is M1  # because M1b and M1 have the same name
    True

As illustrated above, various bases can be introduced on the module created by
``FiniteRankFreeModule``::

    sage: e = M.basis('e') ; e
    Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
    sage: f = M.basis('f', from_family=(-e[0], e[1]-e[2], -2*e[1]+3*e[2])) ; f
    Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring
    sage: M.bases()
    [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
     Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]

Each element of a basis is accessible via its index::

    sage: e[0]
    Element e_0 of the Rank-3 free module M over the Integer Ring
    sage: e[0].parent()
    Rank-3 free module M over the Integer Ring
    sage: f[1]
    Element f_1 of the Rank-3 free module M over the Integer Ring
    sage: f[1].parent()
    Rank-3 free module M over the Integer Ring

while on module ``N``, the element of the (unique) basis is accessible
directly from the module symbol::

    sage: N.0
    (1, 0, 0)
    sage: N.1
    (0, 1, 0)
    sage: N.0.parent()
    Ambient free module of rank 3 over the principal ideal domain Integer Ring

The arithmetic of elements is similar; the difference lies in the display:
a basis has to be specified for elements of ``M``, while elements of ``N`` are
displayed directly as elements of `\ZZ^3`::

    sage: u = 2*e[0] - 3*e[2] ; u
    Element of the Rank-3 free module M over the Integer Ring
    sage: u.display(e)
    2 e_0 - 3 e_2
    sage: u.display(f)
    -2 f_0 - 6 f_1 - 3 f_2
    sage: u[e,:]
    [2, 0, -3]
    sage: u[f,:]
    [-2, -6, -3]
    sage: v = 2*N.0 - 3*N.2 ; v
    (2, 0, -3)

For the case of ``M``, in order to avoid to specify the basis if the user is
always working with the same basis (e.g. only one basis has been defined),
the concept of *default basis* has been introduced::

    sage: M.default_basis()
    Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
    sage: M.print_bases()
    Bases defined on the Rank-3 free module M over the Integer Ring:
     - (e_0,e_1,e_2) (default basis)
     - (f_0,f_1,f_2)

This is different from the *distinguished basis* of ``N``: it simply means that
the mention of the basis can be omitted in function arguments::

    sage: u.display()  # equivalent to u.display(e)
    2 e_0 - 3 e_2
    sage: u[:]         # equivalent to u[e,:]
    [2, 0, -3]

At any time, the default basis can be changed::

    sage: M.set_default_basis(f)
    sage: u.display()
    -2 f_0 - 6 f_1 - 3 f_2

Another difference between ``FiniteRankFreeModule`` and ``FreeModule`` is that
for the former the range of indices can be specified (by default, it starts
from 0)::

    sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1) ; M
    Rank-3 free module M over the Integer Ring
    sage: e = M.basis('e') ; e  # compare with (e_0,e_1,e_2) above
    Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
    sage: e[1], e[2], e[3]
    (Element e_1 of the Rank-3 free module M over the Integer Ring,
     Element e_2 of the Rank-3 free module M over the Integer Ring,
     Element e_3 of the Rank-3 free module M over the Integer Ring)

All the above holds for ``VectorSpace`` instead of ``FreeModule``: the object
created by ``VectorSpace`` is actually a Cartesian power of the base field::

    sage: V = VectorSpace(QQ,3) ; V
    Vector space of dimension 3 over Rational Field
    sage: V.category()
    Category of finite dimensional vector spaces with basis
     over (number fields and quotient fields and metric spaces)
    sage: V is QQ^3
    True
    sage: V.basis()
    [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]

To create a vector space without any distinguished basis, one has to use
``FiniteRankFreeModule``::

    sage: V = FiniteRankFreeModule(QQ, 3, name='V') ; V
    3-dimensional vector space V over the Rational Field
    sage: V.category()
    Category of finite dimensional vector spaces over Rational Field
    sage: V.bases()
    []
    sage: V.print_bases()
    No basis has been defined on the 3-dimensional vector space V over the
     Rational Field

The class :class:`FiniteRankFreeModule` has been created for the needs
of the `SageManifolds project <http://sagemanifolds.obspm.fr/>`_, where
free modules do not have any distinguished basis. Too kinds of free modules
occur in the context of differentiable manifolds (see
`here <http://sagemanifolds.obspm.fr/tensor_modules.html>`_ for more
details):

- the tangent vector space at any point of the manifold (cf.
  :class:`~sage.manifolds.differentiable.tangent_space.TangentSpace`);
- the set of vector fields on a parallelizable open subset `U` of the manifold,
  which is a free module over the algebra of scalar fields on `U` (cf.
  :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`).

For instance, without any specific coordinate choice, no basis can be
distinguished in a tangent space.

On the other side, the modules created by ``FreeModule`` have much more
algebraic functionalities than those created by ``FiniteRankFreeModule``. In
particular, submodules have not been implemented yet in
:class:`FiniteRankFreeModule`. Moreover, modules resulting from ``FreeModule``
are tailored to the specific kind of their base ring:

- free module over a commutative ring that is not an integral domain
  (`\ZZ/6\ZZ`)::

    sage: R = IntegerModRing(6) ; R
    Ring of integers modulo 6
    sage: FreeModule(R, 3)
    Ambient free module of rank 3 over Ring of integers modulo 6
    sage: type(FreeModule(R, 3))
    <class 'sage.modules.free_module.FreeModule_ambient_with_category'>

- free module over an integral domain that is not principal (`\ZZ[X]`)::

    sage: R.<X> = ZZ[] ; R
    Univariate Polynomial Ring in X over Integer Ring
    sage: FreeModule(R, 3)
    Ambient free module of rank 3 over the integral domain Univariate
     Polynomial Ring in X over Integer Ring
    sage: type(FreeModule(R, 3))
    <class 'sage.modules.free_module.FreeModule_ambient_domain_with_category'>

- free module over a principal ideal domain (`\ZZ`)::

    sage: R = ZZ ; R
    Integer Ring
    sage: FreeModule(R,3)
    Ambient free module of rank 3 over the principal ideal domain Integer Ring
    sage: type(FreeModule(R, 3))
    <class 'sage.modules.free_module.FreeModule_ambient_pid_with_category'>

On the contrary, all objects constructed with ``FiniteRankFreeModule`` belong
to the same class::

    sage: R = IntegerModRing(6)
    sage: type(FiniteRankFreeModule(R, 3))
    <class 'sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule_with_category'>
    sage: R.<X> = ZZ[]
    sage: type(FiniteRankFreeModule(R, 3))
    <class 'sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule_with_category'>
    sage: R = ZZ
    sage: type(FiniteRankFreeModule(R, 3))
    <class 'sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule_with_category'>


.. RUBRIC:: Differences between ``FiniteRankFreeModule`` and
  ``CombinatorialFreeModule``

An alternative to construct free modules in Sage is
:class:`~sage.combinat.free_module.CombinatorialFreeModule`.
However, as ``FreeModule``, it leads to a module with a distinguished basis::

    sage: N = CombinatorialFreeModule(ZZ, [1,2,3]) ; N
    Free module generated by {1, 2, 3} over Integer Ring
    sage: N.category()
    Category of finite dimensional modules with basis over Integer Ring

The distinguished basis is returned by the method ``basis()``::

    sage: b = N.basis() ; b
    Finite family {1: B[1], 2: B[2], 3: B[3]}
    sage: b[1]
    B[1]
    sage: b[1].parent()
    Free module generated by {1, 2, 3} over Integer Ring

For the free module ``M`` created above with ``FiniteRankFreeModule``, the
method ``basis`` has at least one argument: the symbol string that
specifies which basis is required::

    sage: e = M.basis('e') ; e
    Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
    sage: e[1]
    Element e_1 of the Rank-3 free module M over the Integer Ring
    sage: e[1].parent()
    Rank-3 free module M over the Integer Ring

The arithmetic of elements is similar::

    sage: u = 2*e[1] - 5*e[3] ; u
    Element of the Rank-3 free module M over the Integer Ring
    sage: v = 2*b[1] - 5*b[3] ; v
    2*B[1] - 5*B[3]

One notices that elements of ``N`` are displayed directly in terms of their
expansions on the distinguished basis. For elements of ``M``, one has to use
the method
:meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.display`
in order to specify the basis::

    sage: u.display(e)
    2 e_1 - 5 e_3

The components on the basis are returned by the square bracket operator for
``M`` and by the method ``coefficient`` for ``N``::

    sage: [u[e,i] for i in {1,2,3}]
    [2, 0, -5]
    sage: u[e,:]  # a shortcut for the above
    [2, 0, -5]
    sage: [v.coefficient(i) for i in {1,2,3}]
    [2, 0, -5]
"""
# ******************************************************************************
#       Copyright (C) 2014-2021 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#                     2014-2016 Travis Scrimshaw <tscrimsh@umn.edu>
#                     2015      Michal Bejger <bejger@camk.edu.pl>
#                     2016      Frédéric Chapoton
#                     2020      Michael Jung
#                     2020-2022 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************
from __future__ import annotations

from typing import Generator, Optional

from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.categories.modules import Modules
from sage.categories.morphism import SetIsomorphism
from sage.categories.rings import Rings
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.sets.family import Family, TrivialFamily
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.reflexive_module import (
    ReflexiveModule_abstract,
    ReflexiveModule_base,
    ReflexiveModule_dual,
)


class FiniteRankFreeModule_abstract(UniqueRepresentation, ReflexiveModule_abstract):
    r"""
    Abstract base class for free modules of finite rank over a commutative ring.
    """

    def __init__(
        self,
        ring,
        rank,
        name=None,
        latex_name=None,
        category=None,
        ambient=None,
    ):
        r"""
        See :class:`FiniteRankFreeModule` for documentation and examples.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: TestSuite(M).run()
            sage: e = M.basis('e')
            sage: TestSuite(M).run()
            sage: f = M.basis('f')
            sage: TestSuite(M).run()
        """
        # This duplicates the normalization done in __classcall_private__,
        # but it is needed for various subclasses.
        if ring not in Rings().Commutative():
            raise TypeError("the module base ring must be commutative")
        category = Modules(ring).FiniteDimensional().or_subcategory(category)
        Parent.__init__(self, base=ring, category=category)
        self._ring = ring # same as self._base
        if ambient is None:
            self._ambient_module = self
        else:
            self._ambient_module = ambient
        self._rank = rank
        self._name = name
        # This duplicates the normalization done in __classcall_private__,
        # but it is needed for various subclasses.
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M._latex_()
            'M'
            sage: latex(M)
            M
            sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', latex_name=r'\mathcal{M}')
            sage: M1._latex_()
            '\\mathcal{M}'
            sage: latex(M1)
            \mathcal{M}
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
            return self._latex_name

    def rank(self) -> int:
        r"""
        Return the rank of the free module ``self``.

        Since the ring over which ``self`` is built is assumed to be
        commutative (and hence has the invariant basis number property), the
        rank is defined uniquely, as the cardinality of any basis of ``self``.

        EXAMPLES:

        Rank of free modules over `\ZZ`::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: M.rank()
            3
            sage: M.tensor_module(0,1).rank()
            3
            sage: M.tensor_module(0,2).rank()
            9
            sage: M.tensor_module(1,0).rank()
            3
            sage: M.tensor_module(1,1).rank()
            9
            sage: M.tensor_module(1,2).rank()
            27
            sage: M.tensor_module(2,2).rank()
            81
        """
        return self._rank

    @cached_method
    def zero(self):
        r"""
        Return the zero element of ``self``.

        EXAMPLES:

        Zero elements of free modules over `\ZZ`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.zero()
            Element zero of the Rank-3 free module M over the Integer Ring
            sage: M.zero().parent() is M
            True
            sage: M.zero() is M(0)
            True
            sage: T = M.tensor_module(1,1)
            sage: T.zero()
            Type-(1,1) tensor zero on the Rank-3 free module M over the Integer Ring
            sage: T.zero().parent() is T
            True
            sage: T.zero() is T(0)
            True

        Components of the zero element with respect to some basis::

            sage: e = M.basis('e')
            sage: M.zero()[e,:]
            [0, 0, 0]
            sage: all(M.zero()[e,i] == M.base_ring().zero() for i in M.irange())
            True
            sage: T.zero()[e,:]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: M.tensor_module(1,2).zero()[e,:]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, 0, 0], [0, 0, 0], [0, 0, 0]]]
        """
        resu = self._element_constructor_(name='zero', latex_name='0')
        for basis in self._known_bases:
            resu._add_comp_unsafe(basis)
            # (since new components are initialized to zero)
        resu._is_zero = True # This element is certainly zero
        resu.set_immutable()
        return resu

    def ambient_module(self): # compatible with sage.modules.free_module.FreeModule_generic
        """
        Return the ambient module associated to this module.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.ambient_module() is M
            True

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym0123x45M = M.tensor_module(6, 0, sym=((0, 1, 2, 3), (4, 5)))
            sage: T60M = M.tensor_module(6, 0)
            sage: Sym0123x45M.ambient_module() is T60M
            True
        """
        return self._ambient_module

    ambient = ambient_module # compatible with sage.modules.with_basis.subquotient.SubmoduleWithBasis

    def is_submodule(self, other):
        """
        Return ``True`` if ``self`` is a submodule of ``other``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 4, name='N')
            sage: M.is_submodule(M)
            True
            sage: M.is_submodule(N)
            False
        """
        return self == other or self.ambient_module() == other

    def isomorphism_with_fixed_basis(self, basis=None, codomain=None):
        r"""
        Construct the canonical isomorphism from the free module ``self``
        to a free module in which ``basis`` of ``self`` is mapped to the
        distinguished basis of ``codomain``.

        INPUT:

        - ``basis`` -- (default: ``None``) the basis of ``self`` which
          should be mapped to the distinguished basis on ``codomain``;
          if ``None``, the default basis is assumed.
        - ``codomain`` -- (default: ``None``) the codomain of the
          isomorphism represented by a free module within the category
          :class:`~sage.categories.modules_with_basis.ModulesWithBasis` with
          the same rank and base ring as ``self``; if ``None``, a free module
          represented by
          :class:`~sage.combinat.free_module.CombinatorialFreeModule` is
          constructed

        OUTPUT:

        - a module isomorphism represented by
          :class:`~sage.categories.morphism.SetIsomorphism`

        EXAMPLES::

            sage: V = FiniteRankFreeModule(QQ, 3, start_index=1); V
            3-dimensional vector space over the Rational Field
            sage: basis = e = V.basis("e"); basis
            Basis (e_1,e_2,e_3) on the 3-dimensional vector space over the
             Rational Field
            sage: phi_e = V.isomorphism_with_fixed_basis(basis); phi_e
            Generic morphism:
              From: 3-dimensional vector space over the Rational Field
              To:   Free module generated by {1, 2, 3} over Rational Field
            sage: phi_e.codomain().category()
            Category of finite dimensional vector spaces with basis over
             Rational Field
            sage: phi_e(e[1] + 2 * e[2])
            e[1] + 2*e[2]

            sage: abc = V.basis(['a', 'b', 'c'], symbol_dual=['d', 'e', 'f']); abc
            Basis (a,b,c) on the 3-dimensional vector space over the Rational Field
            sage: phi_abc = V.isomorphism_with_fixed_basis(abc); phi_abc
            Generic morphism:
            From: 3-dimensional vector space over the Rational Field
            To:   Free module generated by {1, 2, 3} over Rational Field
            sage: phi_abc(abc[1] + 2 * abc[2])
            B[1] + 2*B[2]

        Providing a codomain::

            sage: W = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: phi_eW = V.isomorphism_with_fixed_basis(basis, codomain=W); phi_eW
            Generic morphism:
            From: 3-dimensional vector space over the Rational Field
            To:   Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: phi_eW(e[1] + 2 * e[2])
            B['a'] + 2*B['b']
            sage: ~phi_eW
            Generic morphism:
            From: Free module generated by {'a', 'b', 'c'} over Rational Field
            To:   3-dimensional vector space over the Rational Field
            sage: (~phi_eW)(W.basis()['b']).display()
            e_2

        Providing a :class:`~sage.modules.free_module.Module_free_ambient` as the codomain::

            sage: W = QQ^3
            sage: phi_eW = V.isomorphism_with_fixed_basis(basis, codomain=W); phi_eW
            Generic morphism:
            From: 3-dimensional vector space over the Rational Field
            To:   Vector space of dimension 3 over Rational Field
            sage: phi_eW(e[1] + 2 * e[2])
            (1, 2, 0)

        Sending (1,1)-tensors to matrices::

            sage: T11 = V.tensor_module(1, 1); T11
            Free module of type-(1,1) tensors on the
             3-dimensional vector space over the Rational Field
            sage: e_T11 = T11.basis("e"); e_T11
            Standard basis on the
             Free module of type-(1,1) tensors on the 3-dimensional vector space over the Rational Field
             induced by Basis (e_1,e_2,e_3) on the 3-dimensional vector space over the Rational Field
            sage: W = MatrixSpace(QQ, 3)
            sage: phi_e_T11 = T11.isomorphism_with_fixed_basis(e_T11, codomain=W); phi_e_T11
            Generic morphism:
            From: Free module of type-(1,1) tensors on the
                  3-dimensional vector space over the Rational Field
            To:   Full MatrixSpace of 3 by 3 dense matrices over Rational Field
            sage: t = T11.an_element(); t.display()
            1/2 e_1⊗e^1
            sage: phi_e_T11(t)
            [1/2   0   0]
            [  0   0   0]
            [  0   0   0]
            sage: ~phi_e_T11
            Generic morphism:
            From: Full MatrixSpace of 3 by 3 dense matrices over Rational Field
            To:   Free module of type-(1,1) tensors on the
                  3-dimensional vector space over the Rational Field
            sage: (~phi_e_T11)(W([[0, 1/2, 1/3],
            ....:                 [-1/2, 0, 0],
            ....:                 [-1/3, 0, 0]])).display()
            1/2 e_1⊗e^2 + 1/3 e_1⊗e^3 - 1/2 e_2⊗e^1 - 1/3 e_3⊗e^1

        Sending symmetric bilinear forms to matrices (note that they are currently elements
        of `T^{(0,2)}(M)`, not the symmetric power of `M`)::

            sage: T02 = V.tensor_module(0, 2); T02
            Free module of type-(0,2) tensors on the
             3-dimensional vector space over the Rational Field
            sage: e_T02 = T02.basis("e"); e_T02
            Standard basis on the
             Free module of type-(0,2) tensors on the
              3-dimensional vector space over the Rational Field
             induced by Basis (e_1,e_2,e_3) on the
              3-dimensional vector space over the Rational Field
            sage: W = MatrixSpace(QQ, 3)
            sage: phi_e_T02 = T02.isomorphism_with_fixed_basis(e_T02, codomain=W); phi_e_T02
            Generic morphism:
            From: Free module of type-(0,2) tensors on the
                  3-dimensional vector space over the Rational Field
            To:   Full MatrixSpace of 3 by 3 dense matrices over Rational Field

            sage: a = V.sym_bilinear_form()
            sage: a[1,1], a[1,2], a[1,3] = 1, 2, 3
            sage: a[2,2], a[2,3] = 4, 5
            sage: a[3,3] = 6
            sage: a.display()
            e^1⊗e^1 + 2 e^1⊗e^2 + 3 e^1⊗e^3 + 2 e^2⊗e^1 + 4 e^2⊗e^2 + 5 e^2⊗e^3
             + 3 e^3⊗e^1 + 5 e^3⊗e^2 + 6 e^3⊗e^3
            sage: phi_e_T02(a)
            [1 2 3]
            [2 4 5]
            [3 5 6]

        Same but explicitly in the subspace of symmetric bilinear forms::

            sage: Sym2Vdual = V.dual_symmetric_power(2); Sym2Vdual
            Free module of fully symmetric type-(0,2) tensors on the
             3-dimensional vector space over the Rational Field
            sage: Sym2Vdual.is_submodule(T02)
            True
            sage: Sym2Vdual.rank()
            6
            sage: e_Sym2Vdual = Sym2Vdual.basis("e"); e_Sym2Vdual
            Standard basis on the
             Free module of fully symmetric type-(0,2) tensors on the
              3-dimensional vector space over the Rational Field
             induced by Basis (e_1,e_2,e_3) on the
              3-dimensional vector space over the Rational Field
            sage: W_basis = [phi_e_T02(b) for b in e_Sym2Vdual]; W_basis
            [
            [1 0 0]  [0 1 0]  [0 0 1]  [0 0 0]  [0 0 0]  [0 0 0]
            [0 0 0]  [1 0 0]  [0 0 0]  [0 1 0]  [0 0 1]  [0 0 0]
            [0 0 0], [0 0 0], [1 0 0], [0 0 0], [0 1 0], [0 0 1]
            ]
            sage: W = MatrixSpace(QQ, 3).submodule(W_basis); W
            Free module generated by {0, 1, 2, 3, 4, 5} over Rational Field
            sage: phi_e_Sym2Vdual = Sym2Vdual.isomorphism_with_fixed_basis(e_Sym2Vdual,
            ....:                                                          codomain=W)
            sage: phi_e_Sym2Vdual
            Generic morphism:
              From: Free module of fully symmetric type-(0,2) tensors on the
                    3-dimensional vector space over the Rational Field
              To:   Free module generated by {0, 1, 2, 3, 4, 5} over Rational Field

        Sending tensors to elements of the tensor square of :class:`CombinatorialFreeModule`::

            sage: T20 = V.tensor_module(2, 0); T20
            Free module of type-(2,0) tensors on the
             3-dimensional vector space over the Rational Field
            sage: e_T20 = T02.basis("e"); e_T20
            Standard basis on the
             Free module of type-(0,2) tensors on the
              3-dimensional vector space over the Rational Field
             induced by Basis (e_1,e_2,e_3) on the
              3-dimensional vector space over the Rational Field
            sage: W = CombinatorialFreeModule(QQ, [1, 2, 3]).tensor_square(); W
            Free module generated by {1, 2, 3} over Rational Field
             # Free module generated by {1, 2, 3} over Rational Field
            sage: phi_e_T20 = T20.isomorphism_with_fixed_basis(e_T20, codomain=W); phi_e_T20
            Generic morphism:
            From: Free module of type-(2,0) tensors on the
                  3-dimensional vector space over the Rational Field
            To:   Free module generated by {1, 2, 3} over Rational Field
                  # Free module generated by {1, 2, 3} over Rational Field
            sage: t = T20.an_element(); t.display()
            1/2 e_1⊗e_1
            sage: phi_e_T20(t)
            1/2*B[1] # B[1]

        TESTS::

            sage: V = FiniteRankFreeModule(QQ, 3); V
            3-dimensional vector space over the Rational Field
            sage: e = V.basis("e")
            sage: V.isomorphism_with_fixed_basis(e, codomain=QQ^42)
            Traceback (most recent call last):
            ...
            ValueError: domain and codomain must have the same rank
            sage: V.isomorphism_with_fixed_basis(e, codomain=RR^3)
            Traceback (most recent call last):
            ...
            ValueError: domain and codomain must have the same base ring
        """
        base_ring = self.base_ring()
        if basis is None:
            basis = self.default_basis()
        if codomain is None:
            from sage.combinat.free_module import CombinatorialFreeModule
            if isinstance(basis._symbol, str):
                prefix = basis._symbol
            else:
                prefix = None
            codomain = CombinatorialFreeModule(base_ring, basis.keys(),
                                               prefix=prefix)
        else:
            try:
                codomain_rank = codomain.rank()
            except AttributeError:
                # https://github.com/sagemath/sage/issues/34445: MatrixSpace does not have rank
                codomain_rank = codomain.dimension()
            if codomain_rank != self.rank():
                raise ValueError("domain and codomain must have the same rank")
            if codomain.base_ring() != base_ring:
                raise ValueError("domain and codomain must have the same "
                                 "base ring")

        codomain_basis = Family(codomain.basis())
        if isinstance(codomain_basis, TrivialFamily):
            # assume that codomain basis keys are to be ignored;
            # need them several times, can't keep as generators
            key_pairs = tuple(enumerate(basis.keys()))
            basis_by_codomain_key = basis
        else:
            # assume that the keys of the codomain should be used
            # need them several times, can't keep as generators
            key_pairs = tuple(zip(codomain_basis.keys(), basis.keys()))
            basis_by_codomain_key = {}
            for codomain_key, domain_key in key_pairs:
                basis_by_codomain_key[codomain_key] = basis[domain_key]

        def _isomorphism(x):
            r"""
            Concrete isomorphism from ``self`` to ``codomain``.
            """
            return codomain.sum(x[basis, domain_key] * codomain_basis[codomain_key]
                                for codomain_key, domain_key in key_pairs)

        def _inverse(y):
            r"""
            Concrete isomorphism from ``codomain`` to ``self``.
            """
            return self.linear_combination(
                (basis_by_codomain_key[codomain_key], coefficient)
                for codomain_key, coefficient in y.monomial_coefficients().items())

        category = Modules(self.base_ring())
        homset = Hom(self, codomain, category)
        isomorphism = SetIsomorphism(homset, _isomorphism)
        inverse = SetIsomorphism(homset.reversed(), _inverse)
        isomorphism._set_inverse(inverse)
        inverse._set_inverse(isomorphism)
        return isomorphism

    def _test_isomorphism_with_fixed_basis(self, **options):
        r"""
        Test that the method ``isomorphism_with_fixed_basis`` works correctly.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M._test_isomorphism_with_fixed_basis()
        """
        tester = self._tester(**options)
        try:
            basis = self.basis('test')
        except AttributeError:
            return
        morphism = self.isomorphism_with_fixed_basis(basis)
        tester.assertEqual(morphism.codomain().rank(), self.rank())


class FiniteRankFreeModule(ReflexiveModule_base, FiniteRankFreeModule_abstract):
    r"""
    Free module of finite rank over a commutative ring.

    A *free module of finite rank* over a commutative ring `R` is a module `M`
    over `R` that admits a *finite basis*, i.e. a finite family of linearly
    independent generators. Since `R` is commutative, it has the invariant
    basis number property, so that the rank of the free module `M` is defined
    uniquely, as the cardinality of any basis of `M`.

    No distinguished basis of `M` is assumed. On the contrary, many bases can be
    introduced on the free module along with change-of-basis rules (as module
    automorphisms). Each
    module element has then various representations over the various bases.

    .. NOTE::

        The class :class:`FiniteRankFreeModule` does not inherit from
        class :class:`~sage.modules.free_module.FreeModule_generic`
        nor from class
        :class:`~sage.combinat.free_module.CombinatorialFreeModule`, since
        both classes deal with modules with a *distinguished basis* (see
        details :ref:`above <diff-FreeModule>`).
        Moreover, following the recommendation exposed in :issue:`16427`
        the class :class:`FiniteRankFreeModule` inherits directly from
        :class:`~sage.structure.parent.Parent` (with the category set to
        :class:`~sage.categories.modules.Modules`) and not from the Cython
        class :class:`~sage.modules.module.Module`.

    The class :class:`FiniteRankFreeModule` is a Sage *parent* class,
    the corresponding *element* class being
    :class:`~sage.tensor.modules.free_module_element.FiniteRankFreeModuleElement`.

    INPUT:

    - ``ring`` -- commutative ring `R` over which the free module is
      constructed
    - ``rank`` -- positive integer; rank of the free module
    - ``name`` -- (default: ``None``) string; name given to the free module
    - ``latex_name`` -- (default: ``None``) string;  LaTeX symbol to denote
      the freemodule; if none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower bound of the range of
      indices in bases defined on the free module
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the tensor components;
      ``output_formatter`` must take 1 or 2 arguments: the first argument
      must be an element of the ring `R` and  the second one, if any, some
      format specification

    EXAMPLES:

    Free module of rank 3 over `\ZZ`::

        sage: FiniteRankFreeModule._clear_cache_() # for doctests only
        sage: M = FiniteRankFreeModule(ZZ, 3) ; M
        Rank-3 free module over the Integer Ring
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M') ; M  # declaration with a name
        Rank-3 free module M over the Integer Ring
        sage: M.category()
        Category of finite dimensional modules over Integer Ring
        sage: M.base_ring()
        Integer Ring
        sage: M.rank()
        3

    If the base ring is a field, the free module is in the category of vector
    spaces::

        sage: V = FiniteRankFreeModule(QQ, 3, name='V') ; V
        3-dimensional vector space V over the Rational Field
        sage: V.category()
        Category of finite dimensional vector spaces over Rational Field

    The LaTeX output is adjusted via the parameter ``latex_name``::

        sage: latex(M)  # the default is the symbol provided in the string ``name``
        M
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', latex_name=r'\mathcal{M}')
        sage: latex(M)
        \mathcal{M}

    The free module M has no distinguished basis::

        sage: M in ModulesWithBasis(ZZ)
        False
        sage: M in Modules(ZZ)
        True

    In particular, no basis is initialized at the module construction::

        sage: M.print_bases()
        No basis has been defined on the Rank-3 free module M over the Integer Ring
        sage: M.bases()
        []

    Bases have to be introduced by means of the method :meth:`basis`,
    the first defined basis being considered as the *default basis*, meaning
    it can be skipped in function arguments required a basis (this can
    be changed by means of the method :meth:`set_default_basis`)::

        sage: e = M.basis('e') ; e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: M.default_basis()
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring

    A second basis can be created from a family of linearly independent
    elements expressed in terms of basis ``e``::

        sage: f = M.basis('f', from_family=(-e[0], e[1]+e[2], 2*e[1]+3*e[2]))
        sage: f
        Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring
        sage: M.print_bases()
        Bases defined on the Rank-3 free module M over the Integer Ring:
         - (e_0,e_1,e_2) (default basis)
         - (f_0,f_1,f_2)
        sage: M.bases()
        [Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring,
         Basis (f_0,f_1,f_2) on the Rank-3 free module M over the Integer Ring]

    M is a *parent* object, whose elements are instances of
    :class:`~sage.tensor.modules.free_module_element.FiniteRankFreeModuleElement`
    (actually a dynamically generated subclass of it)::

        sage: v = M.an_element() ; v
        Element of the Rank-3 free module M over the Integer Ring
        sage: from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement
        sage: isinstance(v, FiniteRankFreeModuleElement)
        True
        sage: v in M
        True
        sage: M.is_parent_of(v)
        True
        sage: v.display() # expansion w.r.t. the default basis (e)
        e_0 + e_1 + e_2
        sage: v.display(f)
        -f_0 + f_1

    The test suite of the category of modules is passed::

        sage: TestSuite(M).run()

    Constructing an element of ``M`` from (the integer) 0 yields
    the zero element of ``M``::

        sage: M(0)
        Element zero of the Rank-3 free module M over the Integer Ring
        sage: M(0) is M.zero()
        True

    Non-zero elements are constructed by providing their components in
    a given basis::

        sage: v = M([-1,0,3]) ; v  # components in the default basis (e)
        Element of the Rank-3 free module M over the Integer Ring
        sage: v.display() # expansion w.r.t. the default basis (e)
        -e_0 + 3 e_2
        sage: v.display(f)
        f_0 - 6 f_1 + 3 f_2
        sage: v = M([-1,0,3], basis=f) ; v  # components in a specific basis
        Element of the Rank-3 free module M over the Integer Ring
        sage: v.display(f)
        -f_0 + 3 f_2
        sage: v.display()
        e_0 + 6 e_1 + 9 e_2
        sage: v = M([-1,0,3], basis=f, name='v') ; v
        Element v of the Rank-3 free module M over the Integer Ring
        sage: v.display(f)
        v = -f_0 + 3 f_2
        sage: v.display()
        v = e_0 + 6 e_1 + 9 e_2

    An alternative is to construct the element from an empty list of
    componentsand to set the nonzero components afterwards::

        sage: v = M([], name='v')
        sage: v[e,0] = -1
        sage: v[e,2] = 3
        sage: v.display(e)
        v = -e_0 + 3 e_2

    Indices on the free module, such as indices labelling the element of a
    basis, are provided by the generator method :meth:`irange`. By default,
    they range from 0 to the module's rank minus one::

        sage: list(M.irange())
        [0, 1, 2]

    This can be changed via the parameter ``start_index`` in the module
    construction::

        sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: list(M1.irange())
        [1, 2, 3]

    The parameter ``output_formatter`` in the constructor of the free module
    is used to set the output format of tensor components::

        sage: N = FiniteRankFreeModule(QQ, 3, output_formatter=Rational.numerical_approx)
        sage: e = N.basis('e')
        sage: v = N([1/3, 0, -2], basis=e)
        sage: v[e,:]
        [0.333333333333333, 0.000000000000000, -2.00000000000000]
        sage: v.display(e)  # default format (53 bits of precision)
        0.333333333333333 e_0 - 2.00000000000000 e_2
        sage: v.display(e, format_spec=10)  # 10 bits of precision
        0.33 e_0 - 2.0 e_2
    """

    Element = FiniteRankFreeModuleElement
    _sindex: int

    @staticmethod
    def __classcall_private__(cls, ring, rank, name=None, latex_name=None, start_index=0,
                              output_formatter=None, category=None, ambient=None):
        r"""
        Normalize init arguments for ``UniqueRepresentation``.

        TESTS::

            sage: FiniteRankFreeModule(QQ, 3) is FiniteRankFreeModule(QQ, 3, name=None)
            True
            sage: FiniteRankFreeModule(QQ, 3, name='M') is FiniteRankFreeModule(QQ, 3, name='M', latex_name='M')
            True
            sage: FiniteRankFreeModule(QQ, 3) is FiniteRankFreeModule(QQ, 3, start_index=0)
            True
            sage: FiniteRankFreeModule(QQ, 3) is FiniteRankFreeModule(QQ, 3, output_formatter=None)
            True
            sage: FiniteRankFreeModule(QQ, 3) is FiniteRankFreeModule(QQ, 3, category=Modules(QQ).FiniteDimensional())
            True
        """
        if ring not in Rings().Commutative():
            raise TypeError("the module base ring must be commutative")
        category = Modules(ring).FiniteDimensional().or_subcategory(category)
        if latex_name is None:
            latex_name = name
        return super().__classcall__(
            cls, ring, rank, name, latex_name, start_index, output_formatter, category, ambient)

    def __init__(
        self,
        ring,
        rank,
        name=None,
        latex_name=None,
        start_index: int = 0,
        output_formatter=None,
        category=None,
        ambient=None,
    ):
        r"""
        See :class:`FiniteRankFreeModule` for documentation and examples.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: TestSuite(M).run()
            sage: e = M.basis('e')
            sage: TestSuite(M).run()
            sage: f = M.basis('f')
            sage: TestSuite(M).run()
        """
        super().__init__(ring, rank, name=name, latex_name=latex_name,
                         category=category, ambient=ambient)
        self._sindex = start_index
        self._output_formatter = output_formatter
        # Dictionary of the tensor modules built on self
        #   (keys = (k,l) --the tensor type)
        # This dictionary is to be extended on need by the method tensor_module
        self._tensor_modules = {}
        # Dictionaries of exterior powers of self and of its dual
        #   (keys = p --the power degree)
        # These dictionaries are to be extended on need by the methods
        # exterior_power and dual_exterior_power
        self._exterior_powers = {}
        self._dual_exterior_powers = {}
        # Set of all modules (tensor powers, exterior powers)
        # that depend on self's bases:
        self._all_modules = set([self])
        # List of known bases on the free module:
        self._known_bases = []
        self._def_basis = None # default basis
        self._basis_changes = {} # Dictionary of the changes of bases
        # Identity automorphism:
        self._identity_map = None # to be set by self.identity_map()
        # General linear group:
        self._general_linear_group = None # to be set by
                                          # self.general_linear_group()

    def construction(self):
        """
        The construction functor and base ring for ``self``.

        EXAMPLES::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.construction()
            (VectorFunctor, Integer Ring)
            sage: N = FiniteRankFreeModule(ZZ, 3, name='N', start_index=17)
            sage: N.construction()
            (VectorFunctor, Integer Ring)
        """
        # Try to take it from the category
        c = super().construction()
        if c is not None:
            return c
        # Implementation restrictions:
        if self._output_formatter:
            return None
        from sage.categories.pushout import VectorFunctor
        kwds = dict(is_sparse=False,
                    inner_product_matrix=None,
                    with_basis=None,
                    name_mapping={self.base_ring(): self._name} if self._name else None,
                    latex_name_mapping={self.base_ring(): self._latex_name} if self._latex_name else None)
        if self._sindex:
            return (VectorFunctor(basis_keys=list(self.irange()), **kwds),
                    self.base_ring())
        return (VectorFunctor(n=self.rank(), **kwds),
                self.base_ring())

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: v = M._element_constructor_(comp=[1,0,-2], basis=e, name='v') ; v
            Element v of the Rank-3 free module M over the Integer Ring
            sage: v.display()
            v = e_0 - 2 e_2
            sage: v == M([1,0,-2])
            True
            sage: v = M._element_constructor_(0) ; v
            Element zero of the Rank-3 free module M over the Integer Ring
            sage: v = M._element_constructor_() ; v
            Element of the Rank-3 free module M over the Integer Ring
        """
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) element of ``self``.

        EXAMPLES::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: v = M._an_element_(); v
            Element of the Rank-3 free module M over the Integer Ring
            sage: v.display()
            e_0 + e_1 + e_2
            sage: v == M.an_element()
            True
            sage: v.parent()
            Rank-3 free module M over the Integer Ring
        """
        if self._def_basis is None:
            self.basis('e')
        resu = self.element_class(self)
        resu.set_comp()[:] = [self._ring.an_element() for i in range(self._rank)]
        return resu

    #### End of parent methods

    #### Methods to be redefined by derived classes ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FiniteRankFreeModule(ZZ, 3, name='M')
            Rank-3 free module M over the Integer Ring
        """
        if self._ring in Fields():
            description = "{}-dimensional vector space ".format(self._rank)
        else:
            description = "Rank-{} free module ".format(self._rank)
        if self._name is not None:
            description += self._name + " "
        description += "over the {}".format(self._ring)
        return description

    def _Hom_(self, other, category=None):
        r"""
        Construct the set of homomorphisms ``self`` --> ``other``.

        INPUT:

        - ``other`` -- another free module of finite rank over the same ring
          as ``self``
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the hom-set Hom(M,N), where M is ``self`` and N is ``other``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: H = M._Hom_(N) ; H
            Set of Morphisms from Rank-3 free module M over the Integer Ring
             to Rank-2 free module N over the Integer Ring
             in Category of finite dimensional modules over Integer Ring
            sage: H = Hom(M,N) ; H  # indirect doctest
            Set of Morphisms from Rank-3 free module M over the Integer Ring
             to Rank-2 free module N over the Integer Ring
             in Category of finite dimensional modules over Integer Ring
        """
        from .free_module_homset import FreeModuleHomset
        return FreeModuleHomset(self, other)

    def tensor_module(self, k, l, *, sym=None, antisym=None):
        r"""
        Return the free module of all tensors of type `(k, l)` defined on
        ``self``.

        INPUT:

        - ``k`` -- nonnegative integer; the contravariant rank, the tensor
          type being `(k, l)`
        - ``l`` -- nonnegative integer; the covariant rank, the tensor type
          being `(k, l)`
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention ``position = 0`` for the first argument. For instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or list of
          antisymmetries among the arguments, with the same convention
          as for ``sym``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule`
          representing the free module
          `T^{(k,l)}(M)` of type-`(k,l)` tensors on the free module ``self``

        EXAMPLES:

        Tensor modules over a free module over `\ZZ`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T = M.tensor_module(1,2) ; T
            Free module of type-(1,2) tensors on the Rank-3 free module M
             over the Integer Ring
            sage: T.an_element()
            Type-(1,2) tensor on the Rank-3 free module M over the Integer Ring

        Tensor modules are unique::

            sage: M.tensor_module(1,2) is T
            True

        The module of type-`(1,0)` tensors is the base module itself::

            sage: M.tensor_module(1,0) is M
            True

        while the module of type-`(0,1)` tensors is the dual of the base module::

            sage: M.tensor_module(0, 1) is M.dual()
            True

        By using the arguments ``sym`` and ``antisym``, submodules of a full tensor
        module can be constructed::

            sage: T = M.tensor_module(4, 4, sym=((0, 1)), antisym=((4, 5))); T
            Free module of type-(4,4) tensors on the Rank-3 free module M over the Integer Ring,
             with symmetry on the index positions (0, 1),
             with antisymmetry on the index positions (4, 5)
            sage: T._name
            'T^{2,3}(M)⊗T^{6,7}(M*)⊗Sym^{0,1}(M)⊗ASym^{4,5}(M*)'
            sage: latex(T)
            T^{\{2,3\}}(M) \otimes T^{\{6,7\}}(M^*) \otimes \mathrm{Sym}^{\{0,1\}}(M) \otimes \mathrm{ASym}^{\{4,5\}}(M^*)

        See :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule`
        and :class:`~sage.tensor.modules.tensor_free_module.TensorFreeSubmodule_sym`
        for more documentation.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 2)
            sage: M.tensor_module(2, 0, sym=(0,1)) is M.symmetric_power(2)
            True
        """
        from .comp import CompWithSym

        sym, antisym = CompWithSym._canonicalize_sym_antisym(k + l, sym, antisym)
        if sym or antisym:
            key = (k, l, sym, antisym)
        else:
            key = (k, l)
        try:
            return self._tensor_modules[key]
        except KeyError:
            if key == (1, 0):
                T = self
            elif key == (0, 1):
                T = self.dual()
            elif sym or antisym:
                from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_sym
                T = TensorFreeSubmodule_sym(self, (k, l), sym=sym, antisym=antisym)
            else:
                from sage.tensor.modules.tensor_free_module import TensorFreeModule
                T = TensorFreeModule(self, (k, l))
            self._tensor_modules[key] = T
            return T

    def symmetric_power(self, p):
        r"""
        Return the `p`-th symmetric power of ``self``.

        EXAMPLES:

        Symmetric powers of a free `\ZZ`-module of rank 3::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.symmetric_power(0)
            Free module of type-(0,0) tensors on the Rank-3 free module M over the Integer Ring
            sage: M.symmetric_power(1)  # return the module itself
            Rank-3 free module M over the Integer Ring
            sage: M.symmetric_power(1) is M
            True
            sage: M.symmetric_power(2)
            Free module of fully symmetric type-(2,0) tensors
             on the Rank-3 free module M over the Integer Ring
            sage: M.symmetric_power(2).an_element()
            Type-(2,0) tensor on the Rank-3 free module M over the Integer Ring
            sage: M.symmetric_power(2).an_element().display()
            e_0⊗e_0
            sage: M.symmetric_power(3)
            Free module of fully symmetric type-(3,0) tensors
             on the Rank-3 free module M over the Integer Ring
            sage: M.symmetric_power(3).an_element()
            Type-(3,0) tensor on the Rank-3 free module M over the Integer Ring
            sage: M.symmetric_power(3).an_element().display()
            e_0⊗e_0⊗e_0
        """
        if p <= 1:
            return self.tensor_module(p, 0)
        return self.tensor_module(p, 0, sym=(tuple(range(p)),))

    def dual_symmetric_power(self, p):
        r"""
        Return the `p`-th symmetric power of the dual of ``self``.

        EXAMPLES:

        Symmetric powers of the dual of a free `\ZZ`-module of rank 3::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.dual_symmetric_power(0)
            Free module of type-(0,0) tensors on the Rank-3 free module M over the Integer Ring
            sage: M.dual_symmetric_power(1)  # return the dual module
            Dual of the Rank-3 free module M over the Integer Ring
            sage: M.dual_symmetric_power(2)
            Free module of fully symmetric type-(0,2) tensors
             on the Rank-3 free module M over the Integer Ring
            sage: M.dual_symmetric_power(2).an_element()
            Symmetric bilinear form  on the Rank-3 free module M over the Integer Ring
            sage: M.dual_symmetric_power(2).an_element().display()
            e^0⊗e^0
            sage: M.dual_symmetric_power(3)
            Free module of fully symmetric type-(0,3) tensors
             on the Rank-3 free module M over the Integer Ring
            sage: M.dual_symmetric_power(3).an_element()
            Type-(0,3) tensor on the Rank-3 free module M over the Integer Ring
            sage: M.dual_symmetric_power(3).an_element().display()
            e^0⊗e^0⊗e^0
        """
        if p <= 1:
            return self.tensor_module(0, p)
        return self.tensor_module(0, p, sym=(tuple(range(p)),))

    def exterior_power(self, p):
        r"""
        Return the `p`-th exterior power of ``self``.

        If `M` stands for the free module ``self``, the *p-th exterior
        power of* `M` is the set `\Lambda^p(M)` of all *alternating
        contravariant tensors* of rank `p`, i.e. of all multilinear maps

        .. MATH::

            \underbrace{M^*\times\cdots\times M^*}_{p\ \; \mbox{times}}
            \longrightarrow R

        that vanish whenever any of two of their arguments are equal.
        `\Lambda^p(M)` is a free module of rank `\binom{n}{p}`
        over the same ring as `M`, where `n` is the rank of `M`.

        INPUT:

        - ``p`` -- nonnegative integer

        OUTPUT:

        - for `p=0`, the base ring `R`
        - for `p=1`, the free module `M`, since `\Lambda^1(M)=M`
        - for `p\geq 2`, instance of
          :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`
          representing the free module `\Lambda^p(M)`

        EXAMPLES:

        Exterior powers of a free `\ZZ`-module of rank 3::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.exterior_power(0)  # return the base ring
            Integer Ring
            sage: M.exterior_power(1)  # return the module itself
            Rank-3 free module M over the Integer Ring
            sage: M.exterior_power(1) is M
            True
            sage: M.exterior_power(2)
            2nd exterior power of the Rank-3 free module M over the Integer Ring
            sage: M.exterior_power(2).an_element()
            Alternating contravariant tensor of degree 2 on the Rank-3
             free module M over the Integer Ring
            sage: M.exterior_power(2).an_element().display()
            e_0∧e_1
            sage: M.exterior_power(3)
            3rd exterior power of the Rank-3 free module M over the Integer Ring
            sage: M.exterior_power(3).an_element()
            Alternating contravariant tensor of degree 3 on the Rank-3
             free module M over the Integer Ring
            sage: M.exterior_power(3).an_element().display()
            e_0∧e_1∧e_2

        See
        :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`
        for more documentation.
        """
        try:
            return self._exterior_powers[p]
        except KeyError:
            if p == 0:
                L = self._ring
            elif p == 1:
                L = self
            else:
                from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
                L = ExtPowerFreeModule(self, p)
            self._exterior_powers[p] = L
            return L

    def dual_exterior_power(self, p):
        r"""
        Return the `p`-th exterior power of the dual of ``self``.

        If `M` stands for the free module ``self``, the *p-th exterior
        power of the dual of* `M` is the set `\Lambda^p(M^*)` of all
        *alternating forms of degree* `p` on `M`, i.e. of all
        multilinear maps

        .. MATH::

            \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
            \longrightarrow R

        that vanish whenever any of two of their arguments are equal.
        `\Lambda^p(M^*)` is a free module of rank `\binom{n}{p}`
        over the same ring as `M`, where `n` is the rank of `M`.

        INPUT:

        - ``p`` -- nonnegative integer

        OUTPUT:

        - for `p=0`, the base ring `R`
        - for `p=1`, instance of
          :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankDualFreeModule`
          representing the dual `M^*`
        - for `p\geq 1`, instance of
          :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerDualFreeModule`
          representing the free module `\Lambda^p(M^*)`

        EXAMPLES:

        Exterior powers of the dual of a free `\ZZ`-module of rank 3::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.dual_exterior_power(0)  # return the base ring
            Integer Ring
            sage: M.dual_exterior_power(1)  # return the dual module
            Dual of the Rank-3 free module M over the Integer Ring
            sage: M.dual_exterior_power(1) is M.dual()
            True
            sage: M.dual_exterior_power(2)
            2nd exterior power of the dual of the Rank-3 free module M over the Integer Ring
            sage: M.dual_exterior_power(2).an_element()
            Alternating form of degree 2 on the Rank-3 free module M over the Integer Ring
            sage: M.dual_exterior_power(2).an_element().display()
            e^0∧e^1
            sage: M.dual_exterior_power(3)
            3rd exterior power of the dual of the Rank-3 free module M over the Integer Ring
            sage: M.dual_exterior_power(3).an_element()
            Alternating form of degree 3 on the Rank-3 free module M over the Integer Ring
            sage: M.dual_exterior_power(3).an_element().display()
            e^0∧e^1∧e^2

        See
        :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerDualFreeModule`
        for more documentation.
        """
        try:
            return self._dual_exterior_powers[p]
        except KeyError:
            if p == 0:
                L = self._ring
            elif p == 1:
                L = FiniteRankDualFreeModule(self)
            else:
                from sage.tensor.modules.ext_pow_free_module import ExtPowerDualFreeModule
                L = ExtPowerDualFreeModule(self, p)
            self._dual_exterior_powers[p] = L
            return L

    def general_linear_group(self):
        r"""
        Return the general linear group of ``self``.

        If ``self`` is the free module `M`, the *general linear group* is the
        group `\mathrm{GL}(M)` of automorphisms of `M`.

        OUTPUT:

        - instance of class
          :class:`~sage.tensor.modules.free_module_linear_group.FreeModuleLinearGroup`
          representing `\mathrm{GL}(M)`

        EXAMPLES:

        The general linear group of a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: GL = M.general_linear_group() ; GL
            General linear group of the Rank-3 free module M over the Integer Ring
            sage: GL.category()
            Category of groups
            sage: type(GL)
            <class 'sage.tensor.modules.free_module_linear_group.FreeModuleLinearGroup_with_category'>

        There is a unique instance of the general linear group::

            sage: M.general_linear_group() is GL
            True

        The group identity element::

            sage: GL.one()
            Identity map of the Rank-3 free module M over the Integer Ring
            sage: GL.one().matrix(e)
            [1 0 0]
            [0 1 0]
            [0 0 1]

        An element::

            sage: GL.an_element()
            Automorphism of the Rank-3 free module M over the Integer Ring
            sage: GL.an_element().matrix(e)
            [ 1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

        See
        :class:`~sage.tensor.modules.free_module_linear_group.FreeModuleLinearGroup`
        for more documentation.
        """
        from sage.tensor.modules.free_module_linear_group import \
                                                          FreeModuleLinearGroup
        if self._general_linear_group is None:
            self._general_linear_group = FreeModuleLinearGroup(self)
        return self._general_linear_group

    def basis(self, symbol, latex_symbol=None, from_family=None,
              indices=None, latex_indices=None, symbol_dual=None,
              latex_symbol_dual=None):
        r"""
        Define or return a basis of the free module ``self``.

        Let `M` denotes the free module ``self`` and `n` its rank.

        The basis can be defined from a set of `n` linearly independent
        elements of `M` by means of the argument ``from_family``.
        If ``from_family`` is not specified, the basis is created from
        scratch and, at this stage, is unrelated to bases that could have been
        defined previously on `M`. It can be related afterwards by means of
        the method :meth:`set_change_of_basis`.

        If the basis specified by the given symbol already exists, it is
        simply returned, whatever the value of the arguments ``latex_symbol``
        or ``from_family``.

        Note that another way to construct a basis of ``self`` is to use
        the method
        :meth:`~sage.tensor.modules.free_module_basis.FreeModuleBasis.new_basis`
        on an existing basis, with the automorphism relating the two bases as
        an argument.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the elements of the basis, or a list/tuple of strings,
          representing the individual symbols of the elements of the basis
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the elements of the basis,
          or a list/tuple of strings, representing the individual LaTeX symbols
          of the elements of the basis; if ``None``, ``symbol`` is used in
          place of ``latex_symbol``
        - ``from_family`` -- (default: ``None``) tuple or list of `n` linearly
          independent elements of the free module ``self`` (`n` being the
          rank of ``self``)
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) list/tuple of strings representing the indices
          labelling the elements of the basis; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) list/tuple of strings
          representing the indices for the LaTeX symbols of the elements of
          the basis; if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual basis; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual basis
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual basis

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing a basis on ``self``

        EXAMPLES:

        Bases on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e') ; e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: e[0]
            Element e_0 of the Rank-3 free module M over the Integer Ring
            sage: latex(e)
            \left(e_{0},e_{1},e_{2}\right)

        The LaTeX symbol can be set explicitly::

            sage: eps = M.basis('eps', latex_symbol=r'\epsilon') ; eps
            Basis (eps_0,eps_1,eps_2) on the Rank-3 free module M
             over the Integer Ring
            sage: latex(eps)
            \left(\epsilon_{0},\epsilon_{1},\epsilon_{2}\right)

        The indices can be customized::

            sage: f = M.basis('f', indices=('x', 'y', 'z')); f
            Basis (f_x,f_y,f_z) on the Rank-3 free module M over the Integer Ring
            sage: latex(f[1])
            f_{y}

        By providing a list or a tuple for the argument ``symbol``, one can
        have a different symbol for each element of the basis; it is then
        mandatory to specify some symbols for the dual basis::

            sage: g = M.basis(('a', 'b', 'c'), symbol_dual=('A', 'B', 'C')); g
            Basis (a,b,c) on the Rank-3 free module M over the Integer Ring
            sage: g.dual_basis()
            Dual basis (A,B,C) on the Rank-3 free module M over the Integer Ring

        If the provided symbol and indices are that of an already defined
        basis, the latter is returned (no new basis is created)::

            sage: M.basis('e') is e
            True
            sage: M.basis('eps') is eps
            True
            sage: M.basis('e', indices=['x', 'y', 'z']) is e
            False
            sage: M.basis('e', indices=['x', 'y', 'z']) is \
            ....:  M.basis('e', indices=['x', 'y', 'z'])
            True

        The individual elements of the basis are labelled according the
        parameter ``start_index`` provided at the free module construction::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e') ; e
            Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
            sage: e[1]
            Element e_1 of the Rank-3 free module M over the Integer Ring

        Construction of a basis from a spanning family of linearly independent
        module elements::

            sage: f1 = -e[2]
            sage: f2 = 4*e[1] + 3*e[3]
            sage: f3 = 7*e[1] + 5*e[3]
            sage: f = M.basis('f', from_family=(f1,f2,f3))
            sage: f[1].display()
            f_1 = -e_2
            sage: f[2].display()
            f_2 = 4 e_1 + 3 e_3
            sage: f[3].display()
            f_3 = 7 e_1 + 5 e_3

        The change-of-basis automorphisms have been registered::

            sage: M.change_of_basis(e,f).matrix(e)
            [ 0  4  7]
            [-1  0  0]
            [ 0  3  5]
            sage: M.change_of_basis(f,e).matrix(e)
            [ 0 -1  0]
            [-5  0  7]
            [ 3  0 -4]
            sage: M.change_of_basis(f,e) == M.change_of_basis(e,f).inverse()
            True

        Check of the change-of-basis e --> f::

            sage: a = M.change_of_basis(e,f) ; a
            Automorphism of the Rank-3 free module M over the Integer Ring
            sage: all( f[i] == a(e[i]) for i in M.irange() )
            True

        Providing a family of module elements that are not linearly independent
        raise an error::

            sage: g = M.basis('g', from_family=(f1, f2, f1+f2))
            Traceback (most recent call last):
            ...
            ValueError: the provided module elements are not linearly
             independent

        For more documentation on bases see
        :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`.
        """
        from .free_module_basis import FreeModuleBasis
        for other in self._known_bases:
            if symbol == other._symbol and indices == other._indices:
                return other
        resu = FreeModuleBasis(self, symbol, latex_symbol=latex_symbol,
                               indices=indices, latex_indices=latex_indices,
                               symbol_dual=symbol_dual,
                               latex_symbol_dual=latex_symbol_dual)
        if from_family:
            try:
                resu._init_from_family(from_family)
            except ZeroDivisionError:
                raise ValueError("the provided module elements are not "
                                 "linearly independent")
        return resu

    def _test_basis(self, tester=None, **options):
        r"""
        Test that the ``basis`` method works correctly.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M._test_basis(verbose=True)
            <BLANKLINE>
              Running the test suite of self.basis('test')
              running ._test_an_element() . . . pass
              running ._test_cardinality() . . . pass
              running ._test_category() . . . pass
              running ._test_construction() . . . pass
              running ._test_elements() . . .
              Running the test suite of self.an_element()
                running ._test_category() . . . pass
                running ._test_eq() . . . pass
                running ._test_new() . . . pass
                running ._test_nonzero_equal() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                pass
              running ._test_elements_eq_reflexive() . . . pass
              running ._test_elements_eq_symmetric() . . . pass
              running ._test_elements_eq_transitive() . . . pass
              running ._test_elements_neq() . . . pass
              running ._test_enumerated_set_contains() . . . pass
              running ._test_enumerated_set_iter_cardinality() . . . pass
              running ._test_enumerated_set_iter_list() . . . pass
              running ._test_eq() . . . pass
              running ._test_iter_len() . . . pass
              running ._test_new() . . . pass
              running ._test_not_implemented_methods() . . . pass
              running ._test_pickling() . . . pass
              running ._test_some_elements() . . . pass
        """
        from sage.misc.sage_unittest import TestSuite
        # The intention is to raise an exception only if this is
        # run as a sub-testsuite of a larger testsuite.
        # (from _test_elements)
        is_sub_testsuite = (tester is not None)
        tester = self._tester(tester=tester, **options)
        try:
            b = self.basis('test')
        except NotImplementedError:
            return
        # Test uniqueness
        b_again = self.basis('test')
        tester.assertTrue(b is b_again)
        # Test rank
        tester.assertEqual(len(b), self.rank())
        indices = list(self.irange())
        tester.assertEqual(len(b), len(indices))
        # Test basis indexing
        for index, element in zip(indices, b):
            tester.assertTrue(element is b[index])
        # Run test suite of the basis object (similar to _test_elements)
        tester.info("\n  Running the test suite of self.basis('test')")
        TestSuite(b).run(verbose=tester._verbose, prefix=tester._prefix + "  ",
                         raise_on_failure=is_sub_testsuite)

    def _tensor(self, tensor_type, name=None, latex_name=None, sym=None,
               antisym=None):
        r"""
        Construct a tensor on the free module ``self``.

        INPUT:

        - ``tensor_type`` -- pair ``(k, l)`` with ``k`` being the
          contravariant rank and ``l`` the covariant rank

        - ``name`` -- (default: ``None``) string; name given to the tensor

        - ``latex_name`` -- (default: ``None``) string;  LaTeX symbol to
          denote the tensor; if none is provided, the LaTeX symbol is set
          to ``name``

        - ``sym`` -- (default: ``None``) a symmetry or an iterable of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention ``position = 0`` for the first argument. For instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or iterable of
          antisymmetries among the arguments, with the same convention
          as for ``sym``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor defined on ``self`` with the provided
          characteristics

        EXAMPLES:

        Tensors on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M._tensor((1,0), name='t') ; t
            Element t of the Rank-3 free module M over the Integer Ring
        """
        from .comp import CompWithSym
        sym, antisym = CompWithSym._canonicalize_sym_antisym(
            tensor_type[0] + tensor_type[1], sym, antisym)
        # Special cases:
        if tensor_type == (1,0):
            return self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type == (0,1):
            return self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type[0] == 0 and tensor_type[1] > 1 and antisym:
            if len(antisym[0]) == tensor_type[1]:
                return self.alternating_form(tensor_type[1], name=name,
                                             latex_name=latex_name)
        elif tensor_type[0] > 1 and tensor_type[1] == 0 and antisym:
            if len(antisym[0]) == tensor_type[0]:
                return self.alternating_contravariant_tensor(tensor_type[0],
                                           name=name, latex_name=latex_name)
        # Generic case:
        return self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name,
                                 sym=sym, antisym=antisym)

    def tensor(self, *args, **kwds):
        r"""
        Construct a tensor on the free module ``self`` or a tensor product with other modules.

        If ``args`` consist of other parents, just delegate to :meth:`tensor_product`.

        Otherwise, construct a tensor from the following input.

        INPUT:

        - ``tensor_type`` -- pair ``(k, l)`` with ``k`` being the
          contravariant rank and ``l`` the covariant rank

        - ``name`` -- (default: ``None``) string; name given to the tensor

        - ``latex_name`` -- (default: ``None``) string;  LaTeX symbol to
          denote the tensor; if none is provided, the LaTeX symbol is set
          to ``name``

        - ``sym`` -- (default: ``None``) a symmetry or an iterable of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention ``position = 0`` for the first argument. For instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or iterable of
          antisymmetries among the arguments, with the same convention
          as for ``sym``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor defined on ``self`` with the provided
          characteristics

        EXAMPLES:

        Tensors on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((1,0), name='t') ; t
            Element t of the Rank-3 free module M over the Integer Ring
            sage: t = M.tensor((0,1), name='t') ; t
            Linear form t on the Rank-3 free module M over the Integer Ring
            sage: t = M.tensor((1,1), name='t') ; t
            Type-(1,1) tensor t on the Rank-3 free module M over the Integer Ring
            sage: t = M.tensor((0,2), name='t', sym=(0,1)) ; t
            Symmetric bilinear form t on the
             Rank-3 free module M over the Integer Ring
            sage: t = M.tensor((0,2), name='t', antisym=(0,1)) ; t
            Alternating form t of degree 2 on the
             Rank-3 free module M over the Integer Ring
            sage: t = M.tensor((1,2), name='t') ; t
            Type-(1,2) tensor t on the Rank-3 free module M over the Integer Ring

        See :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
        for more examples and documentation.

        TESTS:

        Trivial symmetries in the list of symmetries or antisymmetries are silently
        ignored::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.tensor((3,0), sym=[[1]])
            Type-(3,0) tensor on the Rank-3 free module M over the Integer Ring
            sage: M.tensor((3,0), antisym=[[]])
            Type-(3,0) tensor on the Rank-3 free module M over the Integer Ring
        """
        # Until https://github.com/sagemath/sage/issues/30373 is done,
        # TensorProductFunctor._functor_name is "tensor", so this method
        # also needs to double as the tensor product construction
        if isinstance(args[0], Parent):
            return self.tensor_product(*args, **kwds)
        return self._tensor(*args, **kwds)

    def tensor_from_comp(self, tensor_type, comp, name=None, latex_name=None):
        r"""
        Construct a tensor on ``self`` from a set of components.

        The tensor symmetries are deduced from those of the components.

        INPUT:

        - ``tensor_type`` -- pair ``(k, l)`` with ``k`` being the
          contravariant rank and ``l`` the covariant rank
        - ``comp`` -- instance of :class:`~sage.tensor.modules.comp.Components`
          representing the tensor components in a given basis
        - ``name`` -- (default: ``None``) string; name given to the tensor
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the tensor; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor defined on ``self`` with the provided
          characteristics.

        EXAMPLES:

        Construction of a tensor of rank 1::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, CompFullySym, CompFullyAntiSym
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e') ; e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: c = Components(ZZ, e, 1)
            sage: c[:]
            [0, 0, 0]
            sage: c[:] = [-1,4,2]
            sage: t = M.tensor_from_comp((1,0), c)
            sage: t
            Element of the Rank-3 free module M over the Integer Ring
            sage: t.display(e)
            -e_0 + 4 e_1 + 2 e_2
            sage: t = M.tensor_from_comp((0,1), c) ; t
            Linear form on the Rank-3 free module M over the Integer Ring
            sage: t.display(e)
            -e^0 + 4 e^1 + 2 e^2

        Construction of a tensor of rank 2::

            sage: c = CompFullySym(ZZ, e, 2)
            sage: c[0,0], c[1,2] = 4, 5
            sage: t = M.tensor_from_comp((0,2), c) ; t
            Symmetric bilinear form on the
             Rank-3 free module M over the Integer Ring
            sage: t.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: t.display(e)
            4 e^0⊗e^0 + 5 e^1⊗e^2 + 5 e^2⊗e^1
            sage: c = CompFullyAntiSym(ZZ, e, 2)
            sage: c[0,1], c[1,2] = 4, 5
            sage: t = M.tensor_from_comp((0,2), c) ; t
            Alternating form of degree 2 on the
             Rank-3 free module M over the Integer Ring
            sage: t.display(e)
            4 e^0∧e^1 + 5 e^1∧e^2
        """
        from .comp import CompWithSym, CompFullyAntiSym

        # 0/ Compatibility checks:
        if comp._ring is not self._ring:
            raise TypeError("the components are not defined on the same"
                            " ring as the module")
        if comp._frame not in self._known_bases:
            raise TypeError("the components are not defined on a basis of"
                            " the module")
        if comp._nid != tensor_type[0] + tensor_type[1]:
            raise TypeError("number of component indices not compatible with "
                            " the tensor type")

        # 1/ Construction of the tensor:
        if tensor_type == (1, 0):
            resu = self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type == (0, 1):
            resu = self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type[0] == 0 and tensor_type[1] > 1 and \
                                        isinstance(comp, CompFullyAntiSym):
            resu = self.alternating_form(tensor_type[1], name=name,
                                         latex_name=latex_name)
        elif tensor_type[0] > 1 and tensor_type[1] == 0 and \
                                        isinstance(comp, CompFullyAntiSym):
            resu = self.alternating_contravariant_tensor(tensor_type[0],
                                         name=name,
                                         latex_name=latex_name)
        else:
            resu = self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name)
            # Tensor symmetries deduced from those of comp:
            if isinstance(comp, CompWithSym):
                resu._sym = comp._sym
                resu._antisym = comp._antisym

        # 2/ Tensor components set to comp:
        resu._components[comp._frame] = comp

        return resu

    def alternating_contravariant_tensor(self, degree, name=None,
                                         latex_name=None):
        r"""
        Construct an alternating contravariant tensor on the free module.

        INPUT:

        - ``degree`` -- degree of the alternating contravariant tensor
          (i.e. its tensor rank)
        - ``name`` -- (default: ``None``) string; name given to the
          alternating contravariant tensor
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the alternating contravariant tensor; if none is
          provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`

        EXAMPLES:

        Alternating contravariant tensor on a rank-3 module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_contravariant_tensor(2, 'a') ; a
            Alternating contravariant tensor a of degree 2 on the
             Rank-3 free module M over the Integer Ring

        The nonzero components in a given basis have to be set in a second
        step, thereby fully specifying the alternating form::

            sage: e = M.basis('e') ; e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: a.set_comp(e)[0,1] = 2
            sage: a.set_comp(e)[1,2] = -3
            sage: a.display(e)
            a = 2 e_0∧e_1 - 3 e_1∧e_2

        An alternating contravariant tensor of degree 1 is simply
        an element of the module::

            sage: a = M.alternating_contravariant_tensor(1, 'a') ; a
            Element a of the Rank-3 free module M over the Integer Ring

        See
        :class:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor`
        for more documentation.
        """
        if degree == 1:
            return self.element_class(self, name=name,
                                      latex_name=latex_name)
        return self.exterior_power(degree).element_class(self, degree,
                                       name=name, latex_name=latex_name)

    def alternating_form(self, degree, name=None, latex_name=None):
        r"""
        Construct an alternating form on the free module.

        INPUT:

        - ``degree`` -- the degree of the alternating form (i.e. its
          tensor rank)
        - ``name`` -- (default: ``None``) string; name given to the
          alternating form
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the alternating form; if none is provided, the LaTeX symbol
          is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`

        EXAMPLES:

        Alternating forms on a rank-3 module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.alternating_form(2, 'a') ; a
            Alternating form a of degree 2 on the
             Rank-3 free module M over the Integer Ring

        The nonzero components in a given basis have to be set in a second
        step, thereby fully specifying the alternating form::

            sage: e = M.basis('e') ; e
            Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
            sage: a.set_comp(e)[0,1] = 2
            sage: a.set_comp(e)[1,2] = -3
            sage: a.display(e)
            a = 2 e^0∧e^1 - 3 e^1∧e^2

        An alternating form of degree 1 is a linear form::

            sage: a = M.alternating_form(1, 'a') ; a
            Linear form a on the Rank-3 free module M over the Integer Ring

        To construct such a form, it is preferable to call the method
        :meth:`linear_form` instead::

            sage: a = M.linear_form('a') ; a
            Linear form a on the Rank-3 free module M over the Integer Ring

        See
        :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`
        for more documentation.
        """
        if degree == 0:
            try:
                return self._ring.element_class(self._ring, name=name,
                                                latex_name=latex_name)
            except (KeyError, AttributeError):
                raise NotImplementedError('{} apparently '.format(self._ring) +
                                          'does not provide generic elements')
        return self.dual_exterior_power(degree).element_class(self, degree,
                                              name=name, latex_name=latex_name)

    def linear_form(self, name=None, latex_name=None):
        r"""
        Construct a linear form on the free module ``self``.

        A *linear form* on a free module `M` over a ring `R` is a map
        `M \rightarrow R` that is linear. It can be viewed as a tensor of type
        `(0,1)` on `M`.

        INPUT:

        - ``name`` -- (default: ``None``) string; name given to the linear
          form
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the linear form; if none is provided, the LaTeX symbol
          is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`

        EXAMPLES:

        Linear form on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.linear_form('A') ; a
            Linear form A on the Rank-3 free module M over the Integer Ring
            sage: a[:] = [2,-1,3]  # components w.r.t. the module's default basis (e)
            sage: a.display()
            A = 2 e^0 - e^1 + 3 e^2

        A linear form maps module elements to ring elements::

            sage: v = M([1,1,1])
            sage: a(v)
            4

        Test of linearity::

            sage: u = M([-5,-2,7])
            sage: a(3*u - 4*v) == 3*a(u) - 4*a(v)
            True

        See
        :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`
        for more documentation.
        """
        return self.dual_exterior_power(1).element_class(self, 1, name=name,
                                                         latex_name=latex_name)

    def automorphism(self, matrix=None, basis=None, name=None,
                     latex_name=None):
        r"""
        Construct a module automorphism of ``self``.

        Denoting ``self`` by `M`, an automorphism of ``self`` is an element
        of the general linear group `\mathrm{GL}(M)`.

        INPUT:

        - ``matrix`` -- (default: ``None``) matrix of size rank(M)*rank(M)
          representing the automorphism with respect to ``basis``;
          this entry can actually be any material from which a matrix of
          elements of ``self`` base ring can be constructed; the *columns* of
          ``matrix`` must be the components w.r.t. ``basis`` of
          the images of the elements of ``basis``. If ``matrix`` is ``None``,
          the automorphism has to be initialized afterwards by
          method :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.set_comp`
          or via the operator [].
        - ``basis`` -- (default: ``None``) basis of ``self`` defining the
          matrix representation; if ``None`` the default basis of ``self`` is
          assumed.
        - ``name`` -- (default: ``None``) string; name given to the
          automorphism
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the automorphism; if none is provided, the LaTeX symbol
          is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`

        EXAMPLES:

        Automorphism of a rank-2 free `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism(matrix=[[1,2],[1,3]], basis=e, name='a') ; a
            Automorphism a of the Rank-2 free module M over the Integer Ring
            sage: a.parent()
            General linear group of the Rank-2 free module M over the Integer Ring
            sage: a.matrix(e)
            [1 2]
            [1 3]

        An automorphism is a tensor of type (1,1)::

            sage: a.tensor_type()
            (1, 1)
            sage: a.display(e)
            a = e_0⊗e^0 + 2 e_0⊗e^1 + e_1⊗e^0 + 3 e_1⊗e^1

        The automorphism components can be specified in a second step, as
        components of a type-`(1,1)` tensor::

            sage: a1 = M.automorphism(name='a')
            sage: a1[e,:] = [[1,2],[1,3]]
            sage: a1.matrix(e)
            [1 2]
            [1 3]
            sage: a1 == a
            True

        Component by component specification::

            sage: a2 = M.automorphism(name='a')
            sage: a2[0,0] = 1  # component set in the module's default basis (e)
            sage: a2[0,1] = 2
            sage: a2[1,0] = 1
            sage: a2[1,1] = 3
            sage: a2.matrix(e)
            [1 2]
            [1 3]
            sage: a2 == a
            True

        See
        :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
        for more documentation.
        """
        resu = self.general_linear_group().element_class(self, name=name,
                                                         latex_name=latex_name)
        if matrix:
            if basis is None:
                basis = self.default_basis()
            resu.set_comp(basis)[:] = matrix
        return resu

    def sym_bilinear_form(self, name=None, latex_name=None):
        r"""
        Construct a symmetric bilinear form on the free module ``self``.

        INPUT:

        - ``name`` -- (default: ``None``) string; name given to the symmetric
          bilinear form
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the symmetric bilinear form; if none is provided, the LaTeX
          symbol is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          of tensor type `(0,2)` and symmetric

        EXAMPLES:

        Symmetric bilinear form on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.sym_bilinear_form('A') ; a
            Symmetric bilinear form A on the
             Rank-3 free module M over the Integer Ring

        A symmetric bilinear form is a type-`(0,2)` tensor that is symmetric::

            sage: a.parent()
            Free module of type-(0,2) tensors on the
             Rank-3 free module M over the Integer Ring
            sage: a.tensor_type()
            (0, 2)
            sage: a.tensor_rank()
            2
            sage: a.symmetries()
            symmetry: (0, 1);  no antisymmetry

        Components with respect to a given basis::

            sage: e = M.basis('e')
            sage: a[0,0], a[0,1], a[0,2] = 1, 2, 3
            sage: a[1,1], a[1,2] = 4, 5
            sage: a[2,2] = 6

        Only independent components have been set; the other ones are
        deduced by symmetry::

            sage: a[1,0], a[2,0], a[2,1]
            (2, 3, 5)
            sage: a[:]
            [1 2 3]
            [2 4 5]
            [3 5 6]

        A symmetric bilinear form acts on pairs of module elements::

            sage: u = M([2,-1,3]) ; v = M([-2,4,1])
            sage: a(u,v)
            61
            sage: a(v,u) == a(u,v)
            True

        The sum of two symmetric bilinear forms is another symmetric bilinear
        form::

            sage: b = M.sym_bilinear_form('B')
            sage: b[0,0], b[0,1], b[1,2] = -2, 1, -3
            sage: s = a + b ; s
            Symmetric bilinear form A+B on the
             Rank-3 free module M over the Integer Ring
            sage: a[:], b[:], s[:]
            (
            [1 2 3]  [-2  1  0]  [-1  3  3]
            [2 4 5]  [ 1  0 -3]  [ 3  4  2]
            [3 5 6], [ 0 -3  0], [ 3  2  6]
            )

        Adding a symmetric bilinear from with a non-symmetric one results in a
        generic type-`(0,2)` tensor::

            sage: c = M.tensor((0,2), name='C')
            sage: c[0,1] = 4
            sage: s = a + c ; s
            Type-(0,2) tensor A+C on the Rank-3 free module M over the Integer Ring
            sage: s.symmetries()
            no symmetry;  no antisymmetry
            sage: s[:]
            [1 6 3]
            [2 4 5]
            [3 5 6]

        See :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
        for more documentation.
        """
        return self.tensor_module(0,2).element_class(self, (0,2), name=name,
                                              latex_name=latex_name, sym=(0,1))

    #### End of methods to be redefined by derived classes ####

    def dual(self):
        r"""
        Return the dual module of ``self``.

        EXAMPLES:

        Dual of a free module over `\ZZ`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.dual()
            Dual of the Rank-3 free module M over the Integer Ring
            sage: latex(M.dual())
            M^*

        The dual is a free module of the same rank as M::

            sage: from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule_abstract
            sage: isinstance(M.dual(), FiniteRankFreeModule_abstract)
            True
            sage: M.dual().rank()
            3

        It is formed by alternating forms of degree 1, i.e. linear forms::

            sage: M.dual() is M.dual_exterior_power(1)
            True
            sage: M.dual().an_element()
            Linear form on the Rank-3 free module M over the Integer Ring
            sage: a = M.linear_form()
            sage: a in M.dual()
            True

        The elements of a dual basis belong of course to the dual module::

            sage: e = M.basis('e')
            sage: e.dual_basis()[0] in M.dual()
            True
        """
        return self.dual_exterior_power(1)

    def irange(self, start: Optional[int] = None) -> Generator[int, None, None]:
        r"""
        Single index generator, labelling the elements of a basis of ``self``.

        INPUT:

        - ``start`` -- (default: ``None``) integer; initial value of the
          index; if none is provided, ``self._sindex`` is assumed

        OUTPUT:

        - an iterable index, starting from ``start`` and ending at
          ``self._sindex + self.rank() - 1``

        EXAMPLES:

        Index range on a rank-3 module::

            sage: M = FiniteRankFreeModule(ZZ, 3)
            sage: list(M.irange())
            [0, 1, 2]
            sage: list(M.irange(start=1))
            [1, 2]

        The default starting value corresponds to the parameter ``start_index``
        provided at the module construction (the default value being 0)::

            sage: M1 = FiniteRankFreeModule(ZZ, 3, start_index=1)
            sage: list(M1.irange())
            [1, 2, 3]
            sage: M2 = FiniteRankFreeModule(ZZ, 3, start_index=-4)
            sage: list(M2.irange())
            [-4, -3, -2]
        """
        si = self._sindex
        imax = self._rank + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1

    def default_basis(self):
        r"""
        Return the default basis of the free module ``self``.

        The *default basis* is simply a basis whose name can be skipped in
        methods requiring a basis as an argument. By default, it is the first
        basis introduced on the module. It can be changed by the method
        :meth:`set_default_basis`.

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`

        EXAMPLES:

        At the module construction, no default basis is assumed::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M', start_index=1)
            sage: M.default_basis()
            No default basis has been defined on the
             Rank-2 free module M over the Integer Ring

        The first defined basis becomes the default one::

            sage: e = M.basis('e') ; e
            Basis (e_1,e_2) on the Rank-2 free module M over the Integer Ring
            sage: M.default_basis()
            Basis (e_1,e_2) on the Rank-2 free module M over the Integer Ring
            sage: f =  M.basis('f') ; f
            Basis (f_1,f_2) on the Rank-2 free module M over the Integer Ring
            sage: M.default_basis()
            Basis (e_1,e_2) on the Rank-2 free module M over the Integer Ring
        """
        if self._def_basis is None:
            print("No default basis has been defined on the {}".format(self))
        return self._def_basis

    def set_default_basis(self, basis):
        r"""
        Set the default basis of ``self``.

        The *default basis* is simply a basis whose name can be skipped in
        methods requiring a basis as an argument. By default, it is the first
        basis introduced on the module.

        INPUT:

        - ``basis`` -- instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing a basis on ``self``

        EXAMPLES:

        Changing the default basis on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e') ; e
            Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
            sage: f =  M.basis('f') ; f
            Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring
            sage: M.default_basis()
            Basis (e_1,e_2,e_3) on the Rank-3 free module M over the Integer Ring
            sage: M.set_default_basis(f)
            sage: M.default_basis()
            Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring
        """
        from .free_module_basis import FreeModuleBasis
        if not isinstance(basis, FreeModuleBasis):
            raise TypeError("the argument is not a free module basis")
        if basis._fmodule is not self:
            raise ValueError("the basis is not defined on the current module")
        self._def_basis = basis

    def print_bases(self):
        r"""
        Display the bases that have been defined on the free module ``self``.

        Use the method :meth:`bases` to get the raw list of bases.

        EXAMPLES:

        Bases on a rank-4 free module::

            sage: M = FiniteRankFreeModule(ZZ, 4, name='M', start_index=1)
            sage: M.print_bases()
            No basis has been defined on the
             Rank-4 free module M over the Integer Ring
            sage: e = M.basis('e')
            sage: M.print_bases()
            Bases defined on the Rank-4 free module M over the Integer Ring:
             - (e_1,e_2,e_3,e_4) (default basis)
            sage: f = M.basis('f')
            sage: M.print_bases()
            Bases defined on the Rank-4 free module M over the Integer Ring:
             - (e_1,e_2,e_3,e_4) (default basis)
             - (f_1,f_2,f_3,f_4)
            sage: M.set_default_basis(f)
            sage: M.print_bases()
            Bases defined on the Rank-4 free module M over the Integer Ring:
             - (e_1,e_2,e_3,e_4)
             - (f_1,f_2,f_3,f_4) (default basis)

        """
        if not self._known_bases:
            print("No basis has been defined on the {}".format(self))
        else:
            print("Bases defined on the {}:".format(self))
            for basis in self._known_bases:
                item = " - " + basis._name
                if basis is self._def_basis:
                    item += " (default basis)"
                print(item)

    def bases(self):
        r"""
        Return the list of bases that have been defined on the free module
        ``self``.

        Use the method :meth:`print_bases` to get a formatted output with more
        information.

        OUTPUT:

        - list of instances of class
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`

        EXAMPLES:

        Bases on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M_3', start_index=1)
            sage: M.bases()
            []
            sage: e = M.basis('e')
            sage: M.bases()
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M_3 over the Integer Ring]
            sage: f = M.basis('f')
            sage: M.bases()
            [Basis (e_1,e_2,e_3) on the Rank-3 free module M_3 over the Integer Ring,
             Basis (f_1,f_2,f_3) on the Rank-3 free module M_3 over the Integer Ring]
        """
        return list(self._known_bases)

    def change_of_basis(self, basis1, basis2):
        r"""
        Return a module automorphism linking two bases defined on the free
        module ``self``.

        If the automorphism has not been recorded yet (in the internal
        dictionary ``self._basis_changes``), it is computed by transitivity,
        i.e. by performing products of recorded changes of basis.

        INPUT:

        - ``basis1`` -- a basis of ``self``, denoted `(e_i)` below
        - ``basis2`` -- a basis of ``self``, denoted `(f_i)` below

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the basis `(e_i)` to the
          basis `(f_i)` according to `f_i = P(e_i)`

        EXAMPLES:

        Changes of basis on a rank-2 free module::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(ZZ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: f = M.basis('f', from_family=(e[1]+2*e[2], e[1]+3*e[2]))
            sage: P = M.change_of_basis(e,f) ; P
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: P.matrix(e)
            [1 1]
            [2 3]

        Note that the columns of this matrix contain the components of the
        elements of basis ``f`` w.r.t. to basis ``e``::

            sage: f[1].display(e)
            f_1 = e_1 + 2 e_2
            sage: f[2].display(e)
            f_2 = e_1 + 3 e_2

        The change of basis is cached::

            sage: P is M.change_of_basis(e,f)
            True

        Check of the change-of-basis automorphism::

            sage: f[1] == P(e[1])
            True
            sage: f[2] == P(e[2])
            True

        Check of the reverse change of basis::

            sage: M.change_of_basis(f,e) == P^(-1)
            True

        We have of course::

            sage: M.change_of_basis(e,e)
            Identity map of the Rank-2 free module M over the Integer Ring
            sage: M.change_of_basis(e,e) is M.identity_map()
            True

        Let us introduce a third basis on ``M``::

            sage: h = M.basis('h', from_family=(3*e[1]+4*e[2], 5*e[1]+7*e[2]))

        The change of basis ``e`` --> ``h`` has been recorded directly from the
        definition of ``h``::

            sage: Q = M.change_of_basis(e,h) ; Q.matrix(e)
            [3 5]
            [4 7]

        The change of basis ``f`` --> ``h`` is computed by transitivity, i.e.
        from the changes of basis ``f`` --> ``e`` and ``e`` --> ``h``::

            sage: R = M.change_of_basis(f,h) ; R
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: R.matrix(e)
            [-1  2]
            [-2  3]
            sage: R.matrix(f)
            [ 5  8]
            [-2 -3]

        Let us check that ``R`` is indeed the change of basis ``f`` --> ``h``::

            sage: h[1] == R(f[1])
            True
            sage: h[2] == R(f[2])
            True

        A related check is::

            sage: R == Q*P^(-1)
            True
        """
        if basis1 == basis2:
            return self.identity_map()
        bc = self._basis_changes
        if (basis1, basis2) not in bc:
            if basis1 not in self._known_bases:
                raise TypeError("{} is not a basis of the {}".format(basis1,
                                                                     self))
            if basis2 not in self._known_bases:
                raise TypeError("{} is not a basis of the {}".format(basis2,
                                                                     self))
            # Is the inverse already registered ?
            if (basis2, basis1) in bc:
                inv = bc[(basis2, basis1)].inverse()
                bc[(basis1, basis2)] = inv
                return inv
            # Search for a third basis, basis say, such that either the changes
            # basis1 --> basis and basis --> basis2
            # or
            # basis2 --> basis and basis --> basis1
            # are known:
            for basis in self._known_bases:
                if (basis1, basis) in bc and (basis, basis2) in bc:
                    transf = bc[(basis, basis2)] * bc[(basis1, basis)]
                    bc[(basis1, basis2)] = transf
                    bc[(basis2, basis1)] = transf.inverse()
                    break
                if (basis2, basis) in bc and (basis, basis1) in bc:
                    inv = bc[(basis, basis1)] * bc[(basis2, basis)]
                    bc[(basis2, basis1)] = inv
                    bc[(basis1, basis2)] = inv.inverse()
                    break
            else:
                raise ValueError(("the change of basis from '{!r}' to '{!r}'"
                                + " cannot be computed"
                                ).format(basis1, basis2))
        return bc[(basis1, basis2)]

    def set_change_of_basis(self, basis1, basis2, change_of_basis,
                            compute_inverse=True):
        r"""
        Relates two bases by an automorphism of ``self``.

        This updates the internal dictionary ``self._basis_changes``.

        INPUT:

        - ``basis1`` -- basis 1, denoted `(e_i)` below
        - ``basis2`` -- basis 2, denoted `(f_i)` below
        - ``change_of_basis`` -- instance of class
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the basis `(e_i)` to
          the basis `(f_i)` according to `f_i = P(e_i)`
        - ``compute_inverse`` -- boolean (default: ``True``); if set to ``True``,
          the inverse automorphism is computed and the change from basis `(f_i)`
          to `(e_i)` is set to it in the internal dictionary
          ``self._basis_changes``

        EXAMPLES:

        Defining a change of basis on a rank-2 free module::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: e = M.basis('e')
            sage: f = M.basis('f')
            sage: a = M.automorphism()
            sage: a[:] = [[1, 2], [-1, 3]]
            sage: M.set_change_of_basis(e, f, a)

        The change of basis and its inverse have been recorded::

            sage: M.change_of_basis(e,f).matrix(e)
            [ 1  2]
            [-1  3]
            sage: M.change_of_basis(f,e).matrix(e)
            [ 3/5 -2/5]
            [ 1/5  1/5]

        and are effective::

            sage: f[0].display(e)
            f_0 = e_0 - e_1
            sage: e[0].display(f)
            e_0 = 3/5 f_0 + 1/5 f_1
        """
        if basis1 not in self._known_bases:
            raise TypeError("{} is not a basis of the {}".format(basis1,
                                                                 self))
        if basis2 not in self._known_bases:
            raise TypeError("{} is not a basis of the {}".format(basis2,
                                                                 self))
        if change_of_basis not in self.general_linear_group():
            raise TypeError("{} is not an automorphism of the {}".format(
                                                        change_of_basis, self))
        self._basis_changes[(basis1, basis2)] = change_of_basis
        if compute_inverse:
            self._basis_changes[(basis2, basis1)] = change_of_basis.inverse()

    def hom(self, codomain, matrix_rep, bases=None, name=None,
            latex_name=None):
        r"""
        Homomorphism from ``self`` to a free module.

        Define a module homomorphism

        .. MATH::

            \phi:\ M \longrightarrow N,

        where `M` is ``self`` and  `N` is a free module of finite rank
        over the same ring `R` as ``self``.

        .. NOTE::

            This method is a redefinition of
            :meth:`sage.structure.parent.Parent.hom` because the latter assumes
            that ``self`` has some privileged generators, while an instance of
            :class:`FiniteRankFreeModule` has no privileged basis.

        INPUT:

        - ``codomain`` -- the target module `N`
        - ``matrix_rep`` -- matrix of size rank(N)*rank(M) representing the
          homomorphism with respect to the pair of bases defined by ``bases``;
          this entry can actually be any material from which a matrix of
          elements of `R` can be constructed; the *columns* of
          ``matrix_rep`` must be the components w.r.t. ``basis_N`` of
          the images of the elements of ``basis_M``.
        - ``bases`` -- (default: ``None``) pair ``(basis_M, basis_N)`` defining
          the matrix representation, ``basis_M`` being a basis of ``self`` and
          ``basis_N`` a basis of module `N` ; if None the pair formed by the
          default bases of each module is assumed.
        - ``name`` -- (default: ``None``) string; name given to the
          homomorphism
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the homomorphism. If ``None``, ``name`` will be used.

        OUTPUT:

        - the homomorphism `\phi: M \rightarrow N` corresponding to the given
          specifications, as an instance of
          :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

        EXAMPLES:

        Homomorphism between two free modules over `\ZZ`::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e')
            sage: f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]]) ; phi
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring

        Homomorphism defined by a matrix w.r.t. bases that are not the
        default ones::

            sage: ep = M.basis('ep', latex_symbol=r"e'")
            sage: fp = N.basis('fp', latex_symbol=r"f'")
            sage: phi = M.hom(N, [[3,2,1], [1,2,3]], bases=(ep, fp)) ; phi
            Generic morphism:
              From: Rank-3 free module M over the Integer Ring
              To:   Rank-2 free module N over the Integer Ring

        Call with all arguments specified::

            sage: phi = M.hom(N, [[3,2,1], [1,2,3]], bases=(ep, fp),
            ....:             name='phi', latex_name=r'\phi')

        The parent::

            sage: phi.parent() is Hom(M,N)
            True

        See class
        :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`
        for more documentation.
        """
        homset = Hom(self, codomain)
        return homset(matrix_rep, bases=bases, name=name,
                      latex_name=latex_name)

    def endomorphism(self, matrix_rep, basis=None, name=None, latex_name=None):
        r"""
        Construct an endomorphism of the free module ``self``.

        The returned object is a module morphism `\phi: M \rightarrow M`,
        where `M` is ``self``.

        INPUT:

        - ``matrix_rep`` -- matrix of size rank(M)*rank(M) representing the
          endomorphism with respect to ``basis``;
          this entry can actually be any material from which a matrix of
          elements of ``self`` base ring can be constructed; the *columns* of
          ``matrix_rep`` must be the components w.r.t. ``basis`` of
          the images of the elements of ``basis``.
        - ``basis`` -- (default: ``None``) basis of ``self`` defining the
          matrix representation; if None the default basis of ``self`` is
          assumed.
        - ``name`` -- (default: ``None``) string; name given to the
          endomorphism
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the endomorphism. If none is provided, ``name`` will be used.

        OUTPUT:

        - the endomorphism `\phi: M \rightarrow M` corresponding to the given
          specifications, as an instance of
          :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

        EXAMPLES:

        Construction of an endomorphism with minimal data (module's default
        basis and no name)::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: phi = M.endomorphism([[1,-2], [-3,4]]) ; phi
            Generic endomorphism of Rank-2 free module M over the Integer Ring
            sage: phi.matrix()  # matrix w.r.t the default basis
            [ 1 -2]
            [-3  4]

        Construction with full list of arguments (matrix given a basis
        different from the default one)::

            sage: a = M.automorphism() ; a[0,1], a[1,0] = 1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: phi = M.endomorphism([[1,-2], [-3,4]], basis=ep, name='phi',
            ....:                      latex_name=r'\phi')
            sage: phi
            Generic endomorphism of Rank-2 free module M over the Integer Ring
            sage: phi.matrix(ep)  # the input matrix
            [ 1 -2]
            [-3  4]
            sage: phi.matrix()  # matrix w.r.t the default basis
            [4 3]
            [2 1]

        See :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`
        for more documentation.
        """
        from sage.categories.homset import End
        if basis is None:
            basis = self.default_basis()
        return End(self)(matrix_rep, bases=(basis,basis), name=name,
                         latex_name=latex_name)

    def identity_map(self, name='Id', latex_name=None):
        r"""
        Return the identity map of the free module ``self``.

        INPUT:

        - ``name`` -- (string; default: 'Id') name given to the identity
          identity map
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the identity map; if none is provided, the LaTeX symbol is set to
          '\mathrm{Id}' if ``name`` is 'Id' and to ``name`` otherwise

        OUTPUT:

        - the identity map of ``self`` as an instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`

        EXAMPLES:

        Identity map of a rank-3 `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: Id = M.identity_map() ; Id
            Identity map of the Rank-3 free module M over the Integer Ring
            sage: Id.parent()
            General linear group of the Rank-3 free module M over the Integer Ring
            sage: Id.matrix(e)
            [1 0 0]
            [0 1 0]
            [0 0 1]

        The default LaTeX symbol::

            sage: latex(Id)
            \mathrm{Id}

        It can be changed by means of the method
        :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.set_name`::

            sage: Id.set_name(latex_name=r'\mathrm{1}_M')
            sage: latex(Id)
            \mathrm{1}_M

        The identity map is actually the identity element of GL(M)::

            sage: Id is M.general_linear_group().one()
            True

        It is also a tensor of type-`(1,1)` on M::

            sage: Id.tensor_type()
            (1, 1)
            sage: Id.comp(e)
            Kronecker delta of size 3x3
            sage: Id[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        Example with a LaTeX symbol different from the default one and set
        at the creation of the object::

            sage: N = FiniteRankFreeModule(ZZ, 3, name='N')
            sage: f = N.basis('f')
            sage: Id = N.identity_map(name='Id_N', latex_name=r'\mathrm{Id}_N')
            sage: Id
            Identity map of the Rank-3 free module N over the Integer Ring
            sage: latex(Id)
            \mathrm{Id}_N
        """
        if self._identity_map is None:
            self._identity_map = self.general_linear_group().one()
            if name != 'Id':
                if latex_name is None:
                    latex_name = name
                self._identity_map.set_name(name=name, latex_name=latex_name)
        return self._identity_map


class FiniteRankDualFreeModule(ReflexiveModule_dual, FiniteRankFreeModule_abstract):
    r"""
    Dual of a free module of finite rank over a commutative ring.

    Given a free module `M` of finite rank over a commutative ring `R`,
    the *dual of* `M` is the set `M^*` of all linear forms on `M`,
    i.e., linear maps

    .. MATH::

        M \longrightarrow R

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank, as an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (default: ``None``) string; name given to `M^*`
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote `M^*`

    EXAMPLES::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: A = M.dual(); A
        Dual of the Rank-3 free module M over the Integer Ring

    ``A`` is a module (actually a free module) over `\ZZ`::

        sage: A.category()
        Category of finite dimensional modules over Integer Ring
        sage: A in Modules(ZZ)
        True
        sage: A.rank()
        3
        sage: A.base_ring()
        Integer Ring
        sage: A.base_module()
        Rank-3 free module M over the Integer Ring

    ``A`` is a *parent* object, whose elements are linear forms,
    represented by instances of the class
    :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`::

        sage: a = A.an_element() ; a
        Linear form on the Rank-3 free module M over the Integer Ring
        sage: a.display() # expansion with respect to M's default basis (e)
        e^0
        sage: from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
        sage: isinstance(a, FreeModuleAltForm)
        True
        sage: a in A
        True
        sage: A.is_parent_of(a)
        True

    Elements can be constructed from ``A``. In particular, 0 yields
    the zero element of ``A``::

        sage: A(0)
        Linear form zero on the Rank-3 free module M over the Integer Ring
        sage: A(0) is A.zero()
        True

    while nonzero elements are constructed by providing their components in a
    given basis::

        sage: e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: comp = [0,3,-1]
        sage: a = A(comp, basis=e, name='a') ; a
        Linear form a on the Rank-3 free module M over the Integer Ring
        sage: a.display(e)
        a = 3 e^1 - e^2

    An alternative is to construct the alternating form from an empty list of
    components and to set the nonzero components afterwards::

        sage: a = A([], name='a')
        sage: a.set_comp(e)[0] = 3
        sage: a.set_comp(e)[1] = -1
        sage: a.set_comp(e)[2] = 4
        sage: a.display(e)
        a = 3 e^0 - e^1 + 4 e^2

    The dual is unique::

        sage: A is M.dual()
        True

    The exterior power `\Lambda^1(M^*)` is nothing but `M^*`::

        sage: M.dual_exterior_power(1) is M.dual()
        True

    It also coincides with the module of type-`(0,1)` tensors::

        sage: M.dual_exterior_power(1) is M.tensor_module(0,1)
        True
    """

    Element = FreeModuleAltForm

    def __init__(self, fmodule, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.finite_rank_free_module import FiniteRankDualFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = FiniteRankDualFreeModule(M) ; A
            Dual of the Rank-3 free module M over the Integer Ring
            sage: TestSuite(A).run()
        """
        self._fmodule = fmodule
        rank = fmodule._rank
        if name is None and fmodule._name is not None:
            name = fmodule._name + '*'
        if latex_name is None and fmodule._latex_name is not None:
            latex_name = fmodule._latex_name + r'^*'
        super().__init__(fmodule._ring, rank, name=name,
                         latex_name=latex_name)
        fmodule._all_modules.add(self)

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct a linear form.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.dual()
            sage: a = A._element_constructor_(0) ; a
            Linear form zero on the Rank-3 free module M over the Integer Ring
            sage: a = A._element_constructor_([2,0,-1], name='a') ; a
            Linear form a on the Rank-3 free module M over the Integer Ring
            sage: a.display()
            a = 2 e^0 - e^2
        """
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        if isinstance(comp, FreeModuleTensor):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if tensor.tensor_type() == (0,1) and self._degree == 1 and \
                                         tensor.base_module() is self._fmodule:
                resu = self.element_class(self._fmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for basis, comp in tensor._components.items():
                    resu._components[basis] = comp.copy()
                return resu
            else:
                raise TypeError("cannot coerce the {} ".format(tensor) +
                                "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._fmodule, 1, name=name, latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) alternating form.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 4, name='M')
            sage: e = M.basis('e')
            sage: a = M.dual()._an_element_() ; a
            Linear form on the 4-dimensional vector space M over the Rational
             Field
            sage: a.display()
            1/2 e^0

        TESTS:

        When the base module has no default basis, a default
        basis will be set for it::

            sage: M2 = FiniteRankFreeModule(QQ, 4, name='M2')
            sage: a = M2.dual()._an_element_(); a
            Linear form on the 4-dimensional vector space M2 over the Rational Field
            sage: a + a
            Linear form on the 4-dimensional vector space M2 over the Rational Field
            sage: M2.default_basis()
            Basis (e_0,e_1,e_2,e_3) on the 4-dimensional vector space M2 over the Rational Field
        """
        resu = self.element_class(self._fmodule, 1)
        # Make sure that the base module has a default basis
        self._fmodule.an_element()
        sindex = self._fmodule._sindex
        ind = [sindex + i for i in range(resu._tensor_rank)]
        resu.set_comp()[ind] = self._fmodule._ring.an_element()
        return resu

    #### End of parent methods

    @cached_method
    def zero(self):
        r"""
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: A = M.dual()
            sage: A.zero()
            Linear form zero on the Rank-3 free module M over the Integer Ring
            sage: A(0) is A.zero()
            True
        """
        resu = self._element_constructor_(name='zero', latex_name='0')
        for basis in self._fmodule._known_bases:
            resu._components[basis] = resu._new_comp(basis)
            # (since new components are initialized to zero)
        resu._is_zero = True # This element is certainly zero
        resu.set_immutable()
        return resu

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: M.dual_exterior_power(1)._repr_()
            'Dual of the Rank-5 free module M over the Integer Ring'
        """
        return "Dual of the {}".format(self._fmodule)

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module on which the dual is defined.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 5, name='M')
            sage: A = M.dual()
            sage: A.base_module()
            Rank-5 free module M over the Integer Ring
            sage: A.base_module() is M
            True
        """
        return self._fmodule
