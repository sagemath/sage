# -*- coding: utf-8 -*-
"""
Signed Tensor Product Functorial Construction

AUTHORS:

- Travis Scrimshaw (2019-07): initial version
"""
# ****************************************************************************
#  Copyright (C) 2019 Travis Scrimshaw <tcscrims at gamil.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.covariant_functorial_construction import CovariantFunctorialConstruction, CovariantConstructionCategory
from sage.typeset.unicode_characters import unicode_otimes


class SignedTensorProductFunctor(CovariantFunctorialConstruction):
    """
    A singleton class for the signed tensor functor.

    This functor takes a collection of graded algebras (possibly with
    basis) and constructs the signed tensor product of those algebras.
    If this algebra is in a subcategory, say that of
    ``Algebras(QQ).Graded()``, it is automatically endowed with
    its natural algebra structure, thanks to the category
    ``Algebras(QQ).Graded().SignedTensorProducts()`` of signed tensor
    products of graded algebras. For elements, it constructs the natural
    tensor product element in the corresponding tensor product of their
    parents.

    The signed tensor functor is covariant: if ``A`` is a subcategory
    of ``B``, then ``A.SignedTensorProducts()`` is a subcategory of
    ``B.SignedTensorProducts()`` (see also
    :class:`~sage.categories.covariant_functorial_construction.CovariantFunctorialConstruction`). Hence,
    the role of ``Algebras(QQ).Graded().SignedTensorProducts()`` is solely
    to provide mathematical information and algorithms which are relevant to
    signed tensor product of graded algebras.

    Those are implemented in the nested class
    :class:`~sage.categories.algebras.Algebras.SignedTensorProducts`
    of ``Algebras(QQ).Graded()``. This nested class is itself a subclass of
    :class:`~sage.categories.signed_tensor.SignedTensorProductsCategory`.

    EXAMPLES::

        sage: tensor_signed
        The signed tensor functorial construction

    TESTS::

        sage: TestSuite(tensor_signed).run()
    """
    _functor_name = "tensor"
    _functor_category = "SignedTensorProducts"
    symbol = " # "
    unicode_symbol = f" {unicode_otimes} "

    def _repr_(self):
        """
        EXAMPLES::

            sage: tensor_signed
            The signed tensor functorial construction
        """
        # We need this to distinguish it from tensor(), but we want
        #   it to have the same _functor_name (which is used in the
        #   default _repr_
        return "The signed tensor functorial construction"


tensor_signed = SignedTensorProductFunctor()


class SignedTensorProductsCategory(CovariantConstructionCategory):
    r"""
    An abstract base class for all SignedTensorProducts's categories.

    TESTS::

        sage: C = AlgebrasWithBasis(QQ).Graded().SignedTensorProducts()
        sage: C
        Category of signed tensor products of graded algebras with basis over Rational Field
        sage: C.base_category()
        Category of graded algebras with basis over Rational Field
        sage: latex(C)
        \mathbf{SignedTensorProducts}(\mathbf{GradedAlgebrasWithBasis}(\mathbf{AlgebrasWithBasis}_{\Bold{Q}}))
        sage: TestSuite(C).run()
    """
    _functor_category = "SignedTensorProducts"

    def SignedTensorProducts(self):
        r"""
        Return the category of signed tensor products of objects of ``self``.

        By associativity of signed tensor products, this is ``self`` (a tensor
        product of signed tensor products of `Cat`'s is a tensor product of
        `Cat`'s with the same twisting morphism)

        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Graded().SignedTensorProducts().SignedTensorProducts()
            Category of signed tensor products of graded algebras with basis
             over Rational Field
        """
        return self

    def base(self):
        """
        The base of a signed tensor product is the base (usually a ring)
        of the underlying category.

        EXAMPLES::

            sage: AlgebrasWithBasis(ZZ).Graded().SignedTensorProducts().base()
            Integer Ring
        """
        return self.base_category().base()
