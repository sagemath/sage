# sage_setup: distribution = sagemath-categories
r"""
Bimodules
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category, CategoryWithParameters
from sage.categories.left_modules import LeftModules
from sage.categories.right_modules import RightModules

from sage.categories.rings import Rings
_Rings = Rings()

#?class Bimodules(Category_over_base_rng, Category_over_base_rng):


class Bimodules(CategoryWithParameters):
    """
    The category of `(R,S)`-bimodules.

    For `R` and `S` rings, a `(R,S)`-bimodule `X` is a left `R`-module
    and right `S`-module such that the left and right actions commute:
    `r*(x*s) = (r*x)*s`.

    EXAMPLES::

        sage: Bimodules(QQ, ZZ)
        Category of bimodules over Rational Field on the left and Integer Ring on the right
        sage: Bimodules(QQ, ZZ).super_categories()
        [Category of left modules over Rational Field, Category of right modules over Integer Ring]
    """

    def __init__(self, left_base, right_base, name=None):
        """
        The ``name`` parameter is ignored.

        EXAMPLES::

            sage: C = Bimodules(QQ, ZZ)
            sage: TestSuite(C).run()
        """
        if not (left_base in Rings() or
                (isinstance(left_base, Category)
                 and left_base.is_subcategory(Rings()))):
            raise ValueError("the left base must be a ring or a subcategory of Rings()")
        if not (right_base in Rings() or
                (isinstance(right_base, Category)
                 and right_base.is_subcategory(Rings()))):
            raise ValueError("the right base must be a ring or a subcategory of Rings()")
        self._left_base_ring = left_base
        self._right_base_ring = right_base
        Category.__init__(self)

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        Since :issue:`11935`, the element and parent classes of a
        bimodule only depend on the categories of the left and right
        base ring.

        .. SEEALSO::

            - :meth:`CategoryWithParameters`
            - :meth:`CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Bimodules(QQ,ZZ)._make_named_class_key('parent_class')
            (Join of Category of number fields
                 and Category of quotient fields
                 and Category of metric spaces,
             Join of Category of Dedekind domains
                 and Category of euclidean domains
                 and Category of noetherian rings
                 and Category of infinite enumerated sets
                 and Category of metric spaces)


            sage: Bimodules(Fields(), ZZ)._make_named_class_key('element_class')
            (Category of fields,
             Join of Category of Dedekind domains
             and Category of euclidean domains
             and Category of noetherian rings
             and Category of infinite enumerated sets
             and Category of metric spaces)

            sage: Bimodules(QQ, Rings())._make_named_class_key('element_class')
            (Join of Category of number fields
                 and Category of quotient fields
                 and Category of metric spaces,
             Category of rings)

            sage: Bimodules(Fields(), Rings())._make_named_class_key('element_class')
            (Category of fields, Category of rings)
        """
        return (getattr(self._left_base_ring if isinstance(self._left_base_ring, Category) else self._left_base_ring.category(), name),
                getattr(self._right_base_ring if isinstance(self._right_base_ring, Category) else self._right_base_ring.category(), name))

    @classmethod
    def an_instance(cls):
        """
        Return an instance of this class.

        EXAMPLES::

            sage: Bimodules.an_instance()                                               # needs sage.rings.real_mpfr
            Category of bimodules over Rational Field on the left and Real Field with 53 bits of precision on the right
        """
        from sage.rings.rational_field import QQ
        from sage.rings.real_mpfr import RR
        return cls(QQ, RR)

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Bimodules(QQ, ZZ) # indirect doctest
            Category of bimodules over Rational Field on the left and Integer Ring on the right
        """
        return "bimodules over %s on the left and %s on the right" \
            % (self._left_base_ring, self._right_base_ring)

    def left_base_ring(self):
        """
        Return the left base ring over which elements of this category are
        defined.

        EXAMPLES::

            sage: Bimodules(QQ, ZZ).left_base_ring()
            Rational Field
        """
        return self._left_base_ring

    def right_base_ring(self):
        """
        Return the right base ring over which elements of this category are
        defined.

        EXAMPLES::

            sage: Bimodules(QQ, ZZ).right_base_ring()
            Integer Ring
        """
        return self._right_base_ring

    def _latex_(self) -> str:
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: print(Bimodules(QQ, ZZ)._latex_())
            {\mathbf{Bimodules}}_{\Bold{Q}, \Bold{Z}}
        """
        from sage.misc.latex import latex
        return "{{{}}}_{{{}, {}}}".format(Category._latex_(self),
                                          latex(self._left_base_ring),
                                          latex(self._right_base_ring))

    def super_categories(self):
        """
        EXAMPLES::

            sage: Bimodules(QQ, ZZ).super_categories()
            [Category of left modules over Rational Field, Category of right modules over Integer Ring]
        """
        R = self.left_base_ring()
        S = self.right_base_ring()
        return [LeftModules(R), RightModules(S)]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of bimodules defines no additional
        structure: a left and right module morphism between two
        bimodules is a bimodule morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`CategoryWithAxiom`?

        EXAMPLES::

            sage: Bimodules(QQ, ZZ).additional_structure()
        """
        return None

    class ParentMethods:
        pass

    class ElementMethods:
        pass
