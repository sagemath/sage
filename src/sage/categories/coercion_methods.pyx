# sage_setup: distribution = sagemath-categories
"""
Coercion methods for categories

The purpose of this Cython module is to hold special coercion methods,
which are inserted by their respective categories.
"""

from sage.structure.element cimport Element
cimport cython


@cython.binding
def _mul_parent(self, other):
    r"""
    Return the product of the two elements, calculated using
    the ``product`` method of the parent.

    This is inserted by :meth:`Magmas.ParentMethods.__init_extra__` as
    default implementation of ``Magmas.ElementMethods._mul_`` if
    ``product`` is implemented in the parent.

    INPUT:

    - ``other`` -- an element of the parent of ``self``

    OUTPUT: an element of the parent of ``self``

    EXAMPLES::

        sage: S = Semigroups().example("free")
        sage: x = S('a'); y = S('b')
        sage: x._mul_parent(y)
        'ab'

    .. SEEALSO::

        - :meth:`Magmas.ElementMethods._mul_parent`
        - :meth:`Magmas.ElementMethods.__init_extra__`
        - :meth:`Magmas.ParentMethods.product`

    This is :meth:`Magmas.ElementMethods._mul_parent`, implemented as
    a Cython method in :mod:`sage.categories.coercion_methods`::

        sage: from sage.cpython.getattr import raw_getattr
        sage: x._mul_parent.__func__ is raw_getattr(Magmas.ElementMethods,
        ....:                                       '_mul_parent')
        True
        sage: x._mul_parent.__func__ is sage.categories.coercion_methods._mul_parent
        True
    """
    return (<Element>self)._parent.product(self, other)
