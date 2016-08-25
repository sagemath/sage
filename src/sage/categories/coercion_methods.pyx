"""
Coercion methods for categories

Recall that the Sage categories provide implementations for Python's
special methods for the basic arithmetic operations (e.g. ``__mul__``,
``__add__``), that dispatch to the coercion model or call Sage's special
methods (e.g. ``_mul_``, ``_add_``) for the internal operations.

To reduce the induced overhead, we want those methods to be
Cythonized, while for most of the non speed-critical methods
plain Python is more practical to reduce compilation time and
simplify interactive debugging and introspection.

The purpose of this Cython module is to hold all the coercion methods,
which are inserted by their respective categories.
"""

from sage.structure.element cimport Element, have_same_parent_c, coercion_model, parent_c
import operator
cimport cython

# This could eventually be moved to SageObject
@cython.binding
def __add__(Element self, right):
    r"""
    Return the sum of ``self`` and ``right``.

    This calls the ``_add_`` method of ``self``, if it is
    available and the two elements have the same parent.

    Otherwise, the job is delegated to the coercion model.

    Do not override; instead implement an ``_add_`` method in the
    element class or a ``summation`` method in the parent class.

    .. SEEALSO:: :meth:`AdditiveMagmas.ElementMethods._add_`

    EXAMPLES::

        sage: F = CommutativeAdditiveSemigroups().example()
        sage: (a,b,c,d) = F.additive_semigroup_generators()
        sage: a + b
        a + b
        sage: a.__add__(b)
        a + b

    This is :meth:`AdditiveMagmas.ElementMethods.__add__`, implemented as a
    Cython method in :mod:`sage.categories.coercion_methods`::

        sage: a.__add__.im_func is AdditiveMagmas.ElementMethods.__add__.im_func
        True
        sage: a.__add__.im_func is sage.categories.coercion_methods.__add__
        True
    """
    if have_same_parent_c(self, right) and hasattr(self, "_add_"):
        return self._add_(right)
    return coercion_model.bin_op(self, right, operator.add)

@cython.binding
def __radd__(Element self, left):
    r"""
    Handles the sum of two elements, when the left hand side
    needs to be coerced first.

    EXAMPLES::

        sage: F = CommutativeAdditiveSemigroups().example()
        sage: (a,b,c,d) = F.additive_semigroup_generators()
        sage: a.__radd__(b)
        a + b

    This is :meth:`AdditiveMagmas.ElementMethods.__radd__`, implemented
    as a Cython method in :mod:`sage.categories.coercion_methods`::

        sage: a.__radd__.im_func is AdditiveMagmas.ElementMethods.__radd__.im_func
        True
        sage: a.__radd__.im_func is sage.categories.coercion_methods.__radd__
        True
    """
    if have_same_parent_c(left, self) and hasattr(left, "_add_"):
        return left._add_(self)
    return coercion_model.bin_op(left, self, operator.add)


@cython.binding
def Modules__mul__(Element left, right):
    """
    Return the product of ``left`` and ``right``.

    INPUT:

    - ``left`` -- an element of a :class:`module <Modules>`
    - ``right`` -- any object

    EXAMPLES:

    This is used when multiplying an element of a module on the right
    by something, typically a coefficient::

        sage: F = CombinatorialFreeModule(QQ, ["a", "b"])
        sage: x = F.monomial("a")
        sage: x * int(2)
        2*B['a']

    .. SEEALSO:: :meth:`Modules.ElementMethods.__rmul__`

    This is :meth:`Modules.ElementMethods.__rmul__`, implemented as a
    Cython method in :mod:`sage.categories.magmas_cython`::

        sage: x.__mul__.im_func is Modules.ElementMethods.__mul__.im_func
        True
        sage: x.__mul__.im_func is sage.categories.coercion_methods.Modules__mul__
        True

    .. TODO::

        Make a better unit test once ``Modules().example()`` is implemented.
    """
    return coercion_model.bin_op(left, right, operator.mul)

@cython.binding
def Modules__rmul__(Element right, left):
    """
    Return the product of ``left`` and ``right``.

    INPUT:

    - ``right`` -- an element of a :class:`module <Modules>`
    - ``left`` -- any object

    EXAMPLES:

    This is used when multiplying an element of a module on the left
    by something, typically a coefficient::

        sage: F = CombinatorialFreeModule(QQ, ["a", "b"])
        sage: x = F.monomial("a")
        sage: int(2) * x
        2*B['a']
        sage: x.__rmul__(int(2))
        2*B['a']

    .. SEEALSO:: :meth:`Modules.ElementMethods.__mul__`

    This is :meth:`Modules.ElementMethods.__rmul__`, implemented as a Cython
    method in :mod:`sage.categories.coercion_methods`::

        sage: x.__rmul__.im_func is Modules.ElementMethods.__rmul__.im_func
        True
        sage: x.__rmul__.im_func is sage.categories.coercion_methods.Modules__rmul__
        True

    .. TODO::

        Make a better unit test once ``Modules().example()`` is implemented.
    """
    return coercion_model.bin_op(left, right, operator.mul)

@cython.binding
def __mul__(Element self, right):
    r"""
    Return the product of ``self`` and ``right``.

    INPUT:

    - ``self`` -- an element of a :class:`magma <Magmas>`
    - ``right`` -- an object

    This calls the ``_mul_`` method of ``self``, if it is
    available and the two elements have the same parent
    (see :meth:`Magmas.ElementMethods._mul_`).

    Otherwise, the job is delegated to the coercion model.

    Do not override; instead implement a ``_mul_`` method in the
    element class or a ``product`` method in the parent class.

    EXAMPLES::

        sage: S = Semigroups().example("free")
        sage: x = S('a'); y = S('b')
        sage: x * y
        'ab'

    .. SEEALSO::

        - :meth:`Magmas.ElementMethods._mul_`
        - :meth:`Magmas.ElementMethods._mul_parent`
        - :meth:`Magmas.ParentMethods.product`

    This is :meth:`Magmas.ElementMethods.__mul__`, implemented as a
    Cython method in :mod:`sage.categories.coercion_methods`::

        sage: x.__mul__.im_func is Magmas.ElementMethods.__mul__.im_func
        True
        sage: x.__mul__.im_func is sage.categories.coercion_methods.__mul__
        True
    """
    if have_same_parent_c(self, right):
        try:
            return self._mul_(right)
        except AttributeError:
            pass
    return coercion_model.bin_op(self, right, operator.mul)

@cython.binding
def _mul_parent(Element self, other):
    r"""
    Return the product of the two elements, calculated using
    the ``product`` method of the parent.

    This is the default implementation of ``_mul_`` if
    ``product`` is implemented in the parent.

    INPUT:

    - ``other`` -- an element of the parent of ``self``

    OUTPUT:

    - an element of the parent of ``self``

    EXAMPLES::

        sage: S = Semigroups().example("free")
        sage: x = S('a'); y = S('b')
        sage: x._mul_parent(y)
        'ab'

    .. SEEALSO::

        - :meth:`Magmas.ElementMethods._mul_`
        - :meth:`Magmas.ElementMethods._mul_parent`
        - :meth:`Magmas.ParentMethods.product`

    This is :meth:`Magmas.ElementMethods._mul_parent`, implemented as
    a Cython method in :mod:`sage.categories.coercion_methods`::

        sage: x._mul_parent.im_func is Magmas.ElementMethods._mul_parent.im_func
        True
        sage: x._mul_parent.im_func is sage.categories.coercion_methods._mul_parent
        True
    """
    return parent_c(self).product(self, other)


@cython.binding
def __truediv__(left, right):
    """
    Return the result of the division of ``left`` by ``right``, if possible.

    This top-level implementation delegates the work to
    the ``_div_`` method if ``left`` and ``right`` have
    the same parent and to coercion otherwise. See the
    extensive documentation at the top of
    :ref:`sage.structure.element`.

    INPUT:

    - ``self`` -- an element of a :class:`unital magma <Magmas.Unital>`
    - ``right`` -- an object

    .. SEEALSO:: :meth:`Magmas.Unital.ElementMethods._div_`

    EXAMPLES::

        sage: G = FreeGroup(2)
        sage: x0, x1 = G.group_generators()
        sage: c1 = cartesian_product([x0, x1])
        sage: c2 = cartesian_product([x1, x0])
        sage: c1.__div__(c2)
        (x0*x1^-1, x1*x0^-1)
        sage: c1 / c2
        (x0*x1^-1, x1*x0^-1)

    Division supports coercion::

        sage: C = cartesian_product([G, G])
        sage: H = Hom(G, C)
        sage: phi = H(lambda g: cartesian_product([g, g]))
        sage: phi.register_as_coercion()
        sage: x1 / c1
        (x1*x0^-1, 1)
        sage: c1 / x1
        (x0*x1^-1, 1)

    Depending on how the division itself is implemented in
    ``_div_``, division may fail even when ``right``
    actually divides ``left``::

        sage: x = cartesian_product([2, 1])
        sage: y = cartesian_product([1, 1])
        sage: x / y
        (2, 1)
        sage: x / x
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    This is :meth:`Magmas.Unital.ElementMethods.__truediv__`, implemented
    as a Cython method in :mod:`sage.categories.coercion_methods`::

        sage: x.__truediv__.im_func is Magmas.Unital.ElementMethods.__truediv__.im_func
        True
        sage: x.__truediv__.im_func is sage.categories.coercion_methods.__truediv__
        True
    """
    if have_same_parent_c(left, right):
        return left._div_(right)
    return coercion_model.bin_op(left, right, operator.div)

