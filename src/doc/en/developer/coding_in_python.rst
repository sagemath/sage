.. _chapter-python:

=========================
Coding in Python for Sage
=========================

This chapter discusses some issues with, and advice for, coding in
Sage.

.. _section-python-language-standard:

Python language standard
========================

Sage library code needs to be compatible with all versions of Python
that Sage supports.  The information regarding the supported versions
can be found in the files ``build/pkgs/python3/spkg-configure.m4``
(for the Sage distribution), ``src/pyproject.toml`` (for the Sage
library), and ``m4/pyproject_toml_metadata.m4`` (for most other
distribution packages in ``pkgs/``).

Python 3.9 is the oldest supported version.  Hence,
all language and library features that are available in Python 3.9 can
be used; but features introduced in Python 3.10 cannot be used.  If a
feature is deprecated in a newer supported version, it must be ensured
that deprecation warnings issued by Python do not lead to failures in
doctests.

Some key language and library features have been backported to older Python versions
using one of two mechanisms:

- ``from __future__ import annotations`` (see Python reference for
  `__future__ <https://docs.python.org/3/library/__future__.html>`_)
  modernizes type annotations according to `PEP 563
  <https://www.python.org/dev/peps/pep-0563>`_ (Postponed evaluation
  of annotations).  All Sage library code that uses type annotations
  should include this ``__future__`` import and follow PEP 563.

- Backport packages

  - `importlib_metadata <../reference/spkg/importlib_metadata>`_
    (to be used in place of ``importlib.metadata``),
  - `importlib_resources <../reference/spkg/importlib_resources>`_
    (to be used in place of ``importlib.resources``),
  - `typing_extensions <../reference/spkg/typing_extensions>`_
    (to be used in place of ``typing``).

  The Sage library declares these packages as dependencies and ensures that
  versions that provide features of Python 3.12 are available.

Meta :issue:`29756` keeps track of newer Python features and serves
as a starting point for discussions on how to make use of them in the
Sage library.


Design
======

If you are planning to develop some new code for Sage, design is
important. So think about what your program will do and how that fits
into the structure of Sage. In particular, much of Sage is implemented
in the object-oriented language Python, and there is a hierarchy of
classes that organize code and functionality. For example, if you
implement elements of a ring, your class should derive from
``sage.structure.element.RingElement``, rather than starting from
scratch. Try to figure out how your code should fit in with other Sage
code, and design it accordingly.


Special Sage functions
======================

Functions with leading and trailing double underscores ``__XXX__`` are
all predefined by Python. Functions with leading and trailing single
underscores ``_XXX_`` are defined for Sage. Functions with a single
leading underscore are meant to be semi-private, and those with a
double leading underscore are considered really private. Users can
create functions with leading and trailing underscores.

Just as Python has many standard special methods for objects, Sage
also has special methods. They are typically of the form ``_XXX_``.
In a few cases, the trailing underscore is not included, but this will
eventually be changed so that the trailing underscore is always
included. This section describes these special methods.

All objects in Sage should derive from the Cython extension class
``SageObject``:

.. CODE-BLOCK:: python

    from sage.structure.sage_object import SageObject

    class MyClass(SageObject,...):
        ...

or from some other already existing Sage class:

.. CODE-BLOCK:: python

    from sage.structure.parent import Parent

    class MyFavoriteAlgebra(Parent):
        ...

You should implement the ``_latex_`` and ``_repr_`` method for every
object. The other methods depend on the nature of the object.


LaTeX representation
====================

Every object ``x`` in Sage should support the command ``latex(x)``, so
that any Sage object can be easily and accurately displayed via
LaTeX. Here is how to make a class (and therefore its instances)
support the command ``latex``.

#. Define a method ``_latex_(self)`` that returns a LaTeX
   representation of your object. It should be something that can be
   typeset correctly within math mode. Do not include opening and
   closing $'s.

#. Often objects are built up out of other Sage objects, and these
   components should be typeset using the ``latex`` function. For
   example, if ``c`` is a coefficient of your object, and you want to
   typeset ``c`` using LaTeX, use ``latex(c)`` instead of
   ``c._latex_()``, since ``c`` might not have a ``_latex_`` method,
   and ``latex(c)`` knows how to deal with this.

#. Do not forget to include a docstring and an example that
   illustrates LaTeX generation for your object.

#. You can use any macros included in ``amsmath``, ``amssymb``, or
   ``amsfonts``, or the ones defined in :mod:`sage.misc.latex_macros`.

An example template for a ``_latex_`` method follows. Note that the
``.. skip`` line should not be included in your code; it is here to
prevent doctests from running on this fake example.

.. skip

.. CODE-BLOCK:: python

    class X:
       ...
       def _latex_(self):
           r"""
           Return the LaTeX representation of X.

           EXAMPLES::

               sage: a = X(1,2)
               sage: latex(a)
               '\\frac{1}{2}'
           """
           return '\\frac{%s}{%s}'%(latex(self.numer), latex(self.denom))

As shown in the example, ``latex(a)`` will produce LaTeX code
representing the object ``a``. Calling ``view(a)`` will display the
typeset version of this.


Print representation
====================

The standard Python printing method is ``__repr__(self)``. In Sage,
that is for objects that derive from :class:`SageObject` (which is
everything in Sage), instead define ``_repr_(self)``. This is
preferable because if you only define ``_repr_(self)`` and not
``__repr__(self)``, then users can rename your object to print however
they like. Also, some objects should print differently depending on
the context.

Here is an example of the ``_latex_`` and ``_repr_`` functions for the
``Pi`` class. It is from the file
:sage_root:`src/sage/symbolic/constants.py`:

.. CODE-BLOCK:: python

    class Pi(Constant):
        """
        The ratio of a circle's circumference to its diameter.

        EXAMPLES::

            sage: pi
            pi
            sage: float(pi) # rel tol 1e-10
            3.1415926535897931
        """
        ...
        def _repr_(self):
            return "pi"

        def _latex_(self):
            return "\\pi"


Matrix or vector from object
============================

Provide a ``_matrix_`` method for an object that can be coerced to a
matrix over a ring `R`. Then the Sage function ``matrix`` will work
for this object.

The following is from
:sage_root:`src/sage/graphs/generic_graph.py`:

.. CODE-BLOCK:: python

    class GenericGraph(SageObject):
        ...
        def _matrix_(self, R=None):
            if R is None:
                return self.am()
            else:
                return self.am().change_ring(R)


        def adjacency_matrix(self, sparse=None, boundary_first=False):
            ...

Similarly, provide a ``_vector_`` method for an object that can be
coerced to a vector over a ring `R`. Then the Sage function ``vector``
will work for this object. The following is from the file
:sage_root:`src/sage/modules/free_module_element.pyx`:

.. CODE-BLOCK:: python

    cdef class FreeModuleElement(element_Vector):   # abstract base class
        ...
        def _vector_(self, R):
            return self.change_ring(R)


.. _section-preparsing:

Sage preparsing
===============

To make Python even more usable interactively, there are a number of
tweaks to the syntax made when you use Sage from the commandline or
via the notebook (but not for Python code in the Sage
library). Technically, this is implemented by a ``preparse()``
function that rewrites the input string. Most notably, the following
replacements are made:

- Sage supports a special syntax for generating rings or, more
  generally, parents with named generators::

      sage: R.<x,y> = QQ[]
      sage: preparse('R.<x,y> = QQ[]')
      "R = QQ['x, y']; (x, y,) = R._first_ngens(2)"

- Integer and real literals are Sage integers and Sage floating point
  numbers. For example, in pure Python these would be an attribute
  error::

      sage: 16.sqrt()
      4
      sage: 87.factor()
      3 * 29

- Raw literals are not preparsed, which can be useful from an
  efficiency point of view. In Sage raw integer and floating
  literals are followed by an "r" (or "R") for raw, meaning
  not preparsed. For example::

      sage: a = 393939r
      sage: a
      393939
      sage: type(a)
      <... 'int'>
      sage: b = 393939
      sage: type(b)
      <class 'sage.rings.integer.Integer'>
      sage: a == b
      True

- Raw literals can be very useful in certain cases. For instance,
  Python integers can be more efficient than Sage integers when they
  are very small.  Large Sage integers are much more efficient than
  Python integers since they are implemented using the GMP C library.

Consult the file ``preparser.py`` for more details about Sage
preparsing, more examples involving raw literals, etc.

When a file ``foo.sage`` is loaded or attached in a Sage session, a
preparsed version of ``foo.sage`` is created with the name
``foo.sage.py``. The beginning of the preparsed file states::

    This file was *autogenerated* from the file foo.sage.

You can explicitly preparse a file with the ``--preparse``
command-line option: running ::

    sage --preparse foo.sage

creates the file ``foo.sage.py``.

The following files are relevant to preparsing in Sage:

#. :sage_root:`src/bin/sage`

#. :sage_root:`src/bin/sage-preparse`

#. :sage_root:`src/sage/repl/preparse.py`

In particular, the file ``preparse.py`` contains the Sage preparser
code.


The Sage coercion model
=======================

The primary goal of coercion is to be able to transparently do
arithmetic, comparisons, etc. between elements of distinct sets. For
example, when one writes `3 + 1/2`, one wants to perform arithmetic on
the operands as rational numbers, despite the left term being an
integer.  This makes sense given the obvious and natural inclusion of
the integers into the rational numbers. The goal of the coercion
system is to facilitate this (and more complicated arithmetic) without
having to explicitly map everything over into the same domain, and at
the same time being strict enough to not resolve ambiguity or accept
nonsense.

The coercion model for Sage is described in detail, with examples, in
the Coercion section of the Sage Reference Manual.


Mutability
==========

Parent structures (e.g. rings, fields, matrix spaces, etc.) should be
immutable and globally unique whenever possible. Immutability means,
among other things, that properties like generator labels and default
coercion precision cannot be changed.

Global uniqueness while not wasting memory is best implemented using
the standard Python weakref module, a factory function, and module
scope variable.

.. {Rewrite. Difficult to parse. Make gentler}

.. {Put a tutorial on this here}

Certain objects, e.g. matrices, may start out mutable and become
immutable later. See the file
:sage_root:`src/sage/structure/mutability.py`.


The  __hash__ special method
============================

Here is the definition of ``__hash__`` from the Python reference
manual:

    Called by built-in function ``hash()`` and for operations on members
    of hashed collections including ``set``, ``frozenset``, and
    ``dict``. ``__hash__()`` should return an integer. The only required
    property is that objects which compare equal have the same hash
    value; it is advised to mix together the hash values of the
    components of the object that also play a part in comparison of
    objects by packing them into a tuple and hashing the tuple.

    If a class does not define an ``__eq__()`` method it should not define
    a ``__hash__()`` operation either; if it defines ``__eq__()`` but not
    ``__hash__()``, its instances will not be usable as items in hashable
    collections. If a class defines mutable objects and implements an
    ``__eq__()`` method, it should not implement ``__hash__()``, since the
    implementation of hashable collections requires that a key’s hash
    value is immutable (if the object’s hash value changes, it will be
    in the wrong hash bucket).

See https://docs.python.org/3/reference/datamodel.html#object.__hash__ for more
information on the subject.

Notice the phrase, "The only required property is that objects which
compare equal have the same hash value." This is an assumption made by
the Python language, which in Sage we simply cannot make (!), and
violating it has consequences. Fortunately, the consequences are
pretty clearly defined and reasonably easy to understand, so if you
know about them they do not cause you trouble. The following example
illustrates them pretty well:

::

        sage: v = [Mod(2,7)]
        sage: 9 in v
        True
        sage: v = set([Mod(2,7)])
        sage: 9 in v
        False
        sage: 2 in v
        True
        sage: w = {Mod(2,7):'a'}
        sage: w[2]
        'a'
        sage: w[9]
        Traceback (most recent call last):
        ...
        KeyError: 9

Here is another example:

::

        sage: R = RealField(10000)
        sage: a = R(1) + R(10)^-100
        sage: a == RDF(1)  # because the a gets coerced down to RDF
        True

but ``hash(a)`` should not equal ``hash(1)``.

Unfortunately, in Sage we simply cannot require

.. CODE-BLOCK:: text

           (#)   "a == b ==> hash(a) == hash(b)"

because serious mathematics is simply too complicated for this
rule. For example, the equalities ``z == Mod(z, 2)`` and
``z == Mod(z, 3)`` would force ``hash()`` to be constant on the
integers.

The only way we could "fix" this problem for good would be to abandon
using the ``==`` operator for "Sage equality", and implement Sage
equality as a new method attached to each object. Then we could follow
Python rules for ``==`` and our rules for everything else, and all
Sage code would become completely unreadable (and for that matter
unwritable). So we just have to live with it.

So what is done in Sage is to attempt to satisfy ``(#)`` when it is
reasonably easy to do so, but use judgment and not go overboard.
For example,

::

        sage: hash(Mod(2,7))
        2

The output 2 is better than some random hash that also involves the
moduli, but it is of course not right from the Python point of view,
since ``9 == Mod(2,7)``. The goal is to make a hash function that is
fast, but within reason respects any obvious natural inclusions and
coercions.


Exceptions
==========

Please avoid catch-all code like this:

.. CODE-BLOCK:: python

    try:
        some_code()
    except:               # bad
        more_code()

If you do not have any exceptions explicitly listed (as a tuple), your
code will catch absolutely anything, including ``ctrl-C``, typos in
the code, and alarms, and this will lead to confusion. Also, this
might catch real errors which should be propagated to the user.

To summarize, only catch specific exceptions as in the following
example:

.. CODE-BLOCK:: python

    try:
        return self.__coordinate_ring
    except (AttributeError, OtherExceptions) as msg:           # good
        more_code_to_compute_something()

Note that the syntax in ``except`` is to list all the exceptions that
are caught as a tuple, followed by an error message.

A method or a function accepts input described in the ``INPUT`` block of
:ref:`the docstring <section-docstring-function>`. If the input cannot be
handled by the code, then it may raise an exception. The following aims to
guide you in choosing from the most relevant exceptions to Sage. Raise

- :class:`TypeError`: if the input belongs to a class of objects that is not
  supported by the method. For example, a method works only with monic
  polynomials over a finite field, but a polynomial over rationals was given.

- :class:`ValueError`: if the input has a value not supported by the method.
  For example, the above method was given a non-monic polynomial.

- :class:`ArithmeticError`: if the method performs an arithmetic operation
  (sum, product, quotient, and the like) but the input is not appropriate.

- :class:`ZeroDivisionError`: if the method performs division but the input is
  zero. Note that for non-invertible input values, :class:`ArithmeticError` is
  more appropriate. As derived from :class:`ArithmeticError`,
  :class:`ZeroDivisionError` can be caught as :class:`ArithmeticError`.

- :class:`NotImplementedError`: if the input is for a feature not yet
  implemented by the method. Note that this exception is derived from
  :class:`RuntimeError`.

If no specific error seems to apply for your situation, :class:`RuntimeError`
can be used. In all cases, the string associated with the exception should
describe the details of what went wrong.


Integer return values
=====================

Many functions and methods in Sage return integer values.
Those should usually be returned as Sage integers of class
:class:`Integer <sage.rings.integer.Integer>` rather than
as Python integers of class :class:`int`, as users may want
to explore the resulting integers' number-theoretic properties
such as prime factorization. Exceptions should be made when
there are good reasons such as performance or compatibility
with Python code, for instance in methods such as
``__hash__``, ``__len__``, and ``__int__``.

To return a Python integer ``i`` as a Sage integer, use:

.. CODE-BLOCK:: python

    from sage.rings.integer import Integer
    return Integer(i)

To return a Sage integer ``i`` as a Python ineger, use:

.. CODE-BLOCK:: python

    return int(i)


Importing
=========

We mention two issues with importing: circular imports and importing
large third-party modules. See also :ref:`section_dependencies_distributions`
for a discussion of imports from the viewpoint of modularization.

First, you must avoid circular imports. For example, suppose that the
file :sage_root:`src/sage/algebras/steenrod_algebra.py`
started with a line:

.. CODE-BLOCK:: python

    from sage.sage.algebras.steenrod_algebra_bases import *

and that the file
:sage_root:`src/sage/algebras/steenrod_algebra_bases.py`
started with a line:

.. CODE-BLOCK:: python

    from sage.sage.algebras.steenrod_algebra import SteenrodAlgebra

This sets up a loop: loading one of these files requires the other,
which then requires the first, etc.

With this set-up, running Sage will produce an error:

.. CODE-BLOCK:: text

    Exception exceptions.ImportError: 'cannot import name SteenrodAlgebra'
    in 'sage.rings.polynomial.polynomial_element.
    Polynomial_generic_dense.__normalize' ignored
    -------------------------------------------------------------------
    ImportError                       Traceback (most recent call last)

    ...
    ImportError: cannot import name SteenrodAlgebra

Instead, you might replace the ``import *`` line at the top of the
file by more specific imports where they are needed in the code. For
example, the ``basis`` method for the class ``SteenrodAlgebra`` might
look like this (omitting the documentation string):

.. CODE-BLOCK:: python

    def basis(self, n):
        from steenrod_algebra_bases import steenrod_algebra_basis
        return steenrod_algebra_basis(n, basis=self._basis_name, p=self.prime)

Second, do not import at the top level of your module a third-party
module that will take a long time to initialize (e.g. :mod:`matplotlib`). As
above, you might instead import specific components of the module when
they are needed, rather than at the top level of your file.

It is important to try to make ``from sage.all import *`` as fast as
possible, since this is what dominates the Sage startup time, and
controlling the top-level imports helps to do this. One important
mechanism in Sage are lazy imports, which don't actually perform the
import but delay it until the object is actually used. See
:mod:`sage.misc.lazy_import` for more details of lazy imports, and
:ref:`chapter-directory-structure` for an example using lazy imports
for a new module.

If your module needs to make some precomputed data available at the top level,
you can reduce its load time (and thus startup time, unless your module is
imported using :mod:`sage.misc.lazy_import`) by using the decorator
:func:`sage.misc.cachefunc.cached_function` instead. For example, replace

.. CODE-BLOCK:: python

    big_data = initialize_big_data()  # bad: runs at module load time

by

.. CODE-BLOCK:: python

    from sage.misc.cachefunc import cached_function

    @cached_function                  # good: runs on first use
    def big_data():
        return initialize_big_data()


Static typing
=============

Python libraries are increasingly annotated with static typing information;
see the `Python reference on typing <https://docs.python.org/3/library/typing.html>`_.

For typechecking the Sage library, the project uses :ref:`pyright <section-tools-pyright>`;
it automatically runs in the GitHub Actions CI and can also be run locally.

As of Sage 10.2, the Sage library only contains a minimal set of such type
annotations. Pull requests that add more annotations are generally welcome.

The Sage library makes very extensive use of Cython (see chapter :ref:`chapter-cython`).
Although Cython source code often declares static types for the purpose of
compilation to efficient machine code, this typing information is unfortunately
not visible to static checkers such as Pyright. It is necessary to create `type stub
files (".pyi") <https://github.com/microsoft/pyright/blob/main/docs/type-stubs.md>`_
that provide this information. Although various
`tools for writing and maintaining type stub files
<https://typing.readthedocs.io/en/latest/source/writing_stubs.html#writing-and-maintaining-stub-files>`_
are available, creating stub files for Cython files involves manual work.
There is hope that better tools become available soon, see for example
`cython/cython #5744 <https://github.com/cython/cython/pull/5744>`_.
Contributing to the development and testing of such tools likely will have a
greater impact than writing the typestub files manually.

For Cython modules of the Sage library, these type stub files would be placed
next to the ``.pyx`` and ``.pxd`` files.

When importing from other Python libraries that do not provide sufficient typing
information, it is possible to augment the library's typing information for
the purposes of typechecking the Sage library:

- Create typestub files and place them in the directory :file:`SAGE_ROOT/src/typings`.
  For example, the distribution **pplpy** provides the top-level package :mod:`ppl`,
  which publishes no typing information. We can create a typestub file
  :file:`SAGE_ROOT/src/typings/ppl.pyi` or :file:`SAGE_ROOT/src/typings/ppl/__init__.pyi`.

- When these typestub files are working well, it is preferable from the viewpoint
  of the Sage project that they are "upstreamed", i.e., contributed to the
  project that maintains the library. If a new version of the upstream library
  becomes available that provides the necessary typing information, we can
  update the package in the Sage distribution and remove the typestub files again
  from :file:`SAGE_ROOT/src/typings`.

- As a fallback, when neither adding typing annotations to source files
  nor adding typestub files is welcomed by the upstream project, it is possible
  to `contribute typestubs files instead to the typeshed community project
  <https://github.com/python/typeshed/blob/main/CONTRIBUTING.md>`_.


Deprecation
===========

When making a **backward-incompatible** modification in Sage, the old code should
keep working and a message indicating how the code should be updated/written
in the future should be displayed somewhere. We call this *deprecation*. We explain
how to do the deprecation, the deprecation policy, below.

Any class, function, method, or attribute defined in a file under
:sage_root:`src/sage` is subject to the deprecation policy. If its name starts
with an underscore, then it is considered internal, and exempt from the
deprecation policy.

.. NOTE::

    A deprecated class, function, method, or attribute can only be removed one
    year after the first stable release in which it appeared.

When a deprecated function, method, or attribute is used, a deprecation warning
is issued. The warning message contains the number of the GitHub PR that
implemented the deprecation. We use 12345 in the examples below.

.. NOTE::

    For deprecation tools used in the examples, consult the tool's documentation for more
    information on its behaviour and optional arguments.

* **Rename a keyword:** by decorating a function/method with
  :class:`~sage.misc.decorators.rename_keyword`, any user calling
  ``my_function(my_old_keyword=5)`` will see a warning:

  .. CODE-BLOCK:: python

      from sage.misc.decorators import rename_keyword
      @rename_keyword(deprecation=12345, my_old_keyword='my_new_keyword')
      def my_function(my_new_keyword=True):
          return my_new_keyword

* **Rename a function/method:** call
  :func:`~sage.misc.superseded.deprecated_function_alias` to obtain a copy of a
  function that raises a deprecation warning:

  .. CODE-BLOCK:: python

      from sage.misc.superseded import deprecated_function_alias
      def my_new_function():
          ...

      my_old_function = deprecated_function_alias(12345, my_new_function)

* **Moving an object to a different module:**
  if you rename a source file or move some function (or class) to a
  different file, it should still be possible to import that function
  from the old module. This can be done using a
  :func:`~sage.misc.lazy_import.lazy_import` with deprecation.
  In the old module, you would write:

  .. CODE-BLOCK:: python

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.new.module.name', 'name_of_the_function', deprecation=12345)

  You can also lazily import everything using ``*`` or a few functions
  using a tuple:

  .. CODE-BLOCK:: python

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.new.module.name', '*', deprecation=12345)
    lazy_import('sage.other.module', ('func1', 'func2'), deprecation=12345)

* **Remove a name from a global namespace:** this is when you want to
  remove a name from a global namespace (say, ``sage.all`` or some
  other ``all.py`` file) but you want to keep the functionality
  available with an explicit import.
  This case is similar as the previous one: use a lazy import with
  deprecation. One detail: in this case, you don't want the name
  ``lazy_import`` to be visible in the global namespace, so we add
  a leading underscore:

  .. CODE-BLOCK:: python

    from sage.misc.lazy_import import lazy_import as _lazy_import
    _lazy_import('sage.some.package', 'some_function', deprecation=12345)

* **Any other case:** if none of the cases above apply, call
  :func:`~sage.misc.superseded.deprecation` in the function that you want to
  deprecate. It will display the message of your choice (and interact properly
  with the doctest framework):

  .. CODE-BLOCK:: python

      from sage.misc.superseded import deprecation
      deprecation(12345, "Do not use your computer to compute 1 + 1. Use your brain.")

.. NOTE::

    These decorators only work for Python. There is no implementation
    of decorators in Cython. Hence, when in need to rename a keyword/function/method/...
    in a Cython (.pyx) file and/or to deprecate something, forget about decorators and
    just use :func:`~sage.misc.superseded.deprecation_cython` instead. The usage of
    :func:`~sage.misc.superseded.deprecation_cython` is exactly the same as
    :func:`~sage.misc.superseded.deprecation`.

When a class is renamed or removed, it should be deprecated unless it is
internal. A class is internal if its name starts with an underscore, or experts
(authors and reviewers of a PR making changes to the class) agree that the
class is unlikely to be directly imported by user code. Otherwise, or if experts
disagree, it is public.

As a class is imported rather than run by user code, there are some technical
difficulties in using the above deprecation tools. Instead we follow
the procedure below:

* **Renaming a class:** rename ``OldClass`` to ``NewClass`` and add an
  alias ``OldClass = NewClass``:

  .. CODE-BLOCK:: python

    class NewClass:
        ...

    OldClass = NewClass   # OldClass is deprecated. See Issue 12345.

* **Removing a class:**  add a comment:

  .. CODE-BLOCK:: python

    # OldClass is deprecated. See Issue 12345.

    class OldClass:

In both cases, make it sure to display the change in the "Deprecations"
section of the release notes of the next stable release.


Experimental/unstable code
==========================

You can mark your newly created code (classes/functions/methods) as
experimental/unstable. In this case, no deprecation warning is needed
when changing this code, its functionality or its interface.

This should allow you to put your stuff in Sage early, without worrying about
making (design) changes later.

When satisfied with the code (when stable for some time, say, one
year), you can delete this warning.

As usual, all code has to be fully doctested and go through our
reviewing process.

* **Experimental function/method:** use the decorator
  :class:`~sage.misc.superseded.experimental`. Here is an example:

  .. CODE-BLOCK:: python

      from sage.misc.superseded import experimental
      @experimental(12345)
      def experimental_function():
          # do something

* **Experimental class:** use the decorator
  :class:`~sage.misc.superseded.experimental` for its ``__init__``.
  Here is an example:

  .. CODE-BLOCK:: python

      from sage.misc.superseded import experimental
      class experimental_class(SageObject):
          @experimental(12345)
          def __init__(self, some, arguments):
              # do something

* **Any other case:** if none of the cases above apply, call
  :func:`~sage.misc.superseded.experimental_warning` in the code where
  you want to warn. It will display the message of your choice:

  .. CODE-BLOCK:: python

      from sage.misc.superseded import experimental_warning
      experimental_warning(12345, 'This code is not foolproof.')


Using optional packages
=======================

If a function requires an optional package, that function should fail
gracefully---perhaps using a ``try``-``except`` block---when the
optional package is not available, and should give a hint about how to
install it. For example, typing ``sage -optional`` gives a list of all
optional packages, so it might suggest to the user that they type
that. The command ``optional_packages()`` from within Sage also
returns this list.
