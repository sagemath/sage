.. _chapter-cython:

================
Coding in Cython
================

This chapter discusses Cython, which is a compiled language based on
Python.  The major advantage it has over Python is that code can be
much faster (sometimes orders of magnitude) and can directly call
C and C++ code.  As Cython is essentially a superset of the Python
language, one often doesn’t make a distinction between Cython and
Python code in Sage (e.g. one talks of the “Sage Python Library”
and “Python Coding Conventions”).

Python is an interpreted language and has no declared data types for
variables. These features make it easy to write and debug, but Python
code can sometimes be slow. Cython code can look a lot like Python,
but it gets translated into C code (often very efficient C code) and
then compiled. Thus it offers a language which is familiar to Python
developers, but with the potential for much greater speed. Cython also
allows Sage developers to interface with C and C++ much easier than
using the Python C API directly.

Cython is a compiled version of Python. It was originally based on
Pyrex but has changed based on what Sage's developers needed; Cython
has been developed in concert with Sage. However, it is an independent
project now, which is used beyond the scope of Sage. As such, it is a
young, but developing language, with young, but developing
documentation. See its web page, http://www.cython.org/, for the most
up-to-date information or check out the
`Language Basics <http://docs.cython.org/src/userguide/language_basics.html>`_
to get started immediately.


Writing cython code in Sage
===========================

There are several ways to create and build Cython code in Sage.

#. In the Sage Notebook, begin any cell with ``%cython``. When you
   evaluate that cell,

   #. It is saved to a file.

   #. Cython is run on it with all the standard Sage libraries
      automatically linked if necessary.

   #. The resulting shared library file (``.so`` / ``.dll`` /
      ``.dylib``) is then loaded into your running instance of Sage.

   #. The functionality defined in that cell is now available for you
      to use in the notebook. Also, the output cell has a link to the C
      program that was compiled to create the ``.so`` file.

   #. A ``cpdef`` or ``def`` function, say ``testfunction``, defined in
      a ``%cython`` cell in a worksheet can be imported and made available
      in a different ``%cython`` cell within the same worksheet by
      importing it as shown below::

          %cython
          from __main__ import testfunction

#. Create an ``.spyx`` file and :ref:`attach or load it <section-attach-or-load-spyx-files>`
   from the :ref:`command line <section-command-line>`.
   This is similar to creating a ``%cython``
   cell in the notebook but works completely from the command line
   (and not from the notebook).

#. Use ``%%cython`` cell magic in the command line.
   Refer to :meth:`sage.repl.ipython_extension.SageMagics.cython`.

#. Create a ``.pyx`` file and add it to the Sage library.
   Then run ``sage -b`` to rebuild Sage.

.. _section-attach-or-load-spyx-files:

Attaching or loading .spyx files
================================

The easiest way to try out Cython without having to learn anything
about distutils, etc., is to create a file with the extension
``spyx``, which stands for "Sage Pyrex":

#. Create a file ``power2.spyx``.

#. Put the following in it:

   .. CODE-BLOCK:: cython

       def is2pow(n):
           while n != 0 and n%2 == 0:
               n = n >> 1
           return n == 1

#. Start the Sage command line interpreter and load the ``spyx`` file
   (this will fail if you do not have a C compiler installed).

   .. skip

   ::

       sage: load("power2.spyx")
       Compiling power2.spyx...
       sage: is2pow(12)
       False

Note that you can change ``power2.spyx``, then load it again and it
will be recompiled on the fly. You can also attach ``power2.spyx`` so
it is reloaded whenever you make changes:

.. skip

::

    sage: attach("power2.spyx")

Cython is used for its speed. Here is a timed test on a 2.6 GHz
Opteron:

.. skip

::

    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 0.60 s, sys: 0.00 s, total: 0.60 s
    Wall time: 0.60 s

Now, the code in the file ``power2.spyx`` is valid Python, and if we
copy this to a file ``powerslow.py`` and load that, we get the
following:

.. skip

::

    sage: load("powerslow.py")
    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 1.01 s, sys: 0.04 s, total: 1.05 s
    Wall time: 1.05 s

By the way, we could gain even a little more speed with the Cython
version with a type declaration, by changing ``def is2pow(n):`` to
``def is2pow(unsigned int n):``.


.. _section-interrupt:

Interrupt and signal handling
=============================

When writing Cython code for Sage, special care must be taken to ensure
that the code can be interrupted with ``CTRL-C``.

Sage uses the `cysignals package <https://github.com/sagemath/cysignals>`_
for this, see the `cysignals documentation <http://cysignals.readthedocs.org/>`_
for more information.

Unpickling Cython code
======================

Pickling for Python classes and extension classes, such as Cython, is different.
This is discussed in the `Python pickling documentation`_. For the unpickling of
extension classes you need to write a :meth:`__reduce__` method which typically
returns a tuple ``(f, args, ...)`` such that ``f(*args)`` returns (a copy of) the
original object. As an example, the following code snippet is the
:meth:`~sage.rings.integer.Integer.__reduce__` method from
:class:`sage.rings.integer.Integer`:

.. CODE-BLOCK:: cython

    def __reduce__(self):
        '''
        This is used when pickling integers.

        EXAMPLES::

            sage: n = 5
            sage: t = n.__reduce__(); t
            (<cyfunction make_integer at ...>, ('5',))
            sage: t[0](*t[1])
            5
            sage: loads(dumps(n)) == n
            True
        '''
        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle Cython extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        return sage.rings.integer.make_integer, (self.str(32),)


.. _python pickling documentation: http://docs.python.org/library/pickle.html#pickle-protocol

Deprecation
===========

When making a **backward-incompatible** modification in Sage, the old code should
keep working and display a message indicating how it should be updated/written
in the future. We call this a *deprecation*.

.. NOTE::

    Deprecated code can only be removed one year after the first
    stable release in which it appeared.

Each deprecation warning contains the number of the GitHub PR that defines
it. We use 666 in the example below.

.. CODE-BLOCK:: cython

      from sage.misc.superseded import deprecation_cython
      deprecation_cython(666, "Do not use your computer to compute 1+1. Use your brain.")

