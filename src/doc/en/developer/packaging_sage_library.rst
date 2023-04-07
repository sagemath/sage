
.. _chapter-modularization:

============================
 Packaging the Sage Library
============================


Modules, packages, distribution packages
========================================

The Sage library consists of a large number of Python modules,
organized into a hierarchical set of packages that fill the namespace
:mod:`sage`.  All source files are located in a subdirectory of the
directory ``SAGE_ROOT/src/sage/``.

For example,

- the file ``SAGE_ROOT/src/sage/coding/code_bounds.py`` provides the
  module :mod:`sage.coding.code_bounds`;

- the directory containing this file, ``SAGE_ROOT/src/sage/coding/``,
  thus provides the package :mod:`sage.coding`.

There is another notion of "package" in Python, the **distribution
package** (also known as a "distribution" or a "pip-installable
package").  Currently, the entire Sage library is provided by a
single distribution,
`sagemath-standard <https://pypi.org/project/sagemath-standard/>`_,
which is generated from the directory
``SAGE_ROOT/pkgs/sagemath-standard``.

Note that the distribution name is not required to be a Python
identifier. In fact, using dashes (``-``) is preferred to underscores in
distribution names; **setuptools** and other parts of Python's packaging
infrastructure normalize underscores to dashes. (Using dots in
distribution names, to indicate ownership by organizations, still
mentioned in `PEP 423 <https://www.python.org/dev/peps/pep-0423/>`_, appears to
have largely fallen out of favor, and we will not use it in the SageMath
project.)

A distribution that provides Python modules in the :mod:`sage.*` namespace, say
mainly from :mod:`sage.PAC.KAGE`, should be named **sagemath-DISTRI-BUTION**.
Example:

- The distribution
  `sagemath-categories <https://pypi.org/project/sagemath-categories/>`_
  provides a small subset of the modules of the Sage library, mostly
  from the packages :mod:`sage.structure`, :mod:`sage.categories`, and
  :mod:`sage.misc`.

Other distributions should not use the prefix **sagemath-** in the
distribution name. Example:

- The distribution `sage-sws2rst <https://pypi.org/project/sage-sws2rst/>`_
  provides the Python package :mod:`sage_sws2rst`, so it does not fill
  the :mod:`sage.*` namespace and therefore does not use the prefix
  **sagemath-**.

A distribution that provides functionality that does not need to
import anything from the :mod:`sage` namespace should not use the
:mod:`sage` namespace for its own packages/modules. It should be
positioned as part of the general Python ecosystem instead of as a
Sage-specific distribution.  Examples:

- The distribution `pplpy <https://pypi.org/project/pplpy/>`_ provides the Python
  package :mod:`ppl` and is a much extended version of what used to be
  :mod:`sage.libs.ppl`, a part of the Sage library. The package :mod:`sage.libs.ppl` had
  dependencies on :mod:`sage.rings` to convert to/from Sage number
  types. **pplpy** has no such dependencies and is therefore usable in a
  wider range of Python projects.

- The distribution `memory-allocator <https://pypi.org/project/memory-allocator/>`_
  provides the Python package :mod:`memory_allocator`. This used to be
  :mod:`sage.ext.memory_allocator`, a part of the Sage library.


.. _section_namespace_packages:

Ordinary packages vs. implicit namespace packages
-------------------------------------------------

Each module of the Sage library must be packaged in exactly one distribution
package. However, modules in a package may be included in different
distribution packages. In this regard, there is an important constraint that an
ordinary package (directory with ``__init__.py`` file) cannot be split into
more than one distribution package.

By removing the ``__init__.py`` file, however, we can make the package an
"implicit" (or "native") "namespace" package, following
`PEP 420 <https://www.python.org/dev/peps/pep-0420/>`_. Implicit namespace packages can be
included in more than one distribution package. Hence whenever there are two
distribution packages that provide modules with a common prefix of Python
packages, that prefix needs to be a implicit namespace package, i.e., there
cannot be an ``__init__.py`` file.

For example,

- **sagemath-tdlib** will provide :mod:`sage.graphs.graph_decompositions.tdlib`,

- **sagemath-rw** will provide :mod:`sage.graphs.graph_decompositions.rankwidth`,

- **sagemath-graphs** will provide all of the rest of
  :mod:`sage.graphs.graph_decompositions` (and most of :mod:`sage.graphs`).

Then, none of

- :mod:`sage`,

- :mod:`sage.graphs`,

- :mod:`sage.graphs.graph_decomposition`

can be an ordinary package (with an ``__init__.py`` file), but rather
each of them has to be an implicit namespace package (no
``__init__.py`` file).

For an implicit namespace package, ``__init__.py`` cannot be used any more for
initializing the package.

In the Sage 9.6 development cycle, we still use ordinary packages by
default, but several packages are converted to implicit namespace
packages to support modularization.


Source directories of distribution packages
===========================================

The development of the Sage library uses a monorepo strategy for
all distribution packages that fill the :mod:`sage.*` namespace.  This
means that the source trees of these distributions are included in a
single ``git`` repository, in a subdirectory of ``SAGE_ROOT/pkgs``.

All these distribution packages have matching version numbers.  From
the viewpoint of a single distribution, this means that sometimes
there will be a new release of some distribution where the only thing
changing is the version number.

The source directory of a distribution package, such as
``SAGE_ROOT/pkgs/sagemath-standard``, contains the following files:

- ``sage`` -- a relative symbolic link to the monolithic Sage library
  source tree ``SAGE_ROOT/src/sage/``

- `MANIFEST.in <https://packaging.python.org/guides/using-manifest-in/>`_ --
  controls which files and directories of the
  monolithic Sage library source tree are included in the distribution

- `pyproject.toml <https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/>`_,
  `setup.cfg <https://setuptools.pypa.io/en/latest/userguide/declarative_config.html>`_,
  and `requirements.txt <https://pip.pypa.io/en/stable/user_guide/#requirements-files>`_ --
  standard Python packaging metadata, declaring the distribution name, dependencies,
  etc.

- ``README.rst`` -- a description of the distribution

- ``LICENSE.txt`` -- relative symbolic link to the same files
  in ``SAGE_ROOT/src``

- ``VERSION.txt`` -- package version. This file is updated by the release manager by
  running the ``sage-update-version`` script.

  Sometimes it may be necessary to upload a hotfix for a distribution
  package to PyPI. These should be marked by adding a suffix
  ``.post1``, ``.post2``; see `PEP 440 on post-releases
  <https://peps.python.org/pep-0440/#post-releases>`_. For example, if
  the current development release is ``9.7.beta8``, then such a
  version could be marked ``9.7.beta8.post1``.

  Also sometimes when working on tickets it may be necessary to
  increment the version because a new feature is needed in another
  distribution package. Such versions should be marked by using the
  version number of the anticipated next development release and
  adding a suffix ``.dev1``, ``.dev2`` ...  (see `PEP 440 on
  developmental releases
  <https://peps.python.org/pep-0440/#developmental-releases>`_).
  For example, if the current development release is ``9.7.beta8``,
  use ``9.7.beta9.dev1``. If the current development release is
  the stable release ``9.8``, use ``9.9.beta0.dev1``.

  After the ticket is merged in the next development version, it will
  be synchronized again with the other package versions.

- ``setup.py`` -- a `setuptools <https://pypi.org/project/setuptools/>`_-based
  installation script

- ``tox.ini`` -- configuration for testing with `tox <https://pypi.org/project/tox/>`_

The technique of using symbolic links pointing into ``SAGE_ROOT/src``
has allowed the modularization effort to keep the ``SAGE_ROOT/src``
tree monolithic: Modularization has been happening behind the scenes
and will not change where Sage developers find the source files.
When adding a new distribution package that uses a symbolic link pointing into
``SAGE_ROOT/src``, please update ``search.exclude`` in
``SAGE_ROOT/.vscode/settings.json``.

Some of these files may actually be generated from source files with suffix ``.m4`` by the
``SAGE_ROOT/bootstrap`` script via the ``m4`` macro processor.

For every distribution package, there is also a subdirectory of ``SAGE_ROOT/build/pkgs/``,
which contains the build infrastructure that is specific to Sage-the-distribution.
Note that these subdirectories follows a different naming convention,
using underscores instead of dashes, see :ref:`section-directory-structure`.
Because the distribution packages are included in the source tree, we set them
up as "script packages" instead of "normal packages", see :ref:`section-package-source-types`.


.. _section_dependencies_distributions:

Dependencies and distribution packages
======================================

When preparing a portion of the Sage library as a distribution
package, dependencies matter.


Build-time dependencies
-----------------------

If the portion of the library contains any Cython modules, these
modules are compiled during the wheel-building phase of the
distribution package. If the Cython module uses ``cimport`` to pull in
anything from ``.pxd`` files, these files must be either part of the
portion shipped as the distribution being built, or the distribution
that provides these files must be installed in the build
environment. Also, any C/C++ libraries that the Cython module uses
must be accessible from the build environment.

*Declaring build-time dependencies:* Modern Python packaging provides a
mechanism to declare build-time dependencies on other distribution
packages via the file `pyproject.toml <https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/>`_
(``[build-system] requires``); this
has superseded the older ``setup_requires`` declaration. (There is no
mechanism to declare anything regarding the C/C++ libraries.)

While the namespace :mod:`sage.*` is organized roughly according to
mathematical fields or categories, how we partition the implementation
modules into distribution packages has to respect the hard constraints
that are imposed by the build-time dependencies.

We can define some meaningful small distributions that just consist of
a single or a few Cython modules. For example, **sagemath-tdlib**
(:trac:`29864`) would just package the single
Cython module that must be linked with ``tdlib``,
:mod:`sage.graphs.graph_decompositions.tdlib`. Starting with the Sage
9.6 development cycle, as soon as namespace packages are activated, we
can start to create these distributions. This is quite a mechanical
task.

*Reducing build-time dependencies:* Sometimes it is possible to
replace build-time dependencies of a Cython module on a library by a
runtime dependency.  In other cases, it may be possible to split a
module that simultaneously depends on several libraries into smaller
modules, each of which has narrower dependencies.


Module-level runtime dependencies
---------------------------------

Any ``import`` statements at the top level of a Python or Cython
module are executed when the module is imported. Hence, the imported
modules must be part of the distribution, or provided by another
distribution -- which then must be declared as a run-time dependency.

*Declaring run-time dependencies:* These dependencies are declared in
``setup.cfg`` (generated from ``setup.cfg.m4``) as
`install_requires <https://setuptools.pypa.io/en/latest/userguide/dependency_management.html#declaring-required-dependency>`_.

*Reducing module-level run-time dependencies:*

- Avoid importing from :mod:`sage.PAC.KAGE.all` modules when :mod:`sage.PAC.KAGE` is
  a namespace package. The main purpose of the :mod:`*.all` modules is for
  populating the global interactive environment that is available to users at
  the ``sage:`` prompt. In particular, no Sage library code should import from
  :mod:`sage.rings.all`.

  To audit the Sage library for such imports, use ``sage --tox -e relint``.
  In most cases, the imports can be fixed automatically using the
  tool ``sage --fiximports``.

- Replace module-level imports by method-level imports.  Note that
  this comes with a small runtime overhead, which can become
  noticeable if the method is called in tight inner loops.

- Sage provides the :func:`~sage.misc.lazy_import.lazy_import`
  mechanism. Lazy imports can be
  declared at the module level, but the actual importing is only done
  on demand. It is a runtime error at that time if the imported module
  is not present. This can be convenient compared to local imports in
  methods when the same imports are needed in several methods.

- Avoid the "modularization anti-pattern" of importing a class from
  another module just to run an ``isinstance(object, Class)`` test, in
  particular when the module implementing ``Class`` has heavy
  dependencies.  For example, importing the class
  :class:`~sage.rings.padics.generic_nodes.pAdicField` (or the
  function :class:`~sage.rings.padics.generic_nodes.is_pAdicField`)
  requires the libraries NTL and PARI.

  Instead, provide an abstract base class (ABC) in a module that only
  has light dependencies, make ``Class`` a subclass of ``ABC``, and
  use ``isinstance(object, ABC)``. For example, :mod:`sage.rings.abc`
  provides abstract base classes for many ring (parent) classes,
  including :class:`sage.rings.abc.pAdicField`.  So we can replace::

    from sage.rings.padics.generic_nodes import pAdicFieldGeneric  # heavy dependencies
    isinstance(object, pAdicFieldGeneric)

  and::

    from sage.rings.padics.generic_nodes import is_pAdicField      # heavy dependencies
    is_pAdicField(object)                                          # deprecated

  by::

    import sage.rings.abc                                          # no dependencies
    isinstance(object, sage.rings.abc.pAdicField)

  Note that going through the abstract base class only incurs a small
  performance penalty::

    sage: object = Qp(5)

    sage: from sage.rings.padics.generic_nodes import pAdicFieldGeneric
    sage: %timeit isinstance(object, pAdicFieldGeneric)            # fast                           # not tested
    68.7 ns ± 2.29 ns per loop (...)

    sage: import sage.rings.abc
    sage: %timeit isinstance(object, sage.rings.abc.pAdicField)    # also fast                      # not tested
    122 ns ± 1.9 ns per loop (...)

- If it is not possible or desired to create an abstract base class for
  ``isinstance`` testing (for example, when the class is defined in some
  external package), other solutions need to be used.

  Note that Python caches successful module imports, but repeating an
  unsuccessful module import incurs a cost every time::

    sage: from sage.schemes.generic.scheme import Scheme
    sage: sZZ = Scheme(ZZ)

    sage: def is_Scheme_or_Pluffe(x):
    ....:    if isinstance(x, Scheme):
    ....:        return True
    ....:    try:
    ....:        from xxxx_does_not_exist import Pluffe            # slow on every call
    ....:    except ImportError:
    ....:        return False
    ....:    return isinstance(x, Pluffe)

    sage: %timeit is_Scheme_or_Pluffe(sZZ)                         # fast                           # not tested
    111 ns ± 1.15 ns per loop (...)

    sage: %timeit is_Scheme_or_Pluffe(ZZ)                          # slow                           # not tested
    143 µs ± 2.58 µs per loop (...)

  The :func:`~sage.misc.lazy_import.lazy_import` mechanism can be used to simplify
  this pattern via the :meth:`~sage.misc.lazy_import.LazyImport.__instancecheck__`
  method and has similar performance characteristics::

    sage: lazy_import('xxxx_does_not_exist', 'Pluffe')

    sage: %timeit isinstance(sZZ, (Scheme, Pluffe))                # fast                           # not tested
    95.2 ns ± 0.636 ns per loop (...)

    sage: %timeit isinstance(ZZ, (Scheme, Pluffe))                 # slow                           # not tested
    158 µs ± 654 ns per loop (...)

  It is faster to do the import only once, for example when loading the module,
  and to cache the failure.  We can use the following idiom, which makes
  use of the fact that ``isinstance`` accepts arbitrarily nested lists
  and tuples of types::

    sage: try:
    ....:     from xxxx_does_not_exist import Pluffe               # runs once
    ....: except ImportError:
    ....:     # Set to empty tuple of types for isinstance
    ....:     Pluffe = ()

    sage: %timeit isinstance(sZZ, (Scheme, Pluffe))                # fast                           # not tested
    95.9 ns ± 1.52 ns per loop (...)

    sage: %timeit isinstance(ZZ, (Scheme, Pluffe))                 # fast                           # not tested
    126 ns ± 1.9 ns per loop (...)


Other runtime dependencies
--------------------------

If ``import`` statements are used within a method, the imported module
is loaded the first time that the method is called. Hence the module
defining the method can still be imported even if the module needed by
the method is not present.

It is then a question whether a run-time dependency should be
declared. If the method needing that import provides core
functionality, then probably yes. But if it only provides what can be
considered "optional functionality", then probably not, and in this
case it will be up to the user to install the distribution enabling
this optional functionality.

As an example, let us consider designing a distribution that centers
around the package :mod:`sage.coding`. First, let's see if it uses symbolics::

  (9.5.beta6) $ git grep -E 'sage[.](symbolic|functions|calculus)' src/sage/coding
  src/sage/coding/code_bounds.py:        from sage.functions.other import ceil
  ...
  src/sage/coding/grs_code.py:from sage.symbolic.ring import SR
  ...
  src/sage/coding/guruswami_sudan/utils.py:from sage.functions.other import floor

Apparently it does not in a very substantial way:

- The imports of the symbolic functions :func:`~sage.functions.other.ceil`
  and :func:`~sage.functions.other.floor` can
  likely be replaced by the artithmetic functions
  :func:`~sage.arith.misc.integer_floor` and
  :func:`~sage.arith.misc.integer_ceil`.

- Looking at the import of ``SR`` by :mod:`sage.coding.grs_code`, it
  seems that ``SR`` is used for running some symbolic sum, but the
  doctests do not show symbolic results, so it is likely that this can
  be replaced.

- Note though that the above textual search for the module names is
  merely a heuristic. Looking at the source of "entropy", through
  ``log`` from :mod:`sage.misc.functional`, a runtime dependency on
  symbolics comes in. In fact, for this reason, two doctests there are
  already marked as ``# optional - sage.symbolic``.

So if packaged as **sagemath-coding**, now a domain expert would have
to decide whether these dependencies on symbolics are strong enough to
declare a runtime dependency (``install_requires``) on
**sagemath-symbolics**. This declaration would mean that any user who
installs **sagemath-coding** (``pip install sagemath-coding``) would
pull in **sagemath-symbolics**, which has heavy compile-time
dependencies (ECL/Maxima/FLINT/Singular/...).

The alternative is to consider the use of symbolics by
**sagemath-coding** merely as something that provides some extra
features, which will only be working if the user also has installed
**sagemath-symbolics**.

*Declaring optional run-time dependencies:* It is possible to declare
such optional dependencies as `extras_require <https://setuptools.pypa.io/en/latest/userguide/dependency_management.html#optional-dependencies>`_ in ``setup.cfg``
(generated from ``setup.cfg.m4``).  This is a very limited mechanism
-- in particular it does not affect the build phase of the
distribution in any way. It basically only provides a way to give a
nickname to a distribution that can be installed as an add-on.

In our example, we could declare an ``extras_require`` so that users
could use ``pip install sagemath-coding[symbolics]``.


Doctest-only dependencies
-------------------------

Doctests often use examples constructed using functionality provided
by other portions of the Sage library.  This kind of integration
testing is one of the strengths of Sage; but it also creates extra
dependencies.

Fortunately, these dependencies are very mild, and we can deal with
them using the same mechanism that we use for making doctests
conditional on the presence of optional libraries: using ``# optional -
FEATURE`` directives in the doctests.  Adding these directives will
allow developers to test the distribution separately, without
requiring all of Sage to be present.

*Declaring doctest-only dependencies:* The
`extras_require <https://setuptools.pypa.io/en/latest/userguide/dependency_management.html#optional-dependencies>`_
mechanism mentioned above can also be used for this.


Version constraints of dependencies
-----------------------------------

The version information for dependencies comes from the files
``build/pkgs/*/install-requires.txt`` and
``build/pkgs/*/package-version.txt``.  We use the
`m4 <https://www.gnu.org/software/m4/manual/html_node/index.html>`_
macro processor to insert the version information in the generated files
``pyproject.toml``, ``setup.cfg``, ``requirements.txt``.


Hierarchy of distribution packages
==================================

.. PLOT::

    def node(label, pos):
        return text(label, (3*pos[0],2*pos[1]), background_color='pink', color='black')
    def edge(start, end, **kwds):
        return arrow((3*start[0],2*start[1]),(3*end[0],2*end[1]-.28), arrowsize=2, **kwds)
    def extras_require(start, end):
        return edge(start, end, linestyle='dashed')
    g = Graphics()
    g += (extras_require((0.5,0),(0.5,1)) + node("sage_conf", (0.5,0)))
    g += (edge((1.5,0),(0.75,2)) + edge((1.5,0),(1.5,1))
          + node("sagemath-objects", (1.5,0)))
    g += (edge((0.5,1),(0,2)) + edge((0.5,1),(0.6,2)) + edge((0.5,1),(1.25,2)) + edge((0.5,1),(1.8,2))
          + node("sagemath-environment", (0.5,1)))
    g += (edge((1.5,1),(0.2,2)) + edge((1.5,1),(1.41,2)) + edge((1.5,1),(2,2))
          + node("sagemath-categories", (1.5,1)))
    g += (edge((0,2),(0,3)) + edge((0,2),(0.75,3)) + edge((0.67,2),(1,3)) + edge((1.33,2),(1.25,3)) + edge((2,2),(2,3))
          + node("sagemath-graphs", (0,2)) + node("sagemath-repl", (0.67,2)) + node("sagemath-polyhedra", (1.33,2)) + node("sagemath-singular", (2,2)))
    g += (edge((1,3),(1,4)) + edge((2,3),(1.2,4))
          + node("sagemath-tdlib", (0,3)) + node("sagemath-standard-no-symbolics", (1,3)) + node("sagemath-symbolics", (2,3)))
    g += node("sagemath-standard", (1,4))
    sphinx_plot(g, figsize=(8, 4), axes=False)


Solid arrows indicate ``install_requires``, i.e., a declared runtime dependency.
Dashed arrows indicate ``extras_require``, i.e., a declared optional runtime dependency.
Not shown in the diagram are build dependencies and optional dependencies for testing.

- `sage_conf <https://pypi.org/project/sage-conf/>`_ is a configuration
  module. It provides the configuration variable settings determined by the
  ``configure`` script.

- `sagemath-environment <https://pypi.org/project/sagemath-environment/>`_
  provides the connection to the system and software environment. It includes
  :mod:`sage.env`, :mod:`sage.features`, :mod:`sage.misc.package_dir`, etc.

- `sagemath-objects <https://pypi.org/project/sagemath-objects/>`_
  provides a small fundamental subset of the modules of the Sage library,
  in particular all of :mod:`sage.structure`, a small portion of :mod:`sage.categories`,
  and a portion of :mod:`sage.misc`.

- `sagemath-categories <https://pypi.org/project/sagemath-categories/>`_
  provides a small subset of the modules of the Sage library, building upon sagemath-objects.
  It provides all of :mod:`sage.categories` and a small portion of :mod:`sage.rings`.

- `sagemath-repl <https://pypi.org/project/sagemath-repl/>`_ provides
  the IPython kernel and Sage preparser (:mod:`sage.repl`),
  the Sage doctester (:mod:`sage.doctest`), and some related modules from :mod:`sage.misc`.


Testing distribution packages
=============================

Of course, we need tools for testing modularized distributions of
portions of the Sage library.

- Modularized distributions must be testable separately!

- But we want to keep integration testing with other portions of Sage too!

Preparing doctests
------------------

Whenever an optional package is needed for a particular test, we use the
doctest annotation ``# optional``. This mechanism can also be used for making a
doctest conditional on the presence of a portion of the Sage library.

The available tags take the form of package or module names such as
:mod:`sage.combinat`, :mod:`sage.graphs`, :mod:`sage.plot`, :mod:`sage.rings.number_field`,
:mod:`sage.rings.real_double`, and :mod:`sage.symbolic`.  They are defined via
:class:`~sage.features.Feature` subclasses in the module :mod:`sage.features.sagemath`, which
also provides the mapping from features to the distributions providing them
(actually, to SPKG names).  Using this mapping, Sage can issue installation
hints to the user.

For example, the package :mod:`sage.tensor` is purely algebraic and has
no dependency on symbolics. However, there are a small number of
doctests that depend on :class:`sage.symbolic.ring.SymbolicRing` for integration
testing. Hence, these doctests are marked ``# optional -
sage.symbolic``.

When defining new features for the purpose of doctest annotations, it may be a good
idea to hide implementation details from feature names. For example, all doctests that
use finite fields have to depend on PARI. However, we have defined a feature
:mod:`sage.rings.finite_rings` (which implies the presence of :mod:`sage.libs.pari`).
Annotating the doctests ``# optional - sage.rings.finite_rings`` expresses the
dependency in a clearer way than using ``# optional - sage.libs.pari``, and it
will be a smaller maintenance burden when implementation details change.


Testing the distribution in virtual environments with tox
---------------------------------------------------------

So how to test that this works?

Sure, we could go into the installation directory
``SAGE_VENV/lib/python3.9/site-packages/`` and do ``rm -rf
sage/symbolic`` and test that things still work. But that's not a good
way of testing.

Instead, we use a virtual environment in which we only install the
distribution to be tested (and its Python dependencies).

Let's try it out first with the entire Sage library, represented by
the distribution **sagemath-standard**.  Note that after Sage has been
built normally, a set of wheels for all installed Python packages is
available in ``SAGE_VENV/var/lib/sage/wheels/``::

  $ ls venv/var/lib/sage/wheels
  Babel-2.9.1-py2.py3-none-any.whl
  Cython-0.29.24-cp39-cp39-macosx_11_0_x86_64.whl
  Jinja2-2.11.2-py2.py3-none-any.whl
  ...
  sage_conf-9.5b6-py3-none-any.whl
  ...
  scipy-1.7.2-cp39-cp39-macosx_11_0_x86_64.whl
  setuptools-58.2.0-py3-none-any.whl
  ...
  wheel-0.37.0-py2.py3-none-any.whl
  widgetsnbextension-3.5.1-py2.py3-none-any.whl
  zipp-3.5.0-py3-none-any.whl

Note in particular the wheel for **sage-conf**, which provides
configuration variable settings and the connection to the non-Python
packages installed in ``SAGE_LOCAL``.

We can now set up a separate virtual environment, in which we install
these wheels and our distribution to be tested.  This is where
`tox <https://pypi.org/project/tox/>`_
comes into play: It is the standard Python tool for creating
disposable virtual environments for testing.  Every distribution in
``SAGE_ROOT/pkgs/`` provides a configuration file ``tox.ini``.

Following the comments in the file
``SAGE_ROOT/pkgs/sagemath-standard/tox.ini``, we can try the following
command::

  $ ./bootstrap && ./sage -sh -c '(cd pkgs/sagemath-standard && SAGE_NUM_THREADS=16 tox -v -v -v -e sagepython-sagewheels-nopypi)'

This command does not make any changes to the normal installation of
Sage. The virtual environment is created in a subdirectory of
``SAGE_ROOT/pkgs/sagemath-standard-no-symbolics/.tox/``. After the command
finishes, we can start the separate installation of the Sage library
in its virtual environment::

  $ pkgs/sagemath-standard/.tox/sagepython-sagewheels-nopypi/bin/sage

We can also run parts of the testsuite::

  $ pkgs/sagemath-standard/.tox/sagepython-sagewheels-nopypi/bin/sage -tp 4 src/sage/graphs/

The whole ``.tox`` directory can be safely deleted at any time.

We can do the same with other distributions, for example the large
distribution **sagemath-standard-no-symbolics**
(from :trac:`32601`), which is intended to provide
everything that is currently in the standard Sage library, i.e.,
without depending on optional packages, but without the packages
:mod:`sage.symbolic`, :mod:`sage.functions`, :mod:`sage.calculus`, etc.

Again we can run the test with ``tox`` in a separate virtual environment::

  $ ./bootstrap && ./sage -sh -c '(cd pkgs/sagemath-standard-no-symbolics && SAGE_NUM_THREADS=16 tox -v -v -v -e sagepython-sagewheels-nopypi)'

Some small distributions, for example the ones providing the two
lowest levels, `sagemath-objects <https://pypi.org/project/sagemath-objects/>`_
and `sagemath-categories <https://pypi.org/project/sagemath-categories/>`_
(from :trac:`29865`), can be installed and tested
without relying on the wheels from the Sage build::

  $ ./bootstrap && ./sage -sh -c '(cd pkgs/sagemath-objects && SAGE_NUM_THREADS=16 tox -v -v -v -e sagepython)'

This command finds the declared build-time and run-time dependencies
on PyPI, either as source tarballs or as prebuilt wheels, and builds
and installs the distribution
`sagemath-objects <https://pypi.org/project/sagemath-objects/>`_ in a virtual
environment in a subdirectory of ``pkgs/sagemath-objects/.tox``.

Building these small distributions serves as a valuable regression
testsuite.  However, a current issue with both of these distributions
is that they are not separately testable: The doctests for these
modules depend on a lot of other functionality from higher-level parts
of the Sage library.
