.. nodoctest

.. highlight:: shell-session

.. _chapter-tools:

=============================
Development and Testing Tools
=============================

uv
==

`uv <https://docs.astral.sh/uv/>`_ is a versatile tool for 
managing and synchronizing project dependencies. 

The lockfile `uv.lock` in the root captures the exact package versions for 
all systems and ensures consistent, reproducible installations.
It is automatically updated during ``uv`` operations like ``uv add`` 
and ``uv run``, or explicitly with ``uv lock``. 
Moreover, it is periodically updated by `Renovate <https://docs.renovatebot.com/>`_.

.. _section-tools-tox:

Tox
===

`tox <https://tox.readthedocs.io/en/latest/>`_ is a popular package that is
used by a large number of Python projects as the standard entry point
for testing and linting.

Sage includes tox as a standard package and uses it for three purposes:

- For portability testing of the Sage distribution, as we explain in
  :ref:`chapter-portability_testing`.  This is configured in the file
  :sage_root:`tox.ini`.

- For testing modularized distributions of the Sage library. This is configured
  in ``tox.ini`` files in subdirectories of :sage_root:`pkgs/`, such as
  :sage_root:`pkgs/sagemath-standard/tox.ini`. Each distribution's configuration
  defines tox environments for testing the distribution with different Python
  versions and different ways how the dependencies are provided.
  We explain this in :ref:`chapter-modularization`.

- As an entry point for testing and linting of the Sage library, as we describe below.
  This is configured in the file :sage_root:`src/tox.ini`.

The tox configuration :sage_root:`src/tox.ini` can be invoked by using the command
``./sage --tox``.  (If ``tox`` is available in your system installation,
you can just type ``tox`` instead.)

This configuration provides an entry point for various testing/linting methods,
known as "tox environments".  We can type ``./sage --advanced`` to see what is
available::

  $ ./sage --advanced
  SageMath version 9.2
  ...
  Testing files:
  ...
  --tox [options] <files|dirs> -- general entry point for testing
                                  and linting of the Sage library
     -e <envlist>     -- run specific test environments; default:
                         doctest,coverage,startuptime,pycodestyle-minimal,relint,codespell,rst,ruff-minimal
        doctest                -- run the Sage doctester
                                  (same as "sage -t")
        coverage               -- give information about doctest coverage of files
                                  (same as "sage --coverage[all]")
        startuptime            -- display how long each component of Sage takes to start up
                                  (same as "sage --startuptime")
        pycodestyle-minimal    -- check against Sage's minimal style conventions
        relint                 -- check whether some forbidden patterns appear
        codespell              -- check for misspelled words in source code
        rst                    -- validate Python docstrings markup as reStructuredText
        ruff-minimal           -- check against Sage's minimal style conventions
        coverage.py            -- run the Sage doctester with Coverage.py
        coverage.py-html       -- run the Sage doctester with Coverage.py, generate HTML report
        pyright                -- run the static typing checker pyright
        pycodestyle            -- check against the Python style conventions of PEP8
        cython-lint            -- check Cython files for code style
        ruff                   -- check against Python style conventions
     -p auto          -- run test environments in parallel
     --help           -- show tox help


Doctest
=======

The command ``./sage -tox -e doctest`` runs the Sage doctester. This is
equivalent to using the command ``./sage -t``; see :ref:`chapter-doctesting`.

.. NOTE::

   ``doctest`` is a special tox environment that requires that Sage has
   been built already. A virtual environment is created by tox, but
   Sage is invoked in the normal Sage environment.


.. _section-tools-coverage-py:

Doctest with Coverage.py
========================

The command ``./sage -tox -e coverage.py`` runs the Sage doctester
(:ref:`chapter-doctesting`) in the normal Sage environment, but
under the control of
`Coverage.py <https://coverage.readthedocs.io/en/latest/index.html>`_
for code coverage analysis.

If invoked as ``./sage -tox -e coverage.py-html``, additionally a
detailed HTML report is generated.

*Configuration:* ``[coverage:run]`` block in :sage_root:`src/tox.ini`

*Documentation:* https://coverage.readthedocs.io

.. NOTE::

   ``coverage.py`` is a special tox environment that requires that Sage has
   been built already. A virtual environment is created by tox, but the
   **coverage** package is installed into the normal Sage environment, and
   Sage is invoked from there.


.. _section-tools-coverage:

Coverage
========

The command ``./sage -tox -e coverage`` checks that each function has
at least one doctest (typically in an **EXAMPLES** or **TESTS** block,
see :ref:`section-docstring-function`).

Without additional arguments, this command is equivalent to using the
command ``./sage --coverageall`` and gives a short report with a one-line
summary for each module of the Sage library.

If invoked with arguments, for example ``./sage -tox -e coverage
-- src/sage/geometry src/sage/combinat/tableau.py``, it is equivalent to
using the command ``./sage --coverage``, which includes details on
the modules in the given files or directories.

.. NOTE::

   ``coverage`` is a special tox environment that requires that Sage has been
   built already. A virtual environment is created by tox, but
   Sage is invoked in the normal Sage environment.


.. _section-tools-startuptime:

Startuptime
===========

The command ``./sage -tox -e startuptime`` measures the time for loading
each module that is imported during the start up phase of Sage. It is
equivalent to using the command ``./sage --startuptime``.

Without additional arguments, the command gives a short report that lists
the modules with the longest contributions to the overall startup time,
sorted by time.

If invoked with arguments, for example ``sage -tox -e startuptime -- sage.rings
src/sage/geometry/polyhedron``, it provides details on the given modules, packages,
source files, or directories.

.. NOTE::

   ``startuptime`` is a special tox environment that requires that Sage has been
   built already. A virtual environment is created by tox, but
   Sage is invoked in the normal Sage environment.


.. _section-tools-pycodestyle:

Pycodestyle
===========
`Pycodestyle <https://pycodestyle.pycqa.org/en/latest/>`_ (formerly known as pep8)
checks Python code against the style conventions of `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_.
Tox automatically installs pycodestyle in a separate virtual environment
on the first use.

Sage defines two configurations for pycodestyle.  The command ``./sage -tox -e pycodestyle-minimal`` uses
pycodestyle in a minimal configuration.
As of Sage 9.5, the entire Sage library conforms to this configuration::

  $ ./sage -tox -e pycodestyle-minimal -- src/sage/
  pycodestyle-minimal installed: pycodestyle==2.8.0
  pycodestyle-minimal run-test-pre: PYTHONHASHSEED='28778046'
  pycodestyle-minimal run-test: commands[0] | pycodestyle --select E401,E70,W605,E711,E712,E721 sage
  ___________ summary ____________
    pycodestyle-minimal: commands succeeded
    congratulations :)

When preparing a branch for a Sage ticket, developers should verify that ``./sage -tox -e
pycodestyle-minimal`` passes.

The second configuration is used with the command ``./sage -tox -e pycodestyle`` and runs a
more thorough check::

  $ ./sage -tox -e pycodestyle -- src/sage/quadratic_forms/quadratic_form.py
  pycodestyle installed: pycodestyle==2.8.0
  pycodestyle run-test-pre: PYTHONHASHSEED='2520226550'
  pycodestyle run-test: commands[0] | pycodestyle sage/quadratic_forms/quadratic_form.py
  sage/quadratic_forms/quadratic_form.py:135:9: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:163:64: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:165:52: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:173:42: E228 missing whitespace around modulo operator
  ...
  sage/quadratic_forms/quadratic_form.py:1620:9: E266 too many leading '#' for block comment
  sage/quadratic_forms/quadratic_form.py:1621:9: E266 too many leading '#' for block comment
  25      E111 indentation is not a multiple of 4
  2       E117 over-indented
  129     E127 continuation line over-indented for visual indent
  1       E128 continuation line under-indented for visual indent
  4       E201 whitespace after '['
  4       E202 whitespace before ']'
  2       E222 multiple spaces after operator
  7       E225 missing whitespace around operator
  1       E228 missing whitespace around modulo operator
  25      E231 missing whitespace after ','
  1       E262 inline comment should start with '# '
  3       E265 block comment should start with '# '
  62      E266 too many leading '#' for block comment
  2       E272 multiple spaces before keyword
  2       E301 expected 1 blank line, found 0
  17      E303 too many blank lines (2)
  ERROR: InvocationError for command .../pycodestyle sage/quadratic_forms/quadratic_form.py (exited with code 1)
  ___________ summary ____________
  ERROR:   pycodestyle: commands failed

When preparing a branch for a PR that adds new code,
developers should verify that ``./sage -tox -e pycodestyle`` does not
issue warnings for the added code.  This will avoid later cleanup
PRs as the Sage codebase is moving toward full PEP 8 compliance.

On the other hand, it is usually not advisable to mix coding-style
fixes with productive changes on the same PR because this would
makes it harder for reviewers to evaluate the changes.

By passing the options ``--count -qq`` we can reduce the output to
only show the number of style violation warnings.  This can be helpful
for planning work on coding-style clean-up PRs that focus on one
or a few related issues::

  $ ./sage -tox -e pycodestyle -- --count -qq src/sage
  pycodestyle installed: pycodestyle==2.8.0
  pycodestyle run-test-pre: PYTHONHASHSEED='3166223974'
  pycodestyle run-test: commands[0] | pycodestyle --count -qq sage
  557     E111 indentation is not a multiple of 4
  1       E112 expected an indented block
  194     E114 indentation is not a multiple of 4 (comment)
  ...
  7       E743 ambiguous function definition 'l'
  335     W291 trailing whitespace
  4       W292 no newline at end of file
  229     W293 blank line contains whitespace
  459     W391 blank line at end of file
  97797
  ERROR: InvocationError for command .../pycodestyle --count -qq sage (exited with code 1)
  ___________ summary ____________
  ERROR:   pycodestyle: commands failed

*Installation:* (for manual use:) ``pip install -U pycodestyle --user``

*Usage:*

- With tox: See above.

- Manual: Run ``pycodestyle path/to/the/file.py``.

- VS Code: The minimal version of pycodestyle is activated by default in
  :sage_root:`.vscode/settings.json` (the corresponding setting is
  ``"python.linting.pycodestyleEnabled": true``). Note that the
  ``settings.json`` file is not ignored by Git so be aware to keep it in sync
  with the Sage repo on GitHub. For further details, see the
  `official VS Code documentation <https://code.visualstudio.com/docs/python/linting>`__.

*Configuration:* ``[pycodestyle]`` block in :sage_root:`src/tox.ini`

*Documentation:* https://pycodestyle.pycqa.org/en/latest/index.html


.. _section-tools-cython-lint:

Cython-lint
===========

`Cython-lint <https://pypi.org/project/cython-lint/>`_ checks Cython source files
for coding style.


.. _section-tools-ruff:

Ruff
====

`Ruff <https://pypi.org/project/ruff/>`_ is a powerful and fast linter
for Python code, written in Rust.

It comes with a large choice of possible checks, and has the capacity
to fix some of the warnings it emits.

Sage defines two configurations for ruff.  The command ``./sage -tox -e ruff-minimal`` uses
ruff in a minimal configuration. As of Sage 10.3, the entire Sage library conforms to this
configuration. When preparing a Sage PR, developers should verify that
``./sage -tox -e ruff-minimal`` passes.

The second configuration is used with the command ``./sage -tox -e ruff`` and runs a
more thorough check.  When preparing a PR that adds new code,
developers should verify that ``./sage -tox -e ruff`` does not
issue warnings for the added code.  This will avoid later cleanup
PRs as the Sage codebase is moving toward full PEP 8 compliance.

On the other hand, it is usually not advisable to mix coding-style
fixes with productive changes on the same PR because this would
makes it harder for reviewers to evaluate the changes.

.. _section-tools-relint:

Relint
======

`Relint <https://pypi.org/project/relint/>`_ checks all source files for forbidden
text patterns specified by regular expressions.

Our configuration of relint flags some outdated Python constructions, plain TeX
commands when equivalent LaTeX commands are available, common mistakes in
documentation markup, and modularization anti-patterns.

*Configuration:* :sage_root:`src/.relint.yml`

*Documentation:* https://pypi.org/project/relint/


.. _section-tools-codespell:

Codespell
=========
`Codespell <https://pypi.org/project/codespell/>`_ uses a dictionary to check for
misspelled words in source code.

Sage defines a configuration for codespell::

  $ ./sage -tox -e codespell -- src/sage/homology/
  codespell installed: codespell==2.1.0
  codespell run-test-pre: PYTHONHASHSEED='1285169064'
  codespell run-test: commands[0] | codespell '--skip=*.png,*.jpg,*.JPG,*.inv,*.dia,*.pdf,*.ico,*#*,*~*,*.bak,*.orig,*.log,*.sobj,*.tar,*.gz,*.pyc,*.o,*.sws,*.so,*.a,.DS_Store' --skip=doc/ca,doc/de,doc/es,doc/hu,doc/ja,doc/ru,doc/fr,doc/it,doc/pt,doc/tr --skip=src/doc/ca,src/doc/de,src/doc/es,src/doc/hu,src/doc/ja,src/doc/ru,src/doc/fr,src/doc/it,src/doc/pt,src/doc/tr '--skip=.git,.tox,worktree*,dist,upstream,logs,local,cythonized,scripts-3,autom4te.cache,tmp,lib.*,*.egg-info' --dictionary=- --dictionary=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-dictionary.txt --ignore-words=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-ignore.txt sage/homology
  sage/homology/hochschild_complex.py:271: mone ==> mono, money, none
  sage/homology/hochschild_complex.py:277: mone ==> mono, money, none
  sage/homology/hochschild_complex.py:280: mone ==> mono, money, none
  sage/homology/chain_complex.py:2185: mor ==> more
  sage/homology/chain_complex.py:2204: mor ==> more
  sage/homology/chain_complex.py:2210: mor ==> more
  sage/homology/chain_complex.py:2211: mor ==> more
  sage/homology/chain_complex.py:2214: mor ==> more
  sage/homology/chain_complex.py:2215: mor ==> more
  ERROR: InvocationError for command .../codespell '--skip=*.png,...' --dictionary=- --dictionary=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-dictionary.txt --ignore-words=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-ignore.txt sage/homology (exited with code 65)
  ___________ summary ____________
  ERROR:   codespell: commands failed

*Configuration:*

- ``[testenv:codespell]`` block in :sage_root:`src/tox.ini`

- :sage_root:`src/.codespell-dictionary.txt` and :sage_root:`src/.codespell-ignore.txt`


.. _section-tools-pytest:

Pytest
======
`Pytest <https://docs.pytest.org/en/stable/>`_ is a testing framework.
It is included in the Sage distribution as an optional package.

Currently, Sage only makes very limited use of pytest, for testing the
package :mod:`sage.numerical.backends` and some modules in
:mod:`sage.manifolds`.

*Installation:*

- ``./sage -i pytest pytest_xdist``.

*Usage:*

- Tox, Sage doctester: At the end of ``./sage -t`` (or ``./sage --tox -e doctest``), Pytest is automatically invoked.

- Manual: Run ``./sage -pytest path/to/the/test_file.py`` or ``./sage -pytest``
  to run all tests. The additional argument ``-n`` can be used to
  distribute tests across multiple CPUs to speed up test execution.
  For example, ``./sage -pytest -n 4`` will run 4 tests in parallel, while
  ``./sage -pytest -n auto`` will spawn a number of workers processes equal
  to the number of available CPUs.

- VS Code: Install the `Python extension <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`_ and follow the `official VS Code documentation <https://code.visualstudio.com/docs/python/testing>`__.

*Configuration:* :sage_root:`src/conftest.py`

*Documentation:* https://docs.pytest.org/en/stable/index.html


.. _section-tools-pyright:

Pyright
=======
`Pyright <https://github.com/microsoft/pyright>`_ is static type checker.

*Installation:*

- (for manual use:) ``npm install -g pyright``, see `documentation <https://github.com/microsoft/pyright#installation>`__ for details.

*Usage:*

- Tox: Run ``./sage -tox -e pyright path/to/the/file.py``

- Manual: Run ``pyright path/to/the/file.py``. If you want to check the whole Sage library, you most likely run out of memory with the default settings.
  You can use the following command to check the whole library::

    NODE_OPTIONS="--max-old-space-size=8192" pyright

- VS Code: Install the `Pylance <https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance>`__ extension.

*Configuration:* :sage_root:`pyrightconfig.json`

*Documentation:* https://github.com/microsoft/pyright#documentation


.. _section-tools-pyflakes:

Pyflakes
========
`Pyflakes <https://github.com/PyCQA/pyflakes>`_ checks for common coding errors.


.. _section-act:

Act
===

`act <https://github.com/nektos/act>`_ is a tool, written in Go, and using Docker,
to run GitHub Actions locally; in particular, it speeds up developing Actions.
We recommend using ``gh extension`` facility to install ``act``. ::

    [alice@localhost sage]$ gh extension install https://github.com/nektos/gh-act

Extra steps needed for configuration of Docker to run Actions locally can be found on
`act's GitHub <https://github.com/nektos/act>`_

Here we give a very short sampling of ``act``'s capabilities. If you installed standalone
``act``, it should be invoked as ``act``, not as ``gh act``.
After the set up, one can e.g. list all the available linting actions::

    [alice@localhost sage]$ gh act -l | grep lint
    0      lint-pycodestyle        Code style check with pycodestyle                          Lint                                               lint.yml                push,pull_request
    0      lint-relint             Code style check with relint                               Lint                                               lint.yml                push,pull_request
    0      lint-rst                Validate docstring markup as RST                           Lint                                               lint.yml                push,pull_request
    [alice@localhost sage]$

run a particular action ``lint-rst`` ::

    [alice@localhost sage]$ gh act -j lint-rst
    ...

and so on.

By default, ``act`` pulls all the data needed from the next, but it can also cache it,
speeding up repeated runs quite a lot. The following repeats running of ``lint-rst`` using cached data::

    [alice@localhost sage]$ gh act -p false -r -j lint-rst
    [Lint/Validate docstring markup as RST]   Start image=catthehacker/ubuntu:act-latest
    ...
    | rst: commands[0] /home/alice/work/software/sage/src> flake8 --select=RST
    |   rst: OK (472.60=setup[0.09]+cmd[472.51] seconds)
    |   congratulations :) (474.10 seconds)
    ...
    [Lint/Validate docstring markup as RST]     Success - Main Lint using tox -e rst
    [Lint/Validate docstring markup as RST]  Run Post Set up Python
    [Lint/Validate docstring markup as RST]     docker exec cmd=[node /var/run/act/actions/actions-setup-python@v4/dist/cache-save/index.js] user= workdir=
    [Lint/Validate docstring markup as RST]     Success - Post Set up Python
    [Lint/Validate docstring markup as RST]   Job succeeded

Here ``-p false`` means using already pulled Docker images, and ``-r`` means do not remove Docker images
after a successful run which used them. This, and many more details, can be found by running ``gh act -h``, as well
as reading ``act``'s documentation.

.. This section is a stub.
   More Sage-specfic details for using ``act`` should be added. PRs welcome!

